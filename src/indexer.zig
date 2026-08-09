//! FASTA indexing through one bounded reader for `.zfi` and `.fai` output.
//!
//! The parser skips empty lines for geometry, sets `seq_offset` to the first base-bearing
//! line, invents a separator when the final line has no newline, and rejects names longer
//! than `MAX_INDEX_NAME_LEN`. Proven dense blocks update the same `ChunkParseState`; the
//! first anomaly resumes the general line machine at the exact unconsumed byte.
//!
//! Duplicate-name filtering is first-wins and collision-safe: lookup may use a hash, but
//! identity is always full-string equality (`NameDedup`).

const std = @import("std");
const builtin = @import("builtin");
const index_format = @import("index_format.zig");

const IndexRecord = index_format.IndexRecord;
const ZfiHeader = index_format.ZfiHeader;
const ZFI_MAGIC = index_format.ZFI_MAGIC;

const SIMD_CHUNK_SIZE = 32;
const SimdVec = @Vector(SIMD_CHUNK_SIZE, u8);
const STRIDE_VECTOR_LEN = std.simd.suggestVectorLength(u8) orelse 1;
const STRIDE_VALIDATION_BLOCK_LINES = 256;
pub const INDEX_READ_BUFFER_SIZE = 1 * 1024 * 1024;
const FILE_IO_BUF_SIZE = 8 * 1024;
const INDEX_OUTPUT_BUFFER_SIZE = 64 * 1024;
const INDEX_PATH_BUFFER_SIZE = std.Io.Dir.max_path_bytes;
// Catalogs start modestly and grow with records or names, never sequence payload bytes.
const ZFI_INITIAL_RECORD_CAPACITY = 4096;
const DEDUP_INITIAL_CAPACITY: u32 = 16384;
const FAI_SPOOL_NAME_BUFFER_SIZE = 64;
const FAI_U64_DECIMAL_DIGITS = 20;
const FAI_SUFFIX_BUFFER_SIZE = 4 * (1 + FAI_U64_DECIMAL_DIGITS) + 1;
const FAI_SPOOL_CREATE_ATTEMPTS = 16;
pub const MAX_INDEX_NAME_LEN = std.math.maxInt(u16);

pub const ScanOptions = struct {
    enable_dedup: bool,
    max_name_len: ?usize = null,
    collect_side_tables: bool = false,
    require_initial_header: bool = false,
};

pub const IndexOptions = struct {
    emit_fai: bool = false,
    enable_dedup: bool = true,
};

// Offset into a name blob for `.zfi` catalog dedup (one store for map and output).
const NameBlobRef = struct { offset: u32, len: u32 };

fn NameDedupWith(comptime Context: type) type {
    return struct {
        map: std.HashMap([]const u8, void, Context, std.hash_map.default_max_load_percentage),
        key_arena: ?std.heap.ArenaAllocator = null,
        // When set, `.zfi` dedup stores `NameBlobRef` keys into this blob (not `map`).
        name_blob: ?*std.ArrayList(u8) = null,
        ref_map: ?BlobRefMap = null,
        // After a successful first-seen `observe` in blob mode, the embedded span.
        last_offset: u32 = 0,
        last_len: u16 = 0,

        const BlobHashContext = struct {
            blob: *std.ArrayList(u8),
            pub fn hash(self: @This(), key: NameBlobRef) u64 {
                return std.hash_map.hashString(self.blob.items[key.offset..][0..key.len]);
            }
            pub fn eql(self: @This(), a: NameBlobRef, b: NameBlobRef) bool {
                if (a.len != b.len) return false;
                return std.mem.eql(
                    u8,
                    self.blob.items[a.offset..][0..a.len],
                    self.blob.items[b.offset..][0..b.len],
                );
            }
        };

        const BlobRefMap = std.HashMap(
            NameBlobRef,
            void,
            BlobHashContext,
            std.hash_map.default_max_load_percentage,
        );

        const SliceLookup = struct {
            blob: *std.ArrayList(u8),
            pub fn hash(_: @This(), name: []const u8) u64 {
                return std.hash_map.hashString(name);
            }
            pub fn eql(self: @This(), name: []const u8, ref: NameBlobRef) bool {
                if (name.len != ref.len) return false;
                return std.mem.eql(u8, name, self.blob.items[ref.offset..][0..ref.len]);
            }
        };

        const Self = @This();

        pub fn init(allocator: std.mem.Allocator) Self {
            return .{
                .map = std.HashMap([]const u8, void, Context, std.hash_map.default_max_load_percentage).init(allocator),
                .key_arena = std.heap.ArenaAllocator.init(allocator),
            };
        }

        // `.zfi`: one catalog in `blob` for both dedup and embedded names.
        fn initOwningBlob(allocator: std.mem.Allocator, blob: *std.ArrayList(u8)) Self {
            return .{
                .map = std.HashMap([]const u8, void, Context, std.hash_map.default_max_load_percentage).init(allocator),
                .name_blob = blob,
                .ref_map = BlobRefMap.initContext(allocator, .{ .blob = blob }),
            };
        }

        fn bindsToBlob(self: *const Self) bool {
            return self.name_blob != null;
        }

        fn ensureCapacity(self: *Self, min_cap: u32) !void {
            if (self.ref_map) |*rm| {
                try rm.ensureTotalCapacity(min_cap);
                return;
            }
            try self.map.ensureTotalCapacity(min_cap);
        }

        pub fn deinit(self: *Self) void {
            if (self.ref_map) |*rm| {
                rm.deinit();
                self.ref_map = null;
            }
            // Map entries borrow arena memory; free the table first, then the arena.
            self.map.deinit();
            if (self.key_arena) |*arena| arena.deinit();
            self.key_arena = null;
            self.name_blob = null;
        }

        // Returns true when `name` was already observed (caller should skip the record).
        pub fn observe(self: *Self, name: []const u8) !bool {
            if (self.name_blob) |blob| {
                var rm = &(self.ref_map orelse unreachable);
                const lookup = SliceLookup{ .blob = blob };
                const gop = try rm.getOrPutAdapted(name, lookup);
                if (gop.found_existing) {
                    self.last_offset = gop.key_ptr.offset;
                    self.last_len = @intCast(gop.key_ptr.len);
                    return true;
                }
                const off: u32 = @intCast(blob.items.len);
                blob.appendSlice(rm.allocator, name) catch |err| {
                    _ = rm.removeByPtr(gop.key_ptr);
                    return err;
                };
                gop.key_ptr.* = .{ .offset = off, .len = @intCast(name.len) };
                self.last_offset = off;
                self.last_len = @intCast(name.len);
                return false;
            }
            const arena = &(self.key_arena orelse unreachable);
            const gop = try self.map.getOrPut(name);
            if (gop.found_existing) return true;
            // Replace the temporary lookup key with an arena-owned copy (stable until deinit).
            const owned = arena.allocator().dupe(u8, name) catch |err| {
                _ = self.map.remove(name);
                return err;
            };
            gop.key_ptr.* = owned;
            return false;
        }
    };
}

/// First-wins name set whose decisions use full-string equality despite hash collisions.
/// FAI owns stable keys in a child arena; `.zfi` shares names with its output blob.
pub const NameDedup = NameDedupWith(std.hash_map.StringContext);

/// One parsed FASTA record passed to reader emit callbacks.
pub const FastaRecordEmit = struct {
    record: IndexRecord,
    name: []const u8,
    uses_uniform_formula: bool,
    /// Incremental side-table bytes for `.zfi` (line count plus `SideTableLine` rows).
    side_table: []const u8 = &.{},
    /// When true, `record` already has blob-relative `name_offset`/`name_len` and `NAME_IN_ZFI_FLAG`
    /// (`.zfi` blob-backed `NameDedup`); skip `embedZfiName`.
    name_embedded: bool = false,
};

pub const ZfiIndex = struct {
    records: std.ArrayList(IndexRecord),
    side_tables: std.ArrayList(u8),
    name_blob: std.ArrayList(u8),

    pub fn deinit(self: *ZfiIndex, allocator: std.mem.Allocator) void {
        self.records.deinit(allocator);
        self.side_tables.deinit(allocator);
        self.name_blob.deinit(allocator);
    }
};

const LineMetrics = struct {
    seq_len: u64 = 0,
    line_bases: u32 = 0,
    line_bytes: u32 = 0,
    is_uniform_width: bool = true,
    line_count: u64 = 0,
    // First base line: every content byte was a base (no trailing space/tab) and the
    // on-wire separator was LF or CRLF. Needed so `AAAA \\n` is not treated as CRLF-dense.
    first_content_dense: bool = false,
};

const SequenceLengthInfo = struct {
    seq_len: u64,
    uses_uniform_formula: bool,
};

fn countBases(data: []const u8) u64 {
    var count: u64 = 0;
    var pos: usize = 0;
    while (pos + SIMD_CHUNK_SIZE <= data.len) {
        const chunk: SimdVec = data[pos..][0..SIMD_CHUNK_SIZE].*;
        const space_char: SimdVec = @splat(' ');
        var chunk_count: u32 = 0;
        inline for (0..SIMD_CHUNK_SIZE) |j| {
            if (chunk[j] > space_char[j]) chunk_count += 1;
        }
        count += chunk_count;
        pos += SIMD_CHUNK_SIZE;
    }
    while (pos < data.len) {
        if (data[pos] > ' ') count += 1;
        pos += 1;
    }
    return count;
}

fn hasLineSeparatorAt(data: []const u8, pos: usize, sep_len: u32) bool {
    return switch (sep_len) {
        1 => pos < data.len and data[pos] == '\n',
        2 => pos + 1 < data.len and data[pos] == '\r' and data[pos + 1] == '\n',
        else => false,
    };
}

fn lineSliceHasWhitespace(data: []const u8, start: usize, end: usize) bool {
    var pos = start;
    while (pos + SIMD_CHUNK_SIZE <= end) {
        const chunk: SimdVec = data[pos..][0..SIMD_CHUNK_SIZE].*;
        const space_char: SimdVec = @splat(' ');
        const mask = chunk <= space_char;
        if (@reduce(.Or, mask)) return true;
        pos += SIMD_CHUNK_SIZE;
    }
    while (pos < end) : (pos += 1) {
        if (data[pos] <= ' ') return true;
    }
    return false;
}

// Bases in one sequence line's content slice (no trailing `\n`/`\r`).
// Dense lines use length only; other lines fall back to `countBases`.
fn countBasesInLineSlice(data: []const u8, start: usize, end: usize) u32 {
    if (start >= end) return 0;
    if (!lineSliceHasWhitespace(data, start, end)) {
        return @intCast(end - start);
    }
    return @intCast(countBases(data[start..end]));
}

fn firstLineIsDense(line_bases: u32, line_bytes: u32) bool {
    if (line_bases == 0 or line_bytes <= line_bases) return false;
    const sep_len = line_bytes - line_bases;
    // Count-only check: callers that saw the bytes must also require content density
    // (`LineMetrics.first_content_dense`) so space+LF is not mistaken for CRLF.
    return (sep_len == 1 or sep_len == 2) and line_bases + sep_len == line_bytes;
}

fn countWhitespace(comptime vector_len: usize, data: []const u8) usize {
    const Vec = @Vector(vector_len, u8);
    var count: usize = 0;
    var pos: usize = 0;
    while (pos + vector_len <= data.len) : (pos += vector_len) {
        const chunk: Vec = data[pos..][0..vector_len].*;
        count += std.simd.countTrues(chunk <= @as(Vec, @splat(' ')));
    }
    while (pos < data.len) : (pos += 1) {
        count += @intFromBool(data[pos] <= ' ');
    }
    return count;
}

fn validatedStrideBlock(
    comptime vector_len: usize,
    data: []const u8,
    start: usize,
    line_count: usize,
    line_bases: usize,
    line_bytes: usize,
) usize {
    const sep_len = line_bytes - line_bases;
    const block_end = start + line_count * line_bytes;
    var line_index: usize = 0;
    while (line_index < line_count) : (line_index += 1) {
        const line_start = start + line_index * line_bytes;
        if (data[line_start] == '>') return line_index;
        if (!hasLineSeparatorAt(data, line_start + line_bases, @intCast(sep_len))) {
            return line_index;
        }
    }
    if (countWhitespace(vector_len, data[start..block_end]) == line_count * sep_len) {
        return line_count;
    }

    line_index = 0;
    while (line_index < line_count) : (line_index += 1) {
        const line_start = start + line_index * line_bytes;
        if (countWhitespace(vector_len, data[line_start .. line_start + line_bases]) != 0) {
            return line_index;
        }
    }
    return 0;
}

fn validatedStrideRun(
    comptime vector_len: usize,
    data: []const u8,
    start: usize,
    line_bases_u32: u32,
    line_bytes_u32: u32,
) usize {
    const line_bases: usize = line_bases_u32;
    const line_bytes: usize = line_bytes_u32;
    if (line_bases == 0 or line_bytes <= line_bases) return 0;
    const sep_len = line_bytes - line_bases;
    if (sep_len != 1 and sep_len != 2) return 0;

    var cursor = start;
    var accepted_lines: usize = 0;
    while (data.len - cursor >= line_bytes) {
        const available_lines = (data.len - cursor) / line_bytes;
        const block_lines = @min(available_lines, STRIDE_VALIDATION_BLOCK_LINES);
        const valid_lines = validatedStrideBlock(
            vector_len,
            data,
            cursor,
            block_lines,
            line_bases,
            line_bytes,
        );
        accepted_lines += valid_lines;
        if (valid_lines != block_lines) return accepted_lines;
        cursor += block_lines * line_bytes;
    }
    return accepted_lines;
}

fn findNextNewline(data: []const u8, start: usize, end: usize) usize {
    const relative = std.mem.findScalar(u8, data[start..end], '\n') orelse return end;
    return start + relative;
}

const LineMetricsBuilder = struct {
    metrics: LineMetrics = .{},
    first_actual_bytes: u32 = 0,
    first_content_len: u32 = 0,
    pending_bases: u32 = 0,
    pending_actual_bytes: u32 = 0,
    pending_content_len: u32 = 0,
    have_pending: bool = false,
    // Blank seen after at least one base line; applied only when another base line follows.
    // Trailing blanks before EOF or the next header do not change record geometry.
    blank_after_bases: bool = false,

    fn ingestLine(self: *LineMetricsBuilder, bases: u32, actual_bytes: u32, has_lf: bool, content_len: u32) void {
        if (bases == 0) {
            // Leading blanks are ignored. Trailing blanks remain pending until a later base
            // line proves that they were interior.
            if (self.metrics.line_count > 0) {
                self.blank_after_bases = true;
            }
            return;
        }
        if (self.blank_after_bases) {
            self.metrics.is_uniform_width = false;
            self.blank_after_bases = false;
        }
        self.metrics.seq_len += bases;
        self.metrics.line_count += 1;

        if (self.metrics.line_count == 1) {
            self.first_actual_bytes = actual_bytes;
            self.first_content_len = content_len;
            self.metrics.line_bases = bases;
            self.metrics.line_bytes = if (has_lf) actual_bytes else actual_bytes + 1;
            // Dense FAI geometry: content is all bases; trailing bytes are only LF or CRLF.
            const sep_len = if (has_lf) actual_bytes -| content_len else @as(u32, 1);
            self.metrics.first_content_dense = bases == content_len and (sep_len == 1 or sep_len == 2);
        } else if (self.have_pending) {
            // Interior full wraps must match on-wire bytes. LF vs CRLF is a real break.
            if (self.pending_bases != self.metrics.line_bases or
                self.pending_actual_bytes != self.first_actual_bytes)
            {
                self.metrics.is_uniform_width = false;
            }
        }

        self.pending_bases = bases;
        self.pending_actual_bytes = actual_bytes;
        self.pending_content_len = content_len;
        self.have_pending = true;
    }

    fn finish(self: *LineMetricsBuilder) LineMetrics {
        // Trailing blanks do not break uniformity.
        self.blank_after_bases = false;
        if (self.metrics.line_count > 1 and self.have_pending) {
            // CRLF vs LF on the final wrap does not break uniformity when content width
            // does not grow.
            if (self.pending_bases > self.metrics.line_bases or
                self.pending_content_len > self.first_content_len)
            {
                self.metrics.is_uniform_width = false;
            }
        }
        return self.metrics;
    }
};

fn sequenceLengthInfoFromMetrics(metrics: LineMetrics) SequenceLengthInfo {
    if (metrics.seq_len == 0) return .{ .seq_len = 0, .uses_uniform_formula = false };
    if (!metrics.is_uniform_width or metrics.line_bases == 0 or !metrics.first_content_dense) {
        return .{ .seq_len = metrics.seq_len, .uses_uniform_formula = false };
    }
    return .{
        .seq_len = metrics.seq_len,
        .uses_uniform_formula = firstLineIsDense(metrics.line_bases, metrics.line_bytes),
    };
}

fn finalizeSideTableOffsets(index: *ZfiIndex) !void {
    const side_table_base = @sizeOf(index_format.ZfiHeader) + index.records.items.len * @sizeOf(IndexRecord);
    for (index.records.items) |*rec| {
        if (!rec.isUniformWidth()) {
            const local_offset = rec.sideTableOffset();
            try rec.markNonUniform(side_table_base + local_offset);
        }
    }
}

fn embedZfiName(
    rec: *IndexRecord,
    name: []const u8,
    name_blob: *std.ArrayList(u8),
    allocator: std.mem.Allocator,
) !void {
    rec.name_offset = name_blob.items.len;
    rec.name_len = @intCast(name.len);
    try name_blob.appendSlice(allocator, name);
    rec._pad[0] |= index_format.NAME_IN_ZFI_FLAG;
}

pub fn scanZfiReader(
    reader: *std.Io.Reader,
    read_buf: []u8,
    enable_dedup: bool,
    allocator: std.mem.Allocator,
) !ZfiIndex {
    return scanZfiReaderWithOptions(
        reader,
        read_buf,
        .{ .enable_dedup = enable_dedup },
        allocator,
    );
}

fn scanZfiReaderWithOptions(
    reader: *std.Io.Reader,
    read_buf: []u8,
    options: ScanOptions,
    allocator: std.mem.Allocator,
) !ZfiIndex {
    var index = ZfiIndex{
        .records = .empty,
        .side_tables = .empty,
        .name_blob = .empty,
    };
    errdefer index.deinit(allocator);
    // Bounds fixed catalog startup cost; larger catalogs grow geometrically.
    try index.records.ensureTotalCapacity(allocator, ZFI_INITIAL_RECORD_CAPACITY);

    const Ctx = struct {
        index: *ZfiIndex,
        allocator: std.mem.Allocator,

        fn emit(ctx: *@This(), emit_info: FastaRecordEmit) !void {
            var rec = emit_info.record;
            if (!emit_info.uses_uniform_formula) {
                const local_offset = ctx.index.side_tables.items.len;
                try ctx.index.side_tables.appendSlice(ctx.allocator, emit_info.side_table);
                try rec.markNonUniform(local_offset);
            }
            if (!emit_info.name_embedded) {
                try embedZfiName(&rec, emit_info.name, &ctx.index.name_blob, ctx.allocator);
            }
            try ctx.index.records.append(ctx.allocator, rec);
        }
    };

    var ctx = Ctx{ .index = &index, .allocator = allocator };
    // One catalog avoids a second copy of embedded names.
    var blob_dedup: ?NameDedup = null;
    if (options.enable_dedup) {
        blob_dedup = NameDedup.initOwningBlob(allocator, &index.name_blob);
        try blob_dedup.?.ensureCapacity(DEDUP_INITIAL_CAPACITY);
    }
    defer if (blob_dedup) |*seen| seen.deinit();
    const external_dedup: ?*NameDedup = if (blob_dedup) |*s| s else null;
    _ = try scanFastaReader(
        reader,
        read_buf,
        .{
            .enable_dedup = options.enable_dedup,
            .collect_side_tables = true,
            .require_initial_header = options.require_initial_header,
        },
        allocator,
        external_dedup,
        &ctx,
        Ctx.emit,
    );
    try finalizeSideTableOffsets(&index);
    return index;
}

pub fn scanZfiData(
    data: []const u8,
    enable_dedup: bool,
    allocator: std.mem.Allocator,
) !ZfiIndex {
    var r = std.Io.Reader.fixed(data);
    var read_buf: [INDEX_READ_BUFFER_SIZE]u8 = undefined;
    return scanZfiReader(&r, &read_buf, enable_dedup, allocator);
}

// Single on-disk `.zfi` serialization path:
// header, records, side tables, name blob, `ZFID` source identity, `ZFNM` footer.
fn writeZfiIndex(
    writer: *std.Io.Writer,
    index: *const ZfiIndex,
    source_size: u64,
    source_mtime_ns: u64,
) !void {
    const header = ZfiHeader{
        .magic = ZFI_MAGIC,
        .record_count = @intCast(index.records.items.len),
        .source_size = source_size,
    };
    try writer.writeAll(std.mem.asBytes(&header));
    try writer.writeAll(std.mem.sliceAsBytes(index.records.items));
    try writer.writeAll(index.side_tables.items);
    try writer.writeAll(index.name_blob.items);
    // Identity before the name footer so legacy loaders that seek `ZFNM` at EOF still
    // find the footer; they ignore the orphan `ZFID` gap between blob and footer.
    try writer.writeAll(&index_format.encodeZfiSourceId(source_mtime_ns));
    try writer.writeAll(&index_format.encodeZfiNameFooter(index.name_blob.items.len));
}

pub fn writeZfiIndexFile(
    io: std.Io,
    path: []const u8,
    index: *const ZfiIndex,
    source_size: u64,
    source_mtime_ns: u64,
) !void {
    const file = try std.Io.Dir.cwd().createFile(io, path, .{});
    defer file.close(io);

    var file_buf: [INDEX_OUTPUT_BUFFER_SIZE]u8 = undefined;
    var file_fw = file.writer(io, &file_buf);
    try writeZfiIndex(&file_fw.interface, index, source_size, source_mtime_ns);
    try file_fw.flush();
}

/// Caller owns the returned production-layout bytes.
pub fn zfiIndexToBytes(
    index: *const ZfiIndex,
    source_size: u64,
    source_mtime_ns: u64,
    allocator: std.mem.Allocator,
) ![]u8 {
    var aw: std.Io.Writer.Allocating = .init(allocator);
    errdefer aw.deinit();
    try writeZfiIndex(&aw.writer, index, source_size, source_mtime_ns);
    return try aw.toOwnedSlice();
}

const ChunkParseState = struct {
    allocator: std.mem.Allocator,
    collect_side_tables: bool = false,
    active: bool = false,
    in_header: bool = false,
    parsing_name: bool = false,
    at_line_start: bool = true,
    header_start_offset: u64 = 0,
    name: std.ArrayList(u8) = .empty,
    seq_offset: u64 = 0,
    seq_offset_set: bool = false,
    line_builder: LineMetricsBuilder = .{},
    in_pending_line: bool = false,
    pending_line_bases: u32 = 0,
    pending_line_bytes: u32 = 0,
    // Pre-LF byte count. The LF is counted when the line commits.
    pending_pre_lf_bytes: u32 = 0,
    pending_pre_lf_ends_with_cr: bool = false,
    pending_line_file_offset: u64 = 0,
    side_table: std.ArrayList(u8) = .empty,
    side_base_start: u64 = 0,
    side_table_active: bool = false,
    deferred_side_line: ?struct {
        bases: u32,
        actual_bytes: u32,
        line_file_offset: u64,
    } = null,
    // A short final line remains FAI-uniform. Side-table rows stay deferred until a later
    // base line proves the short line was interior.
    deferred_short_tail: ?struct {
        bases: u32,
        actual_bytes: u32,
        line_file_offset: u64,
    } = null,

    fn ensureSideTableHeader(self: *ChunkParseState) !void {
        if (self.side_table.items.len >= @sizeOf(u64)) return;
        var line_count: u64 = 0;
        try self.side_table.appendSlice(self.allocator, std.mem.asBytes(&line_count));
    }

    // Materialize uniform stride rows only after the record proves non-uniform.
    fn activateSideTable(self: *ChunkParseState, lines_to_materialize: u64) !void {
        if (self.side_table_active) return;
        self.side_table_active = true;
        try self.ensureSideTableHeader();
        if (lines_to_materialize == 0) {
            self.deferred_side_line = null;
            return;
        }
        const first = self.deferred_side_line orelse return error.MissingSideTable;
        var i: u64 = 0;
        while (i < lines_to_materialize) : (i += 1) {
            const off = first.line_file_offset + i * @as(u64, first.actual_bytes);
            try self.appendSideTableLine(first.bases, first.actual_bytes, off);
        }
        self.deferred_side_line = null;
    }

    fn startHeader(self: *ChunkParseState, header_start_offset: u64) !void {
        self.active = true;
        self.in_header = true;
        self.parsing_name = true;
        self.at_line_start = false;
        self.header_start_offset = header_start_offset;
        self.name.clearRetainingCapacity();
        self.seq_offset = 0;
        self.seq_offset_set = false;
        self.line_builder = .{};
        self.in_pending_line = false;
        self.pending_line_bases = 0;
        self.pending_line_bytes = 0;
        self.pending_pre_lf_bytes = 0;
        self.pending_pre_lf_ends_with_cr = false;
        self.pending_line_file_offset = 0;
        self.side_base_start = 0;
        self.side_table.clearRetainingCapacity();
        self.side_table_active = false;
        self.deferred_side_line = null;
        self.deferred_short_tail = null;
    }

    fn appendSideTableLine(
        self: *ChunkParseState,
        bases: u32,
        actual_bytes: u32,
        line_file_offset: u64,
    ) !void {
        const entry = index_format.SideTableLine{
            .base_start = self.side_base_start,
            .byte_offset = line_file_offset,
            .line_bytes = actual_bytes,
            .line_bases = bases,
        };
        try self.side_table.appendSlice(self.allocator, std.mem.asBytes(&entry));
        self.side_base_start += bases;
    }

    // Flush a deferred short tail once the current line proves it was interior.
    fn flushDeferredShortTailThen(
        self: *ChunkParseState,
        metrics_line_count: u64,
        bases: u32,
        actual_bytes: u32,
        line_file_offset: u64,
    ) !void {
        const short = self.deferred_short_tail orelse return error.MissingSideTable;
        self.deferred_short_tail = null;
        // Uniform stride rows before the short line, then short, then current.
        if (metrics_line_count < 2) return error.MissingSideTable;
        try self.activateSideTable(metrics_line_count - 2);
        try self.appendSideTableLine(short.bases, short.actual_bytes, short.line_file_offset);
        try self.appendSideTableLine(bases, actual_bytes, line_file_offset);
    }

    fn recordSequenceLine(
        self: *ChunkParseState,
        bases: u32,
        actual_bytes: u32,
        has_lf: bool,
        content_len: u32,
        line_file_offset: u64,
    ) !void {
        self.line_builder.ingestLine(bases, actual_bytes, has_lf, content_len);

        if (!self.collect_side_tables) return;

        // Blank lines never become side-table rows.
        // Interior blanks flip uniformity when the next base line arrives.
        if (bases == 0) return;

        const metrics = self.line_builder.metrics;
        const matches_stride = metrics.line_bases > 0 and
            bases == metrics.line_bases and
            actual_bytes == self.line_builder.first_actual_bytes;
        const dense = metrics.first_content_dense and firstLineIsDense(metrics.line_bases, metrics.line_bytes);
        const can_use_formula = metrics.is_uniform_width and matches_stride and dense;

        if (can_use_formula) {
            // The first line anchors reconstruction if a later line forces a side table.
            if (metrics.line_count == 1) {
                self.deferred_side_line = .{
                    .bases = bases,
                    .actual_bytes = actual_bytes,
                    .line_file_offset = line_file_offset,
                };
            }
            return;
        }

        // A later base line proves the deferred short line was interior.
        if (self.deferred_short_tail != null) {
            try self.flushDeferredShortTailThen(metrics.line_count, bases, actual_bytes, line_file_offset);
            return;
        }

        // Short-or-equal mismatch: same rule as `LineMetricsBuilder.finish` (only a *longer*
        // final line flips non-uniform). Defer materialization so Genome-scale uniform
        // records with a short final wrap do not allocate O(lines) side-table RAM.
        const short_ok_if_final = metrics.is_uniform_width and
            metrics.line_bases > 0 and
            dense and
            bases <= metrics.line_bases and
            content_len <= self.line_builder.first_content_len and
            actual_bytes <= self.line_builder.first_actual_bytes + 1; // allow final CRLF vs LF (+1)
        if (short_ok_if_final) {
            self.deferred_short_tail = .{
                .bases = bases,
                .actual_bytes = actual_bytes,
                .line_file_offset = line_file_offset,
            };
            return;
        }

        // Wider line or other hard break: reconstruct the uniform prefix, then append.
        try self.activateSideTable(metrics.line_count - 1);
        try self.appendSideTableLine(bases, actual_bytes, line_file_offset);
    }

    fn noteSeqOffset(self: *ChunkParseState, line_file_offset: u64) void {
        if (self.seq_offset_set) return;
        // `seq_offset` is the first base-bearing line start (side-table invariant).
        self.seq_offset = line_file_offset;
        self.seq_offset_set = true;
    }

    // True when further body wraps can use fixed `line_bytes` strides.
    fn strideEligible(self: *const ChunkParseState) bool {
        if (!self.active or self.in_header or self.in_pending_line) return false;
        if (self.side_table_active or self.deferred_short_tail != null) return false;
        // Interior blank must hit `ingestLine` so uniformity flips.
        if (self.line_builder.blank_after_bases) return false;
        const m = self.line_builder.metrics;
        if (m.line_count == 0 or !m.is_uniform_width) return false;
        if (!m.first_content_dense or !firstLineIsDense(m.line_bases, m.line_bytes)) return false;
        // Stride uses the first line's on-wire byte length (includes LF/CRLF).
        if (self.line_builder.first_actual_bytes != m.line_bytes) return false;
        // FAI skips `deferred_short_tail`, so a mid-body CRLF vs LF mismatch lives only in
        // `pending_*` until the next `ingestLine`. Do not stride over that next line or the
        // pending comparison (and uniformity flip) is skipped.
        if (self.line_builder.have_pending and
            (self.line_builder.pending_bases != m.line_bases or
                self.line_builder.pending_actual_bytes != self.line_builder.first_actual_bytes))
        {
            return false;
        }
        return true;
    }

    // Bulk-account full wraps only after `strideEligible` and block validation prove them.
    fn ingestStrideRun(self: *ChunkParseState, n_lines: u64) void {
        std.debug.assert(n_lines >= 1);
        std.debug.assert(self.line_builder.metrics.line_count >= 1);
        const bases = self.line_builder.metrics.line_bases;
        const line_bytes = self.line_builder.metrics.line_bytes;
        self.line_builder.metrics.seq_len += n_lines * @as(u64, bases);
        self.line_builder.metrics.line_count += n_lines;
        self.line_builder.pending_bases = bases;
        self.line_builder.pending_actual_bytes = line_bytes;
        self.line_builder.pending_content_len = self.line_builder.first_content_len;
        self.line_builder.have_pending = true;
        self.at_line_start = true;
    }

    fn ingestSequenceLine(
        self: *ChunkParseState,
        data: []const u8,
        line_start: usize,
        line_end: usize,
        has_lf: bool,
        line_file_offset: u64,
    ) !void {
        var content_end = line_end;
        if (content_end > line_start and data[content_end - 1] == '\r') {
            content_end -= 1;
        }

        const actual_bytes: u32 = @intCast((if (has_lf) line_end + 1 else line_end) - line_start);
        const content_len: u32 = @intCast(content_end - line_start);
        const bases = countBasesInLineSlice(data, line_start, content_end);
        if (bases > 0) self.noteSeqOffset(line_file_offset);

        try self.recordSequenceLine(bases, actual_bytes, has_lf, content_len, line_file_offset);
        self.at_line_start = has_lf;
        self.in_pending_line = false;
        self.pending_line_bases = 0;
        self.pending_line_bytes = 0;
        self.pending_pre_lf_bytes = 0;
        self.pending_pre_lf_ends_with_cr = false;
    }

    fn appendPendingFragment(self: *ChunkParseState, fragment: []const u8, fragment_file_offset: u64) !void {
        if (!self.in_pending_line) {
            self.in_pending_line = true;
            self.pending_line_bases = 0;
            self.pending_line_bytes = 0;
            self.pending_pre_lf_bytes = 0;
            self.pending_pre_lf_ends_with_cr = false;
            self.pending_line_file_offset = fragment_file_offset;
        }
        const fragment_bases = countBasesInLineSlice(fragment, 0, fragment.len);
        const next_bases = std.math.add(u32, self.pending_line_bases, fragment_bases) catch
            return error.SequenceLineTooLong;
        const next_bytes = std.math.add(u32, self.pending_line_bytes, @intCast(fragment.len)) catch
            return error.SequenceLineTooLong;
        self.pending_line_bases = next_bases;
        self.pending_line_bytes = next_bytes;
        self.pending_pre_lf_bytes = next_bytes;
        if (fragment.len > 0) {
            self.pending_pre_lf_ends_with_cr = fragment[fragment.len - 1] == '\r';
        }
    }

    fn commitPendingLine(self: *ChunkParseState, has_lf: bool) !void {
        if (!self.in_pending_line and self.pending_line_bytes == 0) return;
        if (!has_lf and self.pending_line_bytes == std.math.maxInt(u32)) {
            return error.SequenceLineTooLong;
        }
        if (self.pending_line_bases > 0) self.noteSeqOffset(self.pending_line_file_offset);
        var content_len = if (has_lf) self.pending_pre_lf_bytes else self.pending_line_bytes;
        if (content_len > 0 and self.pending_pre_lf_ends_with_cr) {
            content_len -= 1;
        }
        try self.recordSequenceLine(
            self.pending_line_bases,
            self.pending_line_bytes,
            has_lf,
            content_len,
            self.pending_line_file_offset,
        );
        self.at_line_start = has_lf;
        self.in_pending_line = false;
        self.pending_line_bases = 0;
        self.pending_line_bytes = 0;
        self.pending_pre_lf_bytes = 0;
        self.pending_pre_lf_ends_with_cr = false;
    }

    fn finalizeHeader(self: *ChunkParseState) void {
        // Dedup runs at emit time (after seq_len > 0), so empty sequences never claim a
        // name. Blob-backed `.zfi` also embeds only kept records.
        self.in_header = false;
        self.at_line_start = true;
    }

    fn emitRecordIfReady(
        self: *ChunkParseState,
        ctx: anytype,
        seen_names: ?*NameDedup,
        comptime emitRecord: fn (@TypeOf(ctx), FastaRecordEmit) anyerror!void,
        record_count: *u32,
    ) !void {
        try self.commitPendingLine(false);
        const metrics = self.line_builder.finish();
        const seq_info = sequenceLengthInfoFromMetrics(metrics);
        if (seq_info.seq_len == 0) {
            self.deferred_short_tail = null;
            return;
        }

        var name_embedded = false;
        var name_offset: u64 = self.header_start_offset + 1;
        var out_name_len: u16 = @intCast(self.name.items.len);
        var name_pad: [6]u8 = .{0} ** 6;
        if (seen_names) |seen| {
            if (try seen.observe(self.name.items)) {
                self.deferred_short_tail = null;
                return;
            }
            if (seen.bindsToBlob()) {
                name_offset = seen.last_offset;
                out_name_len = seen.last_len;
                name_pad[0] = index_format.NAME_IN_ZFI_FLAG;
                name_embedded = true;
            }
        }

        var side_table_slice: []const u8 = &.{};
        if (seq_info.uses_uniform_formula) {
            // Short final wrap was deferred and never needed; drop it.
            self.deferred_short_tail = null;
        } else if (self.collect_side_tables) {
            // finish() can mark non-uniform after every line looked formula-ready (e.g. last
            // line longer than the stride). Materialize the full uniform-looking prefix then.
            // A deferred short tail becomes the final side-table row if finish rejects it.
            if (self.deferred_short_tail) |short| {
                if (!self.side_table_active) {
                    if (metrics.line_count < 1) return error.MissingSideTable;
                    try self.activateSideTable(metrics.line_count - 1);
                    try self.appendSideTableLine(short.bases, short.actual_bytes, short.line_file_offset);
                }
                self.deferred_short_tail = null;
            } else if (!self.side_table_active) {
                try self.activateSideTable(metrics.line_count);
            }
            if (!self.side_table_active) {
                return error.MissingSideTable;
            }
            @memcpy(self.side_table.items[0..@sizeOf(u64)], std.mem.asBytes(&metrics.line_count));
            side_table_slice = self.side_table.items;
        }
        try emitRecord(ctx, .{
            .record = .{
                .name_offset = name_offset,
                .name_len = out_name_len,
                .seq_offset = self.seq_offset,
                .seq_len = seq_info.seq_len,
                .line_bases = metrics.line_bases,
                .line_bytes = metrics.line_bytes,
                ._pad = name_pad,
            },
            .name = self.name.items,
            .uses_uniform_formula = seq_info.uses_uniform_formula,
            .side_table = side_table_slice,
            .name_embedded = name_embedded,
        });
        record_count.* += 1;
    }
};

fn processChunkBytes(
    state: *ChunkParseState,
    data: []const u8,
    file_offset: u64,
    max_name_len: ?usize,
    seen_names: ?*NameDedup,
    ctx: anytype,
    comptime emitRecord: fn (@TypeOf(ctx), FastaRecordEmit) anyerror!void,
    record_count: *u32,
) !void {
    var i: usize = 0;

    if (state.in_pending_line) {
        const nl = findNextNewline(data, 0, data.len);
        if (nl >= data.len) {
            try state.appendPendingFragment(data, file_offset);
            return;
        }

        try state.appendPendingFragment(data[0..nl], file_offset);
        state.pending_line_bytes = std.math.add(u32, state.pending_line_bytes, 1) catch
            return error.SequenceLineTooLong;
        try state.commitPendingLine(true);
        i = nl + 1;
    }

    while (i < data.len) {
        const byte_offset = file_offset + @as(u64, @intCast(i));
        const byte = data[i];

        if (state.at_line_start and byte == '>') {
            if (state.active) {
                try state.emitRecordIfReady(ctx, seen_names, emitRecord, record_count);
            }
            try state.startHeader(byte_offset);
            i += 1;
            continue;
        }

        if (state.in_header) {
            if (state.parsing_name) {
                var name_end = i;
                while (name_end < data.len and
                    data[name_end] != ' ' and
                    data[name_end] != '\t' and
                    data[name_end] != '\r' and
                    data[name_end] != '\n')
                {
                    name_end += 1;
                }

                const fragment = data[i..name_end];
                if (state.name.items.len > MAX_INDEX_NAME_LEN or
                    fragment.len > MAX_INDEX_NAME_LEN - state.name.items.len)
                {
                    return error.HeaderTooLong;
                }
                if (max_name_len) |limit| {
                    if (state.name.items.len > limit or fragment.len > limit - state.name.items.len) {
                        return error.HeaderTooLong;
                    }
                }
                try state.name.appendSlice(state.allocator, fragment);
                i = name_end;
                if (i >= data.len) return;
                if (data[i] == '\n') {
                    state.finalizeHeader();
                    i += 1;
                    continue;
                }
                state.parsing_name = false;
            }

            const header_end = findNextNewline(data, i, data.len);
            if (header_end >= data.len) return;
            state.finalizeHeader();
            i = header_end + 1;
            continue;
        }

        if (!state.active) {
            state.at_line_start = byte == '\n';
            i += 1;
            continue;
        }

        // After the first dense full-width line, accept blocks only when every separator is
        // exact and every other byte is greater than space. The next header remains in-band.
        if (state.at_line_start and state.strideEligible()) {
            const line_bases = state.line_builder.metrics.line_bases;
            const line_bytes = state.line_builder.metrics.line_bytes;
            const run_lines = validatedStrideRun(
                STRIDE_VECTOR_LEN,
                data,
                i,
                line_bases,
                line_bytes,
            );
            if (run_lines > 0) {
                state.ingestStrideRun(run_lines);
                i += run_lines * @as(usize, line_bytes);
            }
            if (i >= data.len) return;
            // Next header: re-enter so the `>` branch runs. Blank/invalid: fall through
            // to the line machine once (must not re-enter stride on the same offset).
            if (data[i] == '>') continue;
            if (data.len - i < line_bytes) {
                // A partial chunk has no complete line when LF is absent.
                if (findNextNewline(data, i, data.len) >= data.len) {
                    try state.appendPendingFragment(
                        data[i..],
                        file_offset + @as(u64, @intCast(i)),
                    );
                    return;
                }
            }
        }

        const line_start = i;
        const nl = findNextNewline(data, i, data.len);
        if (nl >= data.len) {
            try state.appendPendingFragment(
                data[line_start..],
                file_offset + @as(u64, @intCast(line_start)),
            );
            return;
        }

        try state.ingestSequenceLine(
            data,
            line_start,
            nl,
            true,
            file_offset + @as(u64, @intCast(line_start)),
        );
        i = nl + 1;
    }
}

pub fn scanFastaReader(
    reader: *std.Io.Reader,
    read_buf: []u8,
    options: ScanOptions,
    allocator: std.mem.Allocator,
    external_dedup: ?*NameDedup,
    ctx: anytype,
    comptime emitRecord: fn (@TypeOf(ctx), FastaRecordEmit) anyerror!void,
) !u32 {
    // FAI owns stable arena keys. `.zfi` receives blob-backed dedup so names live once.
    var owned_dedup: ?NameDedup = null;
    const dedup_ptr: ?*NameDedup = blk: {
        if (external_dedup) |p| break :blk p;
        if (options.enable_dedup) {
            owned_dedup = NameDedup.init(allocator);
            try owned_dedup.?.ensureCapacity(DEDUP_INITIAL_CAPACITY);
            break :blk &owned_dedup.?;
        }
        break :blk null;
    };
    defer if (owned_dedup) |*seen| seen.deinit();

    var state = ChunkParseState{
        .allocator = allocator,
        .collect_side_tables = options.collect_side_tables,
    };
    defer state.name.deinit(allocator);
    defer if (options.collect_side_tables) state.side_table.deinit(allocator);
    var record_count: u32 = 0;
    var file_offset: u64 = 0;

    while (true) {
        const n = reader.readSliceShort(read_buf) catch return error.SourceReadFailed;
        if (n == 0) break;
        if (file_offset == 0 and options.require_initial_header and read_buf[0] != '>') {
            return error.NotFasta;
        }
        try processChunkBytes(
            &state,
            read_buf[0..n],
            file_offset,
            options.max_name_len,
            dedup_ptr,
            ctx,
            emitRecord,
            &record_count,
        );
        file_offset += @intCast(n);
    }

    if (state.active) {
        try state.emitRecordIfReady(ctx, dedup_ptr, emitRecord, &record_count);
    }

    return record_count;
}

const FaiEmitBuffer = struct {
    writer: *std.Io.Writer,
    buf: [INDEX_OUTPUT_BUFFER_SIZE]u8 = undefined,
    len: usize = 0,

    fn flush(self: *FaiEmitBuffer) !void {
        if (self.len == 0) return;
        self.writer.writeAll(self.buf[0..self.len]) catch return error.FaiSpoolWriteFailed;
        self.len = 0;
    }

    fn emitRecord(self: *FaiEmitBuffer, emit_info: FastaRecordEmit) !void {
        if (!emit_info.uses_uniform_formula) return error.NonUniformFai;
        const rec = emit_info.record;
        var suffix_buf: [FAI_SUFFIX_BUFFER_SIZE]u8 = undefined;
        const suffix = formatFaiSuffix(&suffix_buf, rec);
        const total = emit_info.name.len + suffix.len;
        if (total > self.buf.len) {
            try self.flush();
            self.writer.writeAll(emit_info.name) catch return error.FaiSpoolWriteFailed;
            self.writer.writeAll(suffix) catch return error.FaiSpoolWriteFailed;
            return;
        }
        if (self.len + total > self.buf.len) try self.flush();
        @memcpy(self.buf[self.len..][0..emit_info.name.len], emit_info.name);
        self.len += emit_info.name.len;
        @memcpy(self.buf[self.len..][0..suffix.len], suffix);
        self.len += suffix.len;
    }
};

fn appendUnsignedDecimal(buffer: []u8, value: u64) usize {
    var reversed: [FAI_U64_DECIMAL_DIGITS]u8 = undefined;
    var remaining = value;
    var digit_count: usize = 0;
    while (true) {
        reversed[digit_count] = @intCast('0' + remaining % 10);
        digit_count += 1;
        remaining /= 10;
        if (remaining == 0) break;
    }

    for (0..digit_count) |index| {
        buffer[index] = reversed[digit_count - index - 1];
    }
    return digit_count;
}

fn formatFaiSuffix(buffer: *[FAI_SUFFIX_BUFFER_SIZE]u8, rec: IndexRecord) []const u8 {
    const values = [_]u64{ rec.seq_len, rec.seq_offset, rec.line_bases, rec.line_bytes };
    var len: usize = 0;
    for (values) |value| {
        buffer[len] = '\t';
        len += 1;
        len += appendUnsignedDecimal(buffer[len..], value);
    }
    buffer[len] = '\n';
    return buffer[0 .. len + 1];
}

const FaiSpool = struct {
    dir: std.Io.Dir,
    close_dir: bool,
    file: std.Io.File,
    name: [FAI_SPOOL_NAME_BUFFER_SIZE]u8,
    name_len: usize,

    fn deinit(self: *FaiSpool, io: std.Io) void {
        self.file.close(io);
        self.dir.deleteFile(io, self.name[0..self.name_len]) catch {};
        if (self.close_dir) self.dir.close(io);
    }
};

fn tryCreateFaiSpool(io: std.Io, dir: std.Io.Dir, close_dir: bool) ?FaiSpool {
    for (0..FAI_SPOOL_CREATE_ATTEMPTS) |_| {
        var nonce: u128 = undefined;
        std.Io.randomSecure(io, std.mem.asBytes(&nonce)) catch return null;

        var name_buf: [FAI_SPOOL_NAME_BUFFER_SIZE]u8 = undefined;
        const name = std.fmt.bufPrint(&name_buf, "z-fasta-{x}.fai.tmp", .{nonce}) catch unreachable;
        const file = dir.createFile(io, name, .{ .read = true, .exclusive = true }) catch |err| switch (err) {
            error.PathAlreadyExists => continue,
            else => return null,
        };

        var spool = FaiSpool{
            .dir = dir,
            .close_dir = close_dir,
            .file = file,
            .name = undefined,
            .name_len = name.len,
        };
        @memcpy(spool.name[0..name.len], name);
        return spool;
    }
    return null;
}

fn tryOpenFaiSpoolDir(io: std.Io, path: []const u8) ?FaiSpool {
    const dir = std.Io.Dir.cwd().openDir(io, path, .{}) catch return null;
    return tryCreateFaiSpool(io, dir, true) orelse {
        dir.close(io);
        return null;
    };
}

fn openFaiSpool(
    io: std.Io,
    environ: std.process.Environ,
    allocator: std.mem.Allocator,
) !FaiSpool {
    const env_names = [_][]const u8{ "TMPDIR", "TEMP", "TMP" };
    if (comptime @hasDecl(std.process.Environ.Block, "view")) {
        for (env_names) |env_name| {
            const path = std.process.Environ.getPosix(environ, env_name) orelse continue;
            if (path.len == 0) continue;
            if (tryOpenFaiSpoolDir(io, path)) |spool| return spool;
        }
    } else {
        var env_map = environ.createMap(allocator) catch return error.OutOfMemory;
        defer env_map.deinit();
        for (env_names) |env_name| {
            const path = env_map.get(env_name) orelse continue;
            if (path.len == 0) continue;
            if (tryOpenFaiSpoolDir(io, path)) |spool| return spool;
        }
    }

    if (builtin.os.tag != .windows) {
        if (tryOpenFaiSpoolDir(io, "/tmp")) |spool| return spool;
    }
    return tryCreateFaiSpool(io, .cwd(), false) orelse error.NoUsableFaiSpool;
}

fn sourceMetadataMatches(before: std.Io.File.Stat, after: std.Io.File.Stat) bool {
    return before.size == after.size and before.mtime.nanoseconds == after.mtime.nanoseconds;
}

fn sourcePathUnchanged(io: std.Io, dir: std.Io.Dir, path: []const u8, before: std.Io.File.Stat) !bool {
    const after = dir.statFile(io, path, .{}) catch |err| switch (err) {
        error.FileNotFound => return false,
        else => return error.SourceStatFailed,
    };
    return sourceMetadataMatches(before, after);
}

fn scanFaiReader(
    reader: *std.Io.Reader,
    read_buf: []u8,
    writer: anytype,
    options: ScanOptions,
    allocator: std.mem.Allocator,
) !u32 {
    var fai_buf = FaiEmitBuffer{ .writer = writer };
    const Ctx = struct {
        buffer: *FaiEmitBuffer,

        fn emit(ctx: *@This(), emit_info: FastaRecordEmit) !void {
            try ctx.buffer.emitRecord(emit_info);
        }
    };

    var ctx = Ctx{ .buffer = &fai_buf };
    const count = try scanFastaReader(reader, read_buf, options, allocator, null, &ctx, Ctx.emit);
    try fai_buf.flush();
    return count;
}

pub fn scanFaiData(
    data: []const u8,
    writer: anytype,
    enable_dedup: bool,
    allocator: std.mem.Allocator,
) !u32 {
    var r = std.Io.Reader.fixed(data);
    var read_buf: [INDEX_READ_BUFFER_SIZE]u8 = undefined;
    return scanFaiReader(&r, &read_buf, writer, .{ .enable_dedup = enable_dedup }, allocator);
}

pub fn runIndex(
    io: std.Io,
    environ: std.process.Environ,
    path: []const u8,
    options: IndexOptions,
) !void {
    const file = std.Io.Dir.cwd().openFile(io, path, .{}) catch |err| switch (err) {
        error.FileNotFound, error.AccessDenied => return err,
        else => return error.SourceOpenFailed,
    };
    defer file.close(io);

    const stat = file.stat(io) catch return error.SourceStatFailed;
    if (stat.size == 0) return error.EmptyFile;

    // Page allocation releases grown catalogs instead of retaining every old buffer.
    // FAI dedup still uses a child arena because its slice keys must remain stable.
    const allocator = std.heap.page_allocator;

    var io_buf: [FILE_IO_BUF_SIZE]u8 = undefined;
    var read_buf: [INDEX_READ_BUFFER_SIZE]u8 = undefined;
    var file_reader = file.reader(io, &io_buf);

    if (options.emit_fai) {
        var spool = try openFaiSpool(io, environ, allocator);
        defer spool.deinit(io);

        var out_buf: [INDEX_OUTPUT_BUFFER_SIZE]u8 = undefined;
        var stdout_fw = std.Io.File.Writer.initStreaming(.stdout(), io, &out_buf);

        var tmp_io_buf: [INDEX_OUTPUT_BUFFER_SIZE]u8 = undefined;
        var tmp_fw = spool.file.writer(io, &tmp_io_buf);

        const record_count = scanFaiReader(
            &file_reader.interface,
            &read_buf,
            &tmp_fw.interface,
            .{ .enable_dedup = options.enable_dedup, .require_initial_header = true },
            allocator,
        ) catch |err| switch (err) {
            error.HeaderTooLong,
            error.SequenceLineTooLong,
            error.NotFasta,
            error.NonUniformFai,
            error.FaiSpoolWriteFailed,
            error.SourceReadFailed,
            error.OutOfMemory,
            => return err,
            else => return error.ProcessingFailed,
        };
        tmp_fw.interface.flush() catch return error.FaiSpoolWriteFailed;
        if (record_count == 0) return error.NoValidSequences;

        if (!try sourcePathUnchanged(io, std.Io.Dir.cwd(), path, stat)) return error.SourceChanged;

        const tmp_size = (spool.file.stat(io) catch return error.FaiSpoolReadFailed).size;
        var offset: u64 = 0;
        var copy_buf: [INDEX_OUTPUT_BUFFER_SIZE]u8 = undefined;
        while (offset < tmp_size) {
            const want: usize = @intCast(@min(copy_buf.len, tmp_size - offset));
            const n = std.Io.File.readPositionalAll(
                spool.file,
                io,
                copy_buf[0..want],
                offset,
            ) catch return error.FaiSpoolReadFailed;
            if (n == 0) return error.FaiSpoolReadFailed;
            stdout_fw.interface.writeAll(copy_buf[0..n]) catch return error.StdoutReplayFailed;
            offset += n;
        }
        stdout_fw.flush() catch return error.StdoutFlushFailed;
        return;
    }

    var zfi_path_buf: [INDEX_PATH_BUFFER_SIZE]u8 = undefined;
    const zfi_path = std.fmt.bufPrint(&zfi_path_buf, "{s}.zfi", .{path}) catch
        return error.OutputPathTooLong;

    const cwd = std.Io.Dir.cwd();

    var zfi_index = scanZfiReaderWithOptions(
        &file_reader.interface,
        &read_buf,
        .{ .enable_dedup = options.enable_dedup, .require_initial_header = true },
        allocator,
    ) catch |err| switch (err) {
        error.HeaderTooLong,
        error.SequenceLineTooLong,
        error.NotFasta,
        error.SourceReadFailed,
        error.OutOfMemory,
        => return err,
        error.MissingSideTable => return error.ProcessingFailed,
        else => return error.ProcessingFailed,
    };
    defer zfi_index.deinit(allocator);

    if (zfi_index.records.items.len == 0) return error.NoValidSequences;

    const source_mtime_ns = index_format.timestampToNs(stat.mtime) catch return error.UnsupportedTimestamp;
    var atomic_file = cwd.createFileAtomic(io, zfi_path, .{ .replace = true }) catch
        return error.ZfiWriteFailed;
    defer atomic_file.deinit(io);

    var out_buf: [INDEX_OUTPUT_BUFFER_SIZE]u8 = undefined;
    var file_fw = atomic_file.file.writer(io, &out_buf);
    writeZfiIndex(&file_fw.interface, &zfi_index, stat.size, source_mtime_ns) catch
        return error.ZfiWriteFailed;
    file_fw.flush() catch return error.ZfiWriteFailed;

    if (!try sourcePathUnchanged(io, cwd, path, stat)) return error.SourceChanged;

    atomic_file.replace(io) catch return error.ZfiFinalizeFailed;

    std.debug.print("wrote {s} ({d} sequences)\n", .{ zfi_path, zfi_index.records.items.len });
}

test "[property] - [stride scanner]: matches scalar structural boundaries" {
    const line_bases: u32 = 4;
    const line_count = 520;
    const target_line = 256;
    const cases = [_][]const u8{
        "AAAA\n",
        "AAAA\r\n",
    };
    const mutations = [_]u8{ ' ', '\t', '\r', '\n', 0, '>' };

    var storage: [line_count * 6]u8 = undefined;
    for (cases) |line| {
        const line_bytes: u32 = @intCast(line.len);
        const block = storage[0 .. line_count * line.len];
        for (0..line_count) |line_index| {
            @memcpy(block[line_index * line.len ..][0..line.len], line);
        }

        for (0..line.len + 1) |line_offset| {
            const mutation_offset = target_line * line.len + line_offset;
            const original = block[mutation_offset];
            for (mutations) |mutation| {
                block[mutation_offset] = mutation;

                const native = validatedStrideRun(
                    STRIDE_VECTOR_LEN,
                    block,
                    0,
                    line_bases,
                    line_bytes,
                );
                const scalar = validatedStrideRun(1, block, 0, line_bases, line_bytes);

                try std.testing.expectEqual(scalar, native);
                const at_next_line = line_offset == line.len;
                const invalid = if (line_offset < line_bases)
                    mutation <= ' ' or (line_offset == 0 and mutation == '>')
                else if (!at_next_line)
                    mutation != line[line_offset]
                else
                    mutation <= ' ' or mutation == '>';
                const stop_line = target_line + @as(usize, @intFromBool(at_next_line));
                try std.testing.expectEqual(if (invalid) stop_line else line_count, native);
            }
            block[mutation_offset] = original;
        }
    }
}

test "[property] - [FAI formatter]: matches standard formatting at integer boundaries" {
    const cases = [_]IndexRecord{
        .{
            .name_offset = 0,
            .name_len = 0,
            .seq_len = 0,
            .seq_offset = 0,
            .line_bases = 0,
            .line_bytes = 0,
        },
        .{
            .name_offset = 0,
            .name_len = 0,
            .seq_len = std.math.maxInt(u64),
            .seq_offset = std.math.maxInt(u64),
            .line_bases = std.math.maxInt(u32),
            .line_bytes = std.math.maxInt(u32),
        },
        .{
            .name_offset = 0,
            .name_len = 0,
            .seq_len = 0,
            .seq_offset = std.math.maxInt(u64),
            .line_bases = 0,
            .line_bytes = std.math.maxInt(u32),
        },
        .{
            .name_offset = 0,
            .name_len = 0,
            .seq_len = std.math.maxInt(u64),
            .seq_offset = 0,
            .line_bases = std.math.maxInt(u32),
            .line_bytes = 0,
        },
    };

    for (cases) |rec| {
        var actual_buf: [FAI_SUFFIX_BUFFER_SIZE]u8 = undefined;
        const actual = formatFaiSuffix(&actual_buf, rec);
        var expected_buf: [128]u8 = undefined;
        const expected = try std.fmt.bufPrint(
            &expected_buf,
            "\t{d}\t{d}\t{d}\t{d}\n",
            .{ rec.seq_len, rec.seq_offset, rec.line_bases, rec.line_bytes },
        );

        try std.testing.expectEqualStrings(expected, actual);
    }
}

test "[unit] - [FASTA scanner]: gates initial-header validation by option" {
    const data = "not-fasta\n";
    var read_buf: [4]u8 = undefined;
    var output = std.Io.Writer.Allocating.init(std.testing.allocator);
    defer output.deinit();

    var unchecked_reader = std.Io.Reader.fixed(data);
    try std.testing.expectEqual(
        @as(u32, 0),
        try scanFaiReader(
            &unchecked_reader,
            &read_buf,
            &output.writer,
            .{ .enable_dedup = true },
            std.testing.allocator,
        ),
    );

    var checked_reader = std.Io.Reader.fixed(data);
    try std.testing.expectError(
        error.NotFasta,
        scanFaiReader(
            &checked_reader,
            &read_buf,
            &output.writer,
            .{ .enable_dedup = true, .require_initial_header = true },
            std.testing.allocator,
        ),
    );
}

test "[failure] - [FASTA scanner]: rejects sequence lines above u32 geometry" {
    var state = ChunkParseState{
        .allocator = std.testing.allocator,
        .in_pending_line = true,
        .pending_line_bases = std.math.maxInt(u32),
    };
    try std.testing.expectError(error.SequenceLineTooLong, state.appendPendingFragment("A", 0));

    state.pending_line_bases = 0;
    state.pending_line_bytes = std.math.maxInt(u32);
    const Ctx = struct {
        fn emit(_: *@This(), _: FastaRecordEmit) !void {
            unreachable;
        }
    };
    var ctx = Ctx{};
    var record_count: u32 = 0;
    try std.testing.expectError(
        error.SequenceLineTooLong,
        processChunkBytes(&state, "\n", 0, null, null, &ctx, Ctx.emit, &record_count),
    );

    state.pending_line_bases = 1;
    try std.testing.expectError(error.SequenceLineTooLong, state.commitPendingLine(false));
}

test "[property] - [name deduplication]: remains exact under hash collisions" {
    const CollisionContext = struct {
        pub fn hash(_: @This(), _: []const u8) u64 {
            return 0;
        }

        pub fn eql(_: @This(), a: []const u8, b: []const u8) bool {
            return std.mem.eql(u8, a, b);
        }
    };
    const CollidingDedup = NameDedupWith(CollisionContext);
    const names = [_][]const u8{ "alpha", "beta", "alpha", "gamma", "beta" };
    var seen = CollidingDedup.init(std.testing.allocator);
    defer seen.deinit();

    var kept: usize = 0;
    for (names) |name| {
        if (!(try seen.observe(name))) kept += 1;
    }

    try std.testing.expectEqual(@as(usize, 3), kept);
    try std.testing.expectEqual(@as(usize, 3), seen.map.count());
}

test "[failure] - [FASTA scanner]: maps reader failures to SourceReadFailed" {
    var reader: std.Io.Reader = .failing;
    var read_buf: [4]u8 = undefined;
    var output = std.Io.Writer.Allocating.init(std.testing.allocator);
    defer output.deinit();

    try std.testing.expectError(
        error.SourceReadFailed,
        scanFaiReader(
            &reader,
            &read_buf,
            &output.writer,
            .{ .enable_dedup = true },
            std.testing.allocator,
        ),
    );
}

test "[failure] - [FAI scanner]: maps writer failures to FaiSpoolWriteFailed" {
    var reader = std.Io.Reader.fixed(">seq\nA\n");
    var read_buf: [4]u8 = undefined;
    var writer: std.Io.Writer = .failing;

    try std.testing.expectError(
        error.FaiSpoolWriteFailed,
        scanFaiReader(
            &reader,
            &read_buf,
            &writer,
            .{ .enable_dedup = true },
            std.testing.allocator,
        ),
    );
}

test "[integration] - [FAI spool]: uses exclusive files and cleans them" {
    var tmp = std.testing.tmpDir(.{ .iterate = true });
    defer tmp.cleanup();

    {
        var first = tryCreateFaiSpool(std.testing.io, tmp.dir, false) orelse
            return error.TestUnexpectedResult;
        defer first.deinit(std.testing.io);
        var second = tryCreateFaiSpool(std.testing.io, tmp.dir, false) orelse
            return error.TestUnexpectedResult;
        defer second.deinit(std.testing.io);

        try std.testing.expect(!std.mem.eql(
            u8,
            first.name[0..first.name_len],
            second.name[0..second.name_len],
        ));
    }

    var entries = tmp.dir.iterate();
    try std.testing.expectEqual(null, try entries.next(std.testing.io));
}

test "[integration] - [source identity]: detects path replacement" {
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();

    const before = blk: {
        const original = try tmp.dir.createFile(std.testing.io, "source.fa", .{});
        defer original.close(std.testing.io);
        try std.Io.File.writeStreamingAll(original, std.testing.io, ">seq\nAAAA\n");
        break :blk try original.stat(std.testing.io);
    };

    try tmp.dir.rename("source.fa", tmp.dir, "old.fa", std.testing.io);
    {
        const replacement = try tmp.dir.createFile(std.testing.io, "source.fa", .{});
        defer replacement.close(std.testing.io);
        try std.Io.File.writeStreamingAll(replacement, std.testing.io, ">seq\nA\n");
    }

    try std.testing.expect(!try sourcePathUnchanged(std.testing.io, tmp.dir, "source.fa", before));
}
