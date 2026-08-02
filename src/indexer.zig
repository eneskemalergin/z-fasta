//! FASTA indexing: mmap slice scan and chunked `--low-mem` reader.
//!
//! `scanFastaRecords` walks a byte slice and invokes a per-record callback. Offsets in
//! `IndexRecord` and side-table lines are absolute positions in the FASTA file.
//!
//! Mmap and streaming share line-metrics rules (`countBasesInLineSlice`,
//! `LineMetricsBuilder`): skip empty lines for geometry, set `seq_offset` to the first
//! base-bearing line, invent a phantom separator when the final line has no newline,
//! and reject names longer than `max_index_name_len`.
//!
//! Mmap length uses the fast `measureSequenceLength` / `countFixedWidthBases` stride on
//! the region from the first base-bearing line. Streaming and side-table builds still
//! walk lines via `LineMetricsBuilder` (interior blank-after-bases → non-uniform;
//! trailing blanks match mmap body trim and stay uniform).
//!
//! Duplicate-name filtering is first-wins and collision-safe on both paths: lookup may
//! use a hash, but identity is always full-string equality (`NameDedup`).

const std = @import("std");
const index_format = @import("index_format.zig");

pub const IndexRecord = index_format.IndexRecord;
pub const ZfiHeader = index_format.ZfiHeader;
pub const ZFI_MAGIC = index_format.ZFI_MAGIC;

const SIMD_CHUNK_SIZE = 32;
const SimdVec = @Vector(SIMD_CHUNK_SIZE, u8);
const STREAM_VECTOR_LEN = std.simd.suggestVectorLength(u8) orelse 1;
// Measured block size from the retained experiment; validation uses no extra storage.
const STRIDE_VALIDATION_BLOCK_LINES = 256;
pub const low_mem_chunk_size = 1 * 1024 * 1024;
const FILE_IO_BUF_SIZE = 8 * 1024;
pub const max_index_name_len = std.math.maxInt(u16);

/// Offset into a name blob for stream `.zfi` catalog dedup (one store for map + embed).
pub const NameBlobRef = struct { off: u32, len: u32 };

/// First-wins name set for indexing. A hash accelerates the map; discarding a record
/// requires full-string equality (`Context.eql`), never hash equality alone.
/// Production uses `NameDedup` (`StringContext`). Tests may inject a colliding hasher via
/// `NameDedupWith` to prove identity stays byte equality under forced collisions.
///
/// Modes:
/// - Borrow (`owns_keys=false`): mmap; keys borrow the source buffer.
/// - Arena (`owns_keys=true`, no blob): stream FAI; keys in a child arena.
/// - Blob (`initOwningBlob`): stream `.zfi`; keys are offsets into `name_blob` so embed
///   does not copy each name a second time (Tx RSS was ~2× catalog without this).
pub fn NameDedupWith(comptime Context: type) type {
    return struct {
        map: std.HashMap([]const u8, void, Context, std.hash_map.default_max_load_percentage),
        /// When true, keys are allocated from `key_arena` and freed only when the arena dies.
        /// Mmap scans set this false and borrow name slices from the source buffer.
        owns_keys: bool,
        key_arena: ?std.heap.ArenaAllocator = null,
        /// When set, stream `.zfi` dedup stores `NameBlobRef` keys into this blob (not `map`).
        name_blob: ?*std.ArrayList(u8) = null,
        ref_map: ?BlobRefMap = null,
        /// After a successful first-seen `observe` in blob mode, the embedded span.
        last_offset: u32 = 0,
        last_len: u16 = 0,

        const BlobHashContext = struct {
            blob: *std.ArrayList(u8),
            pub fn hash(self: @This(), key: NameBlobRef) u64 {
                return std.hash_map.hashString(self.blob.items[key.off..][0..key.len]);
            }
            pub fn eql(self: @This(), a: NameBlobRef, b: NameBlobRef) bool {
                if (a.len != b.len) return false;
                return std.mem.eql(
                    u8,
                    self.blob.items[a.off..][0..a.len],
                    self.blob.items[b.off..][0..b.len],
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
                return std.mem.eql(u8, name, self.blob.items[ref.off..][0..ref.len]);
            }
        };

        const Self = @This();

        /// Soft start capacity for stream catalogs; map grows geometrically after this.
        /// Keep modest so Genome `--low-mem` RSS stays in the few-MB class.
        pub const stream_dedup_pre_cap: u32 = 16384;

        pub fn init(allocator: std.mem.Allocator, owns_keys: bool) Self {
            var self: Self = .{
                .map = std.HashMap([]const u8, void, Context, std.hash_map.default_max_load_percentage).init(allocator),
                .owns_keys = owns_keys,
                .key_arena = null,
            };
            if (owns_keys) {
                self.key_arena = std.heap.ArenaAllocator.init(allocator);
            }
            return self;
        }

        /// Stream `.zfi`: one catalog in `blob` for both dedup and embedded names.
        pub fn initOwningBlob(allocator: std.mem.Allocator, blob: *std.ArrayList(u8)) Self {
            return .{
                .map = std.HashMap([]const u8, void, Context, std.hash_map.default_max_load_percentage).init(allocator),
                .owns_keys = true,
                .name_blob = blob,
                .ref_map = BlobRefMap.initContext(allocator, .{ .blob = blob }),
            };
        }

        pub fn bindsToBlob(self: *const Self) bool {
            return self.name_blob != null;
        }

        pub fn ensureStreamCapacity(self: *Self, min_cap: u32) !void {
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

        /// Returns true when `name` was already observed (caller should skip the record).
        pub fn observe(self: *Self, name: []const u8) !bool {
            if (self.name_blob) |blob| {
                var rm = &(self.ref_map orelse return error.OutOfMemory);
                const lookup = SliceLookup{ .blob = blob };
                const gop = try rm.getOrPutAdapted(name, lookup);
                if (gop.found_existing) {
                    self.last_offset = gop.key_ptr.off;
                    self.last_len = @intCast(gop.key_ptr.len);
                    return true;
                }
                const off: u32 = @intCast(blob.items.len);
                blob.appendSlice(rm.allocator, name) catch |err| {
                    _ = rm.removeByPtr(gop.key_ptr);
                    return err;
                };
                gop.key_ptr.* = .{ .off = off, .len = @intCast(name.len) };
                self.last_offset = off;
                self.last_len = @intCast(name.len);
                return false;
            }
            if (!self.owns_keys) {
                const gop = try self.map.getOrPut(name);
                return gop.found_existing;
            }
            const arena = &(self.key_arena orelse return error.OutOfMemory);
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

pub const NameDedup = NameDedupWith(std.hash_map.StringContext);

pub const OutputMode = enum { fai, zfi };

/// One parsed FASTA record passed to `scanFastaRecords` emit callbacks.
pub const FastaRecordEmit = struct {
    record: IndexRecord,
    name: []const u8,
    /// Sequence region in the scan buffer (byte after header newline through next header or EOF).
    /// May include leading blank lines before `record.seq_offset`.
    seq_data: []const u8,
    /// Buffer index where `seq_data` begins (byte after header newline).
    seq_region_start: usize = 0,
    uses_uniform_formula: bool,
    /// Incremental side-table bytes for streaming `.zfi` (line count + `SideTableLine` rows).
    streaming_side_table: []const u8 = &.{},
    /// When true, `record` already has blob-relative `name_offset`/`name_len` and `name_in_zfi_flag`
    /// (stream `.zfi` blob-backed `NameDedup`); skip `embedZfiName`.
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
    /// First base line: every content byte was a base (no trailing space/tab) and the
    /// on-wire separator was LF or CRLF. Needed so `AAAA \\n` is not treated as CRLF-dense.
    first_content_dense: bool = false,
};

const SequenceLengthInfo = struct {
    seq_len: u64,
    uses_uniform_formula: bool,
};

fn findNextGt(data: []const u8, start: usize) usize {
    var pos = start;
    while (pos + SIMD_CHUNK_SIZE <= data.len) {
        const chunk: SimdVec = data[pos..][0..SIMD_CHUNK_SIZE].*;
        const mask = chunk == @as(SimdVec, @splat('>'));
        if (@reduce(.Or, mask)) {
            inline for (0..SIMD_CHUNK_SIZE) |j| {
                if (mask[j]) return pos + j;
            }
        }
        pos += SIMD_CHUNK_SIZE;
    }
    while (pos < data.len) {
        if (data[pos] == '>') return pos;
        pos += 1;
    }
    return data.len;
}

fn findNextHeaderStart(data: []const u8, start: usize) usize {
    var pos = start;
    while (pos < data.len) {
        const gt_pos = findNextGt(data, pos);
        if (gt_pos >= data.len) return data.len;
        if (gt_pos == 0 or data[gt_pos - 1] == '\n') return gt_pos;
        pos = gt_pos + 1;
    }
    return data.len;
}

pub fn countBases(data: []const u8) u64 {
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

/// Bases in one sequence line's content slice (no trailing `\n`/`\r`).
/// Dense lines (no bytes `<= ' '`) use length only; otherwise falls back to `countBases`.
fn countBasesInLineSlice(data: []const u8, start: usize, end: usize) u32 {
    if (start >= end) return 0;
    if (!lineSliceHasWhitespace(data, start, end)) {
        return @intCast(end - start);
    }
    return @intCast(countBases(data[start..end]));
}

fn countFixedWidthBases(data: []const u8, line_bases: u32, line_bytes: u32) ?u64 {
    if (data.len == 0) return 0;
    if (line_bases == 0 or line_bytes <= line_bases) return null;

    const sep_len = line_bytes - line_bases;
    if (sep_len != 1 and sep_len != 2) return null;

    var cursor: usize = 0;
    var bases: u64 = 0;
    while (cursor < data.len) {
        const remaining = data.len - cursor;
        if (remaining >= line_bytes) {
            if (!hasLineSeparatorAt(data, cursor + line_bases, sep_len)) return null;
            if (data[cursor + line_bases - 1] <= ' ') return null;
            bases += line_bases;
            cursor += line_bytes;
        } else {
            if (remaining <= line_bases) {
                const tail = data[cursor..];
                if (tail.len > 1) {
                    const interior = tail[0 .. tail.len - 1];
                    if (std.mem.indexOfScalar(u8, interior, '\n') != null) return null;
                    if (std.mem.indexOfScalar(u8, interior, '\r') != null) return null;
                }
                return bases + countBasesInLineSlice(data, cursor, data.len);
            }

            const tail_bases = remaining - sep_len;
            if (tail_bases > line_bases) return null;
            if (!hasLineSeparatorAt(data, cursor + tail_bases, sep_len)) return null;
            if (tail_bases > 0 and data[cursor + tail_bases - 1] <= ' ') return null;
            return bases + tail_bases;
        }
    }
    return bases;
}

fn trimTrailingRecordNewlines(data: []const u8) []const u8 {
    var end = data.len;
    while (end > 0) {
        const c = data[end - 1];
        if (c != '\n' and c != '\r') break;
        end -= 1;
    }
    return data[0..end];
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

fn measureSequenceLength(data: []const u8, line_bases: u32, line_bytes: u32) SequenceLengthInfo {
    const body = trimTrailingRecordNewlines(data);
    const fixed_width = countFixedWidthBases(body, line_bases, line_bytes) orelse return .{
        .seq_len = countBases(body),
        .uses_uniform_formula = false,
    };
    if (!firstLineIsDense(line_bases, line_bytes)) {
        const fallback = countBases(body);
        return .{
            .seq_len = if (fixed_width == fallback) fixed_width else fallback,
            .uses_uniform_formula = fixed_width == fallback,
        };
    }
    return .{
        .seq_len = fixed_width,
        .uses_uniform_formula = true,
    };
}

fn findNextNewline(data: []const u8, start: usize, end: usize) usize {
    var pos = start;
    while (pos + SIMD_CHUNK_SIZE <= end) {
        const chunk: SimdVec = data[pos..][0..SIMD_CHUNK_SIZE].*;
        const mask = chunk == @as(SimdVec, @splat('\n'));
        if (@reduce(.Or, mask)) {
            inline for (0..SIMD_CHUNK_SIZE) |j| {
                if (mask[j]) return pos + j;
            }
        }
        pos += SIMD_CHUNK_SIZE;
    }
    while (pos < end) : (pos += 1) {
        if (data[pos] == '\n') return pos;
    }
    return end;
}

fn findNextStreamingNewline(data: []const u8, start: usize, end: usize) usize {
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
    /// Blank seen after at least one base line; applied when another base line follows
    /// (trailing blanks before EOF/next header match mmap `trimTrailingRecordNewlines`).
    blank_after_bases: bool = false,

    fn ingestLine(self: *LineMetricsBuilder, bases: u32, actual_bytes: u32, has_lf: bool, content_len: u32) void {
        if (bases == 0) {
            // Leading blanks ignored. Trailing blanks ignored until a later base line proves
            // an interior blank (mmap trims trailing newlines from the sequence body).
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
            // Interior full wraps must match on-wire bytes (LF vs CRLF is a real break;
            // mmap `countFixedWidthBases` rejects mixed separators in the body).
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
        // Trailing blanks do not break uniformity (same as mmap body trim).
        self.blank_after_bases = false;
        if (self.metrics.line_count > 1 and self.have_pending) {
            // Final wrap only: mmap trims trailing `\r`/`\n`, so CRLF vs LF on the last
            // line must not flip when content width does not grow.
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

const SequenceRegionScan = struct {
    metrics: LineMetrics,
    /// Buffer index of the first base-bearing line start, if any.
    first_base_line_start: ?usize,
};

// Walk a sequence region with the same rules as streaming `LineMetricsBuilder`:
// blank lines do not set geometry; `first_base_line_start` is the index `seq_offset` uses.
fn scanSequenceRegion(data: []const u8, seq_region_start: usize, seq_end: usize) SequenceRegionScan {
    var builder = LineMetricsBuilder{};
    var first_base_line_start: ?usize = null;
    var pos = seq_region_start;
    while (pos < seq_end) {
        const line_start = pos;
        const line_end = findNextNewline(data, pos, seq_end);
        const has_lf = line_end < seq_end;
        var content_end = line_end;
        if (content_end > line_start and data[content_end - 1] == '\r') {
            content_end -= 1;
        }

        const bases = countBasesInLineSlice(data, line_start, content_end);
        const actual_bytes: u32 = @intCast((if (has_lf) line_end + 1 else line_end) - line_start);
        const content_len: u32 = @intCast(content_end - line_start);
        builder.ingestLine(bases, actual_bytes, has_lf, content_len);
        if (bases > 0 and first_base_line_start == null) {
            first_base_line_start = line_start;
        }

        if (!has_lf) break;
        pos = line_end + 1;
    }
    return .{ .metrics = builder.finish(), .first_base_line_start = first_base_line_start };
}

fn scanLineMetrics(data: []const u8, seq_offset: usize, seq_end: usize) LineMetrics {
    return scanSequenceRegion(data, seq_offset, seq_end).metrics;
}

/// Skip leading blank lines; return geometry of the first base-bearing line.
fn peekFirstBaseLineGeometry(data: []const u8, seq_region_start: usize, seq_end: usize) ?struct {
    first_base_line: usize,
    line_bases: u32,
    line_bytes: u32,
} {
    var pos = seq_region_start;
    while (pos < seq_end) {
        const line_start = pos;
        const line_end = findNextNewline(data, pos, seq_end);
        const has_lf = line_end < seq_end;
        var content_end = line_end;
        if (content_end > line_start and data[content_end - 1] == '\r') {
            content_end -= 1;
        }
        const bases = countBasesInLineSlice(data, line_start, content_end);
        if (bases > 0) {
            const line_bytes: u32 = if (has_lf)
                @intCast((line_end + 1) - line_start)
            else
                @intCast((line_end - line_start) + 1);
            return .{
                .first_base_line = line_start,
                .line_bases = bases,
                .line_bytes = line_bytes,
            };
        }
        if (!has_lf) break;
        pos = line_end + 1;
    }
    return null;
}

fn appendSideTable(
    data: []const u8,
    seq_offset: usize,
    seq_end: usize,
    file_base_offset: u64,
    table: *std.ArrayList(u8),
    allocator: std.mem.Allocator,
) !u64 {
    const local_offset = table.items.len;
    var line_count: u64 = 0;
    try table.appendSlice(allocator, std.mem.asBytes(&line_count));

    var builder = LineMetricsBuilder{};
    var pos = seq_offset;
    var base_start: u64 = 0;
    while (pos < seq_end) {
        const line_start = pos;
        const line_end = findNextNewline(data, pos, seq_end);
        const has_lf = line_end < seq_end;
        var content_end = line_end;
        if (content_end > line_start and data[content_end - 1] == '\r') {
            content_end -= 1;
        }

        const bases = countBasesInLineSlice(data, line_start, content_end);
        const actual_bytes: u32 = @intCast((if (has_lf) line_end + 1 else line_end) - line_start);
        const content_len: u32 = @intCast(content_end - line_start);
        builder.ingestLine(bases, actual_bytes, has_lf, content_len);

        if (bases > 0) {
            const entry = index_format.SideTableLine{
                .base_start = base_start,
                .byte_offset = @intCast(file_base_offset + line_start),
                .line_bytes = actual_bytes,
                .line_bases = bases,
            };
            try table.appendSlice(allocator, std.mem.asBytes(&entry));
            base_start += bases;
        }

        if (!has_lf) break;
        pos = line_end + 1;
    }

    const metrics = builder.finish();
    @memcpy(table.items[local_offset..][0..@sizeOf(u64)], std.mem.asBytes(&metrics.line_count));
    return @intCast(local_offset);
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

pub fn scanFastaRecords(
    data: []const u8,
    file_base_offset: u64,
    enable_dedup: bool,
    dedup_names: ?*NameDedup,
    max_name_len: ?usize,
    allocator: std.mem.Allocator,
    ctx: anytype,
    comptime emitRecord: fn (@TypeOf(ctx), FastaRecordEmit) anyerror!void,
) !u32 {
    // Mmap keys borrow `data`; do not own copies.
    var local_dedup: ?NameDedup = null;
    if (enable_dedup and dedup_names == null) {
        local_dedup = NameDedup.init(allocator, false);
    }
    defer if (local_dedup) |*s| s.deinit();

    const seen_names: ?*NameDedup = if (enable_dedup)
        dedup_names orelse &local_dedup.?
    else
        null;

    var record_count: u32 = 0;
    var pos: usize = 0;

    while (pos < data.len) {
        pos = findNextHeaderStart(data, pos);
        if (pos >= data.len) break;

        const name_local = pos + 1;
        var name_end = name_local;
        while (name_end < data.len and
            data[name_end] != ' ' and
            data[name_end] != '\t' and
            data[name_end] != '\n' and
            data[name_end] != '\r')
        {
            name_end += 1;
        }
        const name_len_usize = name_end - name_local;
        if (name_len_usize > max_index_name_len) return error.HeaderTooLong;
        if (max_name_len) |limit| {
            if (name_len_usize > limit) return error.HeaderTooLong;
        }
        const name_len: u16 = @intCast(name_len_usize);

        var header_end = name_end;
        while (header_end < data.len and data[header_end] != '\n') {
            header_end += 1;
        }

        const seq_region_start: usize = if (header_end < data.len) header_end + 1 else header_end;
        const seq_end = findNextHeaderStart(data, seq_region_start);
        // Fast path: only walk leading blanks, then stride-count with measureSequenceLength.
        // Irregular / blank-after-bases regions fall back via uses_uniform_formula=false → side table.
        const first = peekFirstBaseLineGeometry(data, seq_region_start, seq_end) orelse {
            pos = seq_end;
            continue;
        };
        const seq_info = measureSequenceLength(
            data[first.first_base_line..seq_end],
            first.line_bases,
            first.line_bytes,
        );

        pos = seq_end;
        if (seq_info.seq_len == 0) continue;

        const name = data[name_local..][0..name_len];
        if (seen_names) |seen| {
            if (try seen.observe(name)) continue;
        }

        const seq_data = data[seq_region_start..seq_end];
        const rec = IndexRecord{
            .name_offset = file_base_offset + @as(u64, @intCast(name_local)),
            .name_len = name_len,
            .seq_offset = file_base_offset + @as(u64, @intCast(first.first_base_line)),
            .seq_len = seq_info.seq_len,
            .line_bases = first.line_bases,
            .line_bytes = first.line_bytes,
        };
        try emitRecord(ctx, .{
            .record = rec,
            .name = name,
            .seq_data = seq_data,
            .seq_region_start = seq_region_start,
            .uses_uniform_formula = seq_info.uses_uniform_formula,
        });
        record_count += 1;
    }
    return record_count;
}

pub fn streamingScan(
    data: []const u8,
    writer: anytype,
    mode: OutputMode,
    enable_dedup: bool,
    allocator: std.mem.Allocator,
) !u32 {
    const Ctx = struct {
        writer: @TypeOf(writer),
        mode: OutputMode,

        fn emit(ctx: *@This(), emit_info: FastaRecordEmit) !void {
            const rec = emit_info.record;
            switch (ctx.mode) {
                .fai => {
                    // `.fai` can only represent one fixed line geometry per record.
                    if (!emit_info.uses_uniform_formula) return error.NonUniformFai;
                    try ctx.writer.print("{s}\t{d}\t{d}\t{d}\t{d}\n", .{
                        emit_info.name, rec.seq_len, rec.seq_offset, rec.line_bases, rec.line_bytes,
                    });
                },
                .zfi => {
                    try ctx.writer.writeAll(std.mem.asBytes(&rec));
                },
            }
        }
    };

    var ctx = Ctx{ .writer = writer, .mode = mode };
    return scanFastaRecords(data, 0, enable_dedup, null, null, allocator, &ctx, Ctx.emit);
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
    rec._pad[0] |= index_format.name_in_zfi_flag;
}

pub fn scanZfiIndex(data: []const u8, enable_dedup: bool, allocator: std.mem.Allocator) !ZfiIndex {
    var index = ZfiIndex{
        .records = .empty,
        .side_tables = .empty,
        .name_blob = .empty,
    };
    errdefer index.deinit(allocator);

    const Ctx = struct {
        index: *ZfiIndex,
        data: []const u8,
        file_base_offset: u64,
        allocator: std.mem.Allocator,

        fn emit(ctx: *@This(), emit_info: FastaRecordEmit) !void {
            var rec = emit_info.record;
            if (!emit_info.uses_uniform_formula) {
                const local_offset = try appendSideTable(
                    ctx.data,
                    emit_info.seq_region_start,
                    emit_info.seq_region_start + emit_info.seq_data.len,
                    ctx.file_base_offset,
                    &ctx.index.side_tables,
                    ctx.allocator,
                );
                try rec.markNonUniform(local_offset);
            }
            try embedZfiName(&rec, emit_info.name, &ctx.index.name_blob, ctx.allocator);
            try ctx.index.records.append(ctx.allocator, rec);
        }
    };

    var ctx = Ctx{ .index = &index, .data = data, .file_base_offset = 0, .allocator = allocator };
    _ = try scanFastaRecords(data, 0, enable_dedup, null, null, allocator, &ctx, Ctx.emit);
    try finalizeSideTableOffsets(&index);
    return index;
}

pub fn scanZfiIndexStreaming(
    reader: *std.Io.Reader,
    read_buf: []u8,
    enable_dedup: bool,
    allocator: std.mem.Allocator,
) !ZfiIndex {
    var index = ZfiIndex{
        .records = .empty,
        .side_tables = .empty,
        .name_blob = .empty,
    };
    errdefer index.deinit(allocator);
    // Modest pre-size: enough for Proteome/small catalogs; Tx grows geometrically.
    try index.records.ensureTotalCapacity(allocator, 4096);

    const Ctx = struct {
        index: *ZfiIndex,
        allocator: std.mem.Allocator,

        fn emit(ctx: *@This(), emit_info: FastaRecordEmit) !void {
            var rec = emit_info.record;
            if (!emit_info.uses_uniform_formula) {
                const local_offset = ctx.index.side_tables.items.len;
                try ctx.index.side_tables.appendSlice(ctx.allocator, emit_info.streaming_side_table);
                try rec.markNonUniform(local_offset);
            }
            if (!emit_info.name_embedded) {
                try embedZfiName(&rec, emit_info.name, &ctx.index.name_blob, ctx.allocator);
            }
            try ctx.index.records.append(ctx.allocator, rec);
        }
    };

    var ctx = Ctx{ .index = &index, .allocator = allocator };
    // One catalog: dedup refs into `name_blob` (no arena + embed double copy).
    var blob_dedup: ?NameDedup = null;
    if (enable_dedup) {
        blob_dedup = NameDedup.initOwningBlob(allocator, &index.name_blob);
        try blob_dedup.?.ensureStreamCapacity(NameDedup.stream_dedup_pre_cap);
    }
    defer if (blob_dedup) |*seen| seen.deinit();
    const external_dedup: ?*NameDedup = if (blob_dedup) |*s| s else null;
    _ = try scanFastaReader(
        reader,
        read_buf,
        enable_dedup,
        null,
        allocator,
        true,
        external_dedup,
        &ctx,
        Ctx.emit,
    );
    try finalizeSideTableOffsets(&index);
    return index;
}

pub fn scanZfiIndexStreamingData(
    data: []const u8,
    enable_dedup: bool,
    allocator: std.mem.Allocator,
) !ZfiIndex {
    var r = std.Io.Reader.fixed(data);
    var read_buf: [low_mem_chunk_size]u8 = undefined;
    return scanZfiIndexStreaming(&r, &read_buf, enable_dedup, allocator);
}

/// Wrap prebuilt records as a production `ZfiIndex` with empty side tables and
/// an empty name blob. Prefer `scanZfiIndex` for real FASTA fixtures so embedded
/// names and side tables match the CLI. Use this for deliberately crafted records
/// (including corrupt fixtures) that still must go through `writeZfiIndex`.
pub fn zfiIndexFromRecords(records: []const IndexRecord, allocator: std.mem.Allocator) !ZfiIndex {
    var index = ZfiIndex{
        .records = .empty,
        .side_tables = .empty,
        .name_blob = .empty,
    };
    errdefer index.deinit(allocator);
    try index.records.appendSlice(allocator, records);
    return index;
}

/// Single on-disk `.zfi` serialization path (`plan/zfi-format.md`):
/// header, records, side tables, name blob, `ZFID` source identity, `ZFNM` footer.
pub fn writeZfiIndex(
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

/// Writes a `ZfiIndex` to a `.zfi` file via `writeZfiIndex`.
pub fn writeZfiIndexFile(
    io: std.Io,
    path: []const u8,
    index: *const ZfiIndex,
    source_size: u64,
    source_mtime_ns: u64,
) !void {
    const file = try std.Io.Dir.cwd().createFile(io, path, .{});
    defer file.close(io);

    var file_buf: [65536]u8 = undefined;
    var file_fw = file.writer(io, &file_buf);
    try writeZfiIndex(&file_fw.interface, index, source_size, source_mtime_ns);
    try file_fw.flush();
}

/// Serializes a `ZfiIndex` with the same bytes as `writeZfiIndexFile`.
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

/// Scans FASTA data and returns index records as ArrayList (for testing).
pub fn scanHeaders(data: []const u8, allocator: std.mem.Allocator) !std.ArrayList(IndexRecord) {
    return scanHeadersAt(data, 0, allocator);
}

/// Like `scanHeaders` but `file_base_offset` is added to every stored file offset.
pub fn scanHeadersAt(
    data: []const u8,
    file_base_offset: u64,
    allocator: std.mem.Allocator,
) !std.ArrayList(IndexRecord) {
    var records: std.ArrayList(IndexRecord) = .empty;
    errdefer records.deinit(allocator);

    const Ctx = struct {
        list: *std.ArrayList(IndexRecord),
        record_allocator: std.mem.Allocator,

        fn emit(ctx: *@This(), emit_info: FastaRecordEmit) !void {
            try ctx.list.append(ctx.record_allocator, emit_info.record);
        }
    };

    var ctx = Ctx{ .list = &records, .record_allocator = allocator };
    _ = try scanFastaRecords(data, file_base_offset, true, null, null, allocator, &ctx, Ctx.emit);
    return records;
}

const ChunkParseState = struct {
    allocator: std.mem.Allocator,
    zfi_mode: bool = false,
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
    /// Pre-LF byte count (updated in appendPendingFragment; LF +1 happens after).
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
    /// Short-or-equal stride mismatch that `finish()` still treats as FAI-uniform when it is
    /// the final base line. Do not materialize O(lines) side-table rows until a later base
    /// line proves the short line was interior (Genome `--low-mem` RSS bug).
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

    /// Materialize `lines_to_materialize` uniform stride rows from the deferred first line.
    /// Call with `line_count - 1` before appending the line that broke uniformity, or with
    /// full `line_count` when `finish()` flips non-uniform after all lines looked formula-ready.
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

    /// Flush a deferred short tail plus the current base line into the side table.
    /// `metrics.line_count` includes the current line; the short tail was the previous base line.
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

        if (!self.zfi_mode) return;

        // Blank lines never become side-table rows (matches mmap `appendSideTable`).
        // Interior blanks flip uniformity when the next base line arrives.
        if (bases == 0) return;

        const metrics = self.line_builder.metrics;
        const matches_stride = metrics.line_bases > 0 and
            bases == metrics.line_bases and
            actual_bytes == self.line_builder.first_actual_bytes;
        const dense = metrics.first_content_dense and firstLineIsDense(metrics.line_bases, metrics.line_bytes);
        const can_use_formula = metrics.is_uniform_width and matches_stride and dense;

        if (can_use_formula) {
            // Keep the first line deferred for the whole uniform prefix. Clearing it on
            // line 2+ dropped prior rows when a later line forced a side table.
            if (metrics.line_count == 1) {
                self.deferred_side_line = .{
                    .bases = bases,
                    .actual_bytes = actual_bytes,
                    .line_file_offset = line_file_offset,
                };
            }
            return;
        }

        // A later base line after a deferred short last-line candidate: short was interior.
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
        // Match mmap: `seq_offset` is the first base-bearing line start (side-table invariant).
        self.seq_offset = line_file_offset;
        self.seq_offset_set = true;
    }

    /// True when further body wraps can use fixed `line_bytes` strides (Track A2).
    fn strideEligible(self: *const ChunkParseState) bool {
        if (!self.active or self.in_header or self.in_pending_line) return false;
        if (self.side_table_active or self.deferred_short_tail != null) return false;
        // Interior blank must hit `ingestLine` so uniformity flips like mmap.
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

    /// Bulk-account `n_lines` consecutive full wraps already validated as one stride block.
    /// Caller must have ingested the first dense line via the line machine (`strideEligible`).
    /// Avoids per-wrap `recordSequenceLine` on long Genome bodies (Track A2 wall).
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

    fn appendPendingFragment(self: *ChunkParseState, fragment: []const u8, fragment_file_offset: u64) void {
        if (!self.in_pending_line) {
            self.in_pending_line = true;
            self.pending_line_bases = 0;
            self.pending_line_bytes = 0;
            self.pending_pre_lf_bytes = 0;
            self.pending_pre_lf_ends_with_cr = false;
            self.pending_line_file_offset = fragment_file_offset;
        }
        self.pending_line_bases += countBasesInLineSlice(fragment, 0, fragment.len);
        self.pending_line_bytes += @intCast(fragment.len);
        self.pending_pre_lf_bytes = self.pending_line_bytes;
        if (fragment.len > 0) {
            self.pending_pre_lf_ends_with_cr = fragment[fragment.len - 1] == '\r';
        }
    }

    fn commitPendingLine(self: *ChunkParseState, has_lf: bool) !void {
        if (!self.in_pending_line and self.pending_line_bytes == 0) return;
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
        // Dedup runs at emit time (after seq_len > 0), matching mmap: empty sequences
        // never claim a name. Blob-backed stream `.zfi` also embeds only kept records.
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
                name_pad[0] = index_format.name_in_zfi_flag;
                name_embedded = true;
            }
        }

        var side_table_slice: []const u8 = &.{};
        if (seq_info.uses_uniform_formula) {
            // Short final wrap was deferred and never needed; drop it.
            self.deferred_short_tail = null;
        } else if (self.zfi_mode) {
            // finish() can mark non-uniform after every line looked formula-ready (e.g. last
            // line longer than the stride). Materialize the full uniform-looking prefix then.
            // If a short tail was deferred then finish still went non-uniform (rare denseness
            // / policy edge), flush short as the final side-table row.
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
            .seq_data = &.{},
            .uses_uniform_formula = seq_info.uses_uniform_formula,
            .streaming_side_table = side_table_slice,
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
        const nl = findNextStreamingNewline(data, 0, data.len);
        if (nl >= data.len) {
            state.appendPendingFragment(data, file_offset);
            return;
        }

        state.appendPendingFragment(data[0..nl], file_offset);
        state.pending_line_bytes += 1;
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
                if (state.name.items.len > max_index_name_len or
                    fragment.len > max_index_name_len - state.name.items.len)
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

            const header_end = findNextStreamingNewline(data, i, data.len);
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
                STREAM_VECTOR_LEN,
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
                if (findNextStreamingNewline(data, i, data.len) >= data.len) {
                    state.appendPendingFragment(
                        data[i..],
                        file_offset + @as(u64, @intCast(i)),
                    );
                    return;
                }
            }
        }

        const line_start = i;
        const nl = findNextStreamingNewline(data, i, data.len);
        if (nl >= data.len) {
            state.appendPendingFragment(
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
    enable_dedup: bool,
    max_name_len: ?usize,
    allocator: std.mem.Allocator,
    zfi_mode: bool,
    external_dedup: ?*NameDedup,
    ctx: anytype,
    comptime emitRecord: fn (@TypeOf(ctx), FastaRecordEmit) anyerror!void,
) !u32 {
    // Stream FAI: arena-backed NameDedup. Stream `.zfi`: caller passes blob-backed dedup
    // via `external_dedup` so catalog lives once in `name_blob`.
    var owned_dedup: ?NameDedup = null;
    const dedup_ptr: ?*NameDedup = blk: {
        if (external_dedup) |p| break :blk p;
        if (enable_dedup) {
            owned_dedup = NameDedup.init(allocator, true);
            try owned_dedup.?.ensureStreamCapacity(NameDedup.stream_dedup_pre_cap);
            break :blk &owned_dedup.?;
        }
        break :blk null;
    };
    defer if (owned_dedup) |*seen| seen.deinit();

    var state = ChunkParseState{ .allocator = allocator, .zfi_mode = zfi_mode };
    defer state.name.deinit(allocator);
    defer if (zfi_mode) state.side_table.deinit(allocator);
    var record_count: u32 = 0;
    var file_offset: u64 = 0;

    while (true) {
        const n = reader.readSliceShort(read_buf) catch |err| return err;
        if (n == 0) break;
        try processChunkBytes(
            &state,
            read_buf[0..n],
            file_offset,
            max_name_len,
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
    buf: [65536]u8 = undefined,
    len: usize = 0,

    fn flush(self: *FaiEmitBuffer) !void {
        if (self.len == 0) return;
        try self.writer.writeAll(self.buf[0..self.len]);
        self.len = 0;
    }

    fn emitRecord(self: *FaiEmitBuffer, emit_info: FastaRecordEmit) !void {
        if (!emit_info.uses_uniform_formula) return error.NonUniformFai;
        const rec = emit_info.record;
        var suffix_buf: [128]u8 = undefined;
        const suffix = try std.fmt.bufPrint(
            &suffix_buf,
            "\t{d}\t{d}\t{d}\t{d}\n",
            .{ rec.seq_len, rec.seq_offset, rec.line_bases, rec.line_bytes },
        );
        const total = emit_info.name.len + suffix.len;
        if (total > self.buf.len) {
            try self.flush();
            try self.writer.writeAll(emit_info.name);
            try self.writer.writeAll(suffix);
            return;
        }
        if (self.len + total > self.buf.len) try self.flush();
        @memcpy(self.buf[self.len..][0..emit_info.name.len], emit_info.name);
        self.len += emit_info.name.len;
        @memcpy(self.buf[self.len..][0..suffix.len], suffix);
        self.len += suffix.len;
    }
};

fn scanChunkedReader(
    reader: *std.Io.Reader,
    read_buf: []u8,
    writer: anytype,
    enable_dedup: bool,
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
    const count = try scanFastaReader(reader, read_buf, enable_dedup, null, allocator, false, null, &ctx, Ctx.emit);
    try fai_buf.flush();
    return count;
}

pub fn scanChunkedData(
    data: []const u8,
    writer: anytype,
    enable_dedup: bool,
    allocator: std.mem.Allocator,
) !u32 {
    var r = std.Io.Reader.fixed(data);
    var read_buf: [low_mem_chunk_size]u8 = undefined;
    return scanChunkedReader(&r, &read_buf, writer, enable_dedup, allocator);
}

pub fn runChunkedMode(io: std.Io, path: []const u8, enable_dedup: bool) void {
    runIndexLowMem(io, path, true, enable_dedup);
}

pub fn runIndexLowMem(
    io: std.Io,
    path: []const u8,
    emit_fai: bool,
    enable_dedup: bool,
) void {
    const err_exit = index_format.printErrorAndExit;

    const file = std.Io.Dir.cwd().openFile(io, path, .{}) catch |err| {
        switch (err) {
            error.FileNotFound => err_exit("error: file not found: {s}\n", .{path}),
            error.AccessDenied => err_exit("error: access denied: {s}\n", .{path}),
            else => err_exit("error: failed to open file: {s}\n", .{path}),
        }
    };
    defer file.close(io);

    const stat = file.stat(io) catch {
        err_exit("error: failed to stat file: {s}\n", .{path});
    };
    if (stat.size == 0) {
        err_exit("error: file is empty: {s}\n", .{path});
    }

    var first_byte: [1]u8 = undefined;
    const first_read = file.readPositional(io, &.{&first_byte}, 0) catch {
        err_exit("error: failed to read file: {s}\n", .{path});
    };
    if (first_read == 0 or first_byte[0] != '>') {
        err_exit("error: not a FASTA file: {s}\n", .{path});
    }

    // Must free on growth: ArrayList/HashMap realloc under an arena keeps every old
    // buffer until process end (Tx `--low-mem` `.zfi` was ~3× payload RSS). Stream FAI
    // still uses a child arena inside `NameDedup` for stable name keys.
    const allocator = std.heap.page_allocator;

    var io_buf: [FILE_IO_BUF_SIZE]u8 = undefined;
    var read_buf: [low_mem_chunk_size]u8 = undefined;
    var file_reader = file.reader(io, &io_buf);

    if (emit_fai) {
        var out_buf: [65536]u8 = undefined;
        var stdout_fw = std.Io.File.Writer.initStreaming(.stdout(), io, &out_buf);

        // Spill to a sibling temp file so a mid-file NonUniformFai does not leave a partial
        // `.fai` on stdout and so Tx-scale catalogs do not keep the whole FAI text in RSS.
        var fai_tmp_buf: [4096]u8 = undefined;
        const fai_tmp_path = std.fmt.bufPrint(&fai_tmp_buf, "{s}.fai.stdout.tmp", .{path}) catch {
            err_exit("error: path too long\n", .{});
        };
        const cwd = std.Io.Dir.cwd();
        cwd.deleteFile(io, fai_tmp_path) catch {};

        const tmp_file = cwd.createFile(io, fai_tmp_path, .{}) catch {
            err_exit("error: failed to create temp fai: {s}\n", .{fai_tmp_path});
        };
        var tmp_open = true;
        defer {
            if (tmp_open) tmp_file.close(io);
            cwd.deleteFile(io, fai_tmp_path) catch {};
        }

        var tmp_io_buf: [65536]u8 = undefined;
        var tmp_fw = tmp_file.writer(io, &tmp_io_buf);

        const record_count = scanChunkedReader(&file_reader.interface, &read_buf, &tmp_fw.interface, enable_dedup, allocator) catch |err| switch (err) {
            error.HeaderTooLong => err_exit("error: sequence name exceeds {d} bytes: {s}\n", .{ max_index_name_len, path }),
            error.NonUniformFai => err_exit(
                "error: cannot emit .fai for non-uniform sequence layout; run 'z-fasta index' (default) to write .zfi\n",
                .{},
            ),
            else => err_exit("error: processing failed\n", .{}),
        };
        tmp_fw.interface.flush() catch {
            err_exit("error: write failed\n", .{});
        };
        // Streaming writer mode on `tmp_file` can block positional reads; reopen for replay.
        tmp_file.close(io);
        tmp_open = false;

        if (record_count == 0) {
            err_exit("error: no valid sequences found in: {s}\n", .{path});
        }

        const read_file = cwd.openFile(io, fai_tmp_path, .{}) catch {
            err_exit("error: failed to reopen temp fai\n", .{});
        };
        defer read_file.close(io);

        const tmp_size = (read_file.stat(io) catch {
            err_exit("error: failed to stat temp fai\n", .{});
        }).size;
        var offset: u64 = 0;
        var copy_buf: [65536]u8 = undefined;
        while (offset < tmp_size) {
            const want: usize = @intCast(@min(copy_buf.len, tmp_size - offset));
            const n = std.Io.File.readPositionalAll(read_file, io, copy_buf[0..want], offset) catch {
                err_exit("error: failed to read temp fai\n", .{});
            };
            if (n == 0) break;
            stdout_fw.interface.writeAll(copy_buf[0..n]) catch {
                err_exit("error: write failed\n", .{});
            };
            offset += n;
        }
        stdout_fw.flush() catch {};
        return;
    }

    var zfi_path_buf: [4096]u8 = undefined;
    const zfi_path = std.fmt.bufPrint(&zfi_path_buf, "{s}.zfi", .{path}) catch {
        err_exit("error: path too long\n", .{});
    };

    var zfi_tmp_buf: [4096]u8 = undefined;
    const zfi_tmp_path = std.fmt.bufPrint(&zfi_tmp_buf, "{s}.zfi.tmp", .{path}) catch {
        err_exit("error: path too long\n", .{});
    };

    const cwd = std.Io.Dir.cwd();
    cwd.deleteFile(io, zfi_tmp_path) catch {};

    var zfi_index = scanZfiIndexStreaming(&file_reader.interface, &read_buf, enable_dedup, allocator) catch |err| switch (err) {
        error.HeaderTooLong => err_exit("error: sequence name exceeds {d} bytes: {s}\n", .{ max_index_name_len, path }),
        error.MissingSideTable => err_exit("error: scan failed\n", .{}),
        else => err_exit("error: scan failed\n", .{}),
    };
    defer zfi_index.deinit(allocator);

    if (zfi_index.records.items.len == 0) {
        err_exit("error: no valid sequences found in: {s}\n", .{path});
    }

    writeZfiIndexFile(io, zfi_tmp_path, &zfi_index, stat.size, index_format.timestampToNs(stat.mtime)) catch {
        cwd.deleteFile(io, zfi_tmp_path) catch {};
        err_exit("error: write failed\n", .{});
    };

    cwd.rename(zfi_tmp_path, cwd, zfi_path, io) catch {
        cwd.deleteFile(io, zfi_tmp_path) catch {};
        err_exit("error: failed to finalize index: {s}\n", .{zfi_path});
    };

    std.debug.print("wrote {s} ({d} sequences)\n", .{ zfi_path, zfi_index.records.items.len });
}

/// Validates that the data is a FASTA file (starts with '>')
pub fn validateFasta(data: []const u8) bool {
    return data.len > 0 and data[0] == '>';
}

test "stride block validation matches scalar for structural mutations" {
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
                    STREAM_VECTOR_LEN,
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
