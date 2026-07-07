//! FASTA indexing: mmap slice scan and chunked `--low-mem` reader.
//!
//! `scanFastaRecords` walks a byte slice and invokes a per-record callback. Offsets in
//! `IndexRecord` and side-table lines are absolute positions in the FASTA file.
//!
//! Shared line-metrics rules: `countBasesInLineSlice`, `LineMetricsBuilder`, and
//! `scanLineMetrics`. Mmap record length uses `measureSequenceLength` on the full
//! sequence slice; streaming FAI sums per-line metrics from `ChunkParseState`.

const std = @import("std");
const index_format = @import("index_format.zig");

pub const IndexRecord = index_format.IndexRecord;
pub const ZfiHeader = index_format.ZfiHeader;
pub const ZFI_MAGIC = index_format.ZFI_MAGIC;
pub const writeZfi = index_format.writeZfi;

const SIMD_CHUNK_SIZE = 32;
const SimdVec = @Vector(SIMD_CHUNK_SIZE, u8);
pub const low_mem_chunk_size = 1 * 1024 * 1024;
const FILE_IO_BUF_SIZE = 8 * 1024;
pub const max_index_name_len = std.math.maxInt(u16);

const DedupSeen = std.HashMap(
    u64,
    void,
    std.hash_map.AutoContext(u64),
    std.hash_map.default_max_load_percentage,
);

fn hashSequenceName(name: []const u8) u64 {
    return std.hash.Wyhash.hash(0, name);
}

pub const OutputMode = enum { fai, zfi };

/// One parsed FASTA record passed to `scanFastaRecords` emit callbacks.
pub const FastaRecordEmit = struct {
    record: IndexRecord,
    name: []const u8,
    /// Sequence region in the scan buffer (header newline through next header or EOF).
    seq_data: []const u8,
    uses_uniform_formula: bool,
    /// Incremental side-table bytes for streaming `.zfi` (line count + `SideTableLine` rows).
    streaming_side_table: []const u8 = &.{},
};

pub const ZfiIndex = struct {
    records: std.ArrayList(IndexRecord),
    side_tables: std.ArrayList(u8),

    pub fn deinit(self: *ZfiIndex, allocator: std.mem.Allocator) void {
        self.records.deinit(allocator);
        self.side_tables.deinit(allocator);
    }
};

const LineMetrics = struct {
    seq_len: u64 = 0,
    line_bases: u32 = 0,
    line_bytes: u32 = 0,
    is_uniform_width: bool = true,
    line_count: u64 = 0,
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
    return (sep_len == 1 or sep_len == 2) and line_bases + sep_len == line_bytes;
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

const LineMetricsBuilder = struct {
    metrics: LineMetrics = .{},
    first_actual_bytes: u32 = 0,
    pending_bases: u32 = 0,
    pending_actual_bytes: u32 = 0,
    have_pending: bool = false,

    fn ingestLine(self: *LineMetricsBuilder, bases: u32, actual_bytes: u32, has_lf: bool) void {
        if (bases == 0) return;
        self.metrics.seq_len += bases;
        self.metrics.line_count += 1;

        if (self.metrics.line_count == 1) {
            self.first_actual_bytes = actual_bytes;
            self.metrics.line_bases = bases;
            self.metrics.line_bytes = if (has_lf) actual_bytes else actual_bytes + 1;
        } else if (self.have_pending) {
            if (self.pending_bases != self.metrics.line_bases or
                self.pending_actual_bytes != self.first_actual_bytes)
            {
                self.metrics.is_uniform_width = false;
            }
        }

        self.pending_bases = bases;
        self.pending_actual_bytes = actual_bytes;
        self.have_pending = true;
    }

    fn finish(self: *LineMetricsBuilder) LineMetrics {
        if (self.metrics.line_count > 1 and self.have_pending) {
            if (self.pending_bases > self.metrics.line_bases or
                self.pending_actual_bytes > self.first_actual_bytes)
            {
                self.metrics.is_uniform_width = false;
            }
        }
        return self.metrics;
    }
};

fn sequenceLengthInfoFromMetrics(metrics: LineMetrics) SequenceLengthInfo {
    if (metrics.seq_len == 0) return .{ .seq_len = 0, .uses_uniform_formula = false };
    if (!metrics.is_uniform_width or metrics.line_bases == 0) {
        return .{ .seq_len = metrics.seq_len, .uses_uniform_formula = false };
    }
    return .{
        .seq_len = metrics.seq_len,
        .uses_uniform_formula = firstLineIsDense(metrics.line_bases, metrics.line_bytes),
    };
}

fn scanLineMetrics(data: []const u8, seq_offset: usize, seq_end: usize) LineMetrics {
    var builder = LineMetricsBuilder{};
    var pos = seq_offset;

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
        builder.ingestLine(bases, actual_bytes, has_lf);

        if (!has_lf) break;
        pos = line_end + 1;
    }

    return builder.finish();
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
        builder.ingestLine(bases, actual_bytes, has_lf);

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
    dedup_names: ?*std.StringHashMap(void),
    max_name_len: ?usize,
    allocator: std.mem.Allocator,
    ctx: anytype,
    comptime emitRecord: fn (@TypeOf(ctx), FastaRecordEmit) anyerror!void,
) !u32 {
    var local_dedup: ?std.StringHashMap(void) = null;
    if (enable_dedup and dedup_names == null) {
        local_dedup = std.StringHashMap(void).init(allocator);
    }
    defer if (local_dedup) |*s| s.deinit();

    const seen_names: ?*std.StringHashMap(void) = if (enable_dedup)
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
        if (max_name_len) |limit| {
            if (name_end - name_local > limit) return error.HeaderTooLong;
        }
        const name_len: u16 = @intCast(name_end - name_local);

        var header_end = name_end;
        while (header_end < data.len and data[header_end] != '\n') {
            header_end += 1;
        }

        const seq_local: usize = if (header_end < data.len) header_end + 1 else header_end;
        const seq_end = findNextHeaderStart(data, seq_local);

        var line_bases: u32 = 0;
        var line_bytes: u32 = 0;
        if (seq_local < seq_end) {
            const first_line_end = findNextNewline(data, seq_local, seq_end);
            const has_lf = first_line_end < seq_end;
            var content_end = first_line_end;
            if (content_end > seq_local and data[content_end - 1] == '\r') {
                content_end -= 1;
            }
            line_bases = countBasesInLineSlice(data, seq_local, content_end);
            line_bytes = if (has_lf)
                @intCast((first_line_end + 1) - seq_local)
            else
                @intCast((first_line_end - seq_local) + 1);
        }

        const seq_data = data[seq_local..seq_end];
        const seq_info = measureSequenceLength(seq_data, line_bases, line_bytes);
        const seq_len = seq_info.seq_len;

        pos = seq_end;
        if (seq_len == 0) continue;

        const name = data[name_local..][0..name_len];
        if (seen_names) |seen| {
            const gop = try seen.getOrPut(name);
            if (gop.found_existing) continue;
        }

        const rec = IndexRecord{
            .name_offset = file_base_offset + @as(u64, @intCast(name_local)),
            .name_len = name_len,
            .seq_offset = file_base_offset + @as(u64, @intCast(seq_local)),
            .seq_len = seq_len,
            .line_bases = line_bases,
            .line_bytes = line_bytes,
        };
        try emitRecord(ctx, .{
            .record = rec,
            .name = name,
            .seq_data = seq_data,
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

pub fn scanZfiIndex(data: []const u8, enable_dedup: bool, allocator: std.mem.Allocator) !ZfiIndex {
    var index = ZfiIndex{
        .records = .empty,
        .side_tables = .empty,
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
                const seq_local: usize = @intCast(rec.seq_offset - ctx.file_base_offset);
                const local_offset = try appendSideTable(
                    ctx.data,
                    seq_local,
                    seq_local + emit_info.seq_data.len,
                    ctx.file_base_offset,
                    &ctx.index.side_tables,
                    ctx.allocator,
                );
                try rec.markNonUniform(local_offset);
            }
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
    };
    errdefer index.deinit(allocator);

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
            try ctx.index.records.append(ctx.allocator, rec);
        }
    };

    var ctx = Ctx{ .index = &index, .allocator = allocator };
    _ = try scanFastaReader(reader, read_buf, enable_dedup, null, allocator, true, &ctx, Ctx.emit);
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

/// Writes a `ZfiIndex` to a `.zfi` file (header, records, side tables).
pub fn writeZfiIndexFile(
    io: std.Io,
    path: []const u8,
    index: *const ZfiIndex,
    source_size: u64,
) !void {
    const file = try std.Io.Dir.cwd().createFile(io, path, .{});
    defer file.close(io);

    var file_buf: [65536]u8 = undefined;
    var file_fw = file.writer(io, &file_buf);
    const writer = &file_fw.interface;

    const header = ZfiHeader{
        .magic = ZFI_MAGIC,
        .record_count = @intCast(index.records.items.len),
        .source_size = source_size,
    };
    try writer.writeAll(std.mem.asBytes(&header));
    try writer.writeAll(std.mem.sliceAsBytes(index.records.items));
    try writer.writeAll(index.side_tables.items);
    try file_fw.flush();
}

/// Serializes a `ZfiIndex` to the on-disk `.zfi` byte layout (header, records, side tables).
pub fn zfiIndexToBytes(index: *const ZfiIndex, source_size: u64, allocator: std.mem.Allocator) ![]u8 {
    const header = ZfiHeader{
        .magic = ZFI_MAGIC,
        .record_count = @intCast(index.records.items.len),
        .source_size = source_size,
    };
    const records_bytes = std.mem.sliceAsBytes(index.records.items);
    const out = try allocator.alloc(u8, @sizeOf(ZfiHeader) + records_bytes.len + index.side_tables.items.len);
    @memcpy(out[0..@sizeOf(ZfiHeader)], std.mem.asBytes(&header));
    @memcpy(out[@sizeOf(ZfiHeader) .. @sizeOf(ZfiHeader) + records_bytes.len], records_bytes);
    @memcpy(out[@sizeOf(ZfiHeader) + records_bytes.len ..], index.side_tables.items);
    return out;
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
    is_duplicate: bool = false,
    at_line_start: bool = true,
    header_start_offset: u64 = 0,
    name: std.ArrayList(u8) = .empty,
    seq_offset: u64 = 0,
    seq_offset_set: bool = false,
    line_builder: LineMetricsBuilder = .{},
    in_pending_line: bool = false,
    pending_line_bases: u32 = 0,
    pending_line_bytes: u32 = 0,
    pending_line_file_offset: u64 = 0,
    side_table: std.ArrayList(u8) = .empty,
    side_base_start: u64 = 0,
    side_table_active: bool = false,
    deferred_side_line: ?struct {
        bases: u32,
        actual_bytes: u32,
        line_file_offset: u64,
    } = null,

    fn ensureSideTableHeader(self: *ChunkParseState) !void {
        if (self.side_table.items.len >= @sizeOf(u64)) return;
        var line_count: u64 = 0;
        try self.side_table.appendSlice(self.allocator, std.mem.asBytes(&line_count));
    }

    fn activateSideTable(self: *ChunkParseState) !void {
        if (self.side_table_active) return;
        self.side_table_active = true;
        try self.ensureSideTableHeader();
        if (self.deferred_side_line) |deferred| {
            try self.appendSideTableLine(deferred.bases, deferred.actual_bytes, deferred.line_file_offset);
            self.deferred_side_line = null;
        }
    }

    fn startHeader(self: *ChunkParseState, header_start_offset: u64) !void {
        self.active = true;
        self.in_header = true;
        self.parsing_name = true;
        self.is_duplicate = false;
        self.at_line_start = false;
        self.header_start_offset = header_start_offset;
        self.name.clearRetainingCapacity();
        self.seq_offset = 0;
        self.seq_offset_set = false;
        self.line_builder = .{};
        self.in_pending_line = false;
        self.pending_line_bases = 0;
        self.pending_line_bytes = 0;
        self.pending_line_file_offset = 0;
        self.side_base_start = 0;
        self.side_table.clearRetainingCapacity();
        self.side_table_active = false;
        self.deferred_side_line = null;
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

    fn recordSequenceLine(
        self: *ChunkParseState,
        bases: u32,
        actual_bytes: u32,
        has_lf: bool,
        line_file_offset: u64,
    ) !void {
        self.line_builder.ingestLine(bases, actual_bytes, has_lf);

        if (!self.zfi_mode or bases == 0) return;

        const metrics = self.line_builder.metrics;
        const matches_stride = metrics.line_bases > 0 and
            bases == metrics.line_bases and
            actual_bytes == self.line_builder.first_actual_bytes;
        const can_use_formula = metrics.is_uniform_width and matches_stride and
            firstLineIsDense(metrics.line_bases, metrics.line_bytes);

        if (can_use_formula) {
            if (metrics.line_count == 1) {
                self.deferred_side_line = .{
                    .bases = bases,
                    .actual_bytes = actual_bytes,
                    .line_file_offset = line_file_offset,
                };
            } else {
                self.deferred_side_line = null;
            }
            return;
        }

        try self.activateSideTable();
        try self.appendSideTableLine(bases, actual_bytes, line_file_offset);
    }

    fn noteSeqOffset(self: *ChunkParseState, data: []const u8, line_file_offset: u64) void {
        if (self.seq_offset_set) return;
        var pos: usize = 0;
        while (pos < data.len and data[pos] <= ' ') pos += 1;
        if (pos < data.len) {
            self.seq_offset = line_file_offset + @as(u64, @intCast(pos));
            self.seq_offset_set = true;
        }
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
        const bases = countBasesInLineSlice(data, line_start, content_end);
        if (bases > 0 and !self.seq_offset_set) {
            self.noteSeqOffset(data[line_start..content_end], line_file_offset);
        }

        try self.recordSequenceLine(bases, actual_bytes, has_lf, line_file_offset);
        self.at_line_start = has_lf;
        self.in_pending_line = false;
        self.pending_line_bases = 0;
        self.pending_line_bytes = 0;
    }

    fn appendPendingFragment(self: *ChunkParseState, fragment: []const u8, fragment_file_offset: u64) void {
        if (!self.in_pending_line) {
            self.in_pending_line = true;
            self.pending_line_bases = 0;
            self.pending_line_bytes = 0;
            self.pending_line_file_offset = fragment_file_offset;
        }
        self.noteSeqOffset(fragment, fragment_file_offset);
        self.pending_line_bases += countBasesInLineSlice(fragment, 0, fragment.len);
        self.pending_line_bytes += @intCast(fragment.len);
    }

    fn commitPendingLine(self: *ChunkParseState, has_lf: bool) !void {
        if (!self.in_pending_line and self.pending_line_bytes == 0) return;
        try self.recordSequenceLine(
            self.pending_line_bases,
            self.pending_line_bytes,
            has_lf,
            self.pending_line_file_offset,
        );
        self.at_line_start = has_lf;
        self.in_pending_line = false;
        self.pending_line_bases = 0;
        self.pending_line_bytes = 0;
    }

    fn finalizeHeader(self: *ChunkParseState, seen_names: ?*DedupSeen) !void {
        self.in_header = false;
        const name = self.name.items;
        if (seen_names) |seen| {
            const key = hashSequenceName(name);
            const gop = try seen.getOrPut(key);
            if (gop.found_existing) self.is_duplicate = true;
        }
        self.at_line_start = true;
    }

    fn emitRecordIfReady(
        self: *ChunkParseState,
        ctx: anytype,
        comptime emitRecord: fn (@TypeOf(ctx), FastaRecordEmit) anyerror!void,
        record_count: *u32,
    ) !void {
        try self.commitPendingLine(false);
        const metrics = self.line_builder.finish();
        const seq_info = sequenceLengthInfoFromMetrics(metrics);
        if (seq_info.seq_len == 0 or self.is_duplicate) return;

        const name_len: u16 = @intCast(self.name.items.len);
        var side_table_slice: []const u8 = &.{};
        if (!seq_info.uses_uniform_formula and self.zfi_mode) {
            if (!self.side_table_active) {
                return error.MissingSideTable;
            }
            @memcpy(self.side_table.items[0..@sizeOf(u64)], std.mem.asBytes(&metrics.line_count));
            side_table_slice = self.side_table.items;
        }
        try emitRecord(ctx, .{
            .record = .{
                .name_offset = self.header_start_offset + 1,
                .name_len = name_len,
                .seq_offset = self.seq_offset,
                .seq_len = seq_info.seq_len,
                .line_bases = metrics.line_bases,
                .line_bytes = metrics.line_bytes,
            },
            .name = self.name.items,
            .seq_data = &.{},
            .uses_uniform_formula = seq_info.uses_uniform_formula,
            .streaming_side_table = side_table_slice,
        });
        record_count.* += 1;
    }
};

fn processChunkBytes(
    state: *ChunkParseState,
    data: []const u8,
    file_offset: u64,
    max_name_len: ?usize,
    seen_names: ?*DedupSeen,
    ctx: anytype,
    comptime emitRecord: fn (@TypeOf(ctx), FastaRecordEmit) anyerror!void,
    record_count: *u32,
) !void {
    var i: usize = 0;

    if (state.in_pending_line) {
        const nl = findNextNewline(data, 0, data.len);
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
                try state.emitRecordIfReady(ctx, emitRecord, record_count);
            }
            try state.startHeader(byte_offset);
            i += 1;
            continue;
        }

        if (state.in_header) {
            if (byte == '\n') {
                try state.finalizeHeader(seen_names);
                i += 1;
                continue;
            }
            if (state.parsing_name) {
                if (byte == ' ' or byte == '\t' or byte == '\r') {
                    state.parsing_name = false;
                } else {
                    if (max_name_len) |limit| {
                        if (state.name.items.len >= limit) return error.HeaderTooLong;
                    }
                    try state.name.append(state.allocator, byte);
                    if (state.name.items.len > max_index_name_len) return error.HeaderTooLong;
                }
            }
            i += 1;
            continue;
        }

        if (!state.active) {
            state.at_line_start = byte == '\n';
            i += 1;
            continue;
        }

        const line_start = i;
        const nl = findNextNewline(data, i, data.len);
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
    ctx: anytype,
    comptime emitRecord: fn (@TypeOf(ctx), FastaRecordEmit) anyerror!void,
) !u32 {
    var dedup_names: ?DedupSeen = null;
    if (enable_dedup) {
        dedup_names = DedupSeen.init(allocator);
    }
    defer if (dedup_names) |*seen| seen.deinit();

    var state = ChunkParseState{ .allocator = allocator, .zfi_mode = zfi_mode };
    defer state.name.deinit(allocator);
    defer if (zfi_mode) state.side_table.deinit(allocator);
    var record_count: u32 = 0;
    var file_offset: u64 = 0;
    const dedup_ptr: ?*DedupSeen = if (dedup_names) |*s| s else null;

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
        try state.emitRecordIfReady(ctx, emitRecord, &record_count);
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
    const count = try scanFastaReader(reader, read_buf, enable_dedup, null, allocator, false, &ctx, Ctx.emit);
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

    var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var io_buf: [FILE_IO_BUF_SIZE]u8 = undefined;
    var read_buf: [low_mem_chunk_size]u8 = undefined;
    var file_reader = file.reader(io, &io_buf);

    if (emit_fai) {
        var out_buf: [65536]u8 = undefined;
        var stdout_fw = std.Io.File.Writer.initStreaming(.stdout(), io, &out_buf);

        const record_count = scanChunkedReader(&file_reader.interface, &read_buf, &stdout_fw.interface, enable_dedup, allocator) catch |err| switch (err) {
            error.HeaderTooLong => {
                stdout_fw.flush() catch {};
                err_exit("error: sequence name exceeds {d} bytes: {s}\n", .{ max_index_name_len, path });
            },
            else => {
                stdout_fw.flush() catch {};
                err_exit("error: processing failed\n", .{});
            },
        };
        stdout_fw.flush() catch {};

        if (record_count == 0) {
            stdout_fw.flush() catch {};
            err_exit("error: no valid sequences found in: {s}\n", .{path});
        }
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

    writeZfiIndexFile(io, zfi_tmp_path, &zfi_index, stat.size) catch {
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
