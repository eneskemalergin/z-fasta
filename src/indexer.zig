const std = @import("std");
const posix = std.posix;
const index_format = @import("index_format.zig");

pub const IndexRecord = index_format.IndexRecord;
pub const ZfiHeader = index_format.ZfiHeader;
pub const ZFI_MAGIC = index_format.ZFI_MAGIC;
pub const writeZfi = index_format.writeZfi;

// ============================================================================
// SIMD constants
// ============================================================================

const SIMD_CHUNK_SIZE = 32;
const SimdVec = @Vector(SIMD_CHUNK_SIZE, u8);
const CHUNK_SIZE = 4 * 1024 * 1024; // 4 MB buffer for --low-mem mode
const MAX_LOW_MEM_NAME_LEN = 4096;

pub const OutputMode = enum { fai, zfi };

// ============================================================================
// SIMD Helpers
// ============================================================================

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

// ============================================================================
// Streaming Mode (mmap, default)
// ============================================================================

pub fn streamingScan(
    data: []const u8,
    writer: anytype,
    mode: OutputMode,
    enable_dedup: bool,
    allocator: std.mem.Allocator,
) !u32 {
    var seen_names: ?std.StringHashMap(void) = null;
    if (enable_dedup) {
        seen_names = std.StringHashMap(void).init(allocator);
    }
    defer if (seen_names) |*s| s.deinit();

    var record_count: u32 = 0;
    var pos: usize = 0;

    while (pos < data.len) {
        pos = findNextHeaderStart(data, pos);
        if (pos >= data.len) break;

        const name_offset = pos + 1;
        var name_end = name_offset;
        while (name_end < data.len and
            data[name_end] != ' ' and
            data[name_end] != '\t' and
            data[name_end] != '\n' and
            data[name_end] != '\r')
        {
            name_end += 1;
        }
        const name_len: u16 = @intCast(name_end - name_offset);

        var header_end = name_end;
        while (header_end < data.len and data[header_end] != '\n') {
            header_end += 1;
        }

        const seq_offset: u64 = if (header_end < data.len) header_end + 1 else header_end;
        const seq_end = findNextHeaderStart(data, @intCast(seq_offset));
        const seq_len = countBases(data[@intCast(seq_offset)..seq_end]);

        if (seq_len == 0) {
            pos = seq_end;
            continue;
        }

        // Line metrics
        var line_bases: u32 = 0;
        var line_bytes: u32 = 0;
        if (seq_offset < seq_end) {
            var first_line_end: usize = @intCast(seq_offset);
            while (first_line_end < seq_end and data[first_line_end] != '\n') {
                first_line_end += 1;
            }
            var base_count: u32 = 0;
            var j: usize = @intCast(seq_offset);
            while (j < first_line_end) : (j += 1) {
                if (data[j] > ' ') base_count += 1;
            }
            line_bases = base_count;
            if (first_line_end < seq_end) {
                line_bytes = @intCast((first_line_end + 1) - seq_offset);
            } else {
                line_bytes = @intCast((first_line_end - seq_offset) + 1);
            }
        }

        // Dedup check
        const name = data[name_offset..][0..name_len];
        if (seen_names) |*seen| {
            const gop = try seen.getOrPut(name);
            if (gop.found_existing) {
                pos = seq_end;
                continue;
            }
        }

        // Write record
        switch (mode) {
            .fai => {
                try writer.print("{s}\t{d}\t{d}\t{d}\t{d}\n", .{
                    name, seq_len, seq_offset, line_bases, line_bytes,
                });
            },
            .zfi => {
                const rec = IndexRecord{
                    .name_offset = @intCast(name_offset),
                    .name_len = name_len,
                    .seq_offset = seq_offset,
                    .seq_len = seq_len,
                    .line_bases = line_bases,
                    .line_bytes = line_bytes,
                };
                try writer.writeAll(std.mem.asBytes(&rec));
            },
        }
        record_count += 1;
        pos = seq_end;
    }
    return record_count;
}

/// Scans FASTA data and returns index records as ArrayList (for testing).
/// Use streamingScan for production (lower memory).
pub fn scanHeaders(data: []const u8, allocator: std.mem.Allocator) !std.ArrayList(IndexRecord) {
    var records: std.ArrayList(IndexRecord) = .empty;
    errdefer records.deinit(allocator);

    var seen_names = std.StringHashMap(void).init(allocator);
    defer seen_names.deinit();

    var pos: usize = 0;
    while (pos < data.len) {
        pos = findNextHeaderStart(data, pos);
        if (pos >= data.len) break;

        const name_offset = pos + 1;
        var name_end = name_offset;
        while (name_end < data.len and
            data[name_end] != ' ' and
            data[name_end] != '\t' and
            data[name_end] != '\n' and
            data[name_end] != '\r')
        {
            name_end += 1;
        }
        const name_len: u16 = @intCast(name_end - name_offset);

        var header_end = name_end;
        while (header_end < data.len and data[header_end] != '\n') {
            header_end += 1;
        }

        const seq_offset: u64 = if (header_end < data.len) header_end + 1 else header_end;
        const seq_end = findNextHeaderStart(data, @intCast(seq_offset));
        const seq_len = countBases(data[@intCast(seq_offset)..seq_end]);

        if (seq_len == 0) {
            pos = seq_end;
            continue;
        }

        // Line metrics
        var line_bases: u32 = 0;
        var line_bytes: u32 = 0;
        if (seq_offset < seq_end) {
            var first_line_end: usize = @intCast(seq_offset);
            while (first_line_end < seq_end and data[first_line_end] != '\n') {
                first_line_end += 1;
            }
            var base_count: u32 = 0;
            var j: usize = @intCast(seq_offset);
            while (j < first_line_end) : (j += 1) {
                if (data[j] > ' ') base_count += 1;
            }
            line_bases = base_count;
            if (first_line_end < seq_end) {
                line_bytes = @intCast((first_line_end + 1) - seq_offset);
            } else {
                line_bytes = @intCast((first_line_end - seq_offset) + 1);
            }
        }

        // Dedup check
        const name = data[name_offset..][0..name_len];
        const gop = try seen_names.getOrPut(name);
        if (gop.found_existing) {
            pos = seq_end;
            continue;
        }

        try records.append(allocator, IndexRecord{
            .name_offset = @intCast(name_offset),
            .name_len = name_len,
            .seq_offset = seq_offset,
            .seq_len = seq_len,
            .line_bases = line_bases,
            .line_bytes = line_bytes,
        });

        pos = seq_end;
    }
    return records;
}

// ============================================================================
// Chunked Mode (--low-mem, 4 MB constant memory + dedup hash set)
// ============================================================================

const ChunkState = struct {
    name_buf: [MAX_LOW_MEM_NAME_LEN]u8 = undefined,
    name_len: u16 = 0,
    seq_offset: u64 = 0,
    seq_len: u64 = 0,
    line_bases: u32 = 0,
    line_bytes: u32 = 0,
    first_line_seen: bool = false,
    has_header: bool = false,
    is_duplicate: bool = false,
    file_offset: u64 = 0,
    record_count: u32 = 0,
};

const HeaderParseState = struct {
    line_len: usize = 0,
    parsing_name: bool = true,
};

const SequenceLineState = struct {
    line_len: usize = 0,
    base_count: u32 = 0,
};

fn finalizeChunkedRecord(state: *ChunkState, writer: anytype) !void {
    if (!state.has_header or state.seq_len == 0 or state.is_duplicate) return;

    try writer.print("{s}\t{d}\t{d}\t{d}\t{d}\n", .{
        state.name_buf[0..state.name_len],
        state.seq_len,
        state.seq_offset,
        state.line_bases,
        state.line_bytes,
    });
    state.record_count += 1;
}

fn resetChunkedSequence(state: *ChunkState) void {
    state.seq_offset = 0;
    state.seq_len = 0;
    state.line_bases = 0;
    state.line_bytes = 0;
    state.first_line_seen = false;
}

fn startChunkedHeader(state: *ChunkState) void {
    state.has_header = true;
    state.name_len = 0;
    state.is_duplicate = false;
    resetChunkedSequence(state);
}

fn parseChunkedHeaderByte(
    state: *ChunkState,
    parse_state: *HeaderParseState,
    byte: u8,
) !void {
    parse_state.line_len += 1;
    if (!parse_state.parsing_name) return;

    if (byte == ' ' or byte == '\t' or byte == '\r') {
        parse_state.parsing_name = false;
        return;
    }

    if (state.name_len >= MAX_LOW_MEM_NAME_LEN) {
        return error.HeaderTooLong;
    }
    state.name_buf[state.name_len] = byte;
    state.name_len += 1;
}

fn finalizeChunkedHeader(
    state: *ChunkState,
    seen_names: *std.StringHashMap(void),
    allocator: std.mem.Allocator,
) !void {
    const name_slice = state.name_buf[0..state.name_len];
    const name_copy = try allocator.dupe(u8, name_slice);
    const gop = try seen_names.getOrPut(name_copy);
    if (gop.found_existing) {
        allocator.free(name_copy);
        state.is_duplicate = true;
    } else {
        state.is_duplicate = false;
    }
}

fn seedSequenceLine(line_state: *SequenceLineState, byte: u8) void {
    line_state.line_len = 1;
    line_state.base_count = if (byte > ' ') 1 else 0;
}

fn parseSequenceLineByte(line_state: *SequenceLineState, byte: u8) void {
    line_state.line_len += 1;
    if (byte > ' ') line_state.base_count += 1;
}

fn finalizeSequenceLine(state: *ChunkState, line_state: SequenceLineState) void {
    state.seq_len += line_state.base_count;
    if (!state.first_line_seen and line_state.base_count > 0) {
        state.line_bases = line_state.base_count;
        state.line_bytes = @intCast(line_state.line_len + 1);
        state.first_line_seen = true;
    }
}

fn scanChunkedReader(
    reader: *std.Io.Reader,
    writer: anytype,
    seen_names: *std.StringHashMap(void),
    allocator: std.mem.Allocator,
) !u32 {
    var state = ChunkState{};
    var line_start = true;

    while (true) {
        const byte = reader.takeByte() catch |err| switch (err) {
            error.EndOfStream => break,
            else => return err,
        };

        if (line_start and byte == '>') {
            try finalizeChunkedRecord(&state, writer);
            startChunkedHeader(&state);

            var header_state = HeaderParseState{};
            while (true) {
                const header_byte = reader.takeByte() catch |err| switch (err) {
                    error.EndOfStream => {
                        try finalizeChunkedHeader(&state, seen_names, allocator);
                        state.file_offset += header_state.line_len + 1;
                        return state.record_count;
                    },
                    else => return err,
                };

                if (header_byte == '\n') {
                    try finalizeChunkedHeader(&state, seen_names, allocator);
                    state.file_offset += header_state.line_len + 2;
                    line_start = true;
                    break;
                }
                try parseChunkedHeaderByte(&state, &header_state, header_byte);
            }
            continue;
        }

        if (!state.has_header) {
            line_start = byte == '\n';
            state.file_offset += 1;
            continue;
        }

        var seq_line = SequenceLineState{};
        seedSequenceLine(&seq_line, byte);

        while (true) {
            const seq_byte = reader.takeByte() catch |err| switch (err) {
                error.EndOfStream => {
                    if (seq_line.base_count > 0 and state.seq_len == 0) {
                        state.seq_offset = state.file_offset;
                    }
                    finalizeSequenceLine(&state, seq_line);
                    state.file_offset += seq_line.line_len;
                    try finalizeChunkedRecord(&state, writer);
                    return state.record_count;
                },
                else => return err,
            };

            if (seq_byte == '\n') {
                if (seq_line.base_count > 0 and state.seq_len == 0) {
                    state.seq_offset = state.file_offset;
                }
                finalizeSequenceLine(&state, seq_line);
                state.file_offset += seq_line.line_len + 1;
                line_start = true;
                break;
            }

            parseSequenceLineByte(&seq_line, seq_byte);
        }
    }

    try finalizeChunkedRecord(&state, writer);
    return state.record_count;
}

pub fn scanChunkedData(
    data: []const u8,
    writer: anytype,
    allocator: std.mem.Allocator,
) !u32 {
    var r = std.Io.Reader.fixed(data);
    var seen_names = std.StringHashMap(void).init(allocator);
    defer seen_names.deinit();

    return scanChunkedReader(&r, writer, &seen_names, allocator);
}

pub fn runChunkedMode(io: std.Io, path: []const u8) void {
    const err_exit = index_format.printErrorAndExit;

    const file = std.Io.Dir.cwd().openFile(io, path, .{}) catch |err| {
        switch (err) {
            error.FileNotFound => err_exit("error: file not found: {s}\n", .{path}),
            error.AccessDenied => err_exit("error: access denied: {s}\n", .{path}),
            else => err_exit("error: failed to open file: {s}\n", .{path}),
        }
    };
    defer file.close(io);

    var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var seen_names = std.StringHashMap(void).init(allocator);
    defer seen_names.deinit();

    var read_buf: [CHUNK_SIZE]u8 = undefined;
    var file_reader = file.reader(io, &read_buf);
    var out_buf: [65536]u8 = undefined;
    var stdout_fw = std.Io.File.Writer.initStreaming(.stdout(), io, &out_buf);

    const record_count = scanChunkedReader(&file_reader.interface, &stdout_fw.interface, &seen_names, allocator) catch |err| switch (err) {
        error.HeaderTooLong => err_exit("error: sequence name exceeds {d} bytes in --low-mem mode: {s}\n", .{ MAX_LOW_MEM_NAME_LEN, path }),
        else => err_exit("error: processing failed\n", .{}),
    };
    stdout_fw.flush() catch {};

    if (record_count == 0) {
        err_exit("error: no valid sequences found in: {s}\n", .{path});
    }
}

/// Validates that the data is a FASTA file (starts with '>')
pub fn validateFasta(data: []const u8) bool {
    return data.len > 0 and data[0] == '>';
}
