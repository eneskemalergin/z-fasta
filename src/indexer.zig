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
            for (0..SIMD_CHUNK_SIZE) |j| {
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
        for (0..SIMD_CHUNK_SIZE) |j| {
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
    var records = std.ArrayList(IndexRecord).init(allocator);
    errdefer records.deinit();

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

        try records.append(IndexRecord{
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
    name_buf: [4096]u8 = undefined,
    name_len: u16 = 0,
    seq_offset: u64 = 0,
    seq_len: u64 = 0,
    line_bases: u32 = 0,
    line_bytes: u32 = 0,
    first_line_seen: bool = false,
    in_sequence: bool = false,
    is_duplicate: bool = false,
    file_offset: u64 = 0,
    record_count: u32 = 0,
};

fn countLineBases(line: []const u8) u32 {
    var count: u32 = 0;
    for (line) |c| {
        if (c > ' ') count += 1;
    }
    return count;
}

fn processChunk(
    buffer: []const u8,
    state: *ChunkState,
    writer: anytype,
    is_last_chunk: bool,
    seen_names: *std.StringHashMap(void),
    allocator: std.mem.Allocator,
) !void {
    var pos: usize = 0;

    while (pos < buffer.len) {
        if (buffer[pos] == '>') {
            if (state.in_sequence and state.seq_len > 0 and !state.is_duplicate) {
                try writer.print("{s}\t{d}\t{d}\t{d}\t{d}\n", .{
                    state.name_buf[0..state.name_len],
                    state.seq_len,
                    state.seq_offset,
                    state.line_bases,
                    state.line_bytes,
                });
                state.record_count += 1;
            }

            const header_start = pos + 1;
            var name_end = header_start;
            while (name_end < buffer.len and
                buffer[name_end] != ' ' and
                buffer[name_end] != '\t' and
                buffer[name_end] != '\n' and
                buffer[name_end] != '\r')
            {
                name_end += 1;
            }

            const name_len = @min(name_end - header_start, state.name_buf.len);
            @memcpy(state.name_buf[0..name_len], buffer[header_start..][0..name_len]);
            state.name_len = @intCast(name_len);

            const name_slice = state.name_buf[0..name_len];
            const name_copy = try allocator.dupe(u8, name_slice);
            const gop = try seen_names.getOrPut(name_copy);
            if (gop.found_existing) {
                allocator.free(name_copy);
                state.is_duplicate = true;
            } else {
                state.is_duplicate = false;
            }

            while (pos < buffer.len and buffer[pos] != '\n') pos += 1;
            if (pos < buffer.len) pos += 1;

            state.seq_offset = state.file_offset + pos;
            state.seq_len = 0;
            state.line_bases = 0;
            state.line_bytes = 0;
            state.first_line_seen = false;
            state.in_sequence = true;
        } else if (state.in_sequence) {
            const line_start = pos;
            while (pos < buffer.len and buffer[pos] != '\n') pos += 1;
            const line_end = pos;
            const line = buffer[line_start..line_end];
            const bases = countLineBases(line);
            state.seq_len += bases;

            if (!state.first_line_seen and bases > 0) {
                state.line_bases = bases;
                state.line_bytes = @intCast(line_end - line_start + 1);
                state.first_line_seen = true;
            }
            if (pos < buffer.len) pos += 1;
        } else {
            pos += 1;
        }
    }

    state.file_offset += buffer.len;

    if (is_last_chunk and state.in_sequence and state.seq_len > 0 and !state.is_duplicate) {
        try writer.print("{s}\t{d}\t{d}\t{d}\t{d}\n", .{
            state.name_buf[0..state.name_len],
            state.seq_len,
            state.seq_offset,
            state.line_bases,
            state.line_bytes,
        });
        state.record_count += 1;
    }
}

pub fn runChunkedMode(path: []const u8) void {
    const err_exit = index_format.printErrorAndExit;

    const file = std.fs.cwd().openFile(path, .{}) catch |err| {
        switch (err) {
            error.FileNotFound => err_exit("error: file not found: {s}\n", .{path}),
            error.AccessDenied => err_exit("error: access denied: {s}\n", .{path}),
            else => err_exit("error: failed to open file: {s}\n", .{path}),
        }
    };
    defer file.close();

    var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var seen_names = std.StringHashMap(void).init(allocator);
    defer seen_names.deinit();

    var buffer: [CHUNK_SIZE]u8 = undefined;
    var state = ChunkState{};
    var stdout_buffered = std.io.bufferedWriter(std.io.getStdOut().writer());
    const writer = stdout_buffered.writer();

    while (true) {
        const bytes_read = file.read(&buffer) catch {
            err_exit("error: read failed\n", .{});
        };
        if (bytes_read == 0) break;

        const is_last = bytes_read < CHUNK_SIZE;
        processChunk(buffer[0..bytes_read], &state, writer, is_last, &seen_names, allocator) catch {
            err_exit("error: processing failed\n", .{});
        };
        if (is_last) break;
    }
    stdout_buffered.flush() catch {};

    if (state.record_count == 0) {
        err_exit("error: no valid sequences found in: {s}\n", .{path});
    }
}

/// Validates that the data is a FASTA file (starts with '>')
pub fn validateFasta(data: []const u8) bool {
    return data.len > 0 and data[0] == '>';
}
