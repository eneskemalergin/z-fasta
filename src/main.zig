const std = @import("std");
const posix = std.posix;

// z-fasta: Ultra-fast FASTA indexer with minimal memory footprint.
//
// Default mode: mmap + SIMD scan, streaming output (no ArrayList)
// --dedup: Enable duplicate name filtering via HashMap
// --low-mem: Use chunked read() instead of mmap (4 MB constant memory)

const SIMD_CHUNK_SIZE = 32;
const SimdVec = @Vector(SIMD_CHUNK_SIZE, u8);
const CHUNK_SIZE = 4 * 1024 * 1024; // 4 MB buffer for --low-mem mode

/// ZFI binary format
pub const ZFI_MAGIC: [4]u8 = .{ 'Z', 'F', 'I', 0x01 };
pub const ZfiHeader = extern struct {
    magic: [4]u8,
    record_count: u32,
    source_size: u64,
};

/// Index record for ZFI output (40 bytes padded)
pub const IndexRecord = extern struct {
    name_offset: u64,
    name_len: u16,
    _pad: [6]u8 = .{0} ** 6,
    seq_offset: u64,
    seq_len: u64,
    line_bases: u32,
    line_bytes: u32,

    pub fn getName(self: IndexRecord, data: []const u8) []const u8 {
        return data[self.name_offset..][0..self.name_len];
    }
};

const OutputMode = enum { fai, zfi };

fn printErrorAndExit(comptime fmt: []const u8, args: anytype) noreturn {
    std.io.getStdErr().writer().print(fmt, args) catch {};
    std.process.exit(1);
}

/// Validates that the data is a FASTA file (starts with '>')
pub fn validateFasta(data: []const u8) bool {
    return data.len > 0 and data[0] == '>';
}

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

fn countBases(data: []const u8) u64 {
    var count: u64 = 0;
    var pos: usize = 0;
    // samtools counts only printable chars > ' ' (0x20) as bases
    // Skip: control chars (< 0x20), space (0x20), newlines, carriage returns
    while (pos + SIMD_CHUNK_SIZE <= data.len) {
        const chunk: SimdVec = data[pos..][0..SIMD_CHUNK_SIZE].*;
        const space_char: SimdVec = @splat(' ');
        var chunk_count: u32 = 0;
        for (0..SIMD_CHUNK_SIZE) |j| {
            // Count only characters > ' ' (0x20)
            if (chunk[j] > space_char[j]) chunk_count += 1;
        }
        count += chunk_count;
        pos += SIMD_CHUNK_SIZE;
    }
    while (pos < data.len) {
        // Count only characters > ' ' (0x20)
        if (data[pos] > ' ') count += 1;
        pos += 1;
    }
    return count;
}

// ============================================================================
// Streaming Mode (mmap, default)
// ============================================================================

fn streamingScan(
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
                // Count only characters > ' ' (0x20) - matches samtools behavior
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
                // Count only characters > ' ' (0x20) - matches samtools behavior
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

/// Writes the .zfi binary index file.
pub fn writeZfi(
    path: []const u8,
    records: []const IndexRecord,
    source_size: u64,
) !void {
    const file = try std.fs.cwd().createFile(path, .{});
    errdefer std.fs.cwd().deleteFile(path) catch {};
    defer file.close();

    var writer = file.writer();

    const header = ZfiHeader{
        .magic = ZFI_MAGIC,
        .record_count = @intCast(records.len),
        .source_size = source_size,
    };
    try writer.writeAll(std.mem.asBytes(&header));
    try writer.writeAll(std.mem.sliceAsBytes(records));
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
    is_duplicate: bool = false, // Track if current sequence is a duplicate
    file_offset: u64 = 0,
    record_count: u32 = 0,
};

fn countLineBases(line: []const u8) u32 {
    var count: u32 = 0;
    for (line) |c| {
        // Count only characters > ' ' (0x20) - matches samtools behavior
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
            // Flush previous sequence (if not duplicate)
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

            // Parse header
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

            // Check for duplicate
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

fn runChunkedMode(path: []const u8) void {
    const file = std.fs.cwd().openFile(path, .{}) catch |err| {
        switch (err) {
            error.FileNotFound => printErrorAndExit("error: file not found: {s}\n", .{path}),
            error.AccessDenied => printErrorAndExit("error: access denied: {s}\n", .{path}),
            else => printErrorAndExit("error: failed to open file: {s}\n", .{path}),
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
            printErrorAndExit("error: read failed\n", .{});
        };
        if (bytes_read == 0) break;

        const is_last = bytes_read < CHUNK_SIZE;
        processChunk(buffer[0..bytes_read], &state, writer, is_last, &seen_names, allocator) catch {
            printErrorAndExit("error: processing failed\n", .{});
        };
        if (is_last) break;
    }
    stdout_buffered.flush() catch {};

    // Exit with error if no valid sequences found (matches samtools behavior)
    if (state.record_count == 0) {
        printErrorAndExit("error: no valid sequences found in: {s}\n", .{path});
    }
}

// ============================================================================
// Main
// ============================================================================

const VERSION = "0.1.0";

const USAGE =
    \\usage: z-fasta <command> [options]
    \\
    \\Commands:
    \\  index    Build a FASTA index (.zfi binary or .fai text)
    \\
    \\General options:
    \\  --help       Show this help message
    \\  --version    Print version
    \\
    \\Index options:
    \\  --emit-fai   Output FAI format to stdout (default: create .zfi file)
    \\  --no-dedup   Disable duplicate name filtering (default: dedup ON)
    \\  --low-mem    Use chunked reader instead of mmap (4 MB constant memory)
    \\
    \\Examples:
    \\  z-fasta index genome.fa                  Create .zfi binary index
    \\  z-fasta index --emit-fai genome.fa       Output FAI to stdout
    \\  z-fasta index --low-mem genome.fa        Low memory mode (4 MB)
    \\  z-fasta index --no-dedup genome.fa       Allow duplicate headers
    \\
;

fn printUsageAndExit() noreturn {
    std.io.getStdErr().writer().writeAll(USAGE) catch {};
    std.process.exit(1);
}

fn printHelpAndExit() noreturn {
    std.io.getStdOut().writer().writeAll(USAGE) catch {};
    std.process.exit(0);
}

fn printVersionAndExit() noreturn {
    std.io.getStdOut().writer().writeAll("z-fasta " ++ VERSION ++ "\n") catch {};
    std.process.exit(0);
}

pub fn main() void {
    var args = std.process.args();
    _ = args.skip();

    const cmd = args.next() orelse {
        printUsageAndExit();
    };

    if (std.mem.eql(u8, cmd, "--help") or std.mem.eql(u8, cmd, "-h")) {
        printHelpAndExit();
    }
    if (std.mem.eql(u8, cmd, "--version") or std.mem.eql(u8, cmd, "-V")) {
        printVersionAndExit();
    }

    if (!std.mem.eql(u8, cmd, "index")) {
        printUsageAndExit();
    }

    var emit_fai = false;
    var enable_dedup = true; // Default: dedup enabled (matches samtools behavior)
    var low_mem = false;
    var fasta_path: ?[]const u8 = null;

    while (args.next()) |arg| {
        if (std.mem.eql(u8, arg, "--help") or std.mem.eql(u8, arg, "-h")) {
            printHelpAndExit();
        } else if (std.mem.eql(u8, arg, "--version") or std.mem.eql(u8, arg, "-V")) {
            printVersionAndExit();
        } else if (std.mem.eql(u8, arg, "--emit-fai")) {
            emit_fai = true;
        } else if (std.mem.eql(u8, arg, "--no-dedup")) {
            enable_dedup = false;
        } else if (std.mem.eql(u8, arg, "--dedup")) {
            enable_dedup = true; // Explicit dedup (for compatibility)
        } else if (std.mem.eql(u8, arg, "--low-mem")) {
            low_mem = true;
            emit_fai = true; // --low-mem only supports FAI output
        } else {
            fasta_path = arg;
        }
    }

    const path = fasta_path orelse {
        printUsageAndExit();
    };

    // Low-memory chunked mode
    if (low_mem) {
        runChunkedMode(path);
        return;
    }

    // Standard mmap mode
    const file = std.fs.cwd().openFile(path, .{}) catch |err| {
        switch (err) {
            error.FileNotFound => printErrorAndExit("error: file not found: {s}\n", .{path}),
            error.AccessDenied => printErrorAndExit("error: access denied: {s}\n", .{path}),
            else => printErrorAndExit("error: failed to open file: {s}\n", .{path}),
        }
    };
    defer file.close();

    const stat = file.stat() catch {
        printErrorAndExit("error: failed to stat file: {s}\n", .{path});
    };

    if (stat.size == 0) {
        printErrorAndExit("error: file is empty: {s}\n", .{path});
    }

    const data = posix.mmap(
        null,
        stat.size,
        posix.PROT.READ,
        .{ .TYPE = .PRIVATE },
        file.handle,
        0,
    ) catch {
        printErrorAndExit("error: failed to mmap file: {s}\n", .{path});
    };
    defer posix.munmap(data);

    posix.madvise(data.ptr, data.len, posix.MADV.SEQUENTIAL) catch {};

    if (data.len == 0 or data[0] != '>') {
        printErrorAndExit("error: not a FASTA file: {s}\n", .{path});
    }

    var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer arena.deinit();

    if (emit_fai) {
        var buffered = std.io.bufferedWriter(std.io.getStdOut().writer());
        const record_count = streamingScan(data, buffered.writer(), .fai, enable_dedup, arena.allocator()) catch {
            printErrorAndExit("error: failed to scan/write\n", .{});
        };
        buffered.flush() catch {};

        // Exit with error if no valid sequences found (matches samtools behavior)
        if (record_count == 0) {
            printErrorAndExit("error: no valid sequences found in: {s}\n", .{path});
        }
    } else {
        var zfi_path_buf: [4096]u8 = undefined;
        const zfi_path = std.fmt.bufPrint(&zfi_path_buf, "{s}.zfi", .{path}) catch {
            printErrorAndExit("error: path too long\n", .{});
        };

        const out_file = std.fs.cwd().createFile(zfi_path, .{}) catch {
            printErrorAndExit("error: cannot create: {s}\n", .{zfi_path});
        };
        defer out_file.close();

        var buffered = std.io.bufferedWriter(out_file.writer());
        const writer = buffered.writer();

        // Write dummy header (fix record_count later)
        const dummy_header = ZfiHeader{
            .magic = ZFI_MAGIC,
            .record_count = 0,
            .source_size = data.len,
        };
        writer.writeAll(std.mem.asBytes(&dummy_header)) catch {
            printErrorAndExit("error: write failed\n", .{});
        };

        const record_count = streamingScan(data, writer, .zfi, enable_dedup, arena.allocator()) catch {
            printErrorAndExit("error: scan failed\n", .{});
        };

        buffered.flush() catch {};

        // Exit with error if no valid sequences found (matches samtools behavior)
        if (record_count == 0) {
            // Clean up the empty .zfi file
            std.fs.cwd().deleteFile(zfi_path) catch {};
            printErrorAndExit("error: no valid sequences found in: {s}\n", .{path});
        }

        // Fix record_count in header
        out_file.seekTo(4) catch {};
        out_file.writer().writeInt(u32, record_count, .little) catch {};

        std.debug.print("wrote {s} ({d} sequences)\n", .{ zfi_path, record_count });
    }
}
