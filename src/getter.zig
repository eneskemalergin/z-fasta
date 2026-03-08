const std = @import("std");
const index_format = @import("index_format.zig");

const IndexRecord = index_format.IndexRecord;
const LoadedIndex = index_format.LoadedIndex;
const printErrorAndExit = index_format.printErrorAndExit;

// ============================================================================
// Region parsing
// ============================================================================

pub const Region = struct {
    name: []const u8,
    start: u64, // 1-based inclusive
    end: ?u64, // 1-based inclusive, null = to end of sequence
    is_full: bool, // true if no :START-END specified
};

/// Parse a region string. Handles the Ensembl colon trap by parsing from the right.
/// Accepted formats:
///   NAME              — full sequence
///   NAME:START-END    — 1-based, inclusive
///   NAME:START-       — from START to end
pub fn parseRegion(input: []const u8) Region {
    if (input.len == 0) {
        printErrorAndExit("error: empty region string\n", .{});
    }

    // Try to find a region suffix by scanning from the right.
    // Look for the last '-' that separates two valid integers (or START-),
    // then the ':' immediately before the START integer.

    // Find the last ':' in the string
    var colon_pos: ?usize = null;
    {
        var i: usize = input.len;
        while (i > 0) {
            i -= 1;
            if (input[i] == ':') {
                colon_pos = i;
                break;
            }
        }
    }

    if (colon_pos) |cp| {
        const suffix = input[cp + 1 ..];

        // Try to parse as START-END or START-
        if (parseRangeSuffix(suffix)) |range| {
            return Region{
                .name = input[0..cp],
                .start = range.start,
                .end = range.end,
                .is_full = false,
            };
        }

        // If the suffix didn't parse as a valid range, the ':' is part of the name.
        // But there might be other colons further left — keep trying.
        var search_end = cp;
        while (search_end > 0) {
            var i: usize = search_end;
            var found_colon: ?usize = null;
            while (i > 0) {
                i -= 1;
                if (input[i] == ':') {
                    found_colon = i;
                    break;
                }
            }

            if (found_colon) |cp2| {
                const suffix2 = input[cp2 + 1 ..];
                if (parseRangeSuffix(suffix2)) |range| {
                    return Region{
                        .name = input[0..cp2],
                        .start = range.start,
                        .end = range.end,
                        .is_full = false,
                    };
                }
                search_end = cp2;
            } else {
                break;
            }
        }
    }

    // No valid region suffix found — treat entire input as sequence name
    return Region{
        .name = input,
        .start = 1,
        .end = null,
        .is_full = true,
    };
}

const RangeParsed = struct {
    start: u64,
    end: ?u64,
};

fn parseRangeSuffix(suffix: []const u8) ?RangeParsed {
    // Expect: START-END or START-
    // Find the '-' separator
    const dash_pos = std.mem.indexOfScalar(u8, suffix, '-') orelse return null;

    const start_str = suffix[0..dash_pos];
    const end_str = suffix[dash_pos + 1 ..];

    // START must be a valid positive integer
    if (start_str.len == 0) return null;
    const start = std.fmt.parseInt(u64, start_str, 10) catch return null;

    if (end_str.len == 0) {
        // NAME:START- form (to end of sequence)
        return RangeParsed{ .start = start, .end = null };
    }

    const end = std.fmt.parseInt(u64, end_str, 10) catch return null;
    return RangeParsed{ .start = start, .end = end };
}

// ============================================================================
// Sequence extraction
// ============================================================================

/// Run the get command: extract a region from a FASTA file using its index.
pub fn runGet(fasta_path: []const u8, region_str: []const u8) void {
    var idx = index_format.loadIndex(fasta_path);
    defer idx.deinit();

    const region = parseRegion(region_str);

    // Look up the sequence
    const rec_idx = idx.lookupName(region.name) orelse {
        printErrorAndExit("error: sequence not found: {s}\n", .{region.name});
    };
    const rec = idx.records[rec_idx];

    // Resolve coordinates
    var start = region.start;
    var end = region.end orelse rec.seq_len;
    const display_end = end; // Preserve original END for header (samtools keeps it)

    // For full sequence, use the full range
    if (region.is_full) {
        start = 1;
        end = rec.seq_len;
    }

    // Validate coordinates
    if (start < 1) {
        printErrorAndExit("error: start position must be >= 1\n", .{});
    }
    if (start > rec.seq_len) {
        printErrorAndExit("error: start position {d} exceeds sequence length {d}\n", .{ start, rec.seq_len });
    }
    if (end < start) {
        printErrorAndExit("error: end position must be >= start position\n", .{});
    }

    // Clamp end to seq_len silently (samtools behavior)
    if (end > rec.seq_len) {
        end = rec.seq_len;
    }

    const num_bases = end - start + 1;

    // Compute byte offset for start position using O(1) formula
    const base_index = start - 1;
    const line_number = base_index / rec.line_bases;
    const column = base_index % rec.line_bases;
    const start_byte = rec.seq_offset + (line_number * rec.line_bytes) + column;

    // Write output
    var buffered = std.io.bufferedWriter(std.io.getStdOut().writer());
    const writer = buffered.writer();

    // Header
    if (region.is_full) {
        writer.print(">{s}\n", .{region.name}) catch {
            printErrorAndExit("error: write failed\n", .{});
        };
    } else {
        writer.print(">{s}:{d}-{d}\n", .{ region.name, start, display_end }) catch {
            printErrorAndExit("error: write failed\n", .{});
        };
    }

    // Extract bases from mmap'd FASTA, skipping newlines
    const fasta = idx.fasta_data;
    var pos: usize = @intCast(start_byte);
    var bases_written: u64 = 0;
    var line_pos: u32 = 0;
    const wrap_width: u32 = 60; // samtools default output line width

    while (bases_written < num_bases and pos < fasta.len) {
        const byte = fasta[pos];
        pos += 1;

        // Skip newline characters (CRLF-safe)
        if (byte == '\n' or byte == '\r') continue;

        writer.writeByte(byte) catch {
            printErrorAndExit("error: write failed\n", .{});
        };
        bases_written += 1;
        line_pos += 1;

        // Wrap at 60 chars (samtools default)
        if (line_pos >= wrap_width) {
            writer.writeByte('\n') catch {
                printErrorAndExit("error: write failed\n", .{});
            };
            line_pos = 0;
        }
    }

    // Trailing newline if the last line wasn't complete
    if (line_pos > 0) {
        writer.writeByte('\n') catch {
            printErrorAndExit("error: write failed\n", .{});
        };
    }

    buffered.flush() catch {
        printErrorAndExit("error: write failed\n", .{});
    };
}
