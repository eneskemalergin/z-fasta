const std = @import("std");
const posix = std.posix;
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
// Region resolution
// ============================================================================

/// A fully resolved, validated extraction request.
/// Byte offset and length are pre-computed so extraction is a single mmap
/// slice walk — no further index lookups required.
pub const ResolvedRegion = struct {
    name: []const u8, // sequence name for FASTA header
    start: u64, // 1-based inclusive (validated)
    display_end: u64, // end value for header (pre-clamp, samtools convention)
    is_full: bool, // true → emit ">NAME", false → emit ">NAME:start-display_end"
    start_byte: u64, // absolute byte offset into fasta_data for first base
    num_bases: u64, // number of bases to extract
    original_index: usize, // position in CLI argument list (preserves output order)
};

/// Resolve one region string against a loaded index.
/// Validates coordinates and pre-computes the O(1) byte offset.
/// Calls printErrorAndExit on any error (sequence not found, bad coordinates).
pub fn resolveRegion(idx: *const LoadedIndex, region_str: []const u8, original_index: usize) ResolvedRegion {
    const region = parseRegion(region_str);

    const rec_idx = idx.lookupName(region.name) orelse {
        printErrorAndExit("error: sequence not found: {s}\n", .{region.name});
    };

    return resolveParsedRegion(idx, region, rec_idx, original_index);
}

fn resolveParsedRegion(idx: *const LoadedIndex, region: Region, rec_idx: usize, original_index: usize) ResolvedRegion {
    const rec = idx.records[rec_idx];

    var start = region.start;
    var end = region.end orelse rec.seq_len;
    const display_end = end; // capture before clamping (samtools keeps unclamped in header)

    if (region.is_full) {
        start = 1;
        end = rec.seq_len;
    }

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

    // O(1) byte offset formula
    const base_index = start - 1;
    const line_number = base_index / rec.line_bases;
    const column = base_index % rec.line_bases;
    const start_byte = rec.seq_offset + (line_number * rec.line_bytes) + column;

    return ResolvedRegion{
        .name = region.name,
        .start = start,
        .display_end = display_end,
        .is_full = region.is_full,
        .start_byte = start_byte,
        .num_bases = num_bases,
        .original_index = original_index,
    };
}

fn resolveRegionsByRecordScan(idx: *const LoadedIndex, region_strs: []const []const u8, resolved: []ResolvedRegion) void {
    var regions_buf: [1024]Region = undefined;
    var rec_indices_buf: [1024]?usize = undefined;
    const regions = regions_buf[0..region_strs.len];
    const rec_indices = rec_indices_buf[0..region_strs.len];

    for (region_strs, 0..) |rs, i| {
        regions[i] = parseRegion(rs);
        rec_indices[i] = null;
    }

    for (idx.records, 0..) |rec, rec_idx| {
        const rec_name = rec.getName(idx.fasta_data);
        for (regions, 0..) |region, region_idx| {
            if (std.mem.eql(u8, rec_name, region.name)) {
                rec_indices[region_idx] = rec_idx;
            }
        }
    }

    for (regions, 0..) |region, i| {
        const rec_idx = rec_indices[i] orelse {
            printErrorAndExit("error: sequence not found: {s}\n", .{region.name});
        };
        resolved[i] = resolveParsedRegion(idx, region, rec_idx, i);
    }
}

// ============================================================================
// Sequence emission
// ============================================================================

/// Write FASTA output for one resolved region to `writer`.
/// Output is wrapped at 60 bases per line (samtools default).
fn emitRegion(resolved: ResolvedRegion, fasta: []const u8, writer: anytype) void {
    if (resolved.is_full) {
        writer.print(">{s}\n", .{resolved.name}) catch {
            printErrorAndExit("error: write failed\n", .{});
        };
    } else {
        writer.print(">{s}:{d}-{d}\n", .{ resolved.name, resolved.start, resolved.display_end }) catch {
            printErrorAndExit("error: write failed\n", .{});
        };
    }

    const wrap_width: usize = 60;
    var pos: usize = @intCast(resolved.start_byte);
    var bases_written: u64 = 0;
    var line_pos: usize = 0;
    var line_buf: [wrap_width]u8 = undefined;

    while (bases_written < resolved.num_bases and pos < fasta.len) {
        const byte = fasta[pos];
        pos += 1;

        // Skip newline characters (CRLF-safe)
        if (byte == '\n' or byte == '\r') continue;

        line_buf[line_pos] = byte;
        bases_written += 1;
        line_pos += 1;

        if (line_pos >= wrap_width) {
            writer.writeAll(&line_buf) catch {
                printErrorAndExit("error: write failed\n", .{});
            };
            writer.writeByte('\n') catch {
                printErrorAndExit("error: write failed\n", .{});
            };
            line_pos = 0;
        }
    }

    if (line_pos > 0) {
        writer.writeAll(line_buf[0..line_pos]) catch {
            printErrorAndExit("error: write failed\n", .{});
        };
        writer.writeByte('\n') catch {
            printErrorAndExit("error: write failed\n", .{});
        };
    }
}

// ============================================================================
// Public entry point
// ============================================================================

/// Run the get command: extract one or more regions from a FASTA file.
///
/// Regions are emitted in CLI argument order regardless of their position in
/// the file.
///
/// For >= 16 regions the extractions are sorted by file offset before reading
/// to improve sequential page access. Output is buffered per-region and
/// flushed in original CLI order.
pub fn runGet(io: std.Io, fasta_path: []const u8, region_strs: []const []const u8) void {
    const load_mode: index_format.LoadMode = if (region_strs.len < 16) .records_only else .lookup_full_map;
    var idx = index_format.loadIndexWithMode(io, fasta_path, load_mode);
    defer idx.deinit();

    // Point-access pattern: disable kernel readahead.
    posix.madvise(@alignCast(@constCast(idx.fasta_data.ptr)), idx.fasta_data.len, posix.MADV.RANDOM) catch {};

    // Resolve all regions before writing any output — fail-fast on bad names or
    // coordinates so we never emit partial results.
    var resolved_buf: [1024]ResolvedRegion = undefined;
    if (region_strs.len > resolved_buf.len) {
        printErrorAndExit("error: too many regions (max 1024)\n", .{});
    }
    const resolved = resolved_buf[0..region_strs.len];
    if (region_strs.len > 1 and region_strs.len < 16 and idx.source == .zfi and !idx.has_name_map) {
        resolveRegionsByRecordScan(&idx, region_strs, resolved);
    } else {
        for (region_strs, 0..) |rs, i| {
            resolved[i] = resolveRegion(&idx, rs, i);
        }
    }

    var out_buf: [65536]u8 = undefined;
    var stdout_fw = std.Io.File.Writer.initStreaming(.stdout(), io, &out_buf);
    const writer = &stdout_fw.interface;

    if (region_strs.len < 16) {
        // Direct path: small count, emit in CLI order.
        // The mmap page cache handles random access well for a handful of regions.
        for (resolved) |r| {
            emitRegion(r, idx.fasta_data, writer);
        }
    } else {
        // Sorted path: sort by file offset for sequential I/O, write each region into
        // its own buffer, then emit buffers in original CLI order.
        var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        // Copy so we can sort without disturbing original order.
        const sorted = allocator.dupe(ResolvedRegion, resolved) catch {
            printErrorAndExit("error: out of memory\n", .{});
        };
        std.mem.sort(ResolvedRegion, sorted, {}, struct {
            fn lessThan(_: void, a: ResolvedRegion, b: ResolvedRegion) bool {
                return a.start_byte < b.start_byte;
            }
        }.lessThan);

        // One output buffer per region keyed by original_index.
        const output_bufs = allocator.alloc(std.Io.Writer.Allocating, region_strs.len) catch {
            printErrorAndExit("error: out of memory\n", .{});
        };
        for (output_bufs) |*buf| {
            buf.* = std.Io.Writer.Allocating.init(allocator);
        }

        // Extract in file-offset order for sequential page access.
        for (sorted) |r| {
            emitRegion(r, idx.fasta_data, &output_bufs[r.original_index].writer);
        }

        // Emit in original CLI order.
        for (output_bufs) |*buf| {
            const output = buf.toArrayList();
            writer.writeAll(output.items) catch {
                printErrorAndExit("error: write failed\n", .{});
            };
        }
    }

    stdout_fw.flush() catch {
        printErrorAndExit("error: write failed\n", .{});
    };
}
