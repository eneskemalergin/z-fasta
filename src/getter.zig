const std = @import("std");
const posix = std.posix;
const complement = @import("complement.zig");
const bed_parser = @import("bed_parser.zig");
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

pub const GetOptions = struct {
    region_strs: []const []const u8,
    bed_path: ?[]const u8 = null,
    names_path: ?[]const u8 = null,
    honor_strand: bool = false,
    summary: bool = false,
};

const ParsedRequest = struct {
    region: Region,
    reverse_complement: bool,
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
    seq_offset: u64,
    num_bases: u64, // number of bases to extract
    line_bases: u32,
    line_bytes: u32,
    reverse_complement: bool,
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
    return resolveParsedRequest(idx, .{ .region = region, .reverse_complement = false }, rec_idx, original_index);
}

fn resolveParsedRequest(idx: *const LoadedIndex, request: ParsedRequest, rec_idx: usize, original_index: usize) ResolvedRegion {
    const rec = idx.records[rec_idx];
    const region = request.region;

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
        .seq_offset = rec.seq_offset,
        .num_bases = num_bases,
        .line_bases = rec.line_bases,
        .line_bytes = rec.line_bytes,
        .reverse_complement = request.reverse_complement,
        .original_index = original_index,
    };
}

fn resolveParsedRequestsByRecordScan(idx: *const LoadedIndex, requests: []const ParsedRequest, resolved: []ResolvedRegion) void {
    var rec_indices = std.ArrayList(?usize).empty;
    defer rec_indices.deinit(std.heap.page_allocator);
    rec_indices.resize(std.heap.page_allocator, requests.len) catch {
        printErrorAndExit("error: out of memory\n", .{});
    };
    for (rec_indices.items) |*entry| entry.* = null;

    for (idx.records, 0..) |rec, rec_idx| {
        const rec_name = rec.getName(idx.fasta_data);
        for (requests, 0..) |request, request_idx| {
            if (std.mem.eql(u8, rec_name, request.region.name)) {
                rec_indices.items[request_idx] = rec_idx;
            }
        }
    }

    for (requests, 0..) |request, i| {
        const rec_idx = rec_indices.items[i] orelse {
            printErrorAndExit("error: sequence not found: {s}\n", .{request.region.name});
        };
        resolved[i] = resolveParsedRequest(idx, request, rec_idx, i);
    }
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
        writer.print(">{s}{s}\n", .{ resolved.name, if (resolved.reverse_complement) ":rc" else "" }) catch {
            printErrorAndExit("error: write failed\n", .{});
        };
    } else {
        writer.print(">{s}:{d}-{d}{s}\n", .{ resolved.name, resolved.start, resolved.display_end, if (resolved.reverse_complement) ":rc" else "" }) catch {
            printErrorAndExit("error: write failed\n", .{});
        };
    }

    if (resolved.reverse_complement) {
        emitRegionReverseComplement(resolved, fasta, writer);
        return;
    }

    emitRegionForward(resolved, fasta, writer);
}

fn emitRegionForward(resolved: ResolvedRegion, fasta: []const u8, writer: anytype) void {
    const wrap_width: usize = 60;
    var pos: usize = @intCast(resolved.start_byte);
    var bases_written: u64 = 0;
    var line_pos: usize = 0;
    var out_buf: [65536]u8 = undefined;
    var out_len: usize = 0;

    while (bases_written < resolved.num_bases and pos < fasta.len) {
        const byte = fasta[pos];
        pos += 1;

        if (byte == '\n' or byte == '\r') continue;

        if (out_len + 2 > out_buf.len) {
            writer.writeAll(out_buf[0..out_len]) catch {
                printErrorAndExit("error: write failed\n", .{});
            };
            out_len = 0;
        }

        out_buf[out_len] = byte;
        out_len += 1;
        bases_written += 1;
        line_pos += 1;

        if (line_pos >= wrap_width) {
            out_buf[out_len] = '\n';
            out_len += 1;
            line_pos = 0;
        }
    }

    if (line_pos > 0) {
        if (out_len == out_buf.len) {
            writer.writeAll(&out_buf) catch {
                printErrorAndExit("error: write failed\n", .{});
            };
            out_len = 0;
        }
        out_buf[out_len] = '\n';
        out_len += 1;
    }

    if (out_len > 0) {
        writer.writeAll(out_buf[0..out_len]) catch {
            printErrorAndExit("error: write failed\n", .{});
        };
    }
}

fn emitRegionReverseComplement(resolved: ResolvedRegion, fasta: []const u8, writer: anytype) void {
    const wrap_width: usize = 60;
    var bases_remaining = resolved.num_bases;
    var line_pos: usize = 0;
    var out_buf: [65536]u8 = undefined;
    var out_len: usize = 0;

    while (bases_remaining > 0) {
        const region_base_index = resolved.start - 1 + bases_remaining - 1;
        const line_number = region_base_index / resolved.line_bases;
        const column = region_base_index % resolved.line_bases;
        const pos: usize = @intCast(resolved.seq_offset + (line_number * resolved.line_bytes) + column);
        const byte = complement.complement(fasta[pos]);

        if (out_len + 2 > out_buf.len) {
            writer.writeAll(out_buf[0..out_len]) catch {
                printErrorAndExit("error: write failed\n", .{});
            };
            out_len = 0;
        }

        out_buf[out_len] = byte;
        out_len += 1;
        bases_remaining -= 1;
        line_pos += 1;

        if (line_pos >= wrap_width) {
            out_buf[out_len] = '\n';
            out_len += 1;
            line_pos = 0;
        }
    }

    if (line_pos > 0) {
        if (out_len == out_buf.len) {
            writer.writeAll(&out_buf) catch {
                printErrorAndExit("error: write failed\n", .{});
            };
            out_len = 0;
        }
        out_buf[out_len] = '\n';
        out_len += 1;
    }

    if (out_len > 0) {
        writer.writeAll(out_buf[0..out_len]) catch {
            printErrorAndExit("error: write failed\n", .{});
        };
    }
}

fn monotonicNs(io: std.Io) u64 {
    const now = std.Io.Clock.Timestamp.now(io, .awake);
    return @intCast(now.raw.toNanoseconds());
}

fn readAllInput(allocator: std.mem.Allocator, io: std.Io, path: []const u8) []u8 {
    if (std.mem.eql(u8, path, "-")) {
        var stdin_buf: [4096]u8 = undefined;
        var reader = std.Io.File.stdin().reader(io, &stdin_buf);
        return reader.interface.allocRemaining(allocator, .unlimited) catch {
            printErrorAndExit("error: failed to read stdin\n", .{});
        };
    }

    const file = std.Io.Dir.cwd().openFile(io, path, .{}) catch |err| switch (err) {
        error.FileNotFound => printErrorAndExit("error: file not found: {s}\n", .{path}),
        error.AccessDenied => printErrorAndExit("error: access denied: {s}\n", .{path}),
        else => printErrorAndExit("error: failed to open file: {s}\n", .{path}),
    };
    defer file.close(io);

    const stat = file.stat(io) catch {
        printErrorAndExit("error: failed to stat file: {s}\n", .{path});
    };

    const bytes = allocator.alloc(u8, stat.size) catch {
        printErrorAndExit("error: out of memory\n", .{});
    };

    var file_buf: [4096]u8 = undefined;
    var reader = file.reader(io, &file_buf);
    reader.interface.readSliceAll(bytes) catch {
        printErrorAndExit("error: failed to read file: {s}\n", .{path});
    };
    return bytes;
}

fn appendBedRequests(requests: *std.ArrayList(ParsedRequest), bed_data: []const u8, honor_strand: bool, allocator: std.mem.Allocator) void {
    var lines = std.mem.splitScalar(u8, bed_data, '\n');
    var line_number: usize = 0;

    while (lines.next()) |line| {
        line_number += 1;
        const parsed = bed_parser.parseBedLine(line, line_number) catch |err| switch (err) {
            error.MissingChrom => printErrorAndExit("error: invalid BED line {d}: missing chrom\n", .{line_number}),
            error.MissingStart => printErrorAndExit("error: invalid BED line {d}: missing start\n", .{line_number}),
            error.MissingEnd => printErrorAndExit("error: invalid BED line {d}: missing end\n", .{line_number}),
            error.InvalidStart => printErrorAndExit("error: invalid BED line {d}: invalid start\n", .{line_number}),
            error.InvalidEnd => printErrorAndExit("error: invalid BED line {d}: invalid end\n", .{line_number}),
            error.EmptyInterval => printErrorAndExit("error: invalid BED line {d}: end must be greater than start\n", .{line_number}),
        };

        switch (parsed) {
            .skip => continue,
            .region => |region| {
                if (honor_strand and region.strand == .invalid) {
                    printErrorAndExit("error: invalid BED line {d}: invalid strand\n", .{line_number});
                }

                requests.append(allocator, .{
                    .region = .{
                        .name = region.chrom,
                        .start = region.start1Based(),
                        .end = region.end1BasedInclusive(),
                        .is_full = false,
                    },
                    .reverse_complement = honor_strand and region.strand == .minus,
                }) catch {
                    printErrorAndExit("error: out of memory\n", .{});
                };
            },
        }
    }
}

fn appendNamesRequests(requests: *std.ArrayList(ParsedRequest), names_data: []const u8, allocator: std.mem.Allocator) void {
    var lines = std.mem.splitScalar(u8, names_data, '\n');
    while (lines.next()) |line| {
        const trimmed = if (line.len > 0 and line[line.len - 1] == '\r') line[0 .. line.len - 1] else line;
        if (trimmed.len == 0 or trimmed[0] == '#') continue;

        requests.append(allocator, .{
            .region = .{
                .name = trimmed,
                .start = 1,
                .end = null,
                .is_full = true,
            },
            .reverse_complement = false,
        }) catch {
            printErrorAndExit("error: out of memory\n", .{});
        };
    }
}

fn appendCliRequests(requests: *std.ArrayList(ParsedRequest), region_strs: []const []const u8, allocator: std.mem.Allocator) void {
    for (region_strs) |region_str| {
        requests.append(allocator, .{
            .region = parseRegion(region_str),
            .reverse_complement = false,
        }) catch {
            printErrorAndExit("error: out of memory\n", .{});
        };
    }
}

fn writeSummary(io: std.Io, region_count: usize, total_bases: u64, elapsed_ns: u64) void {
    var err_buf: [512]u8 = undefined;
    var stderr_fw = std.Io.File.Writer.initStreaming(.stderr(), io, &err_buf);
    const writer = &stderr_fw.interface;

    const seconds = if (elapsed_ns == 0) 0.0 else @as(f64, @floatFromInt(elapsed_ns)) / @as(f64, @floatFromInt(std.time.ns_per_s));
    const regions_per_second = if (seconds == 0.0) 0.0 else @as(f64, @floatFromInt(region_count)) / seconds;

    writer.print("summary: regions={d} total_bases={d} elapsed_s={d:.6} regions_per_s={d:.1}\n", .{ region_count, total_bases, seconds, regions_per_second }) catch {};
    stderr_fw.flush() catch {};
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
    runGetWithOptions(io, fasta_path, .{ .region_strs = region_strs });
}

pub fn runGetWithOptions(io: std.Io, fasta_path: []const u8, options: GetOptions) void {
    var input_arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer input_arena.deinit();
    const allocator = input_arena.allocator();

    var requests = std.ArrayList(ParsedRequest).empty;
    defer requests.deinit(allocator);

    if (options.bed_path) |bed_path| {
        const bed_data = readAllInput(allocator, io, bed_path);
        appendBedRequests(&requests, bed_data, options.honor_strand, allocator);
    }

    if (options.names_path) |names_path| {
        const names_data = readAllInput(allocator, io, names_path);
        appendNamesRequests(&requests, names_data, allocator);
    }

    appendCliRequests(&requests, options.region_strs, allocator);

    if (requests.items.len == 0) {
        printErrorAndExit("error: no regions provided\n", .{});
    }

    const start_ns = if (options.summary) monotonicNs(io) else 0;

    const load_mode: index_format.LoadMode = if (requests.items.len < 16) .records_only else .lookup_full_map;
    var idx = index_format.loadIndexWithMode(io, fasta_path, load_mode);
    defer idx.deinit();

    // Point-access pattern: disable kernel readahead.
    posix.madvise(@alignCast(@constCast(idx.fasta_data.ptr)), idx.fasta_data.len, posix.MADV.RANDOM) catch {};

    // Resolve all regions before writing any output — fail-fast on bad names or
    // coordinates so we never emit partial results.
    const resolved = allocator.alloc(ResolvedRegion, requests.items.len) catch {
        printErrorAndExit("error: out of memory\n", .{});
    };
    if (requests.items.len > 1 and requests.items.len < 16 and idx.source == .zfi and !idx.has_name_map) {
        resolveParsedRequestsByRecordScan(&idx, requests.items, resolved);
    } else {
        for (requests.items, 0..) |request, i| {
            const rec_idx = idx.lookupName(request.region.name) orelse {
                printErrorAndExit("error: sequence not found: {s}\n", .{request.region.name});
            };
            resolved[i] = resolveParsedRequest(&idx, request, rec_idx, i);
        }
    }

    var out_buf: [65536]u8 = undefined;
    var stdout_fw = std.Io.File.Writer.initStreaming(.stdout(), io, &out_buf);
    const writer = &stdout_fw.interface;

    var total_bases: u64 = 0;

    if (requests.items.len < 16) {
        // Direct path: small count, emit in CLI order.
        // The mmap page cache handles random access well for a handful of regions.
        for (resolved) |r| {
            total_bases += r.num_bases;
            emitRegion(r, idx.fasta_data, writer);
        }
    } else {
        // Sorted path: sort by file offset for sequential I/O, write each region into
        // its own buffer, then emit buffers in original CLI order.
        var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
        defer arena.deinit();
        const sort_allocator = arena.allocator();

        // Copy so we can sort without disturbing original order.
        const sorted = sort_allocator.dupe(ResolvedRegion, resolved) catch {
            printErrorAndExit("error: out of memory\n", .{});
        };
        std.mem.sort(ResolvedRegion, sorted, {}, struct {
            fn lessThan(_: void, a: ResolvedRegion, b: ResolvedRegion) bool {
                return a.start_byte < b.start_byte;
            }
        }.lessThan);

        // One output buffer per region keyed by original_index.
        const output_bufs = sort_allocator.alloc(std.Io.Writer.Allocating, requests.items.len) catch {
            printErrorAndExit("error: out of memory\n", .{});
        };
        for (output_bufs) |*buf| {
            buf.* = std.Io.Writer.Allocating.init(sort_allocator);
        }

        // Extract in file-offset order for sequential page access.
        for (sorted) |r| {
            total_bases += r.num_bases;
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

    if (options.summary) {
        writeSummary(io, resolved.len, total_bases, monotonicNs(io) - start_ns);
    }
}
