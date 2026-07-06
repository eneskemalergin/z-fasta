const std = @import("std");
const posix = std.posix;
const complement = @import("complement.zig");
const bed_parser = @import("bed_parser.zig");
const index_format = @import("index_format.zig");
const stats = @import("stats.zig");

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

pub const Orientation = struct {
    reverse: bool = false,
    complement: bool = false,

    pub fn compose(self: Orientation, other: Orientation) Orientation {
        return .{
            .reverse = self.reverse != other.reverse,
            .complement = self.complement != other.complement,
        };
    }

    pub fn isIdentity(self: Orientation) bool {
        return !self.reverse and !self.complement;
    }

    pub fn reverseComplement() Orientation {
        return .{ .reverse = true, .complement = true };
    }

    pub fn reverseOnly() Orientation {
        return .{ .reverse = true, .complement = false };
    }

    pub fn complementOnly() Orientation {
        return .{ .reverse = false, .complement = true };
    }
};

pub const GetOptions = struct {
    region_strs: []const []const u8,
    bed_path: ?[]const u8 = null,
    names_path: ?[]const u8 = null,
    honor_strand: bool = false,
    summary: bool = false,
    chunk_size: usize = 4_096,
    orientation: Orientation = .{},
    annotate_transform: bool = false,
};

pub const chunk_size_all = std.math.maxInt(usize);

/// Maximum BED/names input size for the all-in-memory path (`--chunk-size -1`).
pub const max_input_file_bytes: usize = 512 * 1024 * 1024;

/// Per-line buffer for chunked BED streaming (`takeDelimiter` limit).
pub const bed_line_reader_buffer_bytes = 4096;

/// Per-region output cap for the multi-region sort buffer path (≥16 regions).
pub const max_sort_path_region_output_bytes: u64 = 64 * 1024 * 1024;

/// Total intermediate output cap for the multi-region sort buffer path.
pub const max_sort_path_total_output_bytes: u64 = 256 * 1024 * 1024;

const ParsedRequest = struct {
    region: Region,
    orientation: Orientation,
};

const BatchStats = struct {
    region_count: usize = 0,
    total_bases: u64 = 0,
};

/// Parse a region string. Handles the Ensembl colon trap by parsing from the right.
/// Accepted formats:
///   NAME              full sequence
///   NAME:START-END    1-based, inclusive
///   NAME:START-       from START to end
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
        // But there might be other colons further left; keep trying.
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

    // No valid region suffix found; treat entire input as sequence name
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
/// slice walk; no further index lookups required.
pub const ResolvedRegion = struct {
    name: []const u8, // sequence name for FASTA header
    start: u64, // 1-based inclusive (validated)
    display_end: u64, // end value for header (pre-clamp, samtools convention)
    is_full: bool, // true -> emit ">NAME", false -> emit ">NAME:start-display_end"
    start_byte: u64, // absolute byte offset into fasta_data for first base
    seq_offset: u64,
    num_bases: u64, // number of bases to extract
    line_bases: u32,
    line_bytes: u32,
    side_table: []const index_format.SideTableLine,
    orientation: Orientation,
    annotate_transform: bool,
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
    return resolveParsedRequest(idx, .{ .region = region, .orientation = .{} }, rec_idx, original_index, false);
}

fn resolveParsedRequest(
    idx: *const LoadedIndex,
    request: ParsedRequest,
    rec_idx: usize,
    original_index: usize,
    annotate_transform: bool,
) ResolvedRegion {
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

    const side_table = idx.sideTableLines(rec);
    const start_byte = byteOffsetForBase(idx.fasta_data, rec, side_table, start - 1);

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
        .side_table = side_table,
        .orientation = request.orientation,
        .annotate_transform = annotate_transform,
        .original_index = original_index,
    };
}

fn byteOffsetForBase(
    fasta: []const u8,
    rec: IndexRecord,
    side_table: []const index_format.SideTableLine,
    base_index: u64,
) u64 {
    if (side_table.len == 0) {
        const line_number = base_index / rec.line_bases;
        const column = base_index % rec.line_bases;
        return rec.seq_offset + (line_number * rec.line_bytes) + column;
    }

    var lo: usize = 0;
    var hi: usize = side_table.len;
    while (lo < hi) {
        const mid = lo + (hi - lo) / 2;
        if (side_table[mid].base_start <= base_index) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    const line_idx = if (lo == 0) 0 else lo - 1;
    const line = side_table[line_idx];
    const column = base_index - line.base_start;

    var seen: u64 = 0;
    var pos: usize = @intCast(line.byte_offset);
    const end = @min(fasta.len, pos + @as(usize, @intCast(line.line_bytes)));
    while (pos < end) : (pos += 1) {
        if (fasta[pos] <= ' ') continue;
        if (seen == column) return @intCast(pos);
        seen += 1;
    }

    printErrorAndExit("error: corrupt non-uniform index side table\n", .{});
}

/// Resolve CLI batches (2..15 regions) loaded with `.records_only` (no name hash map).
/// Scans every index record against every request: O(records x regions). At most 14
/// regions here, so a hash map would add setup and cache pressure without a measurable
/// win. Batches of 16+ load `.lookup_full_map` and use `lookupName` instead. Revisit if
/// the sub-16 threshold rises materially.
fn resolveParsedRequestsByRecordScan(
    idx: *const LoadedIndex,
    requests: []const ParsedRequest,
    resolved: []ResolvedRegion,
    annotate_transform: bool,
) void {
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
        resolved[i] = resolveParsedRequest(idx, request, rec_idx, i, annotate_transform);
    }
}

fn findRecordIndex(idx: *const LoadedIndex, name: []const u8) ?usize {
    if (idx.lookupName(name)) |rec_idx| return rec_idx;

    for (idx.records, 0..) |rec, rec_idx| {
        if (std.mem.eql(u8, rec.getName(idx.fasta_data), name)) return rec_idx;
    }

    return null;
}

/// Classify a record as nucleotide or protein from a prefix of its sequence.
/// `sample_limit` counts sequence bases (non-newline bytes), not raw file bytes.
fn detectRecordType(rec: IndexRecord, fasta: []const u8) stats.SequenceType {
    var counts = [_]u64{0} ** 256;
    var total: u64 = 0;
    var pos: usize = @intCast(rec.seq_offset);
    const sample_limit: u64 = @min(rec.seq_len, 100_000);

    while (pos < fasta.len and total < sample_limit) : (pos += 1) {
        const byte = fasta[pos];
        if (byte <= ' ') continue;
        counts[byte] += 1;
        total += 1;
    }

    return stats.detectType(&counts, total);
}

fn ensureComplementAllowed(idx: *const LoadedIndex, requests: []const ParsedRequest) void {
    var last_name: ?[]const u8 = null;
    var last_rec_idx: usize = 0;
    var last_checked_rec_idx: ?usize = null;
    var last_checked_type: stats.SequenceType = undefined;

    for (requests) |request| {
        if (!request.orientation.complement) continue;

        const rec_idx = if (last_name) |name|
            if (std.mem.eql(u8, name, request.region.name))
                last_rec_idx
            else
                findRecordIndex(idx, request.region.name) orelse {
                    printErrorAndExit("error: sequence not found: {s}\n", .{request.region.name});
                }
        else
            findRecordIndex(idx, request.region.name) orelse {
                printErrorAndExit("error: sequence not found: {s}\n", .{request.region.name});
            };

        last_name = request.region.name;
        last_rec_idx = rec_idx;

        const rec_type = if (last_checked_rec_idx) |cached_rec_idx|
            if (cached_rec_idx == rec_idx)
                last_checked_type
            else blk: {
                const detected = detectRecordType(idx.records[rec_idx], idx.fasta_data);
                last_checked_rec_idx = rec_idx;
                last_checked_type = detected;
                break :blk detected;
            }
        else blk: {
            const detected = detectRecordType(idx.records[rec_idx], idx.fasta_data);
            last_checked_rec_idx = rec_idx;
            last_checked_type = detected;
            break :blk detected;
        };

        if (rec_type == .protein) {
            printErrorAndExit("error: reverse complement is not defined for protein sequences: {s}\n", .{request.region.name});
        }
    }
}

// ============================================================================
// Sequence emission
// ============================================================================

/// Writes one region's FASTA output.
pub fn extractRegion(idx: *const LoadedIndex, region_str: []const u8, writer: anytype) void {
    const resolved = resolveRegion(idx, region_str, 0);
    emitRegion(resolved, idx.fasta_data, writer);
}

/// Write FASTA output for one resolved region to `writer`.
/// Output is wrapped at 60 bases per line (samtools default).
fn emitRegion(resolved: ResolvedRegion, fasta: []const u8, writer: anytype) void {
    const annotation = headerAnnotation(resolved.orientation, resolved.annotate_transform);
    if (resolved.is_full) {
        writer.print(">{s}{s}\n", .{ resolved.name, annotation }) catch {
            printErrorAndExit("error: write failed\n", .{});
        };
    } else {
        writer.print(">{s}:{d}-{d}{s}\n", .{ resolved.name, resolved.start, resolved.display_end, annotation }) catch {
            printErrorAndExit("error: write failed\n", .{});
        };
    }

    if (resolved.orientation.reverse) {
        emitRegionBackward(resolved, fasta, writer);
        return;
    }

    emitRegionForward(resolved, fasta, writer);
}

fn headerAnnotation(orientation: Orientation, annotate_transform: bool) []const u8 {
    if (!annotate_transform or orientation.isIdentity()) return "";
    if (orientation.reverse and orientation.complement) return " (reverse complement)";
    if (orientation.reverse) return " (reverse)";
    if (orientation.complement) return " (complement)";
    return "";
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

        if (byte <= ' ') continue;

        if (out_len + 2 > out_buf.len) {
            writer.writeAll(out_buf[0..out_len]) catch {
                printErrorAndExit("error: write failed\n", .{});
            };
            out_len = 0;
        }

        out_buf[out_len] = if (resolved.orientation.complement) complement.complement(byte) else byte;
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

fn emitRegionBackward(resolved: ResolvedRegion, fasta: []const u8, writer: anytype) void {
    const wrap_width: usize = 60;
    var bases_remaining = resolved.num_bases;
    var line_pos: usize = 0;
    var out_buf: [65536]u8 = undefined;
    var out_len: usize = 0;

    const last_base_index = resolved.start - 1 + resolved.num_bases - 1;
    const rec = IndexRecord{
        .name_offset = 1,
        .name_len = 1,
        .seq_offset = resolved.seq_offset,
        .seq_len = resolved.start - 1 + resolved.num_bases,
        .line_bases = resolved.line_bases,
        .line_bytes = resolved.line_bytes,
    };
    var pos: usize = @intCast(byteOffsetForBase(fasta, rec, resolved.side_table, last_base_index));

    while (bases_remaining > 0) {
        const byte = if (resolved.orientation.complement) complement.complement(fasta[pos]) else fasta[pos];

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

        if (bases_remaining > 0) {
            while (pos > 0) {
                pos -= 1;
                const prev = fasta[pos];
                if (prev > ' ') break;
            }
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

fn estimateRegionOutputBytes(resolved: ResolvedRegion) u64 {
    const wrap_lines = resolved.num_bases / 60 + @as(u64, @intFromBool(resolved.num_bases % 60 != 0));
    const header_len: u64 = if (resolved.is_full)
        @intCast(resolved.name.len + 1)
    else
        @intCast(resolved.name.len + 32);
    const annotation_len: u64 = if (resolved.annotate_transform) 24 else 0;
    return resolved.num_bases + wrap_lines + header_len + annotation_len;
}

fn shouldUseSortBuffers(resolved: []const ResolvedRegion) bool {
    var total: u64 = 0;
    for (resolved) |r| {
        const est = estimateRegionOutputBytes(r);
        if (est > max_sort_path_region_output_bytes) return false;
        total += est;
        if (total > max_sort_path_total_output_bytes) return false;
    }
    return true;
}

fn batchHasReverseReads(resolved: []const ResolvedRegion) bool {
    for (resolved) |r| {
        if (r.orientation.reverse) return true;
    }
    return false;
}

fn shouldAdviseSequentialMmap(
    requests_len: usize,
    already_in_offset_order: bool,
    use_sort_buffers: bool,
    has_reverse_reads: bool,
    allow_sequential_madvise: bool,
) bool {
    if (!allow_sequential_madvise) return false;
    if (requests_len < 16) return false;
    if (has_reverse_reads) return false;
    if (already_in_offset_order) return true;
    return use_sort_buffers;
}

fn adviseFastaMmap(fasta: []const u8, advice: u32) void {
    posix.madvise(@alignCast(@constCast(fasta.ptr)), fasta.len, advice) catch {};
}

fn readAllInput(allocator: std.mem.Allocator, io: std.Io, path: []const u8) []u8 {
    if (std.mem.eql(u8, path, "-")) {
        var stdin_buf: [4096]u8 = undefined;
        var reader = std.Io.File.stdin().reader(io, &stdin_buf);
        return reader.interface.allocRemaining(allocator, .limited(max_input_file_bytes)) catch |err| switch (err) {
            error.StreamTooLong => printErrorAndExit(
                "error: stdin exceeds {d} byte limit; use --chunk-size 4096 for large BED input\n",
                .{max_input_file_bytes},
            ),
            else => printErrorAndExit("error: failed to read stdin\n", .{}),
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

    if (stat.size > max_input_file_bytes) {
        printErrorAndExit(
            "error: input file exceeds {d} byte limit: {s}; use default --chunk-size for large BED files\n",
            .{ max_input_file_bytes, path },
        );
    }

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

fn appendBedRegionRequest(
    requests: *std.ArrayList(ParsedRequest),
    region: bed_parser.BedRegion,
    honor_strand: bool,
    global_orientation: Orientation,
    allocator: std.mem.Allocator,
    name_allocator: std.mem.Allocator,
    duplicate_name: bool,
    last_duplicated_name: ?*?[]const u8,
) void {
    if (honor_strand and region.strand == .invalid) {
        printErrorAndExit("error: invalid BED line {d}: invalid strand\n", .{region.line_number});
    }

    const name = if (duplicate_name) blk: {
        if (last_duplicated_name) |cached_name| {
            if (cached_name.*) |existing| {
                if (std.mem.eql(u8, existing, region.chrom)) break :blk existing;
            }
        }

        const duplicated = name_allocator.dupe(u8, region.chrom) catch {
            printErrorAndExit("error: out of memory\n", .{});
        };
        if (last_duplicated_name) |cached_name| cached_name.* = duplicated;
        break :blk duplicated;
    } else region.chrom;

    requests.append(allocator, .{
        .region = .{
            .name = name,
            .start = region.start1Based(),
            .end = region.end1BasedInclusive(),
            .is_full = false,
        },
        .orientation = (if (honor_strand and region.strand == .minus) Orientation.reverseComplement() else Orientation{}).compose(global_orientation),
    }) catch {
        printErrorAndExit("error: out of memory\n", .{});
    };
}

fn appendBedLineRequest(
    requests: *std.ArrayList(ParsedRequest),
    line: []const u8,
    line_number: usize,
    honor_strand: bool,
    global_orientation: Orientation,
    allocator: std.mem.Allocator,
    name_allocator: std.mem.Allocator,
    duplicate_name: bool,
    last_duplicated_name: ?*?[]const u8,
) void {
    const parsed = bed_parser.parseBedLine(line, line_number) catch |err| switch (err) {
        error.MissingChrom => printErrorAndExit("error: invalid BED line {d}: missing chrom\n", .{line_number}),
        error.MissingStart => printErrorAndExit("error: invalid BED line {d}: missing start\n", .{line_number}),
        error.MissingEnd => printErrorAndExit("error: invalid BED line {d}: missing end\n", .{line_number}),
        error.InvalidStart => printErrorAndExit("error: invalid BED line {d}: invalid start\n", .{line_number}),
        error.InvalidEnd => printErrorAndExit("error: invalid BED line {d}: invalid end\n", .{line_number}),
        error.EmptyInterval => printErrorAndExit("error: invalid BED line {d}: end must be greater than start\n", .{line_number}),
    };

    switch (parsed) {
        .skip => {},
        .region => |region| appendBedRegionRequest(requests, region, honor_strand, global_orientation, allocator, name_allocator, duplicate_name, last_duplicated_name),
    }
}

fn appendBedRequests(requests: *std.ArrayList(ParsedRequest), bed_data: []const u8, honor_strand: bool, global_orientation: Orientation, allocator: std.mem.Allocator) void {
    var lines = std.mem.splitScalar(u8, bed_data, '\n');
    var line_number: usize = 0;

    while (lines.next()) |line| {
        line_number += 1;
        appendBedLineRequest(requests, line, line_number, honor_strand, global_orientation, allocator, allocator, false, null);
    }
}

fn processBedData(
    idx: *const LoadedIndex,
    allocator: std.mem.Allocator,
    bed_data: []const u8,
    honor_strand: bool,
    global_orientation: Orientation,
    annotate_transform: bool,
    writer: anytype,
) BatchStats {
    var requests = std.ArrayList(ParsedRequest).empty;
    defer requests.deinit(allocator);

    appendBedRequests(&requests, bed_data, honor_strand, global_orientation, allocator);
    return processParsedRequests(idx, allocator, requests.items, annotate_transform, writer, false);
}

fn appendNamesRequests(requests: *std.ArrayList(ParsedRequest), names_data: []const u8, orientation: Orientation, allocator: std.mem.Allocator) void {
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
            .orientation = orientation,
        }) catch {
            printErrorAndExit("error: out of memory\n", .{});
        };
    }
}

fn appendCliRequests(requests: *std.ArrayList(ParsedRequest), region_strs: []const []const u8, orientation: Orientation, allocator: std.mem.Allocator) void {
    for (region_strs) |region_str| {
        requests.append(allocator, .{
            .region = parseRegion(region_str),
            .orientation = orientation,
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

fn processParsedRequests(
    idx: *const LoadedIndex,
    allocator: std.mem.Allocator,
    requests: []const ParsedRequest,
    annotate_transform: bool,
    writer: anytype,
    allow_sequential_madvise: bool,
) BatchStats {
    if (requests.len == 0) return .{};

    ensureComplementAllowed(idx, requests);

    const resolved = allocator.alloc(ResolvedRegion, requests.len) catch {
        printErrorAndExit("error: out of memory\n", .{});
    };
    var already_in_offset_order = requests.len >= 16;
    var prev_start_byte: u64 = 0;

    // `.records_only` + 2..15 regions: by-record scan (see resolveParsedRequestsByRecordScan).
    if (requests.len > 1 and requests.len < 16 and idx.source == .zfi and !idx.has_name_map) {
        resolveParsedRequestsByRecordScan(idx, requests, resolved, annotate_transform);
    } else {
        var last_name: ?[]const u8 = null;
        var last_rec_idx: usize = 0;
        for (requests, 0..) |request, i| {
            const rec_idx = if (last_name) |name|
                if (std.mem.eql(u8, name, request.region.name))
                    last_rec_idx
                else
                    idx.lookupName(request.region.name) orelse {
                        printErrorAndExit("error: sequence not found: {s}\n", .{request.region.name});
                    }
            else
                idx.lookupName(request.region.name) orelse {
                    printErrorAndExit("error: sequence not found: {s}\n", .{request.region.name});
                };

            last_name = request.region.name;
            last_rec_idx = rec_idx;

            resolved[i] = resolveParsedRequest(idx, request, rec_idx, i, annotate_transform);
            if (already_in_offset_order) {
                if (i > 0 and resolved[i].start_byte < prev_start_byte) {
                    already_in_offset_order = false;
                }
                prev_start_byte = resolved[i].start_byte;
            }
        }
    }

    const use_sort_buffers = requests.len >= 16 and !already_in_offset_order and shouldUseSortBuffers(resolved);
    const sequential_mmap = shouldAdviseSequentialMmap(
        requests.len,
        already_in_offset_order,
        use_sort_buffers,
        batchHasReverseReads(resolved),
        allow_sequential_madvise,
    );
    const mmap_advice: u32 = if (sequential_mmap) posix.MADV.SEQUENTIAL else posix.MADV.RANDOM;
    adviseFastaMmap(idx.fasta_data, mmap_advice);

    var total_bases: u64 = 0;

    if (requests.len < 16) {
        for (resolved) |r| {
            total_bases += r.num_bases;
            emitRegion(r, idx.fasta_data, writer);
        }
    } else {
        if (already_in_offset_order) {
            for (resolved) |r| {
                total_bases += r.num_bases;
                emitRegion(r, idx.fasta_data, writer);
            }
            return .{
                .region_count = requests.len,
                .total_bases = total_bases,
            };
        }

        if (!use_sort_buffers) {
            for (resolved) |r| {
                total_bases += r.num_bases;
                emitRegion(r, idx.fasta_data, writer);
            }
            return .{
                .region_count = requests.len,
                .total_bases = total_bases,
            };
        }

        var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
        defer arena.deinit();
        const sort_allocator = arena.allocator();

        const sorted = sort_allocator.dupe(ResolvedRegion, resolved) catch {
            printErrorAndExit("error: out of memory\n", .{});
        };
        std.mem.sort(ResolvedRegion, sorted, {}, struct {
            fn lessThan(_: void, a: ResolvedRegion, b: ResolvedRegion) bool {
                return a.start_byte < b.start_byte;
            }
        }.lessThan);

        const output_bufs = sort_allocator.alloc(std.Io.Writer.Allocating, requests.len) catch {
            printErrorAndExit("error: out of memory\n", .{});
        };
        for (output_bufs) |*buf| {
            buf.* = std.Io.Writer.Allocating.init(sort_allocator);
        }

        for (sorted) |r| {
            total_bases += r.num_bases;
            emitRegion(r, idx.fasta_data, &output_bufs[r.original_index].writer);
        }

        for (output_bufs) |*buf| {
            const output = buf.toArrayList();
            writer.writeAll(output.items) catch {
                printErrorAndExit("error: write failed\n", .{});
            };
        }
    }

    return .{
        .region_count = requests.len,
        .total_bases = total_bases,
    };
}

fn processBedReaderChunked(
    idx: *const LoadedIndex,
    reader: *std.Io.Reader,
    honor_strand: bool,
    global_orientation: Orientation,
    chunk_size: usize,
    annotate_transform: bool,
    writer: anytype,
) BatchStats {
    var total = BatchStats{};
    var line_number: usize = 0;
    var reached_end = false;
    var name_arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer name_arena.deinit();
    var last_duplicated_name: ?[]const u8 = null;

    while (true) {
        var chunk_arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
        const chunk_allocator = chunk_arena.allocator();

        var requests = std.ArrayList(ParsedRequest).empty;

        while (requests.items.len < chunk_size) {
            const maybe_line = reader.takeDelimiter('\n') catch |err| switch (err) {
                error.ReadFailed => printErrorAndExit("error: failed to read BED input\n", .{}),
                error.StreamTooLong => printErrorAndExit(
                    "error: BED line {d} exceeds {d}-byte reader buffer (no newline within limit)\n",
                    .{ line_number + 1, bed_line_reader_buffer_bytes },
                ),
            };

            const line = maybe_line orelse {
                reached_end = true;
                break;
            };

            line_number += 1;
            appendBedLineRequest(&requests, line, line_number, honor_strand, global_orientation, chunk_allocator, name_arena.allocator(), true, &last_duplicated_name);
        }

        if (requests.items.len == 0) {
            chunk_arena.deinit();
            break;
        }

        const batch = processParsedRequests(idx, chunk_allocator, requests.items, annotate_transform, writer, false);
        total.region_count += batch.region_count;
        total.total_bases += batch.total_bases;

        chunk_arena.deinit();

        if (reached_end) break;
    }

    return total;
}

fn processBedPathChunked(
    io: std.Io,
    idx: *const LoadedIndex,
    path: []const u8,
    honor_strand: bool,
    global_orientation: Orientation,
    chunk_size: usize,
    annotate_transform: bool,
    writer: anytype,
) BatchStats {
    if (std.mem.eql(u8, path, "-")) {
        var stdin_buf: [bed_line_reader_buffer_bytes]u8 = undefined;
        var stdin_reader = std.Io.File.stdin().reader(io, &stdin_buf);
        return processBedReaderChunked(idx, &stdin_reader.interface, honor_strand, global_orientation, chunk_size, annotate_transform, writer);
    }

    const file = std.Io.Dir.cwd().openFile(io, path, .{}) catch |err| switch (err) {
        error.FileNotFound => printErrorAndExit("error: file not found: {s}\n", .{path}),
        error.AccessDenied => printErrorAndExit("error: access denied: {s}\n", .{path}),
        else => printErrorAndExit("error: failed to open file: {s}\n", .{path}),
    };
    defer file.close(io);

    var file_buf: [bed_line_reader_buffer_bytes]u8 = undefined;
    var file_reader = file.reader(io, &file_buf);
    return processBedReaderChunked(idx, &file_reader.interface, honor_strand, global_orientation, chunk_size, annotate_transform, writer);
}

fn processBedPathAllInMemory(
    io: std.Io,
    idx: *const LoadedIndex,
    allocator: std.mem.Allocator,
    path: []const u8,
    honor_strand: bool,
    global_orientation: Orientation,
    annotate_transform: bool,
    writer: anytype,
) BatchStats {
    const bed_data = readAllInput(allocator, io, path);
    return processBedData(idx, allocator, bed_data, honor_strand, global_orientation, annotate_transform, writer);
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
    if (options.chunk_size == 0) {
        printErrorAndExit("error: --chunk-size must be >= 1\n", .{});
    }

    var input_arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer input_arena.deinit();
    const allocator = input_arena.allocator();

    var requests = std.ArrayList(ParsedRequest).empty;
    defer requests.deinit(allocator);

    if (options.names_path) |names_path| {
        const names_data = readAllInput(allocator, io, names_path);
        appendNamesRequests(&requests, names_data, options.orientation, allocator);
    }

    appendCliRequests(&requests, options.region_strs, options.orientation, allocator);

    const start_ns = if (options.summary) monotonicNs(io) else 0;

    const load_mode: index_format.LoadMode = if (options.bed_path != null or options.names_path != null or requests.items.len >= 16) .lookup_full_map else .records_only;
    var idx = index_format.loadIndexWithMode(io, fasta_path, load_mode);
    defer idx.deinit();

    var out_buf: [65536]u8 = undefined;
    var stdout_fw = std.Io.File.Writer.initStreaming(.stdout(), io, &out_buf);
    const writer = &stdout_fw.interface;

    var totals = BatchStats{};

    if (options.bed_path) |bed_path| {
        const batch = if (options.chunk_size == chunk_size_all)
            processBedPathAllInMemory(io, &idx, allocator, bed_path, options.honor_strand, options.orientation, options.annotate_transform, writer)
        else
            processBedPathChunked(io, &idx, bed_path, options.honor_strand, options.orientation, options.chunk_size, options.annotate_transform, writer);
        totals.region_count += batch.region_count;
        totals.total_bases += batch.total_bases;
    }

    if (requests.items.len > 0) {
        const batch = processParsedRequests(&idx, allocator, requests.items, options.annotate_transform, writer, true);
        totals.region_count += batch.region_count;
        totals.total_bases += batch.total_bases;
    }

    if (totals.region_count == 0) {
        printErrorAndExit("error: no regions provided\n", .{});
    }

    stdout_fw.flush() catch {
        printErrorAndExit("error: write failed\n", .{});
    };

    if (options.summary) {
        writeSummary(io, totals.region_count, totals.total_bases, monotonicNs(io) - start_ns);
    }
}

test "processBedReaderChunked matches non-chunked extraction" {
    const test_io = std.Io.Threaded.global_single_threaded.io();
    var idx = index_format.loadIndex(test_io, "tests/data/simple.fasta");
    defer idx.deinit();

    const bed_data =
        "# comment\n" ++
        "seq2\t0\t4\n" ++
        "seq1\t0\t4\n" ++
        "seq1\t4\t8\n";

    var chunk_reader = std.Io.Reader.fixed(bed_data);
    var chunk_writer = std.Io.Writer.Allocating.init(std.testing.allocator);
    defer chunk_writer.deinit();

    const chunked = processBedReaderChunked(&idx, &chunk_reader, false, .{}, 2, false, &chunk_writer.writer);

    var batch_writer = std.Io.Writer.Allocating.init(std.testing.allocator);
    defer batch_writer.deinit();

    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const batch = processBedData(&idx, arena.allocator(), bed_data, false, .{}, false, &batch_writer.writer);

    try std.testing.expectEqual(batch.region_count, chunked.region_count);
    try std.testing.expectEqual(batch.total_bases, chunked.total_bases);
    try std.testing.expectEqualStrings(batch_writer.written(), chunk_writer.written());
}

test "processBedReaderChunked preserves strand handling across chunk boundaries" {
    const test_io = std.Io.Threaded.global_single_threaded.io();
    var idx = index_format.loadIndex(test_io, "tests/data/simple.fasta");
    defer idx.deinit();

    const bed_data =
        "seq1\t0\t5\tname\t0\t-\n" ++
        "seq2\t0\t4\tname\t0\t+\n";

    var chunk_reader = std.Io.Reader.fixed(bed_data);
    var chunk_writer = std.Io.Writer.Allocating.init(std.testing.allocator);
    defer chunk_writer.deinit();

    const chunked = processBedReaderChunked(&idx, &chunk_reader, true, .{}, 1, false, &chunk_writer.writer);

    try std.testing.expectEqual(@as(usize, 2), chunked.region_count);
    try std.testing.expectEqual(@as(u64, 9), chunked.total_bases);
    try std.testing.expectEqualStrings(
        ">seq1:1-5\nTACGT\n>seq2:1-4\nGGGG\n",
        chunk_writer.written(),
    );
}

test "processBedReaderChunked preserves duplicate chrom cache across chunk boundaries" {
    const test_io = std.Io.Threaded.global_single_threaded.io();
    var idx = index_format.loadIndex(test_io, "tests/data/simple.fasta");
    defer idx.deinit();

    const bed_data =
        "seq1\t0\t4\n" ++
        "seq1\t4\t8\n" ++
        "seq2\t0\t4\n";

    var chunk_reader = std.Io.Reader.fixed(bed_data);
    var chunk_writer = std.Io.Writer.Allocating.init(std.testing.allocator);
    defer chunk_writer.deinit();

    const chunked = processBedReaderChunked(&idx, &chunk_reader, false, .{}, 1, false, &chunk_writer.writer);

    var batch_writer = std.Io.Writer.Allocating.init(std.testing.allocator);
    defer batch_writer.deinit();

    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const batch = processBedData(&idx, arena.allocator(), bed_data, false, .{}, false, &batch_writer.writer);

    try std.testing.expectEqual(batch.region_count, chunked.region_count);
    try std.testing.expectEqual(batch.total_bases, chunked.total_bases);
    try std.testing.expectEqualStrings(batch_writer.written(), chunk_writer.written());
}

test "orientation compose behaves like transform composition" {
    const minus_strand = Orientation.reverseComplement();
    const global_rc = Orientation.reverseComplement();
    const global_reverse = Orientation.reverseOnly();
    const global_complement = Orientation.complementOnly();

    try std.testing.expect(minus_strand.compose(global_rc).isIdentity());
    try std.testing.expectEqual(true, minus_strand.compose(global_reverse).complement);
    try std.testing.expectEqual(false, minus_strand.compose(global_reverse).reverse);
    try std.testing.expectEqual(true, minus_strand.compose(global_complement).reverse);
    try std.testing.expectEqual(false, minus_strand.compose(global_complement).complement);
}

test "processParsedRequests applies complement-only and reverse-only transforms" {
    const test_io = std.Io.Threaded.global_single_threaded.io();
    var idx = index_format.loadIndex(test_io, "tests/data/simple.fasta");
    defer idx.deinit();

    const requests = [_]ParsedRequest{
        .{ .region = parseRegion("seq1:1-5"), .orientation = Orientation.complementOnly() },
        .{ .region = parseRegion("seq1:1-5"), .orientation = Orientation.reverseOnly() },
    };

    var writer = std.Io.Writer.Allocating.init(std.testing.allocator);
    defer writer.deinit();

    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();

    const batch = processParsedRequests(&idx, arena.allocator(), &requests, true, &writer.writer, true);

    try std.testing.expectEqual(@as(usize, 2), batch.region_count);
    try std.testing.expectEqualStrings(
        ">seq1:1-5 (complement)\nTGCAT\n>seq1:1-5 (reverse)\nATGCA\n",
        writer.written(),
    );
}

test "processParsedRequests annotates transforms on by-record-scan path" {
    const test_io = std.Io.Threaded.global_single_threaded.io();
    var idx = index_format.loadIndexWithMode(test_io, "tests/data/simple.fasta", .records_only);
    defer idx.deinit();
    try std.testing.expect(!idx.has_name_map);
    try std.testing.expectEqual(index_format.LoadedIndex.IndexSource.zfi, idx.source);

    const requests = [_]ParsedRequest{
        .{ .region = parseRegion("seq1:1-5"), .orientation = Orientation.complementOnly() },
        .{ .region = parseRegion("seq2:1-4"), .orientation = Orientation.reverseOnly() },
    };

    var writer = std.Io.Writer.Allocating.init(std.testing.allocator);
    defer writer.deinit();

    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();

    const batch = processParsedRequests(&idx, arena.allocator(), &requests, true, &writer.writer, true);

    try std.testing.expectEqual(@as(usize, 2), batch.region_count);
    try std.testing.expectEqualStrings(
        ">seq1:1-5 (complement)\nTGCAT\n>seq2:1-4 (reverse)\nGGGG\n",
        writer.written(),
    );
}

test "detectRecordType classifies nucleotide and protein records" {
    const test_io = std.Io.Threaded.global_single_threaded.io();

    var nucleotide_idx = index_format.loadIndex(test_io, "tests/data/simple.fasta");
    defer nucleotide_idx.deinit();
    try std.testing.expectEqual(
        stats.SequenceType.nucleotide,
        detectRecordType(nucleotide_idx.records[0], nucleotide_idx.fasta_data),
    );

    var protein_idx = index_format.loadIndex(test_io, "tests/data/proteome.fasta");
    defer protein_idx.deinit();
    try std.testing.expectEqual(
        stats.SequenceType.protein,
        detectRecordType(protein_idx.records[0], protein_idx.fasta_data),
    );
}

test "shouldUseSortBuffers accepts typical 20-region batch" {
    var regions: [20]ResolvedRegion = undefined;
    for (&regions, 0..) |*r, i| {
        r.* = .{
            .name = "seq1",
            .start = @intCast(i + 1),
            .display_end = @intCast(i + 1),
            .is_full = false,
            .start_byte = 100 + i,
            .seq_offset = 6,
            .num_bases = 1,
            .line_bases = 24,
            .line_bytes = 25,
            .side_table = &.{},
            .orientation = .{},
            .annotate_transform = false,
            .original_index = i,
        };
    }
    try std.testing.expect(shouldUseSortBuffers(&regions));
}

test "shouldUseSortBuffers rejects oversized single region" {
    const region = ResolvedRegion{
        .name = "chr1",
        .start = 1,
        .display_end = max_sort_path_region_output_bytes,
        .is_full = false,
        .start_byte = 0,
        .seq_offset = 0,
        .num_bases = max_sort_path_region_output_bytes,
        .line_bases = 60,
        .line_bytes = 61,
        .side_table = &.{},
        .orientation = .{},
        .annotate_transform = false,
        .original_index = 0,
    };
    try std.testing.expect(!shouldUseSortBuffers(&.{region}));
}

test "shouldAdviseSequentialMmap for CLI multi-region sort path" {
    try std.testing.expect(shouldAdviseSequentialMmap(100, false, true, false, true));
    try std.testing.expect(shouldAdviseSequentialMmap(100, true, false, false, true));
    try std.testing.expect(!shouldAdviseSequentialMmap(100, false, true, true, true));
    try std.testing.expect(!shouldAdviseSequentialMmap(100, false, true, false, false));
    try std.testing.expect(!shouldAdviseSequentialMmap(10, false, false, false, true));
    try std.testing.expect(!shouldAdviseSequentialMmap(100, false, false, false, true));
}

test "emitRegion handles non-uniform side-table forward extraction" {
    const fasta = ">seq\nAAA\nCCCC \nGG\n";
    const side_table = [_]index_format.SideTableLine{
        .{ .base_start = 0, .byte_offset = 5, .line_bytes = 4, .line_bases = 3 },
        .{ .base_start = 3, .byte_offset = 9, .line_bytes = 6, .line_bases = 4 },
        .{ .base_start = 7, .byte_offset = 15, .line_bytes = 3, .line_bases = 2 },
    };
    const rec = IndexRecord{
        .name_offset = 1,
        .name_len = 3,
        .seq_offset = 5,
        .seq_len = 9,
        .line_bases = 3,
        .line_bytes = 4,
    };

    var writer = std.Io.Writer.Allocating.init(std.testing.allocator);
    defer writer.deinit();

    emitRegion(.{
        .name = "seq",
        .start = 2,
        .display_end = 8,
        .is_full = false,
        .start_byte = byteOffsetForBase(fasta, rec, &side_table, 1),
        .seq_offset = 5,
        .num_bases = 7,
        .line_bases = 3,
        .line_bytes = 4,
        .side_table = &side_table,
        .orientation = .{},
        .annotate_transform = false,
        .original_index = 0,
    }, fasta, &writer.writer);

    try std.testing.expectEqualStrings(">seq:2-8\nAACCCCG\n", writer.written());
}

test "emitRegion handles non-uniform side-table reverse complement extraction" {
    const fasta = ">seq\nAAA\nCCCC \nGG\n";
    const side_table = [_]index_format.SideTableLine{
        .{ .base_start = 0, .byte_offset = 5, .line_bytes = 4, .line_bases = 3 },
        .{ .base_start = 3, .byte_offset = 9, .line_bytes = 6, .line_bases = 4 },
        .{ .base_start = 7, .byte_offset = 15, .line_bytes = 3, .line_bases = 2 },
    };
    const rec = IndexRecord{
        .name_offset = 1,
        .name_len = 3,
        .seq_offset = 5,
        .seq_len = 9,
        .line_bases = 3,
        .line_bytes = 4,
    };

    var writer = std.Io.Writer.Allocating.init(std.testing.allocator);
    defer writer.deinit();

    emitRegion(.{
        .name = "seq",
        .start = 2,
        .display_end = 8,
        .is_full = false,
        .start_byte = byteOffsetForBase(fasta, rec, &side_table, 1),
        .seq_offset = 5,
        .num_bases = 7,
        .line_bases = 3,
        .line_bytes = 4,
        .side_table = &side_table,
        .orientation = Orientation.reverseComplement(),
        .annotate_transform = false,
        .original_index = 0,
    }, fasta, &writer.writer);

    try std.testing.expectEqualStrings(">seq:2-8\nCGGGGTT\n", writer.written());
}
