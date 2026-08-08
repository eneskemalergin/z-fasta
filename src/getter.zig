//! Sequence extraction: positional regions, BED/names batching, strand, and RC transforms.
//!
//! Uniform records use O(1) byte-offset math; non-uniform records use side-table lookup.
//! A present `.zfi` is authoritative; `.fai` is used only when `.zfi` is absent.

const std = @import("std");
const complement = @import("complement.zig");
const bed_parser = @import("bed_parser.zig");
const index_format = @import("index_format.zig");
const stats = @import("stats.zig");

const LoadedIndex = index_format.LoadedIndex;
const printErrorAndExit = index_format.printErrorAndExit;

// --- Region parsing ---

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
    source: RequestSource,
    honor_strand: bool = false,
    summary: bool = false,
    orientation: Orientation = .{},
    annotate_transform: bool = false,
};

pub const RequestSource = union(enum) {
    positional: []const []const u8,
    names: []const u8,
    bed: []const u8,
};

const MAX_REQUEST_NAME_BYTES = std.math.maxInt(u16);
const REQUEST_LINE_READER_BUFFER_BYTES = MAX_REQUEST_NAME_BYTES + 4096;
const ACTIVE_NAME_BYTES = 4 * 1024 * 1024;
const MAX_ACTIVE_NAME_BYTES = ACTIVE_NAME_BYTES + MAX_REQUEST_NAME_BYTES;
const NAMES_REQUEST_BATCH_SIZE = 65536;
const BED_REQUEST_BATCH_SIZE = 4096;

const ParsedRequest = struct {
    region: Region,
    orientation: Orientation,

    fn parsed(self: ParsedRequest, _: []const u8) ParsedRequest {
        return self;
    }
};

const OwnedRequest = struct {
    name_offset: u32,
    name_len: u16,
    start: u64,
    end: u64,
    orientation: Orientation,

    fn parsed(self: OwnedRequest, names: []const u8) ParsedRequest {
        const name_offset: usize = self.name_offset;
        const name_end = name_offset + self.name_len;
        return .{
            .region = .{
                .name = names[name_offset..name_end],
                .start = self.start,
                .end = if (self.end == 0) null else self.end,
                .is_full = self.end == 0,
            },
            .orientation = self.orientation,
        };
    }
};

const RequestWorkspace = struct {
    names: std.ArrayList(u8) = .empty,
    requests: std.ArrayList(OwnedRequest) = .empty,
    scratch: std.heap.ArenaAllocator,

    fn init(allocator: std.mem.Allocator) RequestWorkspace {
        return .{ .scratch = std.heap.ArenaAllocator.init(allocator) };
    }

    fn deinit(self: *RequestWorkspace, allocator: std.mem.Allocator) void {
        self.names.deinit(allocator);
        self.requests.deinit(allocator);
        self.scratch.deinit();
    }

    fn clearRetainingCapacity(self: *RequestWorkspace) void {
        self.names.clearRetainingCapacity();
        self.requests.clearRetainingCapacity();
        _ = self.scratch.reset(.retain_capacity);
    }
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
    // Only the rightmost colon can start a coordinate suffix. Any earlier suffix
    // contains that colon, so its coordinate fields cannot both be integers.
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

        if (parseRangeSuffix(suffix)) |range| {
            return Region{
                .name = input[0..cp],
                .start = range.start,
                .end = range.end,
                .is_full = false,
            };
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

const ParsedRange = struct {
    start: u64,
    end: ?u64,
};

fn parseRangeSuffix(suffix: []const u8) ?ParsedRange {
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
        return .{ .start = start, .end = null };
    }

    const end = std.fmt.parseInt(u64, end_str, 10) catch return null;
    return .{ .start = start, .end = end };
}

// --- Region resolution ---

/// A validated extraction request with the record geometry needed for bounded reads.
pub const ResolvedRegion = struct {
    name: []const u8, // sequence name for FASTA header
    start: u64, // 1-based inclusive (validated)
    display_end: u64, // end value for header (pre-clamp, samtools convention)
    is_full: bool, // true -> emit ">NAME", false -> emit ">NAME:start-display_end"
    seq_offset: u64,
    num_bases: u64, // number of bases to extract
    line_bases: u32,
    line_bytes: u32,
    side_table: []const index_format.SideTableLine,
    record_index: usize,
    orientation: Orientation,
    annotate_transform: bool,
};

/// Resolve one region string against a loaded index.
/// Calls `printErrorAndExit` when the name or coordinates are invalid.
pub fn resolveRegion(idx: *const LoadedIndex, region_str: []const u8) ResolvedRegion {
    const region = parseRegion(region_str);

    const rec_idx = idx.lookupName(region.name) orelse {
        printErrorAndExit("error: sequence not found: {s}\n", .{region.name});
    };

    return resolveParsedRequest(
        idx,
        .{ .region = region, .orientation = .{} },
        rec_idx,
        false,
    );
}

fn resolveParsedRequest(
    idx: *const LoadedIndex,
    request: ParsedRequest,
    rec_idx: usize,
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

    return ResolvedRegion{
        .name = region.name,
        .start = start,
        .display_end = display_end,
        .is_full = region.is_full,
        .seq_offset = rec.seq_offset,
        .num_bases = num_bases,
        .line_bases = rec.line_bases,
        .line_bytes = rec.line_bytes,
        .side_table = side_table,
        .record_index = rec_idx,
        .orientation = request.orientation,
        .annotate_transform = annotate_transform,
    };
}

fn detectRegionType(source: FastaSource, resolved: ResolvedRegion) stats.SequenceType {
    const span = diskSpanForRegion(resolved);
    if (span.end > source.size or span.start >= span.end) {
        printErrorAndExit("error: read past end of FASTA\n", .{});
    }

    var counts = [_]u64{0} ** 256;
    var total: u64 = 0;
    var skip = span.leading_bases;
    var input: [FASTA_INPUT_BUFFER_BYTES]u8 = undefined;
    var offset = span.start;

    while (offset < span.end and total < resolved.num_bases) {
        const wanted: usize = @intCast(@min(@as(u64, input.len), span.end - offset));
        const got = std.Io.File.readPositionalAll(source.file, source.io, input[0..wanted], offset) catch {
            printErrorAndExit("error: failed to read FASTA\n", .{});
        };
        if (got != wanted) printErrorAndExit("error: read past end of FASTA\n", .{});
        for (input[0..got]) |byte| {
            if (byte <= ' ') continue;
            if (skip != 0) {
                skip -= 1;
                continue;
            }
            counts[byte] += 1;
            total += 1;
            if (total == resolved.num_bases) break;
        }
        offset += got;
    }

    if (total != resolved.num_bases) printErrorAndExit("error: read past end of FASTA\n", .{});
    return stats.detectType(&counts, total);
}

// --- Sequence emission ---

const FastaSource = struct {
    io: std.Io,
    file: std.Io.File,
    size: u64,
};

fn fastaSource(idx: *const LoadedIndex) FastaSource {
    return .{ .io = idx.io, .file = idx.fasta_file, .size = idx.fasta_size };
}

/// Writes a resolved region with standard 60-base wrapping.
pub fn extractRegion(idx: *const LoadedIndex, region_str: []const u8, writer: anytype) void {
    const resolved = resolveRegion(idx, region_str);
    var input_buffer: [FASTA_INPUT_BUFFER_BYTES]u8 = undefined;
    var output_buffer: [REGION_OUTPUT_BUFFER_BYTES]u8 = undefined;
    emitRegion(resolved, fastaSource(idx), &input_buffer, &output_buffer, writer);
}

const FASTA_INPUT_BUFFER_BYTES: usize = 256 * 1024;
const REGION_OUTPUT_BUFFER_BYTES: usize = 64 * 1024;

const DiskSpan = struct {
    start: u64,
    end: u64,
    leading_bases: u64,
    trailing_bases: u64,
};

fn spansOverlapOrTouch(first_start: u64, first_end: u64, second_start: u64, second_end: u64) bool {
    return first_start <= second_end and second_start <= first_end;
}

fn sideTableLineForBase(lines: []const index_format.SideTableLine, base_index: u64) usize {
    var lo: usize = 0;
    var hi = lines.len;
    while (lo < hi) {
        const mid = lo + (hi - lo) / 2;
        if (lines[mid].base_start <= base_index) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    return if (lo == 0) 0 else lo - 1;
}

fn uniformByteOffset(resolved: ResolvedRegion, base_index: u64) u64 {
    const line_number = base_index / resolved.line_bases;
    const column = base_index % resolved.line_bases;
    const line_delta = std.math.mul(u64, line_number, resolved.line_bytes) catch {
        printErrorAndExit("error: corrupt index geometry\n", .{});
    };
    const line_start = std.math.add(u64, resolved.seq_offset, line_delta) catch {
        printErrorAndExit("error: corrupt index geometry\n", .{});
    };
    return std.math.add(u64, line_start, column) catch {
        printErrorAndExit("error: corrupt index geometry\n", .{});
    };
}

fn diskSpanForRegion(resolved: ResolvedRegion) DiskSpan {
    const first_base = resolved.start - 1;
    const last_base = first_base + resolved.num_bases - 1;

    if (resolved.side_table.len == 0) {
        const start = uniformByteOffset(resolved, first_base);
        const last = uniformByteOffset(resolved, last_base);
        return .{
            .start = start,
            .end = std.math.add(u64, last, 1) catch printErrorAndExit("error: corrupt index geometry\n", .{}),
            .leading_bases = 0,
            .trailing_bases = 0,
        };
    }

    const first_line_idx = sideTableLineForBase(resolved.side_table, first_base);
    const last_line_idx = sideTableLineForBase(resolved.side_table, last_base);
    const first_line = resolved.side_table[first_line_idx];
    const last_line = resolved.side_table[last_line_idx];
    const last_line_end = std.math.add(u64, last_line.byte_offset, last_line.line_bytes) catch {
        printErrorAndExit("error: corrupt non-uniform index side table\n", .{});
    };
    return .{
        .start = first_line.byte_offset,
        .end = last_line_end,
        .leading_bases = first_base - first_line.base_start,
        .trailing_bases = last_line.base_start + last_line.line_bases - 1 - last_base,
    };
}

const RegionOutput = struct {
    line_pos: usize = 0,
    len: usize = 0,
    bases_written: u64 = 0,
    bytes: []u8,

    inline fn flush(self: *RegionOutput, writer: anytype) void {
        if (self.len == 0) return;
        writer.writeAll(self.bytes[0..self.len]) catch {
            printErrorAndExit("error: write failed\n", .{});
        };
        self.len = 0;
    }

    inline fn appendForward(self: *RegionOutput, src: []const u8, do_complement: bool, writer: anytype) void {
        var start: usize = 0;
        while (start < src.len) {
            const chunk = @min(src.len - start, 60 - self.line_pos);
            if (self.bytes.len - self.len < chunk + 1) self.flush(writer);
            const dst = self.bytes[self.len .. self.len + chunk];
            if (do_complement) {
                complement.complementInto(dst, src[start .. start + chunk]);
            } else {
                @memcpy(dst, src[start .. start + chunk]);
            }
            self.len += chunk;
            self.bases_written += chunk;
            self.line_pos += chunk;
            start += chunk;
            if (self.line_pos == 60) {
                self.bytes[self.len] = '\n';
                self.len += 1;
                self.line_pos = 0;
            }
        }
    }

    inline fn appendReverse(self: *RegionOutput, src: []const u8, do_complement: bool, writer: anytype) void {
        var end = src.len;
        while (end > 0) {
            const chunk = @min(end, 60 - self.line_pos);
            if (self.bytes.len - self.len < chunk + 1) self.flush(writer);
            var j: usize = 0;
            while (j < chunk) : (j += 1) {
                const byte = src[end - 1 - j];
                self.bytes[self.len + j] = if (do_complement) complement.complement(byte) else byte;
            }
            self.len += chunk;
            self.bases_written += chunk;
            self.line_pos += chunk;
            end -= chunk;
            if (self.line_pos == 60) {
                self.bytes[self.len] = '\n';
                self.len += 1;
                self.line_pos = 0;
            }
        }
    }

    fn finish(self: *RegionOutput, writer: anytype) void {
        if (self.line_pos != 0) {
            if (self.len == self.bytes.len) self.flush(writer);
            self.bytes[self.len] = '\n';
            self.len += 1;
        }
        self.flush(writer);
    }
};

fn writeRegionHeader(resolved: ResolvedRegion, writer: anytype) void {
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
}

fn consumeForwardBytes(
    resolved: ResolvedRegion,
    bytes: []const u8,
    skip: *u64,
    output: *RegionOutput,
    writer: anytype,
) void {
    var pos: usize = 0;
    while (pos < bytes.len and output.bases_written < resolved.num_bases) {
        while (pos < bytes.len and bytes[pos] <= ' ') pos += 1;
        const run_start = pos;
        while (pos < bytes.len and bytes[pos] > ' ') pos += 1;
        var run = bytes[run_start..pos];
        if (skip.* != 0) {
            const skipped: usize = @intCast(@min(skip.*, @as(u64, run.len)));
            skip.* -= @intCast(skipped);
            run = run[skipped..];
        }
        const remaining: usize = @intCast(resolved.num_bases - output.bases_written);
        run = run[0..@min(run.len, remaining)];
        output.appendForward(run, resolved.orientation.complement, writer);
    }
}

fn consumeReverseBytes(
    resolved: ResolvedRegion,
    bytes: []const u8,
    skip: *u64,
    output: *RegionOutput,
    writer: anytype,
) void {
    var i = bytes.len;
    while (i > 0 and output.bases_written < resolved.num_bases) {
        while (i > 0 and bytes[i - 1] <= ' ') i -= 1;
        const run_end = i;
        while (i > 0 and bytes[i - 1] > ' ') i -= 1;
        var run = bytes[i..run_end];
        if (skip.* != 0) {
            const skipped: usize = @intCast(@min(skip.*, @as(u64, run.len)));
            skip.* -= @intCast(skipped);
            run = run[0 .. run.len - skipped];
        }
        const remaining: usize = @intCast(resolved.num_bases - output.bases_written);
        if (run.len > remaining) run = run[run.len - remaining ..];
        output.appendReverse(run, resolved.orientation.complement, writer);
    }
}

fn emitUniformFromSpan(
    resolved: ResolvedRegion,
    span_start: u64,
    bytes: []const u8,
    output: *RegionOutput,
    writer: anytype,
) void {
    var remaining = resolved.num_bases;
    var base_index = if (resolved.orientation.reverse)
        resolved.start - 1 + resolved.num_bases - 1
    else
        resolved.start - 1;

    while (remaining != 0) {
        const line_number = base_index / resolved.line_bases;
        const column = base_index % resolved.line_bases;
        const line_start = uniformByteOffset(resolved, line_number * resolved.line_bases);

        const take_u64 = if (resolved.orientation.reverse)
            @min(column + 1, remaining)
        else
            @min(resolved.line_bases - column, remaining);
        const take: usize = @intCast(take_u64);
        const source_start = if (resolved.orientation.reverse)
            line_start + column + 1 - take_u64
        else
            line_start + column;
        if (source_start < span_start) printErrorAndExit("error: corrupt index geometry\n", .{});
        const relative_start: usize = @intCast(source_start - span_start);
        if (relative_start > bytes.len or take > bytes.len - relative_start) {
            printErrorAndExit("error: read past end of FASTA\n", .{});
        }
        const src = bytes[relative_start .. relative_start + take];
        if (resolved.orientation.reverse) {
            output.appendReverse(src, resolved.orientation.complement, writer);
        } else {
            output.appendForward(src, resolved.orientation.complement, writer);
        }

        remaining -= take_u64;
        if (remaining == 0) break;
        if (resolved.orientation.reverse) {
            base_index -= take_u64;
        } else {
            base_index += take_u64;
        }
    }
}

fn consumeUniformForwardChunk(
    resolved: ResolvedRegion,
    chunk_start: u64,
    bytes: []const u8,
    base_index: *u64,
    output: *RegionOutput,
    writer: anytype,
) void {
    const chunk_end = std.math.add(u64, chunk_start, bytes.len) catch {
        printErrorAndExit("error: corrupt index geometry\n", .{});
    };
    while (output.bases_written < resolved.num_bases) {
        const column = base_index.* % resolved.line_bases;
        const source_start = uniformByteOffset(resolved, base_index.*);
        if (source_start >= chunk_end) return;
        if (source_start < chunk_start) printErrorAndExit("error: corrupt index geometry\n", .{});
        const relative_start: usize = @intCast(source_start - chunk_start);
        const available = @min(
            resolved.line_bases - column,
            resolved.num_bases - output.bases_written,
        );
        const take: usize = @intCast(@min(available, chunk_end - source_start));
        output.appendForward(bytes[relative_start .. relative_start + take], resolved.orientation.complement, writer);
        base_index.* += take;
    }
}

fn consumeUniformReverseChunk(
    resolved: ResolvedRegion,
    chunk_start: u64,
    bytes: []const u8,
    base_index: *u64,
    output: *RegionOutput,
    writer: anytype,
) void {
    const chunk_end = std.math.add(u64, chunk_start, bytes.len) catch {
        printErrorAndExit("error: corrupt index geometry\n", .{});
    };
    while (output.bases_written < resolved.num_bases) {
        const column = base_index.* % resolved.line_bases;
        const source_end = std.math.add(u64, uniformByteOffset(resolved, base_index.*), 1) catch {
            printErrorAndExit("error: corrupt index geometry\n", .{});
        };
        if (source_end > chunk_end) return;
        if (source_end <= chunk_start) return;
        const available = @min(
            column + 1,
            resolved.num_bases - output.bases_written,
        );
        const take: usize = @intCast(@min(available, source_end - chunk_start));
        const source_start = source_end - take;
        const relative_start: usize = @intCast(source_start - chunk_start);
        output.appendReverse(bytes[relative_start .. relative_start + take], resolved.orientation.complement, writer);
        if (output.bases_written < resolved.num_bases) base_index.* -= take;
    }
}

fn emitRegionFromSpan(
    resolved: ResolvedRegion,
    span: DiskSpan,
    span_start: u64,
    bytes: []const u8,
    output_buffer: []u8,
    writer: anytype,
) void {
    writeRegionHeader(resolved, writer);
    const relative_start = std.math.cast(usize, span.start - span_start) orelse {
        printErrorAndExit("error: FASTA span is too large\n", .{});
    };
    const relative_end = std.math.cast(usize, span.end - span_start) orelse {
        printErrorAndExit("error: FASTA span is too large\n", .{});
    };
    if (relative_end > bytes.len or relative_start > relative_end) {
        printErrorAndExit("error: read past end of FASTA\n", .{});
    }
    var output = RegionOutput{ .bytes = output_buffer };
    if (resolved.side_table.len == 0) {
        emitUniformFromSpan(resolved, span_start, bytes, &output, writer);
    } else {
        const region_bytes = bytes[relative_start..relative_end];
        var skip = if (resolved.orientation.reverse) span.trailing_bases else span.leading_bases;
        if (resolved.orientation.reverse) {
            consumeReverseBytes(resolved, region_bytes, &skip, &output, writer);
        } else {
            consumeForwardBytes(resolved, region_bytes, &skip, &output, writer);
        }
    }
    if (output.bases_written != resolved.num_bases) {
        printErrorAndExit("error: read past end of FASTA\n", .{});
    }
    output.finish(writer);
}

fn emitRegion(
    resolved: ResolvedRegion,
    source: FastaSource,
    input_buffer: []u8,
    output_buffer: []u8,
    writer: anytype,
) void {
    const span = diskSpanForRegion(resolved);
    if (span.end > source.size or span.start >= span.end) {
        printErrorAndExit("error: read past end of FASTA\n", .{});
    }

    writeRegionHeader(resolved, writer);
    var output = RegionOutput{ .bytes = output_buffer };
    var skip = if (resolved.orientation.reverse) span.trailing_bases else span.leading_bases;
    var base_index = if (resolved.orientation.reverse)
        resolved.start - 1 + resolved.num_bases - 1
    else
        resolved.start - 1;

    if (!resolved.orientation.reverse) {
        var offset = span.start;
        while (offset < span.end and output.bases_written < resolved.num_bases) {
            const wanted: usize = @intCast(@min(@as(u64, input_buffer.len), span.end - offset));
            const got = std.Io.File.readPositionalAll(source.file, source.io, input_buffer[0..wanted], offset) catch {
                printErrorAndExit("error: failed to read FASTA\n", .{});
            };
            if (got != wanted) printErrorAndExit("error: read past end of FASTA\n", .{});
            if (resolved.side_table.len == 0) {
                consumeUniformForwardChunk(resolved, offset, input_buffer[0..got], &base_index, &output, writer);
            } else {
                consumeForwardBytes(resolved, input_buffer[0..got], &skip, &output, writer);
            }
            offset += got;
        }
    } else {
        var end = span.end;
        while (end > span.start and output.bases_written < resolved.num_bases) {
            const wanted: usize = @intCast(@min(@as(u64, input_buffer.len), end - span.start));
            const offset = end - wanted;
            const got = std.Io.File.readPositionalAll(source.file, source.io, input_buffer[0..wanted], offset) catch {
                printErrorAndExit("error: failed to read FASTA\n", .{});
            };
            if (got != wanted) printErrorAndExit("error: read past end of FASTA\n", .{});
            if (resolved.side_table.len == 0) {
                consumeUniformReverseChunk(resolved, offset, input_buffer[0..got], &base_index, &output, writer);
            } else {
                consumeReverseBytes(resolved, input_buffer[0..got], &skip, &output, writer);
            }
            end = offset;
        }
    }

    if (output.bases_written != resolved.num_bases) {
        printErrorAndExit("error: read past end of FASTA\n", .{});
    }
    output.finish(writer);
}

fn emitResolvedBatch(
    resolved: []const ResolvedRegion,
    source: FastaSource,
    writer: anytype,
) u64 {
    var total_bases: u64 = 0;
    var shared: [FASTA_INPUT_BUFFER_BYTES]u8 = undefined;
    var output_buffer: [REGION_OUTPUT_BUFFER_BYTES]u8 = undefined;
    var i: usize = 0;

    while (i < resolved.len) {
        const first_span = diskSpanForRegion(resolved[i]);
        var shared_start = first_span.start;
        var shared_end = first_span.end;
        const shared_record = resolved[i].record_index;
        var shared_first_base = resolved[i].start - 1;
        var shared_last_base = shared_first_base + resolved[i].num_bases - 1;
        var group_end = i + 1;

        while (group_end < resolved.len) : (group_end += 1) {
            const next = diskSpanForRegion(resolved[group_end]);
            const next_first_base = resolved[group_end].start - 1;
            const next_last_base = next_first_base + resolved[group_end].num_bases - 1;
            const merged_start = @min(shared_start, next.start);
            const merged_end = @max(shared_end, next.end);
            const separated_after = next_first_base > shared_last_base and next_first_base - shared_last_base > 1;
            const separated_before = shared_first_base > next_last_base and shared_first_base - next_last_base > 1;
            const logical_neighbor = resolved[group_end].record_index == shared_record and
                !separated_after and !separated_before;
            const disk_neighbor = spansOverlapOrTouch(shared_start, shared_end, next.start, next.end);
            if ((!logical_neighbor and !disk_neighbor) or merged_end - merged_start > shared.len) break;
            shared_start = merged_start;
            shared_end = merged_end;
            shared_first_base = @min(shared_first_base, next_first_base);
            shared_last_base = @max(shared_last_base, next_last_base);
        }

        if (shared_end - shared_start <= shared.len) {
            if (shared_end > source.size) printErrorAndExit("error: read past end of FASTA\n", .{});
            const wanted: usize = @intCast(shared_end - shared_start);
            const got = std.Io.File.readPositionalAll(source.file, source.io, shared[0..wanted], shared_start) catch {
                printErrorAndExit("error: failed to read FASTA\n", .{});
            };
            if (got != wanted) printErrorAndExit("error: read past end of FASTA\n", .{});
            for (resolved[i..group_end]) |region| {
                const span = diskSpanForRegion(region);
                emitRegionFromSpan(region, span, shared_start, shared[0..got], &output_buffer, writer);
                total_bases += region.num_bases;
            }
            i = group_end;
            continue;
        }

        emitRegion(resolved[i], source, &shared, &output_buffer, writer);
        total_bases += resolved[i].num_bases;
        i += 1;
    }

    return total_bases;
}

fn headerAnnotation(orientation: Orientation, annotate_transform: bool) []const u8 {
    if (!annotate_transform or orientation.isIdentity()) return "";
    if (orientation.reverse and orientation.complement) return " (reverse complement)";
    if (orientation.reverse) return " (reverse)";
    if (orientation.complement) return " (complement)";
    return "";
}

fn monotonicNs(io: std.Io) u64 {
    const now = std.Io.Clock.Timestamp.now(io, .awake);
    return @intCast(now.raw.toNanoseconds());
}

const OwnedName = struct {
    offset: u32,
    len: u16,
};

fn ownRequestName(
    allocator: std.mem.Allocator,
    workspace: *RequestWorkspace,
    name: []const u8,
) OwnedName {
    if (workspace.requests.getLastOrNull()) |last| {
        const offset: usize = last.name_offset;
        const end = offset + last.name_len;
        if (std.mem.eql(u8, workspace.names.items[offset..end], name)) {
            return .{ .offset = last.name_offset, .len = last.name_len };
        }
    }

    const offset = workspace.names.items.len;
    const required_capacity = offset + name.len;
    if (required_capacity > workspace.names.capacity) {
        const doubled = @max(workspace.names.capacity, 4096) * 2;
        const new_capacity = @min(MAX_ACTIVE_NAME_BYTES, @max(required_capacity, doubled));
        workspace.names.ensureTotalCapacityPrecise(allocator, new_capacity) catch {
            printErrorAndExit("error: out of memory\n", .{});
        };
    }
    workspace.names.appendSliceAssumeCapacity(name);
    return .{ .offset = @intCast(offset), .len = @intCast(name.len) };
}

fn appendOwnedRequest(
    allocator: std.mem.Allocator,
    workspace: *RequestWorkspace,
    name: []const u8,
    start: u64,
    end: u64,
    orientation: Orientation,
) void {
    const owned_name = ownRequestName(allocator, workspace, name);
    workspace.requests.append(allocator, .{
        .name_offset = owned_name.offset,
        .name_len = owned_name.len,
        .start = start,
        .end = end,
        .orientation = orientation,
    }) catch {
        printErrorAndExit("error: out of memory\n", .{});
    };
}

fn appendBedLineRequest(
    allocator: std.mem.Allocator,
    workspace: *RequestWorkspace,
    line: []const u8,
    line_number: usize,
    honor_strand: bool,
    global_orientation: Orientation,
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
        .region => |region| {
            if (honor_strand and region.strand == .invalid) {
                printErrorAndExit("error: invalid BED line {d}: invalid strand\n", .{line_number});
            }
            if (region.chrom.len > MAX_REQUEST_NAME_BYTES) {
                printErrorAndExit(
                    "error: invalid BED line {d}: chrom exceeds {d} bytes\n",
                    .{ line_number, MAX_REQUEST_NAME_BYTES },
                );
            }
            appendOwnedRequest(
                allocator,
                workspace,
                region.chrom,
                region.start1Based(),
                region.end1BasedInclusive(),
                (if (honor_strand and region.strand == .minus) Orientation.reverseComplement() else Orientation{}).compose(global_orientation),
            );
        },
    }
}

fn appendNamesLine(
    allocator: std.mem.Allocator,
    workspace: *RequestWorkspace,
    line: []const u8,
    line_number: usize,
    orientation: Orientation,
) void {
    const name = if (line.len > 0 and line[line.len - 1] == '\r') line[0 .. line.len - 1] else line;
    if (name.len == 0 or name[0] == '#') return;
    if (name.len > MAX_REQUEST_NAME_BYTES) {
        printErrorAndExit(
            "error: name at line {d} exceeds {d} bytes\n",
            .{ line_number, MAX_REQUEST_NAME_BYTES },
        );
    }
    appendOwnedRequest(allocator, workspace, name, 1, 0, orientation);
}

fn appendCliRequests(
    allocator: std.mem.Allocator,
    requests: *std.ArrayList(ParsedRequest),
    region_strs: []const []const u8,
    orientation: Orientation,
) void {
    for (region_strs) |region_str| {
        requests.append(allocator, .{
            .region = parseRegion(region_str),
            .orientation = orientation,
        }) catch {
            printErrorAndExit("error: out of memory\n", .{});
        };
    }
}

fn writeSummary(writer: *std.Io.Writer, region_count: usize, total_bases: u64, elapsed_ns: u64) !void {
    const seconds = if (elapsed_ns == 0) 0.0 else @as(f64, @floatFromInt(elapsed_ns)) / @as(f64, @floatFromInt(std.time.ns_per_s));
    const regions_per_second = if (seconds == 0.0) 0.0 else @as(f64, @floatFromInt(region_count)) / seconds;

    try writer.print("summary: regions={d} total_bases={d} elapsed_s={d:.6} regions_per_s={d:.1}\n", .{ region_count, total_bases, seconds, regions_per_second });
}

fn processParsedRequests(
    allocator: std.mem.Allocator,
    idx: *index_format.LoadedIndex,
    request_entries: anytype,
    active_names: []const u8,
    annotate_transform: bool,
    writer: anytype,
    source: FastaSource,
) BatchStats {
    if (request_entries.len == 0) return .{};

    const resolved = allocator.alloc(ResolvedRegion, request_entries.len) catch {
        printErrorAndExit("error: out of memory\n", .{});
    };
    defer allocator.free(resolved);

    var last_name: ?[]const u8 = null;
    var last_rec_idx: usize = 0;
    var last_typed_rec_idx: ?usize = null;
    var last_type: stats.SequenceType = undefined;

    for (request_entries, 0..) |entry, i| {
        const request = entry.parsed(active_names);
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
        resolved[i] = resolveParsedRequest(idx, request, rec_idx, annotate_transform);

        if (request.orientation.complement) {
            const rec_type = if (last_typed_rec_idx == rec_idx)
                last_type
            else blk: {
                var sample = resolved[i];
                sample.num_bases = @min(sample.num_bases, stats.GET_TYPE_SAMPLE_BASES);
                const detected = detectRegionType(source, sample);
                last_typed_rec_idx = rec_idx;
                last_type = detected;
                break :blk detected;
            };
            if (rec_type == .protein) {
                printErrorAndExit(
                    "error: reverse complement is not defined for protein sequences: {s} (classified from up to {d} bases)\n",
                    .{ request.region.name, stats.GET_TYPE_SAMPLE_BASES },
                );
            }
        }
    }

    return .{
        .region_count = request_entries.len,
        .total_bases = emitResolvedBatch(resolved, source, writer),
    };
}

const StreamRequestSource = union(enum) {
    names: Orientation,
    bed: struct {
        honor_strand: bool,
        global_orientation: Orientation,
    },
};

fn processRequestReader(
    allocator: std.mem.Allocator,
    idx: *LoadedIndex,
    reader: *std.Io.Reader,
    source: StreamRequestSource,
    name_reservation: usize,
    annotate_transform: bool,
    writer: anytype,
    fasta_source: FastaSource,
) BatchStats {
    var total = BatchStats{};
    var line_number: usize = 0;
    var reached_end = false;
    var workspace = RequestWorkspace.init(allocator);
    defer workspace.deinit(allocator);

    const request_limit: usize = switch (source) {
        .names => NAMES_REQUEST_BATCH_SIZE,
        .bed => BED_REQUEST_BATCH_SIZE,
    };
    workspace.requests.ensureTotalCapacityPrecise(allocator, request_limit) catch {
        printErrorAndExit("error: out of memory\n", .{});
    };
    if (name_reservation > 0) {
        workspace.names.ensureTotalCapacityPrecise(allocator, name_reservation) catch {
            printErrorAndExit("error: out of memory\n", .{});
        };
    }

    while (true) {
        while (workspace.requests.items.len < request_limit and workspace.names.items.len < ACTIVE_NAME_BYTES) {
            const maybe_line = reader.takeDelimiter('\n') catch |err| switch (err) {
                error.ReadFailed => switch (source) {
                    .names => printErrorAndExit("error: failed to read names input\n", .{}),
                    .bed => printErrorAndExit("error: failed to read BED input\n", .{}),
                },
                error.StreamTooLong => switch (source) {
                    .names => printErrorAndExit(
                        "error: name at line {d} exceeds {d} bytes\n",
                        .{ line_number + 1, MAX_REQUEST_NAME_BYTES },
                    ),
                    .bed => printErrorAndExit(
                        "error: BED line {d} exceeds {d}-byte reader buffer (no newline within limit)\n",
                        .{ line_number + 1, REQUEST_LINE_READER_BUFFER_BYTES },
                    ),
                },
            };

            const line = maybe_line orelse {
                reached_end = true;
                break;
            };

            line_number += 1;
            switch (source) {
                .names => |orientation| appendNamesLine(allocator, &workspace, line, line_number, orientation),
                .bed => |bed| appendBedLineRequest(allocator, &workspace, line, line_number, bed.honor_strand, bed.global_orientation),
            }
        }

        if (workspace.requests.items.len == 0) break;

        const batch = processParsedRequests(
            workspace.scratch.allocator(),
            idx,
            workspace.requests.items,
            workspace.names.items,
            annotate_transform,
            writer,
            fasta_source,
        );
        total.region_count += batch.region_count;
        total.total_bases += batch.total_bases;

        workspace.clearRetainingCapacity();

        if (reached_end) break;
    }

    return total;
}

fn processRequestPath(
    allocator: std.mem.Allocator,
    io: std.Io,
    idx: *LoadedIndex,
    path: []const u8,
    source: StreamRequestSource,
    annotate_transform: bool,
    writer: anytype,
    fasta_source: FastaSource,
) BatchStats {
    if (std.mem.eql(u8, path, "-")) {
        var stdin_buf: [REQUEST_LINE_READER_BUFFER_BYTES]u8 = undefined;
        var stdin_reader = std.Io.File.stdin().reader(io, &stdin_buf);
        return processRequestReader(allocator, idx, &stdin_reader.interface, source, 0, annotate_transform, writer, fasta_source);
    }

    const file = std.Io.Dir.cwd().openFile(io, path, .{}) catch |err| switch (err) {
        error.FileNotFound => printErrorAndExit("error: file not found: {s}\n", .{path}),
        error.AccessDenied => printErrorAndExit("error: access denied: {s}\n", .{path}),
        else => printErrorAndExit("error: failed to open file: {s}\n", .{path}),
    };
    defer file.close(io);

    const name_reservation: usize = switch (source) {
        .names => blk: {
            const stat = file.stat(io) catch break :blk 0;
            break :blk @intCast(@min(stat.size, MAX_ACTIVE_NAME_BYTES));
        },
        .bed => 0,
    };

    var file_buf: [REQUEST_LINE_READER_BUFFER_BYTES]u8 = undefined;
    var file_reader = file.reader(io, &file_buf);
    return processRequestReader(allocator, idx, &file_reader.interface, source, name_reservation, annotate_transform, writer, fasta_source);
}

// --- Public entry point ---

pub fn runGetWithOptions(
    backing_allocator: std.mem.Allocator,
    io: std.Io,
    fasta_path: []const u8,
    options: GetOptions,
) void {
    const start_ns = if (options.summary) monotonicNs(io) else 0;

    var input_arena = std.heap.ArenaAllocator.init(backing_allocator);
    defer input_arena.deinit();
    const allocator = input_arena.allocator();

    var requests = std.ArrayList(ParsedRequest).empty;

    // Positional requests are known before loading, so `.fai` retains only their
    // first matching records. Names and BED still require the complete lookup map.
    const load_mode: index_format.LoadMode = switch (options.source) {
        .positional => |region_strs| blk: {
            appendCliRequests(allocator, &requests, region_strs, options.orientation);
            const names = allocator.alloc([]const u8, requests.items.len) catch {
                printErrorAndExit("error: out of memory\n", .{});
            };
            for (requests.items, names) |request, *name| {
                name.* = request.region.name;
            }
            break :blk .{ .positional = names };
        },
        .names, .bed => .lookup_full_map,
    };
    var idx = index_format.loadIndexForGet(backing_allocator, io, fasta_path, load_mode);
    defer idx.deinit();

    const fasta_source = fastaSource(&idx);

    var out_buf: [65536]u8 = undefined;
    var stdout_fw = std.Io.File.Writer.initStreaming(.stdout(), io, &out_buf);
    const writer = &stdout_fw.interface;

    const totals = switch (options.source) {
        .positional => processParsedRequests(allocator, &idx, requests.items, &.{}, options.annotate_transform, writer, fasta_source),
        .names => |path| processRequestPath(
            backing_allocator,
            io,
            &idx,
            path,
            .{ .names = options.orientation },
            options.annotate_transform,
            writer,
            fasta_source,
        ),
        .bed => |path| processRequestPath(
            backing_allocator,
            io,
            &idx,
            path,
            .{ .bed = .{
                .honor_strand = options.honor_strand,
                .global_orientation = options.orientation,
            } },
            options.annotate_transform,
            writer,
            fasta_source,
        ),
    };

    if (totals.region_count == 0) {
        printErrorAndExit("error: no regions provided\n", .{});
    }

    stdout_fw.flush() catch {
        printErrorAndExit("error: write failed\n", .{});
    };

    if (options.summary) {
        var err_buf: [512]u8 = undefined;
        var stderr_fw = std.Io.File.Writer.initStreaming(.stderr(), io, &err_buf);
        writeSummary(&stderr_fw.interface, totals.region_count, totals.total_bases, monotonicNs(io) - start_ns) catch {
            printErrorAndExit("error: write failed\n", .{});
        };
        stderr_fw.flush() catch {
            printErrorAndExit("error: write failed\n", .{});
        };
    }
}

test "summary writing fails when the destination is full" {
    var buffer: [1]u8 = undefined;
    var writer = std.Io.Writer.fixed(&buffer);

    try std.testing.expectError(error.WriteFailed, writeSummary(&writer, 1, 1, 1));
}

test "shared disk spans require overlap or adjacency in either order" {
    try std.testing.expect(spansOverlapOrTouch(10, 20, 20, 30));
    try std.testing.expect(spansOverlapOrTouch(20, 30, 10, 20));
    try std.testing.expect(spansOverlapOrTouch(10, 20, 15, 25));
    try std.testing.expect(spansOverlapOrTouch(15, 25, 10, 20));
    try std.testing.expect(!spansOverlapOrTouch(10, 20, 21, 30));
    try std.testing.expect(!spansOverlapOrTouch(21, 30, 10, 20));
}

test "request workspace reuses consecutive names" {
    var workspace = RequestWorkspace.init(std.testing.allocator);
    defer workspace.deinit(std.testing.allocator);

    appendNamesLine(std.testing.allocator, &workspace, "seq1", 1, .{});
    appendNamesLine(std.testing.allocator, &workspace, "seq1", 2, .{});
    appendNamesLine(std.testing.allocator, &workspace, "seq2", 3, .{});

    try std.testing.expectEqualStrings("seq1seq2", workspace.names.items);
    try std.testing.expectEqual(@as(usize, 3), workspace.requests.items.len);
    try std.testing.expectEqualStrings("seq1", workspace.requests.items[0].parsed(workspace.names.items).region.name);
    try std.testing.expectEqualStrings("seq1", workspace.requests.items[1].parsed(workspace.names.items).region.name);
    try std.testing.expectEqualStrings("seq2", workspace.requests.items[2].parsed(workspace.names.items).region.name);
}

test "request workspace includes the name that crosses its byte target" {
    var workspace = RequestWorkspace.init(std.testing.allocator);
    defer workspace.deinit(std.testing.allocator);
    const name = try std.testing.allocator.alloc(u8, MAX_REQUEST_NAME_BYTES);
    defer std.testing.allocator.free(name);
    @memset(name, 'A');

    for (0..64) |i| {
        name[0] = if (i % 2 == 0) 'A' else 'B';
        appendNamesLine(std.testing.allocator, &workspace, name, i + 1, .{});
    }
    try std.testing.expect(workspace.names.items.len < ACTIVE_NAME_BYTES);

    name[0] = 'A';
    appendNamesLine(std.testing.allocator, &workspace, name, 65, .{});

    try std.testing.expect(workspace.names.items.len >= ACTIVE_NAME_BYTES);
    try std.testing.expect(workspace.names.items.len <= MAX_ACTIVE_NAME_BYTES);
    try std.testing.expect(workspace.names.capacity <= MAX_ACTIVE_NAME_BYTES);
}

test "processRequestReader handles names line endings and final line" {
    const test_io = std.testing.io;
    var idx = index_format.loadIndex(std.testing.allocator, test_io, "tests/data/simple.fasta");
    defer idx.deinit();
    var reader = std.Io.Reader.fixed("# comment\r\n\r\nseq2\r\nseq1\nseq2");
    var writer = std.Io.Writer.Allocating.init(std.testing.allocator);
    defer writer.deinit();

    const result = processRequestReader(
        std.testing.allocator,
        &idx,
        &reader,
        .{ .names = .{} },
        0,
        false,
        &writer.writer,
        fastaSource(&idx),
    );

    try std.testing.expectEqual(@as(usize, 3), result.region_count);
    try std.testing.expectEqual(@as(u64, 48), result.total_bases);
    try std.testing.expectEqualStrings(
        ">seq2\nGGGGCCCCAAAA\n>seq1\nACGTACGTACGTACGTACGTACGT\n>seq2\nGGGGCCCCAAAA\n",
        writer.written(),
    );
}

test "processRequestReader resets names at the request-count boundary" {
    const test_io = std.testing.io;
    var idx = index_format.loadIndex(std.testing.allocator, test_io, "tests/data/simple.fasta");
    defer idx.deinit();
    var names = std.Io.Writer.Allocating.init(std.testing.allocator);
    defer names.deinit();
    const request_count = NAMES_REQUEST_BATCH_SIZE + 1;
    for (0..request_count) |_| try names.writer.writeAll("seq1\n");
    var reader = std.Io.Reader.fixed(names.written());
    var writer = std.Io.Writer.Allocating.init(std.testing.allocator);
    defer writer.deinit();

    const result = processRequestReader(
        std.testing.allocator,
        &idx,
        &reader,
        .{ .names = .{} },
        0,
        false,
        &writer.writer,
        fastaSource(&idx),
    );

    const one_output = ">seq1\nACGTACGTACGTACGTACGTACGT\n";
    try std.testing.expectEqual(request_count, result.region_count);
    try std.testing.expectEqual(@as(u64, request_count * 24), result.total_bases);
    try std.testing.expectEqual(one_output.len * request_count, writer.written().len);
    try std.testing.expectEqualStrings(one_output, writer.written()[0..one_output.len]);
    try std.testing.expectEqualStrings(one_output, writer.written()[writer.written().len - one_output.len ..]);
}

test "processRequestReader streams BED rows in source order" {
    const test_io = std.testing.io;
    var idx = index_format.loadIndex(std.testing.allocator, test_io, "tests/data/simple.fasta");
    defer idx.deinit();

    const bed_data =
        "# comment\r\n" ++
        "seq2\t0\t4\r\n" ++
        "seq1\t0\t4\n" ++
        "seq1\t4\t8";

    var reader = std.Io.Reader.fixed(bed_data);
    var writer = std.Io.Writer.Allocating.init(std.testing.allocator);
    defer writer.deinit();

    const result = processRequestReader(
        std.testing.allocator,
        &idx,
        &reader,
        .{ .bed = .{ .honor_strand = false, .global_orientation = .{} } },
        0,
        false,
        &writer.writer,
        fastaSource(&idx),
    );

    try std.testing.expectEqual(@as(usize, 3), result.region_count);
    try std.testing.expectEqual(@as(u64, 12), result.total_bases);
    try std.testing.expectEqualStrings(
        ">seq2:1-4\nGGGG\n>seq1:1-4\nACGT\n>seq1:5-8\nACGT\n",
        writer.written(),
    );
}

test "processRequestReader preserves output at BED request boundaries" {
    const test_io = std.testing.io;
    var idx = index_format.loadIndex(std.testing.allocator, test_io, "tests/data/simple.fasta");
    defer idx.deinit();

    const counts = [_]usize{
        BED_REQUEST_BATCH_SIZE - 1,
        BED_REQUEST_BATCH_SIZE,
        BED_REQUEST_BATCH_SIZE + 1,
    };
    for (counts) |count| {
        var bed = std.Io.Writer.Allocating.init(std.testing.allocator);
        defer bed.deinit();
        var expected = std.Io.Writer.Allocating.init(std.testing.allocator);
        defer expected.deinit();

        for (0..count) |i| {
            if (count > BED_REQUEST_BATCH_SIZE and i + 1 == count) {
                try bed.writer.writeAll("seq1\t0\t5\tname\t0\t-\n");
                try expected.writer.writeAll(">seq1:1-5\nTACGT\n");
            } else {
                try bed.writer.writeAll("seq1\t0\t1\tname\t0\t+\n");
                try expected.writer.writeAll(">seq1:1-1\nA\n");
            }
        }

        var reader = std.Io.Reader.fixed(bed.written());
        var output = std.Io.Writer.Allocating.init(std.testing.allocator);
        defer output.deinit();
        const result = processRequestReader(
            std.testing.allocator,
            &idx,
            &reader,
            .{ .bed = .{ .honor_strand = true, .global_orientation = .{} } },
            0,
            false,
            &output.writer,
            fastaSource(&idx),
        );

        const expected_bases: u64 = count + @as(usize, @intFromBool(count > BED_REQUEST_BATCH_SIZE)) * 4;
        try std.testing.expectEqual(count, result.region_count);
        try std.testing.expectEqual(expected_bases, result.total_bases);
        try std.testing.expectEqualStrings(expected.written(), output.written());
    }
}

test "parseRegion matches exhaustive backward-scan reference" {
    const Reference = struct {
        fn findLastColon(input: []const u8, end: usize) ?usize {
            var i = end;
            while (i > 0) {
                i -= 1;
                if (input[i] == ':') return i;
            }
            return null;
        }

        fn parse(input: []const u8) Region {
            var search_end = input.len;
            while (findLastColon(input, search_end)) |colon_pos| {
                if (parseRangeSuffix(input[colon_pos + 1 ..])) |range| {
                    return .{
                        .name = input[0..colon_pos],
                        .start = range.start,
                        .end = range.end,
                        .is_full = false,
                    };
                }
                search_end = colon_pos;
            }
            return .{ .name = input, .start = 1, .end = null, .is_full = true };
        }
    };

    const alphabet = "a:0-1";
    var input: [7]u8 = undefined;
    var combinations: usize = 1;
    var tested: usize = 0;
    for (1..input.len + 1) |length| {
        combinations *= alphabet.len;
        for (0..combinations) |encoded| {
            var value = encoded;
            for (input[0..length]) |*byte| {
                byte.* = alphabet[value % alphabet.len];
                value /= alphabet.len;
            }

            const expected = Reference.parse(input[0..length]);
            const actual = parseRegion(input[0..length]);
            try std.testing.expectEqualStrings(expected.name, actual.name);
            try std.testing.expectEqual(expected.start, actual.start);
            try std.testing.expectEqual(expected.end, actual.end);
            try std.testing.expectEqual(expected.is_full, actual.is_full);
            tested += 1;
        }
    }
    try std.testing.expectEqual(@as(usize, 97_655), tested);
}

test "parseRegion handles coordinate overflow at the rightmost colon" {
    const max = parseRegion("name:18446744073709551615-18446744073709551615");
    try std.testing.expectEqualStrings("name", max.name);
    try std.testing.expectEqual(@as(u64, std.math.maxInt(u64)), max.start);
    try std.testing.expectEqual(@as(?u64, std.math.maxInt(u64)), max.end);
    try std.testing.expect(!max.is_full);

    const full_names = [_][]const u8{
        "name:18446744073709551616-1",
        "name:1-18446744073709551616",
        "name:1-2:bad",
        "name:1-2:18446744073709551616-1",
    };
    for (full_names) |name| {
        const region = parseRegion(name);
        try std.testing.expectEqualStrings(name, region.name);
        try std.testing.expect(region.is_full);
    }
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
