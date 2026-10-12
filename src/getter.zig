//! Extracts sequence for `z-fasta get` and for `IndexedReader.fetchInto`.
//! Requests are one-based inclusive; an end past the record is clamped, and a suffix that does
//! not parse as a range stays part of the name. Lookups are first-wins.
//! Byte spans come from the uniform line formula or a side-table search; reads are positional.
//! Output is wrapped at 60 and written per record, so a later failure can leave a valid prefix.

const std = @import("std");
const fasta = @import("fasta.zig");
const bed_parser = @import("bed_parser.zig");
const index_format = @import("index_format.zig");
const diagnostics = @import("diagnostics");

pub const Region = struct {
    name: []const u8,
    start: u64,
    end: ?u64,
    is_full: bool,
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

/// Not validated here; callers pass only combinations the CLI accepts.
pub const Options = struct {
    source: RequestSource,
    honor_strand: bool = false,
    measure_elapsed: bool = false,
    orientation: Orientation = .{},
    annotate_transform: bool = false,
};

/// A names or BED path of `-` reads stdin.
pub const RequestSource = union(enum) {
    positional: []const []const u8,
    names: []const u8,
    bed: []const u8,
};

/// `elapsed_ns` is zero unless `measure_elapsed` was set.
pub const Result = struct { region_count: usize, total_bases: u64, elapsed_ns: u64 };

pub fn parseRegion(input: []const u8) Region {
    if (std.mem.findScalarLast(u8, input, ':')) |cp| {
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

    return Region{
        .name = input,
        .start = 1,
        .end = null,
        .is_full = true,
    };
}

pub const ResolvedRegion = struct {
    name: []const u8,
    start: u64,
    display_end: u64,
    is_full: bool,
    seq_offset: u64,
    num_bases: u64,
    line_bases: u32,
    line_bytes: u32,
    side_table: []const index_format.SideTableLine,
    record_index: usize,
    orientation: Orientation,
    annotate_transform: bool,
};

pub fn resolveRegion(idx: *const LoadedIndex, region_str: []const u8, failure: *diagnostics.Failure) diagnostics.Error!ResolvedRegion {
    failure.reset(.get);
    failure.sources.set(.request_snapshot, .{ .label = .{ .borrowed = region_str }, .kind = .memory });
    var context = Diagnostics{ .failure = failure, .lookup_scope = .retained_matches, .index_source = if (idx.source == .zfi) .zfi else .fai, .at = .{ .source = .request_snapshot, .request = 1 } };
    const result = try resolveRegionText(idx, region_str, &context);
    failure.succeed();
    return result;
}

/// The caller owns `writer` and its flush.
pub fn extractRegion(idx: *const LoadedIndex, region_str: []const u8, writer: anytype, failure: *diagnostics.Failure) diagnostics.Error!void {
    failure.reset(.get);
    failure.sources.set(.output, .{ .label = .{ .borrowed = "caller writer" }, .kind = .stream });
    failure.sources.set(.request_snapshot, .{ .label = .{ .borrowed = region_str }, .kind = .memory });
    var context = Diagnostics{ .failure = failure, .output_source = .output, .lookup_scope = .retained_matches, .index_source = if (idx.source == .zfi) .zfi else .fai, .at = .{ .source = .request_snapshot, .request = 1 } };
    const diag = &context;
    const resolved = try resolveRegionText(idx, region_str, diag);
    var input_buffer: [FASTA_INPUT_BUFFER_BYTES]u8 = undefined;
    var output_buffer: [REGION_OUTPUT_BUFFER_BYTES]u8 = undefined;
    try emitRegion(resolved, fastaSource(idx), &input_buffer, &output_buffer, writer, diag);
    failure.succeed();
}

/// Strict fetch behind `IndexedReader.fetchInto`; `input_buffer` must not alias `dst`.
pub fn fetchSymbols(idx: *const LoadedIndex, ordinal: usize, start: u64, end: u64, dst: []u8, input_buffer: []u8, failure: *diagnostics.Failure) diagnostics.Error![]u8 {
    failure.reset(.get);
    failure.sources.set(.output, .{ .label = .{ .borrowed = "caller buffer" }, .kind = .memory });
    const index_source: diagnostics.SourceId = if (idx.source == .zfi) .zfi else .fai;
    var context = Diagnostics{ .failure = failure, .output_source = .output, .index_source = index_source, .at = .{ .source = index_source } };
    const diag = &context;
    if (ordinal >= idx.records.len) return diag.fail(.{ .@"request.ordinal_out_of_range" = .{ .ordinal = ordinal, .record_count = idx.records.len } }, .resolve, index_source, null, null);
    const rec = idx.records[ordinal];
    if (end < start) return diag.fail(.{ .@"request.end_before_start" = .{ .start = start, .end = end } }, .resolve, index_source, null, null);
    if (end > rec.seq_len) return diag.fail(.{ .@"request.end_after_sequence" = .{ .end = end, .sequence_length = rec.seq_len } }, .resolve, index_source, null, null);
    const count = end - start;
    if (count > dst.len) return diag.fail(.{ .@"request.output_buffer_too_small" = .{ .required_bytes = count, .available_bytes = dst.len } }, .resolve, .output, null, null);
    const symbols = dst[0..@intCast(count)];
    if (count == 0) {
        failure.succeed();
        return symbols;
    }

    const side_table = idx.sideTableLines(rec);
    if (!rec.isUniformWidth() and side_table.len == 0) return diag.consistency(.missing_side_table);
    if (side_table.len == 0 and (rec.line_bases == 0 or rec.line_bytes < rec.line_bases)) return diag.consistency(.traversal_position);
    const resolved = ResolvedRegion{
        .name = idx.recordName(ordinal) orelse "",
        .start = start + 1,
        .display_end = end,
        .is_full = false,
        .seq_offset = rec.seq_offset,
        .num_bases = count,
        .line_bases = rec.line_bases,
        .line_bytes = rec.line_bytes,
        .side_table = side_table,
        .record_index = ordinal,
        .orientation = .{},
        .annotate_transform = false,
    };
    const source = fastaSource(idx);
    const span = try regionSpan(resolved, source, diag);
    var output = RegionOutput{ .bytes = symbols, .diag = diag, .wrap = false };
    var unused: std.Io.Writer = .failing;
    try copyRegion(resolved, span, source, input_buffer, &output, &unused, diag);
    if (side_table.len == 0) {
        const dense = fasta.countBases(symbols);
        if (dense != count) return diag.mismatch(resolved, span, dense);
    }
    failure.succeed();
    return symbols;
}

pub fn runGet(backing_allocator: std.mem.Allocator, io: std.Io, fasta_path: []const u8, options: Options, failure: *diagnostics.Failure) diagnostics.Error!Result {
    failure.reset(.get);
    failure.sources.set(.stdout, .{ .label = .{ .borrowed = "stdout" }, .kind = .stream });
    var context = Diagnostics{ .failure = failure };
    var storage: GetBuffers = undefined;
    const buffers = &storage;
    const start_ns = if (options.measure_elapsed) monotonicNs(io) else 0;
    switch (options.source) {
        .positional => {
            context.at.source = .request_snapshot;
            failure.sources.set(.request_snapshot, .{ .label = .{ .borrowed = "positional requests" }, .kind = .memory });
        },
        .names, .bed => |path| {
            const stdin = std.mem.eql(u8, path, "-");
            context.at.source = if (stdin) .stdin else if (options.source == .bed) .bed else .names;
            failure.sources.set(context.at.source.?, .{ .label = .{ .borrowed = if (stdin) "stdin" else path }, .kind = if (stdin) .stream else .path });
        },
    }
    const totals = switch (options.source) {
        .positional => |regions| try runPositionalGet(backing_allocator, io, fasta_path, options, regions, &context, buffers),
        .names => |path| try runFileGet(backing_allocator, io, fasta_path, path, .{ .names = options.orientation }, options, &context, buffers),
        .bed => |path| try runFileGet(backing_allocator, io, fasta_path, path, .{ .bed = .{ .honor_strand = options.honor_strand, .global_orientation = options.orientation } }, options, &context, buffers),
    };
    failure.succeed();
    return .{ .region_count = totals.region_count, .total_bases = totals.total_bases, .elapsed_ns = if (options.measure_elapsed) monotonicNs(io) - start_ns else 0 };
}

const LoadedIndex = index_format.LoadedIndex;

const Diagnostics = struct {
    failure: *diagnostics.Failure,
    at: diagnostics.Reference = .{},
    index_source: diagnostics.SourceId = .zfi,
    lookup_scope: diagnostics.LookupScope = .complete_catalog,
    request_reader: ?*std.Io.File.Reader = null,
    output_source: diagnostics.SourceId = .stdout,

    fn output(self: *@This()) *diagnostics.OutputState {
        return if (self.output_source == .stdout) &self.failure.outputs.stdout else &self.failure.outputs.file;
    }

    fn fail(self: *@This(), details: diagnostics.Details, stage: diagnostics.Stage, source: diagnostics.SourceId, offset: ?u64, cause: ?anyerror) diagnostics.Error {
        @branchHint(.cold);
        var at = if (self.at.source == source) self.at else diagnostics.Reference{ .source = source, .request = self.at.request };
        if (offset) |value| at.location.offset = value;
        return self.failure.fail(.{
            .details = details,
            .stage = stage,
            .at = at,
            .related = if (self.at.source != null and self.at.source != source) self.at else null,
            .cause = if (cause) |err| .{ .zig = err } else null,
        });
    }

    fn oom(self: *@This(), role: diagnostics.AllocationRole) diagnostics.Error {
        @branchHint(.cold);
        return self.fail(.{ .@"resource.out_of_memory" = .{ .allocation_role = role } }, .allocate, self.at.source orelse .fasta, null, error.OutOfMemory);
    }

    fn missing(self: *@This(), name: []const u8) diagnostics.Error {
        @branchHint(.cold);
        return self.fail(.{ .@"request.sequence_not_found" = .{
            .lookup_scope = self.lookup_scope,
            .index_source = self.index_source,
            .identifier = self.failure.retain(name, 32),
        } }, .resolve, self.at.source orelse .request_snapshot, null, null);
    }

    fn consistency(self: *@This(), check: diagnostics.Check) diagnostics.Error {
        @branchHint(.cold);
        return self.fail(.{ .@"internal.consistency_failed" = .{ .component = .getter, .check = check } }, .resolve, self.index_source, null, null);
    }

    fn read(self: *@This(), offset: u64, wanted: u64, received: ?u64, cause: ?anyerror) diagnostics.Error {
        @branchHint(.cold);
        return self.fail(diagnostics.readFailure(wanted, received), .read, .fasta, offset, cause);
    }

    fn span(self: *@This(), end: u128, size: u64) diagnostics.Error {
        @branchHint(.cold);
        return self.fail(.{ .@"index.sequence_span_out_of_bounds" = .{
            .required_end = end,
            .source_bytes = size,
            .bound_source = .actual,
        } }, .resolve, self.index_source, null, null);
    }

    fn capacity(self: *@This(), role: diagnostics.CapacityRole, required: u128) diagnostics.Error {
        @branchHint(.cold);
        return self.fail(.{ .@"resource.capacity_limit_exceeded" = .{
            .capacity_role = role,
            .limit = std.math.maxInt(usize),
            .required = .{ .exact = required },
        } }, .allocate, self.at.source orelse .fasta, null, null);
    }

    fn add(self: *@This(), counter: diagnostics.Counter, current: u64, increment: u64) diagnostics.Error!u64 {
        return std.math.add(u64, current, increment) catch return self.countOverflow(counter, current, increment);
    }

    fn countOverflow(self: *@This(), counter: diagnostics.Counter, current: u64, increment: u64) diagnostics.Error {
        @branchHint(.cold);
        return self.fail(.{ .@"resource.count_overflow" = .{
            .counter = counter,
            .max_value = std.math.maxInt(u64),
            .current = current,
            .increment = increment,
        } }, if (counter == .source_bytes or counter == .physical_lines) .read else .resolve, self.at.source orelse .fasta, null, null);
    }

    fn mismatch(self: *@This(), resolved: ResolvedRegion, span_value: DiskSpan, observed: u64) diagnostics.Error {
        @branchHint(.cold);
        return self.fail(.{ .@"pair.region_length_mismatch" = .{
            .requested_symbols = resolved.num_bases,
            .observed_symbols = observed,
            .interval = .{ .start = resolved.start, .end = resolved.start + resolved.num_bases - 1 },
            .examined_span = .{ .start = span_value.start, .end = span_value.end },
        } }, .read, .fasta, span_value.start, null);
    }
};

const MAX_REQUEST_NAME_BYTES = std.math.maxInt(u16);
const REQUEST_LINE_READER_BUFFER_BYTES = MAX_REQUEST_NAME_BYTES + 4096;
const REQUEST_NAME_STACK_BYTES = 8 * 1024;
const ACTIVE_NAME_BYTES = 4 * 1024 * 1024;
const MAX_ACTIVE_NAME_BYTES = ACTIVE_NAME_BYTES + MAX_REQUEST_NAME_BYTES;
const NAMES_REQUEST_BATCH_SIZE = 65536;
const BED_REQUEST_BATCH_SIZE = 4096;
const TYPE_SAMPLE_BASES: u64 = 256;

const ParsedRequest = struct {
    region: Region,
    orientation: Orientation,

    fn lineNumber(_: ParsedRequest) ?u64 {
        return null;
    }

    fn parsed(self: ParsedRequest, _: []const u8) ParsedRequest {
        return self;
    }
};

const OwnedRequest = struct {
    line: u64,
    name_offset: u32,
    name_len: u16,
    start: u64,
    end: u64,
    orientation: Orientation,

    fn lineNumber(self: OwnedRequest) ?u64 {
        return self.line;
    }

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

const ParsedRange = struct {
    start: u64,
    end: ?u64,
};

fn parseRangeSuffix(suffix: []const u8) ?ParsedRange {
    const dash_pos = std.mem.findScalar(u8, suffix, '-') orelse return null;

    const start_str = suffix[0..dash_pos];
    const end_str = suffix[dash_pos + 1 ..];

    if (start_str.len == 0) return null;
    const start = std.fmt.parseInt(u64, start_str, 10) catch return null;

    if (end_str.len == 0) return .{ .start = start, .end = null };

    const end = std.fmt.parseInt(u64, end_str, 10) catch return null;
    return .{ .start = start, .end = end };
}

fn resolveRegionText(idx: *const LoadedIndex, region_str: []const u8, diag: *Diagnostics) diagnostics.Error!ResolvedRegion {
    const region = parseRegion(region_str);
    const rec_idx = idx.lookupName(region.name) orelse return diag.missing(region.name);
    return resolveParsedRequest(idx, .{ .region = region, .orientation = .{} }, rec_idx, false, diag);
}

fn resolveParsedRequest(
    idx: *const LoadedIndex,
    request: ParsedRequest,
    rec_idx: usize,
    annotate_transform: bool,
    diag: *Diagnostics,
) diagnostics.Error!ResolvedRegion {
    const rec = idx.records[rec_idx];
    const region = request.region;

    var start = region.start;
    var end = region.end orelse rec.seq_len;
    const display_end = end;

    if (region.is_full) {
        start = 1;
        end = rec.seq_len;
    }

    if (start < 1) {
        return diag.fail(.{ .@"request.start_before_one" = .{ .start = start, .coordinate_system = .one_based_inclusive } }, .resolve, diag.at.source orelse .request_snapshot, null, null);
    }
    if (start > rec.seq_len) {
        return diag.fail(.{ .@"request.start_after_sequence" = .{ .start = start, .sequence_length = rec.seq_len } }, .resolve, diag.at.source orelse .request_snapshot, null, null);
    }
    if (end < start) {
        return diag.fail(.{ .@"request.end_before_start" = .{ .start = start, .end = end } }, .resolve, diag.at.source orelse .request_snapshot, null, null);
    }

    if (end > rec.seq_len) {
        end = rec.seq_len;
    }

    const num_bases = end - start + 1;

    const side_table = idx.sideTableLines(rec);
    if (!rec.isUniformWidth() and side_table.len == 0) return diag.consistency(.missing_side_table);
    if (side_table.len == 0 and (rec.line_bases == 0 or rec.line_bytes < rec.line_bases)) return diag.consistency(.traversal_position);

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

fn detectRegionType(source: FastaSource, resolved: ResolvedRegion, diag: *Diagnostics, buffers: *GetBuffers) diagnostics.Error!fasta.SequenceType {
    const span = try diskSpanForRegion(resolved, diag);
    if (span.end > source.size) return diag.span(span.end, source.size);
    if (span.start >= span.end) return diag.consistency(.traversal_position);

    var counts = [_]u64{0} ** 256;
    var total: u64 = 0;
    var skip = span.leading_bases;
    const input = &buffers.input;
    var offset = span.start;

    while (offset < span.end and total < resolved.num_bases) {
        const wanted: usize = @intCast(@min(@as(u64, input.len), span.end - offset));
        const got = std.Io.File.readPositionalAll(source.file, source.io, input[0..wanted], offset) catch |err| return diag.read(offset, wanted, null, err);
        if (got != wanted) return diag.read(offset, wanted, got, null);
        for (input[0..got]) |byte| {
            if (!fasta.isBase(byte)) continue;
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

    if (total != resolved.num_bases) return diag.mismatch(resolved, span, total);
    return fasta.detectType(&counts, total);
}

const FastaSource = struct {
    io: std.Io,
    file: std.Io.File,
    size: u64,
};

fn fastaSource(idx: *const LoadedIndex) FastaSource {
    return .{ .io = idx.io, .file = idx.fasta_file, .size = idx.fasta_size };
}

const FASTA_INPUT_BUFFER_BYTES: usize = 256 * 1024;
const REGION_OUTPUT_BUFFER_BYTES: usize = 64 * 1024;

const GetBuffers = struct {
    input: [FASTA_INPUT_BUFFER_BYTES]u8,
    output: [REGION_OUTPUT_BUFFER_BYTES]u8,
};

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

fn uniformByteOffset(resolved: ResolvedRegion, base_index: u64, diag: *Diagnostics) diagnostics.Error!u64 {
    const line_number = base_index / resolved.line_bases;
    const column = base_index % resolved.line_bases;
    const line_delta = std.math.mul(u64, line_number, resolved.line_bytes) catch {
        return diag.consistency(.traversal_position);
    };
    const line_start = std.math.add(u64, resolved.seq_offset, line_delta) catch {
        return diag.consistency(.traversal_position);
    };
    return std.math.add(u64, line_start, column) catch {
        return diag.consistency(.traversal_position);
    };
}

fn diskSpanForRegion(resolved: ResolvedRegion, diag: *Diagnostics) diagnostics.Error!DiskSpan {
    const first_base = resolved.start - 1;
    const last_base = first_base + resolved.num_bases - 1;

    if (resolved.side_table.len == 0) {
        const start = try uniformByteOffset(resolved, first_base, diag);
        const last = try uniformByteOffset(resolved, last_base, diag);
        return .{
            .start = start,
            .end = std.math.add(u64, last, 1) catch return diag.consistency(.traversal_position),
            .leading_bases = 0,
            .trailing_bases = 0,
        };
    }

    const first_line_idx = sideTableLineForBase(resolved.side_table, first_base);
    const last_line_idx = sideTableLineForBase(resolved.side_table, last_base);
    const first_line = resolved.side_table[first_line_idx];
    const last_line = resolved.side_table[last_line_idx];
    const last_line_end = std.math.add(u64, last_line.byte_offset, last_line.line_bytes) catch {
        return diag.consistency(.side_table_state);
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
    diag: *Diagnostics,
    wrap: bool = true,

    fn flush(self: *RegionOutput, writer: anytype) diagnostics.Error!void {
        if (self.len == 0) return;
        writer.writeAll(self.bytes[0..self.len]) catch |err| return self.diag.fail(.{ .@"io.write_failed" = {} }, .write, self.diag.output_source, null, err);
        self.len = 0;
    }

    fn appendForward(self: *RegionOutput, src: []const u8, do_complement: bool, writer: anytype) diagnostics.Error!void {
        var start: usize = 0;
        while (start < src.len) {
            const chunk = if (self.wrap) @min(src.len - start, 60 - self.line_pos) else src.len - start;
            if (self.bytes.len - self.len < chunk + @intFromBool(self.wrap)) try self.flush(writer);
            const dst = self.bytes[self.len .. self.len + chunk];
            if (do_complement) {
                fasta.complementInto(dst, src[start .. start + chunk]);
            } else {
                @memcpy(dst, src[start .. start + chunk]);
            }
            self.len += chunk;
            self.bases_written += chunk;
            self.line_pos += chunk;
            start += chunk;
            if (self.wrap and self.line_pos == 60) {
                self.bytes[self.len] = '\n';
                self.len += 1;
                self.line_pos = 0;
            }
        }
    }

    fn appendReverse(self: *RegionOutput, src: []const u8, do_complement: bool, writer: anytype) diagnostics.Error!void {
        var end = src.len;
        while (end > 0) {
            const chunk = if (self.wrap) @min(end, 60 - self.line_pos) else end;
            if (self.bytes.len - self.len < chunk + @intFromBool(self.wrap)) try self.flush(writer);
            var j: usize = 0;
            while (j < chunk) : (j += 1) {
                const byte = src[end - 1 - j];
                self.bytes[self.len + j] = if (do_complement) fasta.complement(byte) else byte;
            }
            self.len += chunk;
            self.bases_written += chunk;
            self.line_pos += chunk;
            end -= chunk;
            if (self.wrap and self.line_pos == 60) {
                self.bytes[self.len] = '\n';
                self.len += 1;
                self.line_pos = 0;
            }
        }
    }

    fn finish(self: *RegionOutput, writer: anytype) diagnostics.Error!void {
        if (self.line_pos != 0) {
            if (self.len == self.bytes.len) try self.flush(writer);
            self.bytes[self.len] = '\n';
            self.len += 1;
        }
        try self.flush(writer);
    }
};

fn writeRegionHeader(resolved: ResolvedRegion, writer: anytype, diag: *Diagnostics) diagnostics.Error!void {
    diag.output().state = .streaming;
    diag.output().may_be_partial = true;
    const annotation = headerAnnotation(resolved.orientation, resolved.annotate_transform);
    if (resolved.is_full) {
        writer.print(">{s}{s}\n", .{ resolved.name, annotation }) catch |err| return diag.fail(.{ .@"io.write_failed" = {} }, .write, diag.output_source, null, err);
    } else {
        writer.print(">{s}:{d}-{d}{s}\n", .{ resolved.name, resolved.start, resolved.display_end, annotation }) catch |err| return diag.fail(.{ .@"io.write_failed" = {} }, .write, diag.output_source, null, err);
    }
}

fn consumeForwardBytes(
    resolved: ResolvedRegion,
    bytes: []const u8,
    skip: *u64,
    output: *RegionOutput,
    writer: anytype,
) diagnostics.Error!void {
    var pos: usize = 0;
    while (pos < bytes.len and output.bases_written < resolved.num_bases) {
        while (pos < bytes.len and !fasta.isBase(bytes[pos])) pos += 1;
        const run_start = pos;
        while (pos < bytes.len and fasta.isBase(bytes[pos])) pos += 1;
        var run = bytes[run_start..pos];
        if (skip.* != 0) {
            const skipped: usize = @intCast(@min(skip.*, @as(u64, run.len)));
            skip.* -= @intCast(skipped);
            run = run[skipped..];
        }
        const remaining: usize = @intCast(resolved.num_bases - output.bases_written);
        run = run[0..@min(run.len, remaining)];
        try output.appendForward(run, resolved.orientation.complement, writer);
    }
}

fn consumeReverseBytes(
    resolved: ResolvedRegion,
    bytes: []const u8,
    skip: *u64,
    output: *RegionOutput,
    writer: anytype,
) diagnostics.Error!void {
    var i = bytes.len;
    while (i > 0 and output.bases_written < resolved.num_bases) {
        while (i > 0 and !fasta.isBase(bytes[i - 1])) i -= 1;
        const run_end = i;
        while (i > 0 and fasta.isBase(bytes[i - 1])) i -= 1;
        var run = bytes[i..run_end];
        if (skip.* != 0) {
            const skipped: usize = @intCast(@min(skip.*, @as(u64, run.len)));
            skip.* -= @intCast(skipped);
            run = run[0 .. run.len - skipped];
        }
        const remaining: usize = @intCast(resolved.num_bases - output.bases_written);
        if (run.len > remaining) run = run[run.len - remaining ..];
        try output.appendReverse(run, resolved.orientation.complement, writer);
    }
}

fn emitUniformFromSpan(
    resolved: ResolvedRegion,
    span_start: u64,
    bytes: []const u8,
    output: *RegionOutput,
    writer: anytype,
    diag: *Diagnostics,
) diagnostics.Error!void {
    var remaining = resolved.num_bases;
    const base_index = if (resolved.orientation.reverse)
        resolved.start - 1 + resolved.num_bases - 1
    else
        resolved.start - 1;
    var column = base_index % resolved.line_bases;
    var line_start = try uniformByteOffset(resolved, base_index, diag) - column;

    while (remaining != 0) {
        const take_u64 = if (resolved.orientation.reverse)
            @min(column + 1, remaining)
        else
            @min(resolved.line_bases - column, remaining);
        const take: usize = @intCast(take_u64);
        const source_start = if (resolved.orientation.reverse)
            line_start + column + 1 - take_u64
        else
            line_start + column;
        if (source_start < span_start) return diag.consistency(.traversal_position);
        const relative_start: usize = @intCast(source_start - span_start);
        if (relative_start > bytes.len or take > bytes.len - relative_start) {
            return diag.consistency(.traversal_position);
        }
        const src = bytes[relative_start .. relative_start + take];
        if (resolved.orientation.reverse) {
            try output.appendReverse(src, resolved.orientation.complement, writer);
        } else {
            try output.appendForward(src, resolved.orientation.complement, writer);
        }

        remaining -= take_u64;
        if (remaining == 0) break;
        if (resolved.orientation.reverse) {
            line_start = std.math.sub(u64, line_start, resolved.line_bytes) catch {
                return diag.consistency(.traversal_position);
            };
            column = resolved.line_bases - 1;
        } else {
            line_start = std.math.add(u64, line_start, resolved.line_bytes) catch {
                return diag.consistency(.traversal_position);
            };
            column = 0;
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
    diag: *Diagnostics,
) diagnostics.Error!void {
    const chunk_end = std.math.add(u64, chunk_start, bytes.len) catch {
        return diag.consistency(.traversal_position);
    };
    while (output.bases_written < resolved.num_bases) {
        const column = base_index.* % resolved.line_bases;
        const source_start = try uniformByteOffset(resolved, base_index.*, diag);
        if (source_start >= chunk_end) return;
        if (source_start < chunk_start) return diag.consistency(.traversal_position);
        const relative_start: usize = @intCast(source_start - chunk_start);
        const available = @min(
            resolved.line_bases - column,
            resolved.num_bases - output.bases_written,
        );
        const take: usize = @intCast(@min(available, chunk_end - source_start));
        try output.appendForward(bytes[relative_start .. relative_start + take], resolved.orientation.complement, writer);
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
    diag: *Diagnostics,
) diagnostics.Error!void {
    const chunk_end = std.math.add(u64, chunk_start, bytes.len) catch {
        return diag.consistency(.traversal_position);
    };
    while (output.bases_written < resolved.num_bases) {
        const column = base_index.* % resolved.line_bases;
        const source_end = std.math.add(u64, try uniformByteOffset(resolved, base_index.*, diag), 1) catch {
            return diag.consistency(.traversal_position);
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
        try output.appendReverse(bytes[relative_start .. relative_start + take], resolved.orientation.complement, writer);
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
    diag: *Diagnostics,
) diagnostics.Error!void {
    try writeRegionHeader(resolved, writer, diag);
    const relative_start = std.math.cast(usize, span.start - span_start) orelse {
        return diag.capacity(.disk_span, span.end - span_start);
    };
    const relative_end = std.math.cast(usize, span.end - span_start) orelse {
        return diag.capacity(.disk_span, span.end - span_start);
    };
    if (relative_end > bytes.len or relative_start > relative_end) {
        return diag.consistency(.traversal_position);
    }
    var output = RegionOutput{ .bytes = output_buffer, .diag = diag };
    if (resolved.side_table.len == 0) {
        try emitUniformFromSpan(resolved, span_start, bytes, &output, writer, diag);
    } else {
        const region_bytes = bytes[relative_start..relative_end];
        var skip = if (resolved.orientation.reverse) span.trailing_bases else span.leading_bases;
        if (resolved.orientation.reverse) {
            try consumeReverseBytes(resolved, region_bytes, &skip, &output, writer);
        } else {
            try consumeForwardBytes(resolved, region_bytes, &skip, &output, writer);
        }
    }
    if (output.bases_written != resolved.num_bases) return diag.mismatch(resolved, span, output.bases_written);
    try output.finish(writer);
    diag.output().records = (diag.output().records orelse 0) + 1;
}

fn emitRegion(
    resolved: ResolvedRegion,
    source: FastaSource,
    input_buffer: []u8,
    output_buffer: []u8,
    writer: anytype,
    diag: *Diagnostics,
) diagnostics.Error!void {
    const span = try regionSpan(resolved, source, diag);
    try writeRegionHeader(resolved, writer, diag);
    var output = RegionOutput{ .bytes = output_buffer, .diag = diag };
    try copyRegion(resolved, span, source, input_buffer, &output, writer, diag);
    try output.finish(writer);
    diag.output().records = (diag.output().records orelse 0) + 1;
}

fn regionSpan(resolved: ResolvedRegion, source: FastaSource, diag: *Diagnostics) diagnostics.Error!DiskSpan {
    const span = try diskSpanForRegion(resolved, diag);
    if (span.end > source.size) return diag.span(span.end, source.size);
    if (span.start >= span.end) return diag.consistency(.traversal_position);
    return span;
}

fn copyRegion(
    resolved: ResolvedRegion,
    span: DiskSpan,
    source: FastaSource,
    input_buffer: []u8,
    output: *RegionOutput,
    writer: anytype,
    diag: *Diagnostics,
) diagnostics.Error!void {
    var skip = if (resolved.orientation.reverse) span.trailing_bases else span.leading_bases;
    var base_index = if (resolved.orientation.reverse)
        resolved.start - 1 + resolved.num_bases - 1
    else
        resolved.start - 1;

    if (!resolved.orientation.reverse) {
        var offset = span.start;
        while (offset < span.end and output.bases_written < resolved.num_bases) {
            const wanted: usize = @intCast(@min(@as(u64, input_buffer.len), span.end - offset));
            const got = std.Io.File.readPositionalAll(source.file, source.io, input_buffer[0..wanted], offset) catch |err| return diag.read(offset, wanted, null, err);
            if (got != wanted) return diag.read(offset, wanted, got, null);
            if (resolved.side_table.len == 0) {
                try consumeUniformForwardChunk(resolved, offset, input_buffer[0..got], &base_index, output, writer, diag);
            } else {
                try consumeForwardBytes(resolved, input_buffer[0..got], &skip, output, writer);
            }
            offset += got;
        }
    } else {
        var end = span.end;
        while (end > span.start and output.bases_written < resolved.num_bases) {
            const wanted: usize = @intCast(@min(@as(u64, input_buffer.len), end - span.start));
            const offset = end - wanted;
            const got = std.Io.File.readPositionalAll(source.file, source.io, input_buffer[0..wanted], offset) catch |err| return diag.read(offset, wanted, null, err);
            if (got != wanted) return diag.read(offset, wanted, got, null);
            if (resolved.side_table.len == 0) {
                try consumeUniformReverseChunk(resolved, offset, input_buffer[0..got], &base_index, output, writer, diag);
            } else {
                try consumeReverseBytes(resolved, input_buffer[0..got], &skip, output, writer);
            }
            end = offset;
        }
    }

    if (output.bases_written != resolved.num_bases) return diag.mismatch(resolved, span, output.bases_written);
}

fn emitResolvedBatch(
    resolved: []const ResolvedRegion,
    source: FastaSource,
    writer: anytype,
    diag: *Diagnostics,
    request_entries: anytype,
    request_base: usize,
    total_base: u64,
    buffers: *GetBuffers,
) diagnostics.Error!u64 {
    var total_bases = total_base;
    const request_source = diag.at.source;
    const shared = &buffers.input;
    const output_buffer = &buffers.output;
    var i: usize = 0;

    while (i < resolved.len) {
        diag.at = requestReference(request_source, request_entries[i].lineNumber(), request_base + i + 1);
        const first_span = try diskSpanForRegion(resolved[i], diag);
        var shared_start = first_span.start;
        var shared_end = first_span.end;
        const shared_record = resolved[i].record_index;
        var shared_first_base = resolved[i].start - 1;
        var shared_last_base = shared_first_base + resolved[i].num_bases - 1;
        var group_end = i + 1;

        while (group_end < resolved.len) : (group_end += 1) {
            diag.at = requestReference(request_source, request_entries[group_end].lineNumber(), request_base + group_end + 1);
            const next = try diskSpanForRegion(resolved[group_end], diag);
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

        diag.at = requestReference(request_source, request_entries[i].lineNumber(), request_base + i + 1);
        if (shared_end - shared_start <= shared.len) {
            if (shared_end > source.size) return diag.span(shared_end, source.size);
            const wanted: usize = @intCast(shared_end - shared_start);
            const got = std.Io.File.readPositionalAll(source.file, source.io, shared[0..wanted], shared_start) catch |err| return diag.read(shared_start, wanted, null, err);
            if (got != wanted) return diag.read(shared_start, wanted, got, null);
            for (resolved[i..group_end], i..) |region, ordinal| {
                diag.at = requestReference(request_source, request_entries[ordinal].lineNumber(), request_base + ordinal + 1);
                const next_total = try diag.add(.total_symbols, total_bases, region.num_bases);
                const span = try diskSpanForRegion(region, diag);
                try emitRegionFromSpan(region, span, shared_start, shared[0..got], output_buffer, writer, diag);
                total_bases = next_total;
            }
            i = group_end;
            continue;
        }

        const next_total = try diag.add(.total_symbols, total_bases, resolved[i].num_bases);
        try emitRegion(resolved[i], source, shared, output_buffer, writer, diag);
        total_bases = next_total;
        i += 1;
    }

    return total_bases - total_base;
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

fn ownRequestName(allocator: std.mem.Allocator, workspace: *RequestWorkspace, name: []const u8, diag: *Diagnostics) diagnostics.Error!OwnedName {
    if (workspace.requests.getLastOrNull()) |last| {
        const offset: usize = last.name_offset;
        if (std.mem.eql(u8, workspace.names.items[offset .. offset + last.name_len], name)) return .{ .offset = last.name_offset, .len = last.name_len };
    }
    const offset = workspace.names.items.len;
    const required_capacity = offset + name.len;
    if (required_capacity > workspace.names.capacity) {
        const doubled = @max(workspace.names.capacity, 4096) * 2;
        const new_capacity = @min(MAX_ACTIVE_NAME_BYTES, @max(required_capacity, doubled));
        workspace.names.ensureTotalCapacityPrecise(allocator, new_capacity) catch return diag.oom(.request_names);
    }
    workspace.names.appendSliceAssumeCapacity(name);
    return .{ .offset = @intCast(offset), .len = @intCast(name.len) };
}

fn appendOwnedRequest(allocator: std.mem.Allocator, workspace: *RequestWorkspace, name: []const u8, start: u64, end: u64, orientation: Orientation, diag: *Diagnostics) diagnostics.Error!void {
    const owned_name = try ownRequestName(allocator, workspace, name, diag);
    workspace.requests.append(allocator, .{
        .line = diag.at.location.line.?,
        .name_offset = owned_name.offset,
        .name_len = owned_name.len,
        .start = start,
        .end = end,
        .orientation = orientation,
    }) catch return diag.oom(.requests);
}

fn appendBedLineRequest(allocator: std.mem.Allocator, workspace: *RequestWorkspace, line: []const u8, line_number: usize, honor_strand: bool, global_orientation: Orientation, diag: *Diagnostics) diagnostics.Error!void {
    diag.at.location.line = line_number;
    const parsed = try bed_parser.parseBedLine(line, honor_strand, diag.at, diag.failure);
    switch (parsed) {
        .skip => {},
        .region => |region| {
            if (region.chrom.len > MAX_REQUEST_NAME_BYTES) return diag.fail(.{ .@"request.identifier_too_long" = .{
                .identifier_bytes = region.chrom.len,
                .limit_bytes = MAX_REQUEST_NAME_BYTES,
            } }, .parse, diag.at.source orelse .bed, null, null);
            try appendOwnedRequest(allocator, workspace, region.chrom, region.start_1based, region.end_1based, (if (honor_strand and region.strand == .reverse) Orientation.reverseComplement() else Orientation{}).compose(global_orientation), diag);
        },
    }
}

fn appendNamesLine(allocator: std.mem.Allocator, workspace: *RequestWorkspace, line: []const u8, line_number: usize, orientation: Orientation, diag: *Diagnostics) diagnostics.Error!void {
    diag.at.location.line = line_number;
    const name = if (line.len > 0 and line[line.len - 1] == '\r') line[0 .. line.len - 1] else line;
    if (name.len == 0 or name[0] == '#') return;
    if (name.len > MAX_REQUEST_NAME_BYTES) return diag.fail(.{ .@"request.identifier_too_long" = .{
        .identifier_bytes = name.len,
        .limit_bytes = MAX_REQUEST_NAME_BYTES,
    } }, .parse, diag.at.source orelse .names, null, null);
    try appendOwnedRequest(allocator, workspace, name, 1, 0, orientation, diag);
}

fn appendCliRequests(allocator: std.mem.Allocator, requests: *std.ArrayList(ParsedRequest), region_strs: []const []const u8, orientation: Orientation, diag: *Diagnostics) diagnostics.Error!void {
    for (region_strs, 0..) |region_str, i| {
        diag.at.request = try diag.add(.requests, i, 1);
        requests.append(allocator, .{ .region = parseRegion(region_str), .orientation = orientation }) catch return diag.oom(.requests);
    }
}

fn requestReference(source: ?diagnostics.SourceId, line: ?u64, ordinal: u64) diagnostics.Reference {
    return .{ .source = source, .location = .{ .line = line }, .request = ordinal };
}

fn processParsedRequests(
    allocator: std.mem.Allocator,
    idx: *LoadedIndex,
    request_entries: anytype,
    active_names: []const u8,
    annotate_transform: bool,
    writer: anytype,
    source: FastaSource,
    diag: *Diagnostics,
    request_base: usize,
    total_base: u64,
    buffers: *GetBuffers,
) diagnostics.Error!BatchStats {
    if (request_entries.len == 0) return .{};
    _ = try diag.add(.requests, request_base, request_entries.len);
    const resolved = allocator.alloc(ResolvedRegion, request_entries.len) catch return diag.oom(.resolved_requests);
    defer allocator.free(resolved);
    var last_name: ?[]const u8 = null;
    var last_rec_idx: usize = 0;
    var last_typed_rec_idx: ?usize = null;
    var last_type: fasta.SequenceType = undefined;
    var last_sample: diagnostics.Sample = undefined;
    const request_source = diag.at.source;
    for (request_entries, 0..) |entry, i| {
        diag.at = requestReference(request_source, entry.lineNumber(), request_base + i + 1);
        const request = entry.parsed(active_names);
        const rec_idx = if (last_name) |name|
            if (std.mem.eql(u8, name, request.region.name)) last_rec_idx else idx.lookupName(request.region.name) orelse return diag.missing(request.region.name)
        else
            idx.lookupName(request.region.name) orelse return diag.missing(request.region.name);
        last_name = request.region.name;
        last_rec_idx = rec_idx;
        resolved[i] = try resolveParsedRequest(idx, request, rec_idx, annotate_transform, diag);
        if (request.orientation.complement) {
            const rec_type = if (last_typed_rec_idx == rec_idx) last_type else blk: {
                var sample = resolved[i];
                sample.num_bases = @min(sample.num_bases, TYPE_SAMPLE_BASES);
                const detected = try detectRegionType(source, sample, diag, buffers);
                last_typed_rec_idx = rec_idx;
                last_type = detected;
                last_sample = .{ .symbols = @intCast(sample.num_bases), .rule = .nucleotide_above_90_percent, .interval = .{ .start = sample.start, .end = sample.start + sample.num_bases - 1 }, .request = diag.at.request };
                break :blk detected;
            };
            if (rec_type == .protein) return diag.fail(.{ .@"transform.complement_on_protein" = .{
                .transform = .{ .reverse = request.orientation.reverse, .complement = request.orientation.complement },
                .classification = .protein,
                .sample = last_sample,
                .identifier = diag.failure.retain(request.region.name, 32),
            } }, .transform, request_source orelse .request_snapshot, null, null);
        }
    }
    return .{ .region_count = request_entries.len, .total_bases = try emitResolvedBatch(resolved, source, writer, diag, request_entries, request_base, total_base, buffers) };
}

const StreamRequestSource = union(enum) {
    names: Orientation,
    bed: struct {
        honor_strand: bool,
        global_orientation: Orientation,
    },
};

const RequestedNames = struct {
    arena: std.heap.ArenaAllocator,
    stack_allocator: std.heap.StackFallbackAllocator(REQUEST_NAME_STACK_BYTES) = undefined,
    names: std.ArrayList([]const u8) = .empty,
    set: std.StringHashMapUnmanaged(void) = .empty,

    fn init(allocator: std.mem.Allocator) RequestedNames {
        return .{ .arena = std.heap.ArenaAllocator.init(allocator) };
    }

    fn deinit(self: *RequestedNames) void {
        self.arena.deinit();
    }

    fn collect(self: *RequestedNames, reader: *std.Io.Reader, source: StreamRequestSource, stable_input: bool) !bool {
        self.stack_allocator = std.heap.stackFallback(REQUEST_NAME_STACK_BYTES, self.arena.allocator());
        const allocator = self.stack_allocator.get();
        const initial_capacity = if (stable_input) @min(128, 1 + std.mem.count(u8, reader.buffered(), "\n")) else 0;
        var previous_name: ?[]const u8 = null;
        while (true) {
            const available = reader.buffered();
            const line = if (std.mem.findScalar(u8, available, '\n')) |end| buffered: {
                reader.toss(end + 1);
                break :buffered available[0..end];
            } else (reader.takeDelimiter('\n') catch return false) orelse return true;
            const name = switch (source) {
                .names => blk: {
                    const name = if (line.len > 0 and line[line.len - 1] == '\r') line[0 .. line.len - 1] else line;
                    if (name.len == 0 or name[0] == '#') continue;
                    break :blk name;
                },
                .bed => bed_parser.requestChrom(line) orelse continue,
            };
            if (name.len > MAX_REQUEST_NAME_BYTES) return false;
            if (previous_name) |previous| {
                if (std.mem.eql(u8, previous, name)) continue;
            }
            if (self.set.getKey(name)) |existing| {
                previous_name = existing;
                continue;
            }
            if (self.names.items.len == 0 and initial_capacity != 0) {
                try self.set.ensureTotalCapacity(allocator, @intCast(initial_capacity));
                try self.names.ensureTotalCapacityPrecise(allocator, initial_capacity);
            }
            const owned = if (stable_input) name else try allocator.dupe(u8, name);
            try self.set.put(allocator, owned, {});
            try self.names.append(allocator, owned);
            previous_name = owned;
        }
    }
};

const SeekableRequest = struct {
    file: std.Io.File,
    size: u64,
    mtime_ns: i128,
};

fn openSeekableRequest(io: std.Io, path: []const u8) ?SeekableRequest {
    if (std.mem.eql(u8, path, "-")) return null;
    const stat = std.Io.Dir.cwd().statFile(io, path, .{}) catch return null;
    if (stat.kind != .file) return null;
    const file = std.Io.Dir.cwd().openFile(io, path, .{}) catch return null;
    const opened_stat = file.stat(io) catch {
        file.close(io);
        return null;
    };
    if (opened_stat.kind != .file) {
        file.close(io);
        return null;
    }
    return .{ .file = file, .size = opened_stat.size, .mtime_ns = opened_stat.mtime.toNanoseconds() };
}

fn checkRequestSource(io: std.Io, opened: SeekableRequest, diag: *Diagnostics) diagnostics.Error!void {
    const stat = opened.file.stat(io) catch |err| return diag.fail(.{ .@"io.stat_failed" = {} }, .compare_source, diag.at.source orelse .names, null, err);
    const mtime = stat.mtime.toNanoseconds();
    if (stat.size != opened.size or mtime != opened.mtime_ns) return diag.fail(.{ .@"input.changed" = .{
        .changed_fields = .{ .size = stat.size != opened.size, .mtime = mtime != opened.mtime_ns },
        .before = .{ .bytes = opened.size, .mtime_ns = opened.mtime_ns },
        .after = .{ .bytes = stat.size, .mtime_ns = mtime },
    } }, .compare_source, diag.at.source orelse .names, null, null);
}

fn processRequestReader(
    allocator: std.mem.Allocator,
    idx: *LoadedIndex,
    reader: *std.Io.Reader,
    source: StreamRequestSource,
    name_reservation: usize,
    annotate_transform: bool,
    writer: anytype,
    fasta_source: FastaSource,
    diag: *Diagnostics,
    replay: ?SeekableRequest,
    buffers: *GetBuffers,
) diagnostics.Error!BatchStats {
    var total = BatchStats{};
    var line_number: usize = 0;
    var line_offset: u64 = 0;
    var reached_end = false;
    var workspace = RequestWorkspace.init(allocator);
    defer workspace.deinit(allocator);
    const request_source = diag.at.source;
    const request_limit: usize = switch (source) {
        .names => NAMES_REQUEST_BATCH_SIZE,
        .bed => BED_REQUEST_BATCH_SIZE,
    };
    workspace.requests.ensureTotalCapacityPrecise(allocator, request_limit) catch return diag.oom(.requests);
    if (name_reservation > 0) workspace.names.ensureTotalCapacityPrecise(allocator, name_reservation) catch return diag.oom(.request_names);
    while (true) {
        const requests_fit = total.region_count <= std.math.maxInt(usize) - request_limit;
        diag.at = .{ .source = request_source, .location = .{ .offset = line_offset, .line = try diag.add(.physical_lines, line_number, 1) } };
        if (replay) |opened| try checkRequestSource(idx.io, opened, diag);
        while (workspace.requests.items.len < request_limit and workspace.names.items.len < ACTIVE_NAME_BYTES) {
            const next_line = try diag.add(.physical_lines, line_number, 1);
            diag.at = .{ .source = request_source, .location = .{ .offset = line_offset, .line = next_line, .column = 1 }, .request = if (requests_fit) total.region_count + workspace.requests.items.len + 1 else try diag.add(.requests, total.region_count, workspace.requests.items.len + 1) };
            const bytes = reader.takeDelimiterInclusive('\n') catch |err| switch (err) {
                error.EndOfStream => remaining: {
                    const buffered = reader.buffered();
                    if (buffered.len == 0) {
                        reached_end = true;
                        break;
                    }
                    reader.toss(buffered.len);
                    break :remaining buffered;
                },
                error.ReadFailed => {
                    const cause: anyerror = if (diag.request_reader) |file_reader| file_reader.err orelse err else err;
                    const offset: ?u64 = if (diag.request_reader) |file_reader| file_reader.pos else null;
                    diag.at.location.column = null;
                    return diag.fail(.{ .@"io.read_failed" = .{} }, .read, request_source orelse .names, offset, cause);
                },
                error.StreamTooLong => return diag.fail(.{ .@"resource.line_buffer_exceeded" = .{
                    .buffer_bytes = reader.buffer.len,
                    .observed_bytes = .{ .at_least = reader.buffer.len },
                } }, .read, request_source orelse .names, null, null),
            };
            const line = if (bytes.len > 0 and bytes[bytes.len - 1] == '\n') bytes[0 .. bytes.len - 1] else bytes;
            line_number = @intCast(next_line);
            switch (source) {
                .names => |orientation| try appendNamesLine(allocator, &workspace, line, line_number, orientation, diag),
                .bed => |bed| try appendBedLineRequest(allocator, &workspace, line, line_number, bed.honor_strand, bed.global_orientation, diag),
            }
            line_offset = try diag.add(.source_bytes, line_offset, bytes.len);
        }
        if (workspace.requests.items.len == 0) break;
        const batch = try processParsedRequests(workspace.scratch.allocator(), idx, workspace.requests.items, workspace.names.items, annotate_transform, writer, fasta_source, diag, total.region_count, total.total_bases, buffers);
        total.region_count = if (requests_fit) total.region_count + batch.region_count else @intCast(try diag.add(.requests, total.region_count, batch.region_count));
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
    opened_file: ?SeekableRequest,
    prepared_reader: ?*std.Io.Reader,
    prepared_file_reader: ?*std.Io.File.Reader,
    source: StreamRequestSource,
    annotate_transform: bool,
    writer: anytype,
    fasta_source: FastaSource,
    diag: *Diagnostics,
    buffers: *GetBuffers,
) diagnostics.Error!BatchStats {
    if (std.mem.eql(u8, path, "-")) {
        var buffer: [REQUEST_LINE_READER_BUFFER_BYTES]u8 = undefined;
        var reader = std.Io.File.stdin().reader(io, &buffer);
        diag.request_reader = &reader;
        return processRequestReader(allocator, idx, &reader.interface, source, 0, annotate_transform, writer, fasta_source, diag, null, buffers);
    }
    const file = if (opened_file) |opened| opened.file else std.Io.Dir.cwd().openFile(io, path, .{}) catch |err| return diag.fail(diagnostics.openFailure(err), .open, diag.at.source orelse .names, null, err);
    defer if (opened_file == null) file.close(io);
    const name_reservation: usize = switch (source) {
        .names => blk: {
            if (opened_file) |opened| break :blk @intCast(@min(opened.size, MAX_ACTIVE_NAME_BYTES));
            const stat = file.stat(io) catch break :blk 0;
            break :blk @intCast(@min(stat.size, MAX_ACTIVE_NAME_BYTES));
        },
        .bed => 0,
    };
    var buffer: [REQUEST_LINE_READER_BUFFER_BYTES]u8 = undefined;
    var reader = file.reader(io, &buffer);
    diag.request_reader = if (prepared_reader != null) prepared_file_reader else &reader;
    return processRequestReader(allocator, idx, prepared_reader orelse &reader.interface, source, name_reservation, annotate_transform, writer, fasta_source, diag, if (prepared_reader != null and prepared_file_reader == null) null else opened_file, buffers);
}

fn runPositionalGet(backing_allocator: std.mem.Allocator, io: std.Io, fasta_path: []const u8, options: Options, regions: []const []const u8, diag: *Diagnostics, buffers: *GetBuffers) diagnostics.Error!BatchStats {
    var input_arena = std.heap.ArenaAllocator.init(backing_allocator);
    defer input_arena.deinit();
    const allocator = input_arena.allocator();
    var requests = std.ArrayList(ParsedRequest).empty;
    try appendCliRequests(allocator, &requests, regions, options.orientation, diag);
    const names = allocator.alloc([]const u8, requests.items.len) catch return diag.oom(.request_names);
    for (requests.items, names) |request, *name| name.* = request.region.name;
    var idx = try index_format.loadIndex(backing_allocator, io, fasta_path, .automatic, .{ .positional = names }, diag.failure);
    defer idx.deinit();
    diag.index_source = if (idx.source == .zfi) .zfi else .fai;
    return emitGet(allocator, io, &idx, options, requests.items, null, null, null, diag, buffers);
}

fn runFileGet(backing_allocator: std.mem.Allocator, io: std.Io, fasta_path: []const u8, request_path: []const u8, source: StreamRequestSource, options: Options, diag: *Diagnostics, buffers: *GetBuffers) diagnostics.Error!BatchStats {
    const request_file = openSeekableRequest(io, request_path);
    defer if (request_file) |opened| opened.file.close(io);
    var request_buffer: [REQUEST_LINE_READER_BUFFER_BYTES]u8 = undefined;
    const buffered_requests = if (request_file) |opened| buffered: {
        if (opened.size == 0 or opened.size > request_buffer.len) break :buffered null;
        const size: usize = @intCast(opened.size);
        const got = std.Io.File.readPositionalAll(opened.file, io, request_buffer[0..size], 0) catch break :buffered null;
        if (got != size) break :buffered null;
        break :buffered request_buffer[0..size];
    } else null;
    var small_reader = std.Io.Reader.fixed(buffered_requests orelse &.{});
    var request_reader: ?std.Io.File.Reader = if (request_file != null and buffered_requests == null) request_file.?.file.reader(io, &request_buffer) else null;
    var idx = load: {
        var requested = RequestedNames.init(backing_allocator);
        defer requested.deinit();
        const load_mode: index_format.LoadMode = if (request_file != null) matched: {
            const reader = if (buffered_requests != null) &small_reader else &request_reader.?.interface;
            const complete = requested.collect(reader, source, buffered_requests != null) catch return diag.oom(.request_collection);
            if (!complete) break :matched .lookup_full_map;
            if (buffered_requests == null) diag.lookup_scope = .retained_matches;
            break :matched .{ .matched = requested.names.items };
        } else .lookup_full_map;
        break :load try index_format.loadIndex(backing_allocator, io, fasta_path, .automatic, load_mode, diag.failure);
    };
    defer idx.deinit();
    diag.index_source = if (idx.source == .zfi) .zfi else .fai;
    if (buffered_requests) |bytes| small_reader = std.Io.Reader.fixed(bytes);
    if (request_reader != null) request_reader = request_file.?.file.reader(io, &request_buffer);
    return emitGet(backing_allocator, io, &idx, options, &.{}, request_file, if (buffered_requests != null) &small_reader else if (request_reader) |*reader| &reader.interface else null, if (request_reader) |*reader| reader else null, diag, buffers);
}

fn emitGet(
    allocator: std.mem.Allocator,
    io: std.Io,
    idx: *LoadedIndex,
    options: Options,
    requests: []const ParsedRequest,
    request_file: ?SeekableRequest,
    prepared_reader: ?*std.Io.Reader,
    prepared_file_reader: ?*std.Io.File.Reader,
    diag: *Diagnostics,
    buffers: *GetBuffers,
) diagnostics.Error!BatchStats {
    const fasta_source = fastaSource(idx);
    var out_buf: [65536]u8 = undefined;
    var stdout_fw = std.Io.File.Writer.initStreaming(.stdout(), io, &out_buf);
    const writer = &stdout_fw.interface;
    errdefer {
        writer.end = 0;
        if (diag.failure.primary()) |primary| {
            if (primary.at.source == .stdout) {
                if (stdout_fw.err) |err| diag.failure.primary_storage.cause = .{ .zig = err };
            }
        }
    }
    const totals = switch (options.source) {
        .positional => try processParsedRequests(allocator, idx, requests, &.{}, options.annotate_transform, writer, fasta_source, diag, 0, 0, buffers),
        .names => |path| try processRequestPath(allocator, io, idx, path, request_file, prepared_reader, prepared_file_reader, .{ .names = options.orientation }, options.annotate_transform, writer, fasta_source, diag, buffers),
        .bed => |path| try processRequestPath(allocator, io, idx, path, request_file, prepared_reader, prepared_file_reader, .{ .bed = .{ .honor_strand = options.honor_strand, .global_orientation = options.orientation } }, options.annotate_transform, writer, fasta_source, diag, buffers),
    };
    if (totals.region_count == 0) return diag.fail(.{ .@"request.no_regions" = .{ .request_source = switch (options.source) {
        .positional => .positional,
        .names => .names,
        .bed => .bed,
    } } }, .resolve, diag.at.source orelse .request_snapshot, null, null);
    diag.failure.work_complete = true;
    if (diagnostics.flushOutput(&stdout_fw, &diag.failure.outputs.stdout)) |cause| return diag.fail(.{ .@"io.flush_failed" = {} }, .flush, .stdout, null, cause);
    return totals;
}

fn loadTestIndex(test_io: std.Io, fasta_path: []const u8) diagnostics.Error!LoadedIndex {
    var failure: diagnostics.Failure = undefined;
    failure.reset(.get);
    return index_format.loadIndex(std.testing.allocator, test_io, fasta_path, .automatic, .lookup_full_map, &failure);
}

fn collectNamesForAllocationCheck(allocator: std.mem.Allocator, stable_input: bool) !void {
    var requested = RequestedNames.init(allocator);
    defer requested.deinit();
    var buffer: [8192]u8 = undefined;
    var writer = std.Io.Writer.fixed(&buffer);
    try writer.writeAll("# skip\na\nb\na\ncolon:name\na\n");
    for (0..300) |i| try writer.print("name{d}\n", .{i});
    var reader = std.Io.Reader.fixed(writer.buffered());
    try std.testing.expect(try requested.collect(&reader, .{ .names = .{} }, stable_input));
    try std.testing.expectEqual(@as(usize, 303), requested.names.items.len);
    try std.testing.expectEqualStrings("name299", requested.names.items[302]);
}

fn checkGetRequestAllocations(allocator: std.mem.Allocator) !void {
    var idx = try loadTestIndex(std.testing.io, "tests/data/simple.fasta");
    defer idx.deinit();
    var output = std.Io.Writer.Allocating.init(std.testing.allocator);
    defer output.deinit();
    var reader = std.Io.Reader.fixed("# skip\nseq1\t0\t2\nseq2\t0\t2\n");
    var failure: diagnostics.Failure = undefined;
    failure.reset(.get);
    var context = Diagnostics{ .failure = &failure, .at = .{ .source = .bed } };
    var storage: GetBuffers = undefined;
    const buffers = &storage;
    const result = processRequestReader(allocator, &idx, &reader, .{ .bed = .{ .honor_strand = false, .global_orientation = .{} } }, 0, false, &output.writer, fastaSource(&idx), &context, null, buffers) catch |err| {
        try std.testing.expectEqual(diagnostics.Code.@"resource.out_of_memory", failure.primary().?.code());
        try std.testing.expectEqual(error.OutOfMemory, failure.primary().?.cause.?.zig);
        return err;
    };
    try std.testing.expectEqual(@as(u64, 4), result.total_bases);
}

test "[failure] - [get requested names]: releases partial collection" {
    for ([_]bool{ false, true }) |stable_input| {
        try std.testing.checkAllAllocationFailures(std.testing.allocator, collectNamesForAllocationCheck, .{stable_input});
    }
}

test "[unit] - [get read grouping]: accepts overlap and adjacency in either order" {
    try std.testing.expect(spansOverlapOrTouch(10, 20, 20, 30));
    try std.testing.expect(spansOverlapOrTouch(20, 30, 10, 20));
    try std.testing.expect(spansOverlapOrTouch(10, 20, 15, 25));
    try std.testing.expect(spansOverlapOrTouch(15, 25, 10, 20));
    try std.testing.expect(!spansOverlapOrTouch(10, 20, 21, 30));
    try std.testing.expect(!spansOverlapOrTouch(21, 30, 10, 20));
}

test "[unit] - [get request workspace]: reuses consecutive names" {
    var failure: diagnostics.Failure = undefined;
    failure.reset(.get);
    var context = Diagnostics{ .failure = &failure, .at = .{ .source = .names } };
    const diag = &context;
    var workspace = RequestWorkspace.init(std.testing.allocator);
    defer workspace.deinit(std.testing.allocator);

    try appendNamesLine(std.testing.allocator, &workspace, "seq1", 1, .{}, diag);
    try appendNamesLine(std.testing.allocator, &workspace, "seq1", 2, .{}, diag);
    try appendNamesLine(std.testing.allocator, &workspace, "seq2", 3, .{}, diag);

    try std.testing.expectEqualStrings("seq1seq2", workspace.names.items);
    try std.testing.expectEqual(@as(usize, 3), workspace.requests.items.len);
    try std.testing.expectEqualStrings("seq1", workspace.requests.items[0].parsed(workspace.names.items).region.name);
    try std.testing.expectEqualStrings("seq1", workspace.requests.items[1].parsed(workspace.names.items).region.name);
    try std.testing.expectEqualStrings("seq2", workspace.requests.items[2].parsed(workspace.names.items).region.name);
}

test "[edge] - [get request workspace]: retains the name crossing its byte target" {
    var failure: diagnostics.Failure = undefined;
    failure.reset(.get);
    var context = Diagnostics{ .failure = &failure, .at = .{ .source = .names } };
    const diag = &context;
    var workspace = RequestWorkspace.init(std.testing.allocator);
    defer workspace.deinit(std.testing.allocator);
    const name = try std.testing.allocator.alloc(u8, MAX_REQUEST_NAME_BYTES);
    defer std.testing.allocator.free(name);
    @memset(name, 'A');

    for (0..64) |i| {
        name[0] = if (i % 2 == 0) 'A' else 'B';
        try appendNamesLine(std.testing.allocator, &workspace, name, i + 1, .{}, diag);
    }
    try std.testing.expect(workspace.names.items.len < ACTIVE_NAME_BYTES);

    name[0] = 'A';
    try appendNamesLine(std.testing.allocator, &workspace, name, 65, .{}, diag);

    try std.testing.expect(workspace.names.items.len >= ACTIVE_NAME_BYTES);
    try std.testing.expect(workspace.names.items.len <= MAX_ACTIVE_NAME_BYTES);
    try std.testing.expect(workspace.names.capacity <= MAX_ACTIVE_NAME_BYTES);
}

test "[unit] - [get names reader]: handles line endings and a final unterminated line" {
    const test_io = std.testing.io;
    var idx = try loadTestIndex(test_io, "tests/data/simple.fasta");
    defer idx.deinit();
    var reader = std.Io.Reader.fixed("# comment\r\n\r\nseq2\r\nseq1\nseq2");
    var writer = std.Io.Writer.Allocating.init(std.testing.allocator);
    defer writer.deinit();

    var failure: diagnostics.Failure = undefined;
    failure.reset(.get);
    var context = Diagnostics{ .failure = &failure, .at = .{ .source = .names } };
    var storage: GetBuffers = undefined;
    const buffers = &storage;
    const result = try processRequestReader(
        std.testing.allocator,
        &idx,
        &reader,
        .{ .names = .{} },
        0,
        false,
        &writer.writer,
        fastaSource(&idx),
        &context,
        null,
        buffers,
    );

    try std.testing.expectEqual(@as(usize, 3), result.region_count);
    try std.testing.expectEqual(@as(u64, 48), result.total_bases);
    try std.testing.expectEqualStrings(
        ">seq2\nGGGGCCCCAAAA\n>seq1\nACGTACGTACGTACGTACGTACGT\n>seq2\nGGGGCCCCAAAA\n",
        writer.written(),
    );
}

test "[edge] - [get names reader]: resets storage at the request-count boundary" {
    const test_io = std.testing.io;
    var idx = try loadTestIndex(test_io, "tests/data/simple.fasta");
    defer idx.deinit();
    var names = std.Io.Writer.Allocating.init(std.testing.allocator);
    defer names.deinit();
    const request_count = NAMES_REQUEST_BATCH_SIZE + 1;
    for (0..request_count) |_| try names.writer.writeAll("seq1\n");
    var reader = std.Io.Reader.fixed(names.written());
    var writer = std.Io.Writer.Allocating.init(std.testing.allocator);
    defer writer.deinit();

    var failure: diagnostics.Failure = undefined;
    failure.reset(.get);
    var context = Diagnostics{ .failure = &failure, .at = .{ .source = .names } };
    var storage: GetBuffers = undefined;
    const buffers = &storage;
    const result = try processRequestReader(
        std.testing.allocator,
        &idx,
        &reader,
        .{ .names = .{} },
        0,
        false,
        &writer.writer,
        fastaSource(&idx),
        &context,
        null,
        buffers,
    );

    const one_output = ">seq1\nACGTACGTACGTACGTACGTACGT\n";
    try std.testing.expectEqual(request_count, result.region_count);
    try std.testing.expectEqual(@as(u64, request_count * 24), result.total_bases);
    try std.testing.expectEqual(one_output.len * request_count, writer.written().len);
    try std.testing.expectEqualStrings(one_output, writer.written()[0..one_output.len]);
    try std.testing.expectEqualStrings(one_output, writer.written()[writer.written().len - one_output.len ..]);
}

test "[unit] - [get BED reader]: streams rows in source order" {
    const test_io = std.testing.io;
    var idx = try loadTestIndex(test_io, "tests/data/simple.fasta");
    defer idx.deinit();

    const bed_data =
        "# comment\r\n" ++
        "seq2\t0\t4\r\n" ++
        "seq1\t0\t4\n" ++
        "seq1\t4\t8";

    var reader = std.Io.Reader.fixed(bed_data);
    var writer = std.Io.Writer.Allocating.init(std.testing.allocator);
    defer writer.deinit();

    var failure: diagnostics.Failure = undefined;
    failure.reset(.get);
    var context = Diagnostics{ .failure = &failure, .at = .{ .source = .names } };
    var storage: GetBuffers = undefined;
    const buffers = &storage;
    const result = try processRequestReader(
        std.testing.allocator,
        &idx,
        &reader,
        .{ .bed = .{ .honor_strand = false, .global_orientation = .{} } },
        0,
        false,
        &writer.writer,
        fastaSource(&idx),
        &context,
        null,
        buffers,
    );

    try std.testing.expectEqual(@as(usize, 3), result.region_count);
    try std.testing.expectEqual(@as(u64, 12), result.total_bases);
    try std.testing.expectEqualStrings(
        ">seq2:1-4\nGGGG\n>seq1:1-4\nACGT\n>seq1:5-8\nACGT\n",
        writer.written(),
    );
}

test "[edge] - [get BED reader]: preserves output across request-batch boundaries" {
    const test_io = std.testing.io;
    var idx = try loadTestIndex(test_io, "tests/data/simple.fasta");
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
        var failure: diagnostics.Failure = undefined;
        failure.reset(.get);
        var context = Diagnostics{ .failure = &failure, .at = .{ .source = .names } };
        var storage: GetBuffers = undefined;
        const buffers = &storage;
        const result = try processRequestReader(
            std.testing.allocator,
            &idx,
            &reader,
            .{ .bed = .{ .honor_strand = true, .global_orientation = .{} } },
            0,
            false,
            &output.writer,
            fastaSource(&idx),
            &context,
            null,
            buffers,
        );

        const expected_bases: u64 = count + @as(usize, @intFromBool(count > BED_REQUEST_BATCH_SIZE)) * 4;
        try std.testing.expectEqual(count, result.region_count);
        try std.testing.expectEqual(expected_bases, result.total_bases);
        try std.testing.expectEqualStrings(expected.written(), output.written());
    }
}

test "[failure] - [GET request facts]: keeps original rows and ordinals across skipped rows and batches" {
    var idx = try loadTestIndex(std.testing.io, "tests/data/simple.fasta");
    defer idx.deinit();
    var rows = std.Io.Writer.Allocating.init(std.testing.allocator);
    defer rows.deinit();
    try rows.writer.writeAll("# skip\n\n");
    for (0..BED_REQUEST_BATCH_SIZE) |_| try rows.writer.writeAll("seq1\t0\t1\n");
    try rows.writer.writeAll("# later skip\nmissing\t0\t1\n");
    var reader = std.Io.Reader.fixed(rows.written());
    var output = std.Io.Writer.Allocating.init(std.testing.allocator);
    defer output.deinit();
    var failure: diagnostics.Failure = undefined;
    failure.reset(.get);
    var context = Diagnostics{ .failure = &failure, .at = .{ .source = .bed } };
    var storage: GetBuffers = undefined;
    const buffers = &storage;
    try std.testing.expectError(error.InvalidRequest, processRequestReader(std.testing.allocator, &idx, &reader, .{ .bed = .{ .honor_strand = false, .global_orientation = .{} } }, 0, false, &output.writer, fastaSource(&idx), &context, null, buffers));
    const primary = failure.primary().?;
    try std.testing.expectEqual(@as(?u64, BED_REQUEST_BATCH_SIZE + 4), primary.at.location.line);
    try std.testing.expectEqual(@as(?u64, BED_REQUEST_BATCH_SIZE + 1), primary.at.request);
    try std.testing.expectEqual(@as(?u64, BED_REQUEST_BATCH_SIZE), failure.outputs.stdout.records);
    try std.testing.expect(primary.at.location.offset == null);
    try std.testing.expectEqualStrings("missing", failure.preview(primary.details.@"request.sequence_not_found".identifier.?));
    try std.testing.expectEqual(@as(usize, 12 * BED_REQUEST_BATCH_SIZE), output.written().len);
    const PreviousRequest = struct { name_offset: u32, name_len: u16, start: u64, end: u64, orientation: Orientation };
    try std.testing.expectEqual(@sizeOf(PreviousRequest) + 8, @sizeOf(OwnedRequest));
}

test "[failure] - [GET totals]: rejects repeated-request overflow before writing another region" {
    var idx = try loadTestIndex(std.testing.io, "tests/data/simple.fasta");
    defer idx.deinit();
    var output = std.Io.Writer.Allocating.init(std.testing.allocator);
    defer output.deinit();
    var failure: diagnostics.Failure = undefined;
    failure.reset(.get);
    var context = Diagnostics{ .failure = &failure, .at = .{ .source = .request_snapshot } };
    var storage: GetBuffers = undefined;
    const buffers = &storage;
    const request = [_]ParsedRequest{.{ .region = parseRegion("seq1:1-2"), .orientation = .{} }};
    try std.testing.expectError(error.CountOverflow, processParsedRequests(std.testing.allocator, &idx, &request, &.{}, false, &output.writer, fastaSource(&idx), &context, 0, std.math.maxInt(u64) - 1, buffers));
    try std.testing.expectEqual(@as(usize, 0), output.written().len);
    const detail = failure.primary().?.details.@"resource.count_overflow";
    try std.testing.expectEqual(diagnostics.Counter.total_symbols, detail.counter);
    try std.testing.expectEqual(@as(?u128, 2), detail.increment);
}

test "[failure] - [GET requests]: releases partial batches on allocation failure" {
    try std.testing.checkAllAllocationFailures(std.testing.allocator, checkGetRequestAllocations, .{});
}

test "[failure] - [GET request replay]: detects metadata changes and limits equal-metadata misses" {
    const io = std.testing.io;
    var directory = std.testing.tmpDir(.{});
    defer directory.cleanup();
    const file = try directory.dir.createFile(io, "requests.txt", .{ .read = true });
    defer file.close(io);
    try file.writePositionalAll(io, "seq1\n", 0);
    const before = try file.stat(io);
    const opened = SeekableRequest{ .file = file, .size = before.size, .mtime_ns = before.mtime.toNanoseconds() };
    var failure: diagnostics.Failure = undefined;
    failure.reset(.get);
    var context = Diagnostics{ .failure = &failure, .at = .{ .source = .names }, .lookup_scope = .retained_matches };
    var storage: GetBuffers = undefined;
    const buffers = &storage;
    try file.writePositionalAll(io, "seq2\nextra\n", 0);
    try std.testing.expectError(error.InputChanged, checkRequestSource(io, opened, &context));
    try std.testing.expect(failure.primary().?.details.@"input.changed".changed_fields.size);

    try file.setLength(io, before.size);
    try file.setTimestamps(io, .{ .modify_timestamp = .{ .new = before.mtime } });
    failure.reset(.get);
    var idx = try index_format.loadIndex(std.testing.allocator, io, "tests/data/simple.fasta", .automatic, .{ .matched = &.{"seq1"} }, &failure);
    defer idx.deinit();
    try checkRequestSource(io, opened, &context);
    var input: [64]u8 = undefined;
    var reader = file.reader(io, &input);
    var output = std.Io.Writer.Allocating.init(std.testing.allocator);
    defer output.deinit();
    try std.testing.expectError(error.InvalidRequest, processRequestReader(std.testing.allocator, &idx, &reader.interface, .{ .names = .{} }, 0, false, &output.writer, fastaSource(&idx), &context, opened, buffers));
    try std.testing.expectEqual(diagnostics.LookupScope.retained_matches, failure.primary().?.details.@"request.sequence_not_found".lookup_scope);
    try std.testing.expectEqualStrings("seq2", failure.preview(failure.primary().?.details.@"request.sequence_not_found".identifier.?));
    try std.testing.expectEqual(@as(usize, 0), output.written().len);
}

test "[failure] - [GET request IO]: retains required read and comparison causes" {
    const io = std.testing.io;
    const Fault = struct {
        fn read(_: ?*anyopaque, _: std.Io.File, _: []const []u8, _: u64) std.Io.File.ReadPositionalError!usize {
            return error.InputOutput;
        }
        fn stat(_: ?*anyopaque, _: std.Io.File) std.Io.File.StatError!std.Io.File.Stat {
            return error.AccessDenied;
        }
    };
    var directory = std.testing.tmpDir(.{});
    defer directory.cleanup();
    const file = try directory.dir.createFile(io, "requests.txt", .{ .read = true });
    defer file.close(io);
    try file.writePositionalAll(io, "seq1\n", 0);
    const stat = try file.stat(io);
    var vtable = io.vtable.*;
    const fault_io = std.Io{ .userdata = io.userdata, .vtable = &vtable };
    var failure: diagnostics.Failure = undefined;
    failure.reset(.get);
    var context = Diagnostics{ .failure = &failure, .at = .{ .source = .names } };
    var storage: GetBuffers = undefined;
    const buffers = &storage;
    vtable.fileStat = Fault.stat;
    try std.testing.expectError(error.Io, checkRequestSource(fault_io, .{
        .file = file,
        .size = stat.size,
        .mtime_ns = stat.mtime.toNanoseconds(),
    }, &context));
    try std.testing.expectEqual(diagnostics.Code.@"io.stat_failed", failure.primary().?.code());
    try std.testing.expectEqual(error.AccessDenied, failure.primary().?.cause.?.zig);
    try std.testing.expectEqual(diagnostics.Stage.compare_source, failure.primary().?.stage);
    failure.reset(.get);
    vtable.fileStat = io.vtable.fileStat;
    vtable.fileReadPositional = Fault.read;
    var buffer: [128]u8 = undefined;
    var reader = file.reader(fault_io, &buffer);
    context.request_reader = &reader;
    var idx = try loadTestIndex(io, "tests/data/simple.fasta");
    defer idx.deinit();
    var writer = std.Io.Writer.Allocating.init(std.testing.allocator);
    defer writer.deinit();
    try std.testing.expectError(error.Io, processRequestReader(std.testing.allocator, &idx, &reader.interface, .{ .names = .{} }, 0, false, &writer.writer, fastaSource(&idx), &context, null, buffers));
    try std.testing.expectEqual(diagnostics.Code.@"io.read_failed", failure.primary().?.code());
    try std.testing.expectEqual(error.InputOutput, failure.primary().?.cause.?.zig);
    try std.testing.expectEqual(@as(?u64, 1), failure.primary().?.at.location.line);
    try std.testing.expectEqual(@as(?u64, 0), failure.primary().?.at.location.offset);
    try std.testing.expect(failure.primary().?.at.location.column == null);
    try std.testing.expectEqual(@as(usize, 0), writer.written().len);
}

test "[regression] - [GET discovery]: retries interrupted discovery with the full catalog" {
    const io = std.testing.io;
    const Fault = struct {
        var scenario: enum { first, partial, persistent } = .first;
        var request_handle: ?std.Io.File.Handle = null;
        var read_count: usize = 0;
        var first_row_bytes: usize = 0;
        var failed: bool = false;
        var output: [128]u8 = undefined;
        var output_len: usize = 0;

        fn read(userdata: ?*anyopaque, file: std.Io.File, data: []const []u8, offset: u64) std.Io.File.ReadPositionalError!usize {
            if (request_handle == null) request_handle = file.handle;
            if (file.handle == request_handle.?) {
                const attempt = read_count;
                read_count += 1;
                if (scenario == .persistent or (!failed and attempt == @as(usize, if (scenario == .partial) 1 else 0))) {
                    failed = true;
                    return error.InputOutput;
                }
                if (scenario == .partial and attempt == 0) {
                    return std.testing.io.vtable.fileReadPositional(userdata, file, &.{data[0][0..first_row_bytes]}, offset);
                }
            }
            return std.testing.io.vtable.fileReadPositional(userdata, file, data, offset);
        }

        fn operate(userdata: ?*anyopaque, operation: std.Io.Operation) std.Io.Cancelable!std.Io.Operation.Result {
            switch (operation) {
                .file_write_streaming => |write| {
                    const start = output_len;
                    @memcpy(output[output_len..][0..write.header.len], write.header);
                    output_len += write.header.len;
                    for (write.data, 0..) |part, i| {
                        const repeats = if (i + 1 == write.data.len) write.splat else 1;
                        for (0..repeats) |_| {
                            @memcpy(output[output_len..][0..part.len], part);
                            output_len += part.len;
                        }
                    }
                    return .{ .file_write_streaming = output_len - start };
                },
                else => return std.testing.io.vtable.operate(userdata, operation),
            }
        }
    };
    var directory = std.testing.tmpDir(.{});
    defer directory.cleanup();
    const path = try directory.dir.realPathFileAlloc(io, ".", std.testing.allocator);
    defer std.testing.allocator.free(path);
    const request_path = try std.fmt.allocPrint(std.testing.allocator, "{s}/requests.txt", .{path});
    defer std.testing.allocator.free(request_path);
    var vtable = io.vtable.*;
    vtable.fileReadPositional = Fault.read;
    vtable.operate = Fault.operate;
    const fault_io = std.Io{ .userdata = io.userdata, .vtable = &vtable };
    for ([_]bool{ false, true }) |bed| {
        const first = if (bed) "seq1\t0\t2\n" else "seq1\n";
        const last = if (bed) "seq2\t0\t2\n" else "seq2\n";
        const padding = "#\n" ** (REQUEST_LINE_READER_BUFFER_BYTES / 2 + 1);
        {
            const file = try directory.dir.createFile(io, "requests.txt", .{});
            defer file.close(io);
            try file.writePositionalAll(io, first, 0);
            try file.writePositionalAll(io, padding, first.len);
            try file.writePositionalAll(io, last, first.len + padding.len);
        }
        const options = Options{ .source = if (bed) .{ .bed = request_path } else .{ .names = request_path } };
        for ([_]@TypeOf(Fault.scenario){ .first, .partial, .persistent }) |scenario| {
            Fault.scenario = scenario;
            Fault.request_handle = null;
            Fault.read_count = 0;
            Fault.first_row_bytes = first.len;
            Fault.failed = false;
            Fault.output_len = 0;
            var failure: diagnostics.Failure = undefined;
            if (scenario == .persistent) {
                try std.testing.expectError(error.Io, runGet(std.testing.allocator, fault_io, "tests/data/simple.fasta", options, &failure));
                const primary = failure.primary().?;
                try std.testing.expectEqual(diagnostics.Code.@"io.read_failed", primary.code());
                try std.testing.expectEqual(error.InputOutput, primary.cause.?.zig);
                try std.testing.expectEqual(if (bed) diagnostics.SourceId.bed else .names, primary.at.source);
                try std.testing.expectEqual(@as(?u64, 1), primary.at.location.line);
                try std.testing.expectEqual(@as(?u64, 0), primary.at.location.offset);
                try std.testing.expectEqual(@as(usize, 0), Fault.output_len);
                try std.testing.expect(!failure.complete and !failure.work_complete);
            } else {
                const result = try runGet(std.testing.allocator, fault_io, "tests/data/simple.fasta", options, &failure);
                try std.testing.expectEqual(@as(usize, 2), result.region_count);
                try std.testing.expectEqual(@as(u64, if (bed) 4 else 36), result.total_bases);
                const expected = if (bed) ">seq1:1-2\nAC\n>seq2:1-2\nGG\n" else ">seq1\nACGTACGTACGTACGTACGTACGT\n>seq2\nGGGGCCCCAAAA\n";
                try std.testing.expectEqualStrings(expected, Fault.output[0..Fault.output_len]);
                try std.testing.expect(failure.primary() == null);
                try std.testing.expect(failure.complete and failure.work_complete);
            }
            try std.testing.expect(Fault.failed);
        }
    }
}
