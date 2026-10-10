//! Assembly and proteome statistics from indexed FASTA files.
//!
//! Length statistics come from index metadata. Type and composition use one
//! descriptor-owned bounded read window over indexed sequence spans.

const std = @import("std");

const index_format = @import("index_format.zig");
const diagnostics = @import("diagnostics");

/// Sequence class shared by stats, GET transformation guards, and validation.
pub const SequenceType = enum { nucleotide, protein };

/// Bases sampled per record for GET transformation guards.
pub const GET_TYPE_SAMPLE_BASES: u64 = 256;

/// Bases sampled across the file for validation alphabet selection.
pub const VALIDATE_TYPE_SAMPLE_BASES: u64 = 100_000;

/// Computes indexed statistics and delivers the report to stdout. Path labels
/// borrow fasta_path until the caller handles the failure or resets storage.
/// Resets once, returns after cleanup, and discards pending output on failure.
pub fn runStats(allocator: std.mem.Allocator, io: std.Io, fasta_path: []const u8, failure: *diagnostics.Failure) diagnostics.Error!void {
    failure.reset(.stats);
    failure.sources.set(.stdout, .{ .label = .{ .borrowed = "stdout" }, .kind = .stream });
    var diag = Diagnostics{ .failure = failure };
    var idx = try index_format.loadIndexCheckedWithMode(allocator, io, fasta_path, .stats_scan, failure);
    defer idx.deinit();
    var computed: ComputedStats = undefined;
    var buffer: [STATS_READ_BUFFER_BYTES]u8 = undefined;
    try computeStats(allocator, &idx, &computed, &buffer, &diag);
    var stdout_fw = std.Io.File.Writer.initStreaming(.stdout(), io, buffer[0..65536]);
    errdefer {
        stdout_fw.interface.end = 0;
        if (stdout_fw.err) |cause| {
            if (failure.primary()) |primary| {
                if (primary.code() == .@"io.write_failed") failure.primary_storage.cause = .{ .zig = cause };
            }
        }
    }
    try writeComputedReport(&stdout_fw.interface, fasta_path, &idx, &computed, &diag);
    stdout_fw.flush() catch |err| {
        failure.outputs.stdout.final_flush = .failed;
        return diag.fail(.{ .@"io.flush_failed" = {} }, .flush, .stdout, null, stdout_fw.err orelse err);
    };
    failure.outputs.stdout.final_flush = .succeeded;
    failure.succeed();
}

/// Detects whether represented sequence counts are nucleotide or protein.
///
/// Shared by stats (full-file counts), GET (per-record sample), and validate
/// (prefix sample). Threshold: nucleotide if IUPAC nuc letters are strictly
/// more than 90% of counted bases (`nuc * 10 > total * 9`). The caller supplies
/// `total` for the same sample represented by `counts`. Empty samples are nucleotide.
pub fn detectType(counts: *const [256]u64, total: u64) SequenceType {
    if (total == 0) return .nucleotide;

    const nuc_count = counts['A'] + counts['a'] +
        counts['C'] + counts['c'] +
        counts['G'] + counts['g'] +
        counts['T'] + counts['t'] +
        counts['R'] + counts['r'] +
        counts['Y'] + counts['y'] +
        counts['S'] + counts['s'] +
        counts['W'] + counts['w'] +
        counts['K'] + counts['k'] +
        counts['M'] + counts['m'] +
        counts['B'] + counts['b'] +
        counts['D'] + counts['d'] +
        counts['H'] + counts['h'] +
        counts['V'] + counts['v'] +
        counts['N'] + counts['n'] +
        counts['U'] + counts['u'];

    if (@as(u128, nuc_count) * 10 > @as(u128, total) * 9) {
        return .nucleotide;
    }
    return .protein;
}

const IndexRecord = index_format.IndexRecord;

const LoadedIndex = index_format.LoadedIndex;

const SideTableLine = index_format.SideTableLine;

const Diagnostics = struct {
    failure: *diagnostics.Failure,
    index_source: diagnostics.SourceId = .zfi,
    record: ?u64 = null,
    output_source: diagnostics.SourceId = .stdout,

    fn fail(self: *@This(), details: diagnostics.Details, stage: diagnostics.Stage, source: diagnostics.SourceId, offset: ?u64, cause: ?anyerror) diagnostics.Error {
        @branchHint(.cold);
        return self.failure.fail(.{
            .details = details,
            .stage = stage,
            .at = .{ .source = source, .location = .{ .offset = offset }, .record = self.record },
            .cause = if (cause) |err| .{ .zig = err } else null,
        });
    }

    fn consistency(self: *@This(), check: diagnostics.Check) diagnostics.Error {
        @branchHint(.cold);
        return self.fail(.{ .@"internal.consistency_failed" = .{ .component = .stats, .check = check } }, .scan, self.index_source, null, null);
    }

    fn overflow(self: *@This(), counter: diagnostics.Counter, current: u128, increment: u128, maximum: u128, source: diagnostics.SourceId) diagnostics.Error {
        @branchHint(.cold);
        return self.fail(.{ .@"resource.count_overflow" = .{
            .counter = counter,
            .current = current,
            .increment = increment,
            .max_value = maximum,
        } }, .scan, source, null, null);
    }

    fn oom(self: *@This(), role: diagnostics.AllocationRole) diagnostics.Error {
        @branchHint(.cold);
        return self.fail(.{ .@"resource.out_of_memory" = .{ .allocation_role = role } }, .allocate, self.index_source, null, error.OutOfMemory);
    }

    fn span(self: *@This(), end: u128, size: u64, entry: ?u64) diagnostics.Error {
        @branchHint(.cold);
        return self.fail(.{ .@"index.sequence_span_out_of_bounds" = .{
            .required_end = end,
            .source_bytes = size,
            .bound_source = .actual,
            .entry = entry,
        } }, .scan, self.index_source, null, null);
    }

    fn read(self: *@This(), offset: u64, wanted: u64, received: ?u64, cause: ?anyerror) diagnostics.Error {
        @branchHint(.cold);
        return self.fail(if (received) |got| .{ .@"io.unexpected_eof" = .{ .exact = .{
            .requested_bytes = wanted,
            .received_bytes = got,
        } } } else .{ .@"io.read_failed" = .{ .requested_bytes = wanted } }, .read, .fasta, offset, cause);
    }
};

const SIMD_CHUNK_SIZE = 32;

const SimdVec = @Vector(SIMD_CHUNK_SIZE, u8);

const TAIL_CHUNK_SIZE = 8;

const TailVec = @Vector(TAIL_CHUNK_SIZE, u8);

const FINAL_CHUNK_SIZE = 4;

const FinalVec = @Vector(FINAL_CHUNK_SIZE, u8);

const CompositionStats = struct {
    counts: [256]u64 = .{0} ** 256,
    total_bases: u64 = 0,
    seq_type: SequenceType = .nucleotide,
    lowercase_count: u64 = 0,
};

const LengthSummary = struct {
    total_symbols: u64,
    shortest_index: usize,
    longest_index: usize,
    mean: u64,
    q1: u64,
    median: u64,
    q3: u64,
    range: u64,
    n50: u64,
    l50: usize,
    n90: u64,
    l90: usize,
    aun_numerator: u128,
};

const NucleotideSummary = struct {
    type_name: []const u8,
    categories: [8]u64,
    canonical: u64,
    gc_total: u64,
};

const ProteinSummary = struct {
    categories: [PROTEIN_FIELDS.len]u64,
    stop: u64,
    invalid: u64,
};

const CompositionSummary = union(enum) {
    nucleotide: NucleotideSummary,
    protein: ProteinSummary,
};

const ProteinField = struct {
    key: []const u8,
    code: u8,
};

const PROTEIN_FIELDS = [_]ProteinField{
    .{ .key = "a_alanine", .code = 'A' },
    .{ .key = "r_arginine", .code = 'R' },
    .{ .key = "n_asparagine", .code = 'N' },
    .{ .key = "d_aspartate", .code = 'D' },
    .{ .key = "c_cysteine", .code = 'C' },
    .{ .key = "e_glutamate", .code = 'E' },
    .{ .key = "q_glutamine", .code = 'Q' },
    .{ .key = "g_glycine", .code = 'G' },
    .{ .key = "h_histidine", .code = 'H' },
    .{ .key = "i_isoleucine", .code = 'I' },
    .{ .key = "l_leucine", .code = 'L' },
    .{ .key = "k_lysine", .code = 'K' },
    .{ .key = "m_methionine", .code = 'M' },
    .{ .key = "f_phenylalanine", .code = 'F' },
    .{ .key = "p_proline", .code = 'P' },
    .{ .key = "s_serine", .code = 'S' },
    .{ .key = "t_threonine", .code = 'T' },
    .{ .key = "w_tryptophan", .code = 'W' },
    .{ .key = "y_tyrosine", .code = 'Y' },
    .{ .key = "v_valine", .code = 'V' },
    .{ .key = "b_asx", .code = 'B' },
    .{ .key = "z_glx", .code = 'Z' },
    .{ .key = "j_xle", .code = 'J' },
    .{ .key = "x_unknown", .code = 'X' },
    .{ .key = "u_selenocysteine", .code = 'U' },
    .{ .key = "o_pyrrolysine", .code = 'O' },
};

fn medianSorted(values: []const u64) u64 {
    const middle = values.len / 2;
    if (values.len % 2 == 1) return values[middle];
    const lower = values[middle - 1];
    return lower + (values[middle] - lower) / 2;
}

fn summarizeLengths(records: []const IndexRecord, lengths: []u64, diag: *Diagnostics) diagnostics.Error!LengthSummary {
    if (records.len == 0 or records.len != lengths.len) return diag.consistency(.record_reservation);
    var total: u64 = 0;
    var shortest_index: usize = 0;
    var longest_index: usize = 0;
    var aun_numerator: u128 = 0;
    for (records, 0..) |record, i| {
        diag.record = i + 1;
        total = std.math.add(u64, total, record.seq_len) catch return diag.overflow(.total_symbols, total, record.seq_len, std.math.maxInt(u64), diag.index_source);
        lengths[i] = record.seq_len;
        if (record.seq_len < records[shortest_index].seq_len) shortest_index = i;
        if (record.seq_len > records[longest_index].seq_len) longest_index = i;
        const square = @as(u128, record.seq_len) * @as(u128, record.seq_len);
        aun_numerator = std.math.add(u128, aun_numerator, square) catch return diag.consistency(.numeric_buffer);
    }

    diag.record = null;
    std.mem.sort(u64, lengths, {}, std.sort.asc(u64));
    const half = lengths.len / 2;
    const q1 = if (lengths.len == 1) lengths[0] else medianSorted(lengths[0..half]);
    const q3 = if (lengths.len == 1) lengths[0] else medianSorted(lengths[lengths.len - half ..]);
    const threshold_50: u128 = (@as(u128, total) * 50 + 99) / 100;
    const threshold_90: u128 = (@as(u128, total) * 90 + 99) / 100;
    var seen: u128 = 0;
    var n50: u64 = 0;
    var l50: usize = 0;
    var n90: u64 = 0;
    var l90: usize = 0;
    var rank: usize = 0;
    var i = lengths.len;
    while (i > 0) {
        i -= 1;
        rank += 1;
        seen += lengths[i];
        if (l50 == 0 and seen >= threshold_50) {
            n50 = lengths[i];
            l50 = rank;
        }
        if (l90 == 0 and seen >= threshold_90) {
            n90 = lengths[i];
            l90 = rank;
            break;
        }
    }

    return .{
        .total_symbols = total,
        .shortest_index = shortest_index,
        .longest_index = longest_index,
        .mean = total / records.len,
        .q1 = q1,
        .median = medianSorted(lengths),
        .q3 = q3,
        .range = records[longest_index].seq_len - records[shortest_index].seq_len,
        .n50 = n50,
        .l50 = l50,
        .n90 = n90,
        .l90 = l90,
        .aun_numerator = aun_numerator,
    };
}

fn writeFixedUnsigned(writer: *std.Io.Writer, numerator: u128, denominator: u128, decimals: u8) std.Io.Writer.Error!void {
    const scale: u128 = switch (decimals) {
        2 => 100,
        3 => 1000,
        else => unreachable,
    };
    var whole = numerator / denominator;
    const remainder = numerator % denominator;
    var fraction = (remainder * scale + denominator / 2) / denominator;
    if (fraction == scale) {
        whole += 1;
        fraction = 0;
    }
    try writer.print("{d}.", .{whole});
    var divisor = scale / 10;
    while (divisor > 0) : (divisor /= 10) {
        try writer.writeByte(@intCast('0' + (fraction / divisor) % 10));
    }
}

fn writeFixedSigned(writer: *std.Io.Writer, numerator: i128, denominator: u128, decimals: u8) std.Io.Writer.Error!void {
    const magnitude: u128 = @intCast(if (numerator < 0) -numerator else numerator);
    const scale: u128 = if (decimals == 3) 1000 else 100;
    const rounded = (magnitude % denominator * scale + denominator / 2) / denominator;
    if (numerator < 0 and (magnitude / denominator != 0 or rounded != 0)) try writer.writeByte('-');
    try writeFixedUnsigned(writer, magnitude, denominator, decimals);
}

fn countLetter(comp: *const CompositionStats, upper: u8) u128 {
    return @as(u128, comp.counts[upper]) + comp.counts[upper + ('a' - 'A')];
}

fn summarizeComposition(comp: *const CompositionStats, diag: *Diagnostics) diagnostics.Error!CompositionSummary {
    if (comp.total_bases == 0 or comp.lowercase_count > comp.total_bases) return diag.consistency(.composition_accounting);
    if (comp.seq_type == .nucleotide) {
        const a = countLetter(comp, 'A');
        const c = countLetter(comp, 'C');
        const g = countLetter(comp, 'G');
        const t = countLetter(comp, 'T');
        const u = countLetter(comp, 'U');
        const n = countLetter(comp, 'N');
        var ambiguous: u128 = 0;
        for ("RYSWKMBDHV") |code| ambiguous += countLetter(comp, code);
        const assigned = @as(u128, a) + c + g + t + u + n + ambiguous;
        if (assigned > comp.total_bases) return diag.consistency(.composition_accounting);
        return .{ .nucleotide = .{
            .type_name = if (t > 0 and u > 0)
                "nucleotide_mixed_tu"
            else if (t > 0)
                "nucleotide_t"
            else if (u > 0)
                "nucleotide_u"
            else
                "nucleotide",
            .categories = .{ @intCast(a), @intCast(c), @intCast(g), @intCast(t), @intCast(u), @intCast(n), @intCast(ambiguous), @intCast(@as(u128, comp.total_bases) - assigned) },
            .canonical = @intCast(a + c + g + t + u),
            .gc_total = @intCast(g + c),
        } };
    }

    var categories: [PROTEIN_FIELDS.len]u64 = undefined;
    var assigned: u128 = 0;
    for (PROTEIN_FIELDS, 0..) |field, i| {
        const count = countLetter(comp, field.code);
        assigned += count;
        if (assigned > comp.total_bases) return diag.consistency(.composition_accounting);
        categories[i] = @intCast(count);
    }
    const stop = comp.counts['*'];
    assigned += stop;
    if (assigned > comp.total_bases) return diag.consistency(.composition_accounting);
    return .{ .protein = .{
        .categories = categories,
        .stop = stop,
        .invalid = @intCast(@as(u128, comp.total_bases) - assigned),
    } };
}

fn writeCountPercent(writer: *std.Io.Writer, key: []const u8, count: u64, total: u64) std.Io.Writer.Error!void {
    try writer.print("  {s}: {d} ", .{ key, count });
    try writeFixedUnsigned(writer, @as(u128, count) * 100, total, 2);
    try writer.writeAll("%\n");
}

fn writeNucleotide(writer: *std.Io.Writer, summary: NucleotideSummary, comp: *const CompositionStats) std.Io.Writer.Error!void {
    try writer.print("  type: {s}\n  percent_denominator: total_symbols\n", .{summary.type_name});
    const keys = [_][]const u8{ "a", "c", "g", "t", "u", "n", "iupac_ambiguous", "invalid" };
    for (keys, summary.categories) |key, count| try writeCountPercent(writer, key, count, comp.total_bases);
    if (summary.canonical == 0) {
        try writer.writeAll("  gc: n/a\n");
    } else {
        try writer.writeAll("  gc: ");
        try writeFixedUnsigned(writer, @as(u128, summary.gc_total) * 100, summary.canonical, 2);
        try writer.writeAll("%\n");
    }
    if (summary.gc_total == 0) {
        try writer.writeAll("  gc_skew: n/a\n");
    } else {
        try writer.writeAll("  gc_skew: ");
        try writeFixedSigned(
            writer,
            @as(i128, summary.categories[2]) - @as(i128, summary.categories[1]),
            summary.gc_total,
            3,
        );
        try writer.writeByte('\n');
    }
    try writeCountPercent(writer, "lowercase", comp.lowercase_count, comp.total_bases);
}

fn writeProtein(writer: *std.Io.Writer, summary: ProteinSummary, comp: *const CompositionStats) std.Io.Writer.Error!void {
    try writer.writeAll("  type: protein\n  percent_denominator: total_symbols\n");
    for (PROTEIN_FIELDS, summary.categories) |field, count| {
        try writeCountPercent(writer, field.key, count, comp.total_bases);
    }
    try writeCountPercent(writer, "stop", summary.stop, comp.total_bases);
    try writeCountPercent(writer, "invalid", summary.invalid, comp.total_bases);
    try writeCountPercent(writer, "lowercase", comp.lowercase_count, comp.total_bases);
}

fn writeReport(
    writer: *std.Io.Writer,
    fasta_path: []const u8,
    index_extension: []const u8,
    fasta_size: u64,
    records: []const IndexRecord,
    summary: LengthSummary,
    shortest_name: []const u8,
    longest_name: []const u8,
    comp: *const CompositionStats,
    composition: CompositionSummary,
) std.Io.Writer.Error!void {
    try writer.print(
        "File:\n  path: {s}\n  index: {s}{s}\n  size_bytes: {d}\n\n",
        .{ fasta_path, fasta_path, index_extension, fasta_size },
    );
    try writer.print(
        "Lengths:\n  indexed_records: {d}\n  total_symbols: {d}\n" ++
            "  shortest_length: {d}\n  shortest_name: {s}\n" ++
            "  longest_length: {d}\n  longest_name: {s}\n" ++
            "  mean: {d}\n  q1: {d}\n  median: {d}\n  q3: {d}\n  range: {d}\n\n",
        .{
            records.len,
            summary.total_symbols,
            records[summary.shortest_index].seq_len,
            shortest_name,
            records[summary.longest_index].seq_len,
            longest_name,
            summary.mean,
            summary.q1,
            summary.median,
            summary.q3,
            summary.range,
        },
    );
    try writer.print(
        "Nx:\n  n50: {d}\n  l50: {d}\n  n90: {d}\n  l90: {d}\n  aun: ",
        .{ summary.n50, summary.l50, summary.n90, summary.l90 },
    );
    try writeFixedUnsigned(writer, summary.aun_numerator, summary.total_symbols, 2);
    try writer.writeAll("\n\nComposition:\n");
    switch (composition) {
        .nucleotide => |nucleotide| try writeNucleotide(writer, nucleotide, comp),
        .protein => |protein| try writeProtein(writer, protein, comp),
    }
}

const ComputedStats = struct {
    summary: LengthSummary,
    comp: CompositionStats,
    composition: CompositionSummary,
    shortest_name: []const u8,
    longest_name: []const u8,
};

fn computeStats(allocator: std.mem.Allocator, idx: *const LoadedIndex, computed: *ComputedStats, buffer: []u8, diag: *Diagnostics) diagnostics.Error!void {
    diag.index_source = if (idx.source == .zfi) .zfi else .fai;
    const records = idx.records;
    if (records.len == 0) return diag.fail(.{ .@"index.no_records" = .{ .format = if (idx.source == .zfi) .zfi else .fai } }, .scan, diag.index_source, null, null);
    const lengths = allocator.alloc(u64, records.len) catch return diag.oom(.lengths);
    defer allocator.free(lengths);
    const summary = try summarizeLengths(records, lengths, diag);
    const shortest_name = try retainedName(idx, summary.shortest_index, .retained_shortest_name, diag);
    const longest_name = try retainedName(idx, summary.longest_index, .retained_longest_name, diag);
    diag.record = null;
    try scanComposition(allocator, idx, buffer, &computed.comp, diag);
    try checkTotal(summary.total_symbols, computed.comp.total_bases, diag);
    computed.composition = try summarizeComposition(&computed.comp, diag);
    computed.summary = summary;
    computed.shortest_name = shortest_name;
    computed.longest_name = longest_name;
    diag.failure.work_complete = true;
}

fn writeComputedReport(writer: *std.Io.Writer, fasta_path: []const u8, idx: *const LoadedIndex, computed: *const ComputedStats, diag: *Diagnostics) diagnostics.Error!void {
    try writeReportChecked(writer, fasta_path, if (idx.source == .zfi) ".zfi" else ".fai", idx.fasta_size, idx.records, computed.summary, computed.shortest_name, computed.longest_name, &computed.comp, computed.composition, diag);
}

fn checkTotal(indexed: u64, observed: u64, diag: *Diagnostics) diagnostics.Error!void {
    if (observed != indexed) return diag.fail(.{ .@"pair.total_length_mismatch" = .{
        .indexed_total = indexed,
        .observed_total = observed,
    } }, .scan, .fasta, null, null);
}

fn retainedName(idx: *const LoadedIndex, record_index: usize, check: diagnostics.Check, diag: *Diagnostics) diagnostics.Error![]const u8 {
    diag.record = record_index + 1;
    return idx.recordName(record_index) orelse diag.consistency(check);
}

fn writeReportChecked(writer: *std.Io.Writer, fasta_path: []const u8, extension: []const u8, size: u64, records: []const IndexRecord, summary: LengthSummary, shortest: []const u8, longest: []const u8, comp: *const CompositionStats, composition: CompositionSummary, diag: *Diagnostics) diagnostics.Error!void {
    const output = if (diag.output_source == .stdout) &diag.failure.outputs.stdout else &diag.failure.outputs.file;
    output.state = .streaming;
    output.may_be_partial = true;
    writeReport(writer, fasta_path, extension, size, records, summary, shortest, longest, comp, composition) catch |err|
        return diag.fail(.{ .@"io.write_failed" = {} }, .write, diag.output_source, null, err);
}

const STATS_READ_BUFFER_BYTES: usize = 256 * 1024;

const BoundedFastaReader = struct {
    io: std.Io,
    file: std.Io.File,
    file_size: u64,
    buffer: []u8,
    start: u64 = 0,
    len: usize = 0,
    diag: *Diagnostics,

    fn read(self: *BoundedFastaReader, offset: u64, limit: u64) diagnostics.Error![]const u8 {
        if (self.buffer.len == 0) return self.diag.consistency(.reader_buffer);
        if (limit == 0) return self.diag.consistency(.traversal_position);
        if (offset >= self.file_size) return self.diag.span(@as(u128, offset) + 1, self.file_size, null);

        const buffered_end = std.math.add(u64, self.start, self.len) catch return self.diag.consistency(.traversal_position);
        if (offset >= self.start and offset < buffered_end) {
            const relative: usize = @intCast(offset - self.start);
            const available = self.len - relative;
            const take: usize = @intCast(@min(@as(u64, available), limit));
            return self.buffer[relative..][0..take];
        }

        const wanted: usize = @intCast(@min(@as(u64, self.buffer.len), self.file_size - offset));
        const got = std.Io.File.readPositionalAll(
            self.file,
            self.io,
            self.buffer[0..wanted],
            offset,
        ) catch |err| return self.diag.read(offset, wanted, null, err);
        if (got != wanted) return self.diag.read(offset, wanted, got, null);
        self.start = offset;
        self.len = got;
        const take: usize = @intCast(@min(@as(u64, got), limit));
        return self.buffer[0..take];
    }
};

fn countCompositionFiltered(data: []const u8, comp: *CompositionStats) void {
    for (data) |byte| {
        if (byte > ' ') {
            comp.counts[byte] += 1;
            comp.total_bases += 1;
            if (byte >= 'a' and byte <= 'z') {
                comp.lowercase_count += 1;
            }
        }
    }
}

fn scanUniformRecord(reader: *BoundedFastaReader, rec: IndexRecord, comp: *CompositionStats) diagnostics.Error!void {
    if (rec.seq_len == 0) return;
    if (rec.line_bases == 0 or rec.line_bytes < rec.line_bases) return reader.diag.consistency(.traversal_position);

    const last_base = rec.seq_len - 1;
    const line = last_base / rec.line_bases;
    const column = last_base % rec.line_bases;
    const line_offset = std.math.mul(u64, line, rec.line_bytes) catch return reader.diag.consistency(.traversal_position);
    const last_offset = std.math.add(u64, line_offset, column) catch return reader.diag.consistency(.traversal_position);
    const span_len = std.math.add(u64, last_offset, 1) catch return reader.diag.consistency(.traversal_position);
    const span_end = std.math.add(u64, rec.seq_offset, span_len) catch return reader.diag.consistency(.traversal_position);
    if (span_end > reader.file_size) return reader.diag.span(span_end, reader.file_size, null);

    var offset = rec.seq_offset;
    var line_pos: usize = 0;
    while (offset < span_end) {
        const data = try reader.read(offset, span_end - offset);

        var pos: usize = 0;
        while (pos < data.len) {
            if (line_pos >= rec.line_bases) {
                const skipped = @min(data.len - pos, @as(usize, rec.line_bytes) - line_pos);
                pos += skipped;
                line_pos += skipped;
                if (line_pos == rec.line_bytes) line_pos = 0;
                continue;
            }
            const bases = @min(data.len - pos, @as(usize, rec.line_bases) - line_pos);
            try countSliceChecked(data[pos..][0..bases], comp, reader.diag);
            pos += bases;
            line_pos += bases;
            if (line_pos == rec.line_bytes) line_pos = 0;
        }
        offset += data.len;
    }
}

fn scanSideTableRecord(
    reader: *BoundedFastaReader,
    lines: []const SideTableLine,
    comp: *CompositionStats,
) diagnostics.Error!void {
    if (lines.len == 0) return reader.diag.consistency(.missing_side_table);

    var line_index: usize = 0;
    while (line_index < lines.len) {
        const start = lines[line_index].byte_offset;
        var end = std.math.add(u64, start, lines[line_index].line_bytes) catch return reader.diag.consistency(.side_table_state);
        line_index += 1;

        while (line_index < lines.len and lines[line_index].byte_offset == end) {
            const next_end = std.math.add(u64, end, lines[line_index].line_bytes) catch return reader.diag.consistency(.side_table_state);
            if (next_end - start > reader.buffer.len) break;
            end = next_end;
            line_index += 1;
        }
        if (end > reader.file_size) return reader.diag.span(end, reader.file_size, line_index);

        var offset = start;
        while (offset < end) {
            const data = try reader.read(offset, end - offset);
            try countFilteredChecked(data, comp, reader.diag);
            offset += data.len;
        }
    }
}

fn scanRecord(
    reader: *BoundedFastaReader,
    rec: IndexRecord,
    side_table: []const SideTableLine,
    comp: *CompositionStats,
) diagnostics.Error!void {
    if (rec.seq_len == 0) return;
    if (rec.isUniformWidth()) {
        try scanUniformRecord(reader, rec, comp);
    } else {
        try scanSideTableRecord(reader, side_table, comp);
    }
}

fn scanComposition(allocator: std.mem.Allocator, idx: *const LoadedIndex, buffer: []u8, comp: *CompositionStats, diag: *Diagnostics) diagnostics.Error!void {
    comp.* = .{};
    var reader = BoundedFastaReader{
        .io = idx.io,
        .file = idx.fasta_file,
        .file_size = idx.fasta_size,
        .buffer = buffer,
        .diag = diag,
    };

    if (recordsInSeqOffsetOrder(idx.records)) {
        for (idx.records, 0..) |rec, i| {
            diag.record = i + 1;
            try scanRecord(&reader, rec, idx.sideTableLines(rec), comp);
        }
    } else {
        const indices = allocator.alloc(usize, idx.records.len) catch return diag.oom(.traversal_order);
        defer allocator.free(indices);
        for (indices, 0..) |*slot, i| slot.* = i;
        std.mem.sort(usize, indices, idx.records, struct {
            fn lessThan(records: []const IndexRecord, a: usize, b: usize) bool {
                return records[a].seq_offset < records[b].seq_offset;
            }
        }.lessThan);
        for (indices) |rec_index| {
            const rec = idx.records[rec_index];
            diag.record = rec_index + 1;
            try scanRecord(&reader, rec, idx.sideTableLines(rec), comp);
        }
    }

    diag.record = null;
    comp.seq_type = detectType(&comp.counts, comp.total_bases);
}

fn recordsInSeqOffsetOrder(records: []const IndexRecord) bool {
    if (records.len < 2) return true;
    var prev = records[0].seq_offset;
    for (records[1..]) |rec| {
        if (rec.seq_offset < prev) return false;
        prev = rec.seq_offset;
    }
    return true;
}

fn countSliceChecked(data: []const u8, comp: *CompositionStats, diag: *Diagnostics) diagnostics.Error!void {
    if (data.len > std.math.maxInt(u64) - comp.total_bases) return diag.overflow(.total_symbols, comp.total_bases, data.len, std.math.maxInt(u64), .fasta);
    countCompositionSlice(data, comp);
}

fn countFilteredChecked(data: []const u8, comp: *CompositionStats, diag: *Diagnostics) diagnostics.Error!void {
    if (data.len <= std.math.maxInt(u64) - comp.total_bases) {
        countCompositionFiltered(data, comp);
        return;
    }
    for (data) |byte| {
        if (byte <= ' ') continue;
        if (comp.total_bases == std.math.maxInt(u64)) return diag.overflow(.total_symbols, comp.total_bases, 1, std.math.maxInt(u64), .fasta);
        comp.counts[byte] += 1;
        comp.total_bases += 1;
        if (byte >= 'a' and byte <= 'z') comp.lowercase_count += 1;
    }
}

fn countCompositionSlice(data: []const u8, comp: *CompositionStats) void {
    const a_vec: SimdVec = @splat('a');
    const z_vec: SimdVec = @splat('z');

    var pos: usize = 0;
    while (pos + SIMD_CHUNK_SIZE <= data.len) {
        const chunk: SimdVec = data[pos..][0..SIMD_CHUNK_SIZE].*;
        inline for (0..SIMD_CHUNK_SIZE) |j| {
            comp.counts[chunk[j]] += 1;
        }
        comp.total_bases += SIMD_CHUNK_SIZE;
        const lower_mask = (chunk >= a_vec) & (chunk <= z_vec);
        comp.lowercase_count += @popCount(@as(u32, @bitCast(lower_mask)));
        pos += SIMD_CHUNK_SIZE;
    }
    while (pos + TAIL_CHUNK_SIZE <= data.len) {
        const chunk: TailVec = data[pos..][0..TAIL_CHUNK_SIZE].*;
        inline for (0..TAIL_CHUNK_SIZE) |j| {
            comp.counts[chunk[j]] += 1;
        }
        comp.total_bases += TAIL_CHUNK_SIZE;
        const lower_mask = (chunk >= @as(TailVec, @splat('a'))) & (chunk <= @as(TailVec, @splat('z')));
        comp.lowercase_count += @popCount(@as(u8, @bitCast(lower_mask)));
        pos += TAIL_CHUNK_SIZE;
    }
    while (pos + FINAL_CHUNK_SIZE <= data.len) {
        const chunk: FinalVec = data[pos..][0..FINAL_CHUNK_SIZE].*;
        inline for (0..FINAL_CHUNK_SIZE) |j| {
            comp.counts[chunk[j]] += 1;
        }
        comp.total_bases += FINAL_CHUNK_SIZE;
        const lower_mask = (chunk >= @as(FinalVec, @splat('a'))) & (chunk <= @as(FinalVec, @splat('z')));
        comp.lowercase_count += @popCount(@as(u4, @bitCast(lower_mask)));
        pos += FINAL_CHUNK_SIZE;
    }
    while (pos < data.len) : (pos += 1) {
        const byte = data[pos];
        comp.counts[byte] += 1;
        comp.total_bases += 1;
        if (byte >= 'a' and byte <= 'z') {
            comp.lowercase_count += 1;
        }
    }
}

fn testIndexRecord(seq_offset: u64, seq_len: u64, line_bases: u32, line_bytes: u32) IndexRecord {
    return .{
        .name_offset = 0,
        .name_len = 0,
        .seq_offset = seq_offset,
        .seq_len = seq_len,
        .line_bases = line_bases,
        .line_bytes = line_bytes,
    };
}

fn checkStatsAllocationFailure(allocator: std.mem.Allocator, path: []const u8, seen: *u64) !void {
    var failure: diagnostics.Failure = undefined;
    failure.reset(.stats);
    var diag = Diagnostics{ .failure = &failure, .output_source = .output };
    var idx = index_format.loadIndexCheckedWithMode(allocator, std.testing.io, path, .stats_scan, &failure) catch |err| {
        try recordStatsAllocationFailure(err, &failure, seen);
        return error.OutOfMemory;
    };
    defer idx.deinit();
    var computed: ComputedStats = undefined;
    var buffer: [STATS_READ_BUFFER_BYTES]u8 = undefined;
    computeStats(allocator, &idx, &computed, &buffer, &diag) catch |err| {
        try recordStatsAllocationFailure(err, &failure, seen);
        return error.OutOfMemory;
    };
    var writer = std.Io.Writer.fixed(buffer[0..8192]);
    try writeComputedReport(&writer, path, &idx, &computed, &diag);
    try std.testing.expect(failure.work_complete and !failure.complete);
    try std.testing.expect(std.mem.indexOf(u8, writer.buffered(), "Composition:\n") != null);
}

fn recordStatsAllocationFailure(err: diagnostics.Error, failure: *const diagnostics.Failure, seen: *u64) !void {
    try std.testing.expectEqual(error.OutOfMemory, err);
    try std.testing.expectEqual(diagnostics.Code.@"resource.out_of_memory", failure.primary().?.code());
    const role = failure.primary().?.details.@"resource.out_of_memory".allocation_role;
    seen.* |= @as(u64, 1) << @intCast(@intFromEnum(role));
    try std.testing.expect(!failure.work_complete);
    try std.testing.expectEqual(.not_started, failure.outputs.file.state);
}

test "[unit] - [stats scan order]: detects sorted and unsorted offsets" {
    const sorted = [_]IndexRecord{
        testIndexRecord(10, 1, 1, 2),
        testIndexRecord(20, 1, 1, 2),
        testIndexRecord(30, 1, 1, 2),
    };
    try std.testing.expect(recordsInSeqOffsetOrder(&sorted));

    const unsorted = [_]IndexRecord{
        testIndexRecord(30, 1, 1, 2),
        testIndexRecord(10, 1, 1, 2),
    };
    try std.testing.expect(!recordsInSeqOffsetOrder(&unsorted));
}

test "[property] - [stats composition]: SIMD matches scalar counts at every tail length" {
    var data: [512]u8 = undefined;
    for (&data, 0..) |*byte, i| byte.* = @truncate(i);

    for (0..data.len + 1) |len| {
        var expected: [256]u64 = .{0} ** 256;
        var expected_lowercase: u64 = 0;
        for (data[0..len]) |byte| {
            expected[byte] += 1;
            if (byte >= 'a' and byte <= 'z') expected_lowercase += 1;
        }
        var comp: CompositionStats = .{};

        countCompositionSlice(data[0..len], &comp);

        try std.testing.expectEqual(@as(u64, @intCast(len)), comp.total_bases);
        try std.testing.expectEqual(expected_lowercase, comp.lowercase_count);
        try std.testing.expectEqualSlices(u64, &expected, &comp.counts);
    }
}

test "[integration] - [stats composition]: scans uniform line-ending variants" {
    var failure: diagnostics.Failure = undefined;
    failure.reset(.stats);
    var diag = Diagnostics{ .failure = &failure };
    const Case = struct {
        bytes: []const u8,
        record: IndexRecord,
    };
    const cases = [_]Case{
        .{
            .bytes = "prefixACGT\nacgt\nNN",
            .record = testIndexRecord(6, 10, 4, 5),
        },
        .{
            .bytes = "xACGT\r\nacgt\r\nNN\r\nsuffix",
            .record = testIndexRecord(1, 10, 4, 6),
        },
    };
    const buffer_sizes = [_]usize{ 1, 2, 3, 4, 5, 6, 7, 15, 16, 17, 32, 33 };

    for (cases) |case| {
        var tmp = std.testing.tmpDir(.{});
        defer tmp.cleanup();
        const file = try tmp.dir.createFile(std.testing.io, "source.fa", .{ .read = true });
        defer file.close(std.testing.io);
        try std.Io.File.writeStreamingAll(file, std.testing.io, case.bytes);
        for (buffer_sizes) |buffer_size| {
            var comp: CompositionStats = .{};
            var storage: [33]u8 = undefined;
            var reader = BoundedFastaReader{
                .io = std.testing.io,
                .file = file,
                .file_size = case.bytes.len,
                .buffer = storage[0..buffer_size],
                .diag = &diag,
            };

            try scanUniformRecord(&reader, case.record, &comp);
            try std.testing.expectEqual(@as(u64, 10), comp.total_bases);
            try std.testing.expectEqual(@as(u64, 4), comp.lowercase_count);
            try std.testing.expectEqual(@as(u64, 1), comp.counts['A']);
            try std.testing.expectEqual(@as(u64, 1), comp.counts['a']);
            try std.testing.expectEqual(@as(u64, 2), comp.counts['N']);
        }
    }
}

test "[integration] - [stats composition]: excludes bytes outside indexed spans" {
    var failure: diagnostics.Failure = undefined;
    failure.reset(.stats);
    var diag = Diagnostics{ .failure = &failure };
    const bytes = ">one\nAC\nGT\n>two\nTT\n";
    const buffer_sizes = [_]usize{ 3, 32 };
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    const file = try tmp.dir.createFile(std.testing.io, "source.fa", .{ .read = true });
    defer file.close(std.testing.io);
    try std.Io.File.writeStreamingAll(file, std.testing.io, bytes);
    for (buffer_sizes) |buffer_size| {
        var storage: [32]u8 = undefined;
        var comp: CompositionStats = .{};
        var reader = BoundedFastaReader{
            .io = std.testing.io,
            .file = file,
            .file_size = bytes.len,
            .buffer = storage[0..buffer_size],
            .diag = &diag,
        };

        try scanUniformRecord(
            &reader,
            testIndexRecord(5, 4, 2, 3),
            &comp,
        );
        try scanUniformRecord(
            &reader,
            testIndexRecord(16, 2, 2, 3),
            &comp,
        );

        try std.testing.expectEqual(@as(u64, 6), comp.total_bases);
        try std.testing.expectEqual(@as(u64, 1), comp.counts['A']);
        try std.testing.expectEqual(@as(u64, 1), comp.counts['C']);
        try std.testing.expectEqual(@as(u64, 1), comp.counts['G']);
        try std.testing.expectEqual(@as(u64, 3), comp.counts['T']);
        try std.testing.expectEqual(@as(u64, 0), comp.counts['>']);
    }
}

test "[integration] - [stats composition]: coalesces only adjacent side-table spans" {
    var failure: diagnostics.Failure = undefined;
    failure.reset(.stats);
    var diag = Diagnostics{ .failure = &failure };
    const bytes = "AC\nGT\nignored\nnn\r\n";
    const lines = [_]SideTableLine{
        .{ .base_start = 0, .byte_offset = 0, .line_bytes = 3, .line_bases = 2 },
        .{ .base_start = 2, .byte_offset = 3, .line_bytes = 3, .line_bases = 2 },
        .{ .base_start = 4, .byte_offset = 14, .line_bytes = 4, .line_bases = 2 },
    };
    const buffer_sizes = [_]usize{ 1, 2, 3, 4, 5, 6, 7 };
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    const file = try tmp.dir.createFile(std.testing.io, "source.fa", .{ .read = true });
    defer file.close(std.testing.io);
    try std.Io.File.writeStreamingAll(file, std.testing.io, bytes);
    for (buffer_sizes) |buffer_size| {
        var comp: CompositionStats = .{};
        var storage: [7]u8 = undefined;
        var reader = BoundedFastaReader{
            .io = std.testing.io,
            .file = file,
            .file_size = bytes.len,
            .buffer = storage[0..buffer_size],
            .diag = &diag,
        };

        try scanSideTableRecord(&reader, &lines, &comp);
        try std.testing.expectEqual(@as(u64, 6), comp.total_bases);
        try std.testing.expectEqual(@as(u64, 2), comp.lowercase_count);
        try std.testing.expectEqual(@as(u64, 0), comp.counts['i']);
        try std.testing.expectEqual(@as(u64, 2), comp.counts['n']);
    }
}

test "[failure] - [stats composition]: rejects invalid uniform geometry" {
    var failure: diagnostics.Failure = undefined;
    failure.reset(.stats);
    var diag = Diagnostics{ .failure = &failure };
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    const file = try tmp.dir.createFile(std.testing.io, "source.fa", .{});
    defer file.close(std.testing.io);
    try std.Io.File.writeStreamingAll(file, std.testing.io, "ACGT");
    var storage: [2]u8 = undefined;
    var comp: CompositionStats = .{};
    var reader = BoundedFastaReader{
        .io = std.testing.io,
        .file = file,
        .file_size = 4,
        .buffer = &storage,
        .diag = &diag,
    };
    const invalid = [_]IndexRecord{
        testIndexRecord(0, 1, 0, 0),
        testIndexRecord(0, 1, 2, 1),
        testIndexRecord(0, 5, 5, 6),
        testIndexRecord(std.math.maxInt(u64), 1, 1, 1),
        testIndexRecord(0, std.math.maxInt(u64), 1, std.math.maxInt(u32)),
    };

    for (invalid) |record| {
        failure.reset(.stats);
        const expected: diagnostics.Error = if (record.seq_offset == 0 and record.seq_len == 5) error.InvalidIndex else error.Internal;
        try std.testing.expectError(expected, scanUniformRecord(&reader, record, &comp));
        const primary = failure.primary().?;
        if (expected == error.InvalidIndex) {
            const detail = primary.details.@"index.sequence_span_out_of_bounds";
            try std.testing.expectEqual(@as(u128, 5), detail.required_end);
            try std.testing.expectEqual(@as(u64, 4), detail.source_bytes);
            try std.testing.expectEqual(diagnostics.BoundSource.actual, detail.bound_source);
        } else {
            try std.testing.expectEqual(diagnostics.Check.traversal_position, primary.details.@"internal.consistency_failed".check);
        }
    }
}

test "[failure] - [stats composition]: propagates positional read failure" {
    var failure: diagnostics.Failure = undefined;
    failure.reset(.stats);
    var diag = Diagnostics{ .failure = &failure };
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    const file = try tmp.dir.createFile(std.testing.io, "source.fa", .{ .read = true });
    std.Io.File.writeStreamingAll(file, std.testing.io, "A") catch |err| {
        file.close(std.testing.io);
        return err;
    };
    file.close(std.testing.io);

    var storage: [1]u8 = undefined;
    var reader = BoundedFastaReader{
        .io = std.testing.io,
        .file = file,
        .file_size = 1,
        .buffer = &storage,
        .diag = &diag,
    };

    try std.testing.expectError(error.Io, reader.read(0, 1));
}

test "[edge] - [stats lengths]: computes a near-maximum median without overflow" {
    const near_max = [_]u64{ std.math.maxInt(u64) - 1, std.math.maxInt(u64) };

    try std.testing.expectEqual(std.math.maxInt(u64) - 1, medianSorted(&near_max));
}

test "[failure] - [stats lengths]: rejects total-symbol overflow" {
    var failure: diagnostics.Failure = undefined;
    failure.reset(.stats);
    var diag = Diagnostics{ .failure = &failure };
    const records = [_]IndexRecord{
        testIndexRecord(0, std.math.maxInt(u64), 1, 2),
        testIndexRecord(0, 1, 1, 2),
    };
    var lengths: [2]u64 = undefined;

    try std.testing.expectError(error.CountOverflow, summarizeLengths(&records, &lengths, &diag));
    const primary = failure.primary().?;
    const detail = primary.details.@"resource.count_overflow";
    try std.testing.expectEqual(diagnostics.Counter.total_symbols, detail.counter);
    try std.testing.expectEqual(@as(?u128, std.math.maxInt(u64)), detail.current);
    try std.testing.expectEqual(@as(?u128, 1), detail.increment);
    try std.testing.expectEqual(@as(?u64, 2), primary.at.record);
    try std.testing.expectEqual(@as(?diagnostics.SourceId, .zfi), primary.at.source);
}

test "[edge] - [stats lengths]: accumulates and formats auN above u64" {
    var failure: diagnostics.Failure = undefined;
    failure.reset(.stats);
    var diag = Diagnostics{ .failure = &failure };
    const length: u64 = 1 << 32;
    const records = [_]IndexRecord{
        testIndexRecord(0, length, 1, 2),
        testIndexRecord(0, length, 1, 2),
    };
    var lengths: [2]u64 = undefined;
    const summary = try summarizeLengths(&records, &lengths, &diag);
    try std.testing.expectEqual(@as(u64, 1) << 33, summary.total_symbols);
    try std.testing.expectEqual(@as(u128, 1) << 65, summary.aun_numerator);

    var buffer: [32]u8 = undefined;
    var writer = std.Io.Writer.fixed(&buffer);
    try writeFixedUnsigned(&writer, summary.aun_numerator, summary.total_symbols, 2);
    try std.testing.expectEqualStrings("4294967296.00", writer.buffered());
}

test "[unit] - [stats formatting]: rounds fixed decimals exactly" {
    var buffer: [64]u8 = undefined;
    var writer = std.Io.Writer.fixed(&buffer);
    try writeFixedUnsigned(&writer, 100, 32, 2);
    try writer.writeByte(' ');
    try writeFixedUnsigned(&writer, 200, 3, 2);
    try writer.writeByte(' ');
    try writeFixedUnsigned(&writer, 226, 16, 2);
    try writer.writeByte(' ');
    try writeFixedSigned(&writer, -30, 32, 3);
    try writer.writeByte(' ');
    try writeFixedSigned(&writer, -1, 10_000, 3);

    try std.testing.expectEqualStrings("3.13 66.67 14.13 -0.938 0.000", writer.buffered());
}

test "[failure] - [stats report]: propagates a full destination" {
    var failure: diagnostics.Failure = undefined;
    failure.reset(.stats);
    var diag = Diagnostics{ .failure = &failure };
    const records = [_]IndexRecord{
        testIndexRecord(0, 4, 4, 5),
    };
    const summary = LengthSummary{
        .total_symbols = 4,
        .shortest_index = 0,
        .longest_index = 0,
        .mean = 4,
        .q1 = 4,
        .median = 4,
        .q3 = 4,
        .range = 0,
        .n50 = 4,
        .l50 = 1,
        .n90 = 4,
        .l90 = 1,
        .aun_numerator = 16,
    };
    var comp: CompositionStats = .{ .total_bases = 4 };
    comp.counts['A'] = 4;
    var buffer: [1]u8 = undefined;
    var writer = std.Io.Writer.fixed(&buffer);

    try std.testing.expectError(
        error.Io,
        writeReportChecked(
            &writer,
            "x.fa",
            ".zfi",
            8,
            &records,
            summary,
            "x",
            "x",
            &comp,
            try summarizeComposition(&comp, &diag),
            &diag,
        ),
    );
    try std.testing.expectEqual(diagnostics.Code.@"io.write_failed", failure.primary().?.code());
    try std.testing.expectEqual(error.WriteFailed, failure.primary().?.cause.?.zig);
    try std.testing.expect(failure.outputs.stdout.may_be_partial);
    try std.testing.expectEqual(.not_attempted, failure.outputs.stdout.final_flush);
}

test "[failure] - [stats composition]: rejects categories above total symbols" {
    var failure: diagnostics.Failure = undefined;
    failure.reset(.stats);
    var diag = Diagnostics{ .failure = &failure };
    const sequence_types = [_]SequenceType{ .nucleotide, .protein };
    for (sequence_types) |sequence_type| {
        var comp: CompositionStats = .{
            .total_bases = 4,
            .seq_type = sequence_type,
        };
        comp.counts['A'] = 5;

        try std.testing.expectError(error.Internal, summarizeComposition(&comp, &diag));
    }
}

test "[failure] - [stats reader]: retains read causes, short counts and original record" {
    const Fault = struct {
        fn read(_: ?*anyopaque, _: std.Io.File, _: []const []u8, _: u64) std.Io.File.ReadPositionalError!usize {
            return error.InputOutput;
        }
        fn short(_: ?*anyopaque, _: std.Io.File, _: []const []u8, _: u64) std.Io.File.ReadPositionalError!usize {
            return 0;
        }
    };
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    const file = try tmp.dir.createFile(std.testing.io, "reader.fa", .{ .read = true });
    defer file.close(std.testing.io);
    for (0..2) |i| {
        var failure: diagnostics.Failure = undefined;
        failure.reset(.stats);
        var diag = Diagnostics{ .failure = &failure, .record = 7 };
        var vtable = std.testing.io.vtable.*;
        vtable.fileReadPositional = if (i == 0) Fault.read else Fault.short;
        var storage: [8]u8 = undefined;
        var reader = BoundedFastaReader{ .io = .{ .userdata = std.testing.io.userdata, .vtable = &vtable }, .file = file, .file_size = 16, .buffer = &storage, .diag = &diag };
        try std.testing.expectError(error.Io, reader.read(3, 2));
        const primary = failure.primary().?;
        try std.testing.expectEqual(@as(?u64, 3), primary.at.location.offset);
        try std.testing.expectEqual(@as(?u64, 7), primary.at.record);
        try std.testing.expectEqual(@as(?diagnostics.SourceId, .fasta), primary.at.source);
        if (i == 0) {
            try std.testing.expectEqual(diagnostics.Code.@"io.read_failed", primary.code());
            try std.testing.expectEqual(error.InputOutput, primary.cause.?.zig);
            try std.testing.expectEqual(@as(?u64, 8), primary.details.@"io.read_failed".requested_bytes);
        } else {
            try std.testing.expectEqual(diagnostics.Code.@"io.unexpected_eof", primary.code());
            const detail = primary.details.@"io.unexpected_eof".exact;
            try std.testing.expectEqual(@as(u64, 8), detail.requested_bytes);
            try std.testing.expectEqual(@as(u64, 0), detail.received_bytes);
            try std.testing.expect(primary.cause == null);
        }
    }
}

test "[failure] - [stats counts]: guards slices and accepts whitespace near the limit" {
    var failure: diagnostics.Failure = undefined;
    failure.reset(.stats);
    var diag = Diagnostics{ .failure = &failure, .record = 3 };
    var comp = CompositionStats{ .total_bases = std.math.maxInt(u64) - 1 };
    comp.counts['A'] = comp.total_bases;
    try std.testing.expectError(error.CountOverflow, countSliceChecked("AA", &comp, &diag));
    try std.testing.expectEqual(@as(u64, std.math.maxInt(u64) - 1), comp.total_bases);
    try std.testing.expectEqual(@as(?u128, 2), failure.primary().?.details.@"resource.count_overflow".increment);
    try std.testing.expectEqual(@as(?diagnostics.SourceId, .fasta), failure.primary().?.at.source);
    failure.reset(.stats);
    try countFilteredChecked(" A ", &comp, &diag);
    try std.testing.expectEqual(std.math.maxInt(u64), comp.total_bases);
    try countFilteredChecked(" \r\n", &comp, &diag);
    try std.testing.expectError(error.CountOverflow, countFilteredChecked(" a", &comp, &diag));
    try std.testing.expectEqual(@as(?u128, 1), failure.primary().?.details.@"resource.count_overflow".increment);
    try std.testing.expectEqual(@as(?u64, 3), failure.primary().?.at.record);
    failure.reset(.stats);
    comp.counts['a'] = 1;
    try std.testing.expectError(error.Internal, summarizeComposition(&comp, &diag));
    try std.testing.expectEqual(diagnostics.Check.composition_accounting, failure.primary().?.details.@"internal.consistency_failed".check);
}

test "[failure] - [stats names and aggregate]: distinguishes extrema and completed totals" {
    var failure: diagnostics.Failure = undefined;
    failure.reset(.stats);
    var idx = try index_format.loadIndexCheckedWithMode(std.testing.allocator, std.testing.io, "tests/data/simple.fasta", .stats_scan, &failure);
    defer idx.deinit();
    var diag = Diagnostics{ .failure = &failure };
    idx.stats_name_indices = .{ std.math.maxInt(usize), std.math.maxInt(usize) };
    for ([_]diagnostics.Check{ .retained_shortest_name, .retained_longest_name }) |check| {
        failure.reset(.stats);
        try std.testing.expectError(error.Internal, retainedName(&idx, 0, check, &diag));
        try std.testing.expectEqual(check, failure.primary().?.details.@"internal.consistency_failed".check);
        try std.testing.expectEqual(@as(?u64, 1), failure.primary().?.at.record);
    }
    failure.reset(.stats);
    diag.record = null;
    try std.testing.expectError(error.InvalidIndex, checkTotal(10, 9, &diag));
    const detail = failure.primary().?.details.@"pair.total_length_mismatch";
    try std.testing.expectEqual(@as(u64, 10), detail.indexed_total);
    try std.testing.expectEqual(@as(u64, 9), detail.observed_total);
    try std.testing.expect(failure.primary().?.at.record == null);
}

test "[failure] - [stats allocations]: releases required lengths and unsorted traversal storage" {
    var seen: u64 = 0;
    try std.testing.checkAllAllocationFailures(std.testing.allocator, checkStatsAllocationFailure, .{ "tests/data/simple.fasta", &seen });
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    const directory = try tmp.dir.realPathFileAlloc(std.testing.io, ".", std.testing.allocator);
    defer std.testing.allocator.free(directory);
    const path = try std.fs.path.join(std.testing.allocator, &.{ directory, "unsorted.fa" });
    defer std.testing.allocator.free(path);
    const file = try tmp.dir.createFile(std.testing.io, "unsorted.fa", .{});
    try std.Io.File.writeStreamingAll(file, std.testing.io, ">one\nA\n>two\nC\n");
    file.close(std.testing.io);
    const fai = try tmp.dir.createFile(std.testing.io, "unsorted.fa.fai", .{});
    try std.Io.File.writeStreamingAll(fai, std.testing.io, "two\t1\t12\t1\t2\none\t1\t5\t1\t2\n");
    fai.close(std.testing.io);
    try std.testing.checkAllAllocationFailures(std.testing.allocator, checkStatsAllocationFailure, .{ path, &seen });
    try std.testing.expect(seen & (@as(u64, 1) << @intFromEnum(diagnostics.AllocationRole.lengths)) != 0);
    try std.testing.expect(seen & (@as(u64, 1) << @intFromEnum(diagnostics.AllocationRole.traversal_order)) != 0);
}
