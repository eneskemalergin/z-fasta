//! Assembly and proteome statistics from indexed FASTA files.

const std = @import("std");
const index_format = @import("index_format.zig");
const platform = @import("platform.zig");

const IndexRecord = index_format.IndexRecord;
const LoadedIndex = index_format.LoadedIndex;
const printErrorAndExit = index_format.printErrorAndExit;

const SIMD_CHUNK_SIZE = 32;
const SimdVec = @Vector(SIMD_CHUNK_SIZE, u8);

// ============================================================================
// Stats computation
// ============================================================================

pub const SequenceType = enum { nucleotide, protein };

const CompositionStats = struct {
    counts: [256]u64,
    total_bases: u64,
    seq_type: SequenceType,
    lowercase_count: u64,
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
    categories: [protein_fields.len]u64,
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

const protein_fields = [_]ProteinField{
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

fn summarizeLengths(records: []const IndexRecord, lengths: []u64) !LengthSummary {
    var total: u64 = 0;
    var shortest_index: usize = 0;
    var longest_index: usize = 0;
    var aun_numerator: u128 = 0;
    for (records, 0..) |record, i| {
        total = try std.math.add(u64, total, record.seq_len);
        lengths[i] = record.seq_len;
        if (record.seq_len < records[shortest_index].seq_len) shortest_index = i;
        if (record.seq_len > records[longest_index].seq_len) longest_index = i;
        const square = @as(u128, record.seq_len) * @as(u128, record.seq_len);
        aun_numerator = try std.math.add(u128, aun_numerator, square);
    }

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

fn writeFixedUnsigned(writer: *std.Io.Writer, numerator: u128, denominator: u128, decimals: u8) !void {
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

fn writeFixedSigned(writer: *std.Io.Writer, numerator: i128, denominator: u128, decimals: u8) !void {
    const magnitude: u128 = @intCast(if (numerator < 0) -numerator else numerator);
    const scale: u128 = if (decimals == 3) 1000 else 100;
    const rounded = (magnitude % denominator * scale + denominator / 2) / denominator;
    if (numerator < 0 and (magnitude / denominator != 0 or rounded != 0)) try writer.writeByte('-');
    try writeFixedUnsigned(writer, magnitude, denominator, decimals);
}

fn countLetter(comp: *const CompositionStats, upper: u8) u64 {
    return comp.counts[upper] + comp.counts[upper + ('a' - 'A')];
}

fn summarizeComposition(comp: *const CompositionStats) !CompositionSummary {
    if (comp.seq_type == .nucleotide) {
        const a = countLetter(comp, 'A');
        const c = countLetter(comp, 'C');
        const g = countLetter(comp, 'G');
        const t = countLetter(comp, 'T');
        const u = countLetter(comp, 'U');
        const n = countLetter(comp, 'N');
        var ambiguous: u64 = 0;
        for ("RYSWKMBDHV") |code| ambiguous += countLetter(comp, code);
        const assigned = @as(u128, a) + c + g + t + u + n + ambiguous;
        if (assigned > comp.total_bases) return error.CompositionMismatch;
        return .{ .nucleotide = .{
            .type_name = if (t > 0 and u > 0)
                "nucleotide_mixed_tu"
            else if (t > 0)
                "nucleotide_t"
            else if (u > 0)
                "nucleotide_u"
            else
                "nucleotide",
            .categories = .{ a, c, g, t, u, n, ambiguous, @intCast(@as(u128, comp.total_bases) - assigned) },
            .canonical = a + c + g + t + u,
            .gc_total = g + c,
        } };
    }

    var categories: [protein_fields.len]u64 = undefined;
    var assigned: u128 = 0;
    for (protein_fields, 0..) |field, i| {
        categories[i] = countLetter(comp, field.code);
        assigned += categories[i];
    }
    const stop = comp.counts['*'];
    assigned += stop;
    if (assigned > comp.total_bases) return error.CompositionMismatch;
    return .{ .protein = .{
        .categories = categories,
        .stop = stop,
        .invalid = @intCast(@as(u128, comp.total_bases) - assigned),
    } };
}

fn writeCountPercent(writer: *std.Io.Writer, key: []const u8, count: u64, total: u64) !void {
    try writer.print("  {s}: {d} ", .{ key, count });
    try writeFixedUnsigned(writer, @as(u128, count) * 100, total, 2);
    try writer.writeAll("%\n");
}

fn writeNucleotide(writer: *std.Io.Writer, summary: NucleotideSummary, comp: *const CompositionStats) !void {
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

fn writeProtein(writer: *std.Io.Writer, summary: ProteinSummary, comp: *const CompositionStats) !void {
    try writer.writeAll("  type: protein\n  percent_denominator: total_symbols\n");
    for (protein_fields, summary.categories) |field, count| {
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
) !void {
    try writer.print(
        "File:\n  path: {s}\n  index: {s}{s}\n  size_bytes: {d}\n\n" ++
            "Lengths:\n  indexed_records: {d}\n  total_symbols: {d}\n" ++
            "  shortest_length: {d}\n  shortest_name: {s}\n" ++
            "  longest_length: {d}\n  longest_name: {s}\n" ++
            "  mean: {d}\n  q1: {d}\n  median: {d}\n  q3: {d}\n  range: {d}\n\n" ++
            "Nx:\n  n50: {d}\n  l50: {d}\n  n90: {d}\n  l90: {d}\n  aun: ",
        .{
            fasta_path,                              fasta_path,    index_extension,                        fasta_size,   records.len,  summary.total_symbols,
            records[summary.shortest_index].seq_len, shortest_name, records[summary.longest_index].seq_len, longest_name, summary.mean, summary.q1,
            summary.median,                          summary.q3,    summary.range,                          summary.n50,  summary.l50,  summary.n90,
            summary.l90,
        },
    );
    try writeFixedUnsigned(writer, summary.aun_numerator, summary.total_symbols, 2);
    try writer.writeAll("\n\nComposition:\n");
    switch (composition) {
        .nucleotide => |nucleotide| try writeNucleotide(writer, nucleotide, comp),
        .protein => |protein| try writeProtein(writer, protein, comp),
    }
}

/// Run the stats command.
pub fn runStats(io: std.Io, fasta_path: []const u8) void {
    var idx = index_format.loadIndexWithMode(io, fasta_path, .stats_scan);
    defer idx.deinit(io);

    const records = idx.records;
    const num_seqs = records.len;

    if (num_seqs == 0) {
        printErrorAndExit("error: no sequences in index\n", .{});
    }

    var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const lengths = allocator.alloc(u64, num_seqs) catch {
        printErrorAndExit("error: out of memory\n", .{});
    };
    const summary = summarizeLengths(records, lengths) catch {
        printErrorAndExit("error: sequence length overflow\n", .{});
    };
    const shortest_name = idx.getRecordNameWithIo(io, summary.shortest_index);
    const longest_name = idx.getRecordNameWithIo(io, summary.longest_index);
    const comp = scanComposition(&idx);
    if (comp.total_bases != summary.total_symbols) {
        printErrorAndExit("error: index sequence lengths do not match FASTA data\n", .{});
    }
    const composition = summarizeComposition(&comp) catch {
        printErrorAndExit("error: composition does not match indexed symbols\n", .{});
    };

    var out_buf: [65536]u8 = undefined;
    var stdout_fw = std.Io.File.Writer.initStreaming(.stdout(), io, &out_buf);
    const writer = &stdout_fw.interface;

    const index_ext: []const u8 = switch (idx.source) {
        .zfi => ".zfi",
        .fai => ".fai",
    };

    writeReport(writer, fasta_path, index_ext, idx.fasta_size, records, summary, shortest_name, longest_name, &comp, composition) catch {
        printErrorAndExit("error: write failed\n", .{});
    };
    stdout_fw.flush() catch {
        printErrorAndExit("error: write failed\n", .{});
    };
}

const fasta_release_batch_bytes: usize = 8 * 1024 * 1024;

/// Release FASTA mmap pages during sequential composition scan; batches `madvise` to limit syscall overhead.
const CompositionReleaseCursor = struct {
    released_end: usize = 0,
    scan_end: usize = 0,

    fn beforeRegion(self: *CompositionReleaseCursor, fasta: []const u8, start_byte: usize) void {
        if (start_byte <= self.released_end) return;
        const pending = start_byte - self.released_end;
        if (pending >= fasta_release_batch_bytes) {
            index_format.dropFastaSpan(fasta, self.released_end, start_byte);
            self.released_end = start_byte;
        }
    }

    fn afterRegion(self: *CompositionReleaseCursor, span_end: usize) void {
        if (span_end > self.scan_end) self.scan_end = span_end;
    }

    fn flush(self: *CompositionReleaseCursor, fasta: []const u8) void {
        const drop_end = @max(self.scan_end, self.released_end);
        if (drop_end > self.released_end) {
            index_format.dropFastaSpan(fasta, self.released_end, drop_end);
            self.released_end = drop_end;
        }
    }
};

fn scanCompositionBytes(
    fasta: []const u8,
    start: usize,
    end: usize,
    counts: *[256]u64,
    total_bases: *u64,
    lowercase_count: *u64,
) void {
    for (fasta[start..end]) |byte| {
        if (byte > ' ') {
            counts[byte] += 1;
            total_bases.* += 1;
            if (byte >= 'a' and byte <= 'z') {
                lowercase_count.* += 1;
            }
        }
    }
}

fn uniformRecordByteSpan(rec: IndexRecord) usize {
    const full_lines = rec.seq_len / rec.line_bases;
    const remainder = rec.seq_len % rec.line_bases;
    return @intCast(full_lines * rec.line_bytes + remainder);
}

fn compositionRegionEnd(fasta: []const u8, rec: IndexRecord, start: usize) usize {
    const full_lines = rec.seq_len / rec.line_bases;
    const remainder = rec.seq_len % rec.line_bases;
    var end: usize = start;
    if (remainder > 0) {
        end = start + (full_lines * rec.line_bytes) + remainder;
        if (end < fasta.len and (fasta[end] == '\n' or fasta[end] == '\r')) {
            end += 1;
            if (end < fasta.len and fasta[end - 1] == '\r' and fasta[end] == '\n') {
                end += 1;
            }
        }
    } else {
        end = start + (full_lines * rec.line_bytes);
    }
    return @min(end, fasta.len);
}

fn scanCompositionRecord(
    idx: *const LoadedIndex,
    rec: IndexRecord,
    counts: *[256]u64,
    total_bases: *u64,
    lowercase_count: *u64,
    release: *CompositionReleaseCursor,
) void {
    const fasta = idx.fasta_data;
    if (rec.seq_len == 0) return;

    if (!rec.isUniformWidth()) {
        const side_table = idx.sideTableLines(rec);
        var record_end: usize = 0;
        for (side_table) |line| {
            const start: usize = @intCast(line.byte_offset);
            const end = @min(fasta.len, start + @as(usize, @intCast(line.line_bytes)));
            release.beforeRegion(fasta, start);
            scanCompositionBytes(fasta, start, end, counts, total_bases, lowercase_count);
            if (end > record_end) record_end = end;
        }
        release.afterRegion(record_end);
        return;
    }

    const start: usize = @intCast(rec.seq_offset);
    release.beforeRegion(fasta, start);

    if (tryScanFixedWidthRecord(fasta, rec, start, counts, total_bases, lowercase_count)) {
        release.afterRegion(compositionRegionEnd(fasta, rec, start));
        return;
    }

    const end = compositionRegionEnd(fasta, rec, start);
    scanCompositionBytes(fasta, start, end, counts, total_bases, lowercase_count);
    release.afterRegion(end);
}

/// Scan all sequence regions in the FASTA for composition.
fn scanComposition(idx: *const LoadedIndex) CompositionStats {
    const fasta = idx.fasta_data;

    // Sequential access hint for full scan (no-op on Windows).
    platform.advise(fasta, .sequential);

    var counts: [256]u64 = .{0} ** 256;
    var lowercase_count: u64 = 0;
    var total_bases: u64 = 0;

    const records = idx.records;

    var release = CompositionReleaseCursor{};

    if (recordsInSeqOffsetOrder(records)) {
        var scan_end: usize = 0;
        for (records) |rec| {
            if (rec.seq_len == 0) continue;

            if (!rec.isUniformWidth()) {
                scanCompositionRecord(idx, rec, &counts, &total_bases, &lowercase_count, &release);
                if (release.scan_end > scan_end) scan_end = release.scan_end;
                continue;
            }

            const start: usize = @intCast(rec.seq_offset);
            if (start >= release.released_end + fasta_release_batch_bytes) {
                release.beforeRegion(fasta, start);
            }

            if (tryScanFixedWidthRecord(fasta, rec, start, &counts, &total_bases, &lowercase_count)) {
                const end = start + uniformRecordByteSpan(rec);
                if (end > scan_end) scan_end = end;
                continue;
            }

            const end = compositionRegionEnd(fasta, rec, start);
            scanCompositionBytes(fasta, start, end, &counts, &total_bases, &lowercase_count);
            if (end > scan_end) scan_end = end;
        }
        release.scan_end = scan_end;
    } else {
        // `.fai` rows are usually name-sorted, not file order; must sort before page release.
        const sorted_indices = std.heap.page_allocator.alloc(usize, records.len) catch {
            printErrorAndExit("error: out of memory\n", .{});
        };
        defer std.heap.page_allocator.free(sorted_indices);
        for (sorted_indices, 0..) |*slot, i| slot.* = i;
        std.mem.sort(usize, sorted_indices, records, struct {
            fn lessThan(ctx: []const IndexRecord, a: usize, b: usize) bool {
                return ctx[a].seq_offset < ctx[b].seq_offset;
            }
        }.lessThan);

        for (sorted_indices) |rec_idx| {
            scanCompositionRecord(
                idx,
                records[rec_idx],
                &counts,
                &total_bases,
                &lowercase_count,
                &release,
            );
        }
    }
    release.flush(fasta);

    const seq_type = detectType(&counts, total_bases);

    return CompositionStats{
        .counts = counts,
        .total_bases = total_bases,
        .seq_type = seq_type,
        .lowercase_count = lowercase_count,
    };
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

fn tryScanFixedWidthRecord(
    fasta: []const u8,
    rec: IndexRecord,
    start: usize,
    counts: *[256]u64,
    total_bases: *u64,
    lowercase_count: *u64,
) bool {
    if (rec.seq_len == 0 or rec.line_bases == 0) return true;
    if (rec.line_bytes < rec.line_bases) return false;

    const newline_bytes = rec.line_bytes - rec.line_bases;
    if (newline_bytes != 1 and newline_bytes != 2) return false;

    var pos = start;
    var bases_remaining = rec.seq_len;
    while (bases_remaining > 0) {
        const bases_this_line: usize = @intCast(@min(bases_remaining, @as(u64, rec.line_bases)));
        const line_end = pos + bases_this_line;
        if (line_end > fasta.len) return false;

        countCompositionSlice(fasta[pos..line_end], counts, total_bases, lowercase_count);

        bases_remaining -= bases_this_line;
        pos += if (bases_remaining > 0) @as(usize, @intCast(rec.line_bytes)) else bases_this_line;
    }

    return true;
}

fn countCompositionSlice(
    data: []const u8,
    counts: *[256]u64,
    total_bases: *u64,
    lowercase_count: *u64,
) void {
    const a_vec: SimdVec = @splat('a');
    const z_vec: SimdVec = @splat('z');

    var pos: usize = 0;
    while (pos + SIMD_CHUNK_SIZE <= data.len) {
        const chunk: SimdVec = data[pos..][0..SIMD_CHUNK_SIZE].*;
        inline for (0..SIMD_CHUNK_SIZE) |j| {
            counts[chunk[j]] += 1;
        }
        total_bases.* += SIMD_CHUNK_SIZE;
        const lower_mask = (chunk >= a_vec) & (chunk <= z_vec);
        lowercase_count.* += @popCount(@as(u32, @bitCast(lower_mask)));
        pos += SIMD_CHUNK_SIZE;
    }
    while (pos < data.len) : (pos += 1) {
        const byte = data[pos];
        counts[byte] += 1;
        total_bases.* += 1;
        if (byte >= 'a' and byte <= 'z') {
            lowercase_count.* += 1;
        }
    }
}

/// Detect if sequences are nucleotide or protein from composition counts.
///
/// Shared by stats (full-file counts), GET (per-record sample), and validate
/// (prefix sample). Threshold: nucleotide if IUPAC nuc letters are strictly
/// more than 90% of counted bases (`nuc * 10 > total * 9`).
pub fn detectType(counts: *const [256]u64, total: u64) SequenceType {
    if (total == 0) return .nucleotide;

    // Sum the full IUPAC nucleotide alphabet (case-insensitive):
    // A C G T U R Y S W K M B D H V N
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

    // u128 avoids overflow on genome-scale totals (total * 9 and nuc * 10).
    if (@as(u128, nuc_count) * 10 > @as(u128, total) * 9) {
        return .nucleotide;
    }
    return .protein;
}

/// Bases sampled per record for GET `--rc` / `--complement-only` protein guard.
pub const get_type_sample_bases: u64 = 256;

/// Bases sampled across the file for validate alphabet selection.
pub const validate_type_sample_bases: u64 = 100_000;

test "recordsInSeqOffsetOrder detects sorted and unsorted offsets" {
    const sorted = [_]IndexRecord{
        .{ .seq_offset = 10, .seq_len = 1, .line_bases = 1, .line_bytes = 2 },
        .{ .seq_offset = 20, .seq_len = 1, .line_bases = 1, .line_bytes = 2 },
        .{ .seq_offset = 30, .seq_len = 1, .line_bases = 1, .line_bytes = 2 },
    };
    try std.testing.expect(recordsInSeqOffsetOrder(&sorted));

    const unsorted = [_]IndexRecord{
        .{ .seq_offset = 30, .seq_len = 1, .line_bases = 1, .line_bytes = 2 },
        .{ .seq_offset = 10, .seq_len = 1, .line_bases = 1, .line_bytes = 2 },
    };
    try std.testing.expect(!recordsInSeqOffsetOrder(&unsorted));
}

test "countCompositionSlice tallies composition and lowercase" {
    var counts: [256]u64 = .{0} ** 256;
    var total: u64 = 0;
    var lowercase: u64 = 0;
    countCompositionSlice("ACGTacgtNN", &counts, &total, &lowercase);
    try std.testing.expectEqual(@as(u64, 10), total);
    try std.testing.expectEqual(@as(u64, 4), lowercase);
    try std.testing.expectEqual(@as(u64, 1), counts['A']);
    try std.testing.expectEqual(@as(u64, 1), counts['a']);
    try std.testing.expectEqual(@as(u64, 2), counts['N']);
}

test "median arithmetic is overflow safe" {
    const near_max = [_]u64{ std.math.maxInt(u64) - 1, std.math.maxInt(u64) };

    try std.testing.expectEqual(std.math.maxInt(u64) - 1, medianSorted(&near_max));
}

test "length total overflow is rejected" {
    const records = [_]IndexRecord{
        .{ .seq_len = std.math.maxInt(u64), .line_bases = 1, .line_bytes = 2 },
        .{ .seq_len = 1, .line_bases = 1, .line_bytes = 2 },
    };
    var lengths: [2]u64 = undefined;

    try std.testing.expectError(error.Overflow, summarizeLengths(&records, &lengths));
}

test "fixed decimals use exact half-away rounding" {
    var buffer: [64]u8 = undefined;
    var writer = std.Io.Writer.fixed(&buffer);
    try writeFixedUnsigned(&writer, 100, 32, 2);
    try writer.writeByte(' ');
    try writeFixedUnsigned(&writer, 226, 16, 2);
    try writer.writeByte(' ');
    try writeFixedSigned(&writer, -30, 32, 3);
    try writer.writeByte(' ');
    try writeFixedSigned(&writer, -1, 10_000, 3);

    try std.testing.expectEqualStrings("3.13 14.13 -0.938 0.000", writer.buffered());
}

test "report writing propagates a full destination" {
    const records = [_]IndexRecord{
        .{ .seq_len = 4, .line_bases = 4, .line_bytes = 5 },
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
    var comp = CompositionStats{
        .counts = .{0} ** 256,
        .total_bases = 4,
        .seq_type = .nucleotide,
        .lowercase_count = 0,
    };
    comp.counts['A'] = 4;
    var buffer: [1]u8 = undefined;
    var writer = std.Io.Writer.fixed(&buffer);

    try std.testing.expectError(
        error.WriteFailed,
        writeReport(
            &writer,
            "x.fa",
            ".zfi",
            8,
            &records,
            summary,
            "x",
            "x",
            &comp,
            try summarizeComposition(&comp),
        ),
    );
}

test "composition categories cannot exceed indexed symbols" {
    var comp = CompositionStats{
        .counts = .{0} ** 256,
        .total_bases = 4,
        .seq_type = .nucleotide,
        .lowercase_count = 0,
    };
    comp.counts['A'] = 5;

    try std.testing.expectError(error.CompositionMismatch, summarizeComposition(&comp));
}
