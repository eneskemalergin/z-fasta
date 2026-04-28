const std = @import("std");
const posix = std.posix;
const index_format = @import("index_format.zig");

const IndexRecord = index_format.IndexRecord;
const LoadedIndex = index_format.LoadedIndex;
const printErrorAndExit = index_format.printErrorAndExit;

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

/// Format a u64 with comma separators into the provided buffer.
pub fn formatComma(buf: []u8, value: u64) []const u8 {
    if (value == 0) {
        buf[0] = '0';
        return buf[0..1];
    }

    // First format without commas (digits in reverse)
    var tmp: [20]u8 = undefined;
    var v = value;
    var tmp_len: usize = 0;
    while (v > 0) {
        tmp[tmp_len] = @intCast((v % 10) + '0');
        v /= 10;
        tmp_len += 1;
    }

    // Write left-to-right with commas every 3 digits from the right
    var pos: usize = 0;
    var d: usize = 0;
    while (d < tmp_len) {
        if (d > 0 and (tmp_len - d) % 3 == 0) {
            buf[pos] = ',';
            pos += 1;
        }
        buf[pos] = tmp[tmp_len - 1 - d];
        pos += 1;
        d += 1;
    }
    return buf[0..pos];
}

/// Format bytes into human-readable size (e.g., "3.1 GB")
pub fn formatSize(buf: []u8, bytes: u64) []const u8 {
    const gb: f64 = @as(f64, @floatFromInt(bytes)) / (1024.0 * 1024.0 * 1024.0);
    const mb: f64 = @as(f64, @floatFromInt(bytes)) / (1024.0 * 1024.0);
    const kb: f64 = @as(f64, @floatFromInt(bytes)) / 1024.0;

    if (gb >= 1.0) {
        return std.fmt.bufPrint(buf, "{d:.1} GB", .{gb}) catch "?";
    } else if (mb >= 1.0) {
        return std.fmt.bufPrint(buf, "{d:.1} MB", .{mb}) catch "?";
    } else if (kb >= 1.0) {
        return std.fmt.bufPrint(buf, "{d:.1} KB", .{kb}) catch "?";
    } else {
        return std.fmt.bufPrint(buf, "{d} B", .{bytes}) catch "?";
    }
}

const AminoAcid = struct {
    code: u8,
    name: []const u8,
};

const amino_acid_names = [_]AminoAcid{
    .{ .code = 'A', .name = "Alanine" },
    .{ .code = 'R', .name = "Arginine" },
    .{ .code = 'N', .name = "Asparagine" },
    .{ .code = 'D', .name = "Aspartate" },
    .{ .code = 'C', .name = "Cysteine" },
    .{ .code = 'E', .name = "Glutamate" },
    .{ .code = 'Q', .name = "Glutamine" },
    .{ .code = 'G', .name = "Glycine" },
    .{ .code = 'H', .name = "Histidine" },
    .{ .code = 'I', .name = "Isoleucine" },
    .{ .code = 'L', .name = "Leucine" },
    .{ .code = 'K', .name = "Lysine" },
    .{ .code = 'M', .name = "Methionine" },
    .{ .code = 'F', .name = "Phenylalanine" },
    .{ .code = 'P', .name = "Proline" },
    .{ .code = 'S', .name = "Serine" },
    .{ .code = 'T', .name = "Threonine" },
    .{ .code = 'W', .name = "Tryptophan" },
    .{ .code = 'Y', .name = "Tyrosine" },
    .{ .code = 'V', .name = "Valine" },
};

fn getAminoAcidName(code: u8) []const u8 {
    const upper = if (code >= 'a' and code <= 'z') code - 32 else code;
    for (amino_acid_names) |aa| {
        if (aa.code == upper) return aa.name;
    }
    return "Unknown";
}

/// Run the stats command.
pub fn runStats(io: std.Io, fasta_path: []const u8, index_only: bool) void {
    var idx = index_format.loadIndexWithMode(io, fasta_path, .records_only);
    defer idx.deinit();

    const records = idx.records;
    const num_seqs = records.len;

    if (num_seqs == 0) {
        printErrorAndExit("error: no sequences in index\n", .{});
    }

    // Compute total bases
    var total_bases: u64 = 0;
    for (records) |rec| {
        total_bases += rec.seq_len;
    }

    // Sort lengths descending for N50/L50/N90/L90/AU
    var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const lengths = allocator.alloc(u64, num_seqs) catch {
        printErrorAndExit("error: out of memory\n", .{});
    };
    for (records, 0..) |rec, i| {
        lengths[i] = rec.seq_len;
    }

    // Sort descending
    std.mem.sort(u64, lengths, {}, struct {
        fn cmp(_: void, a: u64, b: u64) bool {
            return a > b;
        }
    }.cmp);

    // Single pass: N50, L50, N90, L90, AU
    var bases_seen: u64 = 0;
    var au_sum: u128 = 0;
    var n50: u64 = 0;
    var l50: usize = 0;
    var n90: u64 = 0;
    var l90: usize = 0;
    var found_n50 = false;
    var found_n90 = false;

    const threshold_50 = (total_bases + 1) / 2; // ceiling
    const threshold_90 = (total_bases * 9 + 9) / 10; // ceiling

    for (lengths, 0..) |len, i| {
        bases_seen += len;
        au_sum += @as(u128, len) * @as(u128, len);

        if (!found_n50 and bases_seen >= threshold_50) {
            n50 = len;
            l50 = i + 1;
            found_n50 = true;
        }
        if (!found_n90 and bases_seen >= threshold_90) {
            n90 = len;
            l90 = i + 1;
            found_n90 = true;
        }
    }

    const au: u64 = if (total_bases > 0) @intCast(au_sum / @as(u128, total_bases)) else 0;

    // Mean and median
    const mean: u64 = if (num_seqs > 0) total_bases / num_seqs else 0;
    const median: u64 = if (num_seqs > 0) blk: {
        if (num_seqs % 2 == 1) {
            break :blk lengths[num_seqs / 2];
        } else {
            break :blk (lengths[num_seqs / 2 - 1] + lengths[num_seqs / 2]) / 2;
        }
    } else 0;

    // Find shortest and longest with names
    var shortest_idx: usize = 0;
    var longest_idx: usize = 0;
    for (records, 0..) |rec, i| {
        if (rec.seq_len < records[shortest_idx].seq_len) shortest_idx = i;
        if (rec.seq_len > records[longest_idx].seq_len) longest_idx = i;
    }

    // Get names - handle both .zfi and .fai sources
    const shortest_name = getRecordName(&idx, shortest_idx);
    const longest_name = getRecordName(&idx, longest_idx);

    // Count duplicates (check name_map size vs record count)
    // In a deduplicated index, the name_map size == record count
    // For .fai, duplicates were already filtered by samtools
    const duplicates: usize = 0; // Deduplicated during indexing

    // Run composition scan early (if not --index-only) so we can include Type in the header
    const comp: ?CompositionStats = if (!index_only) scanComposition(&idx) else null;

    // Output
    var out_buf: [65536]u8 = undefined;
    var stdout_fw = std.Io.File.Writer.initStreaming(.stdout(), io, &out_buf);
    const writer = &stdout_fw.interface;

    var size_buf: [64]u8 = undefined;
    const size_str = formatSize(&size_buf, idx.fasta_size);

    var cb1: [64]u8 = undefined;
    var cb2: [64]u8 = undefined;
    var cb3: [64]u8 = undefined;
    var cb4: [64]u8 = undefined;
    var cb5: [64]u8 = undefined;
    var cb6: [64]u8 = undefined;
    var cb7: [64]u8 = undefined;
    var cb8: [64]u8 = undefined;
    var cb9: [64]u8 = undefined;
    var cb10: [64]u8 = undefined;
    var cb11: [64]u8 = undefined;
    var cb12: [64]u8 = undefined;

    const index_ext: []const u8 = switch (idx.source) {
        .zfi => ".zfi",
        .fai => ".fai",
    };

    const type_str: []const u8 = if (comp) |c| switch (c.seq_type) {
        .nucleotide => "Nucleotide",
        .protein => "Protein",
    } else "(run without --index-only for composition)";

    writer.print(
        \\File:           {s} ({s} on disk)
        \\Index:          {s}{s}
        \\Type:           {s}
        \\Sequences:      {s}
        \\Total bases:    {s}
        \\Shortest:       {s} ({s})
        \\Longest:        {s} ({s})
        \\Mean:           {s}
        \\Median:         {s}
        \\N50:            {s}
        \\L50:            {s}
        \\N90:            {s}
        \\L90:            {s}
        \\AU:             {s}
        \\Duplicates:     {d}
        \\
    , .{
        fasta_path,
        size_str,
        fasta_path,
        index_ext,
        type_str,
        formatComma(&cb1, @intCast(num_seqs)),
        formatComma(&cb2, total_bases),
        formatComma(&cb3, records[shortest_idx].seq_len),
        shortest_name,
        formatComma(&cb4, records[longest_idx].seq_len),
        longest_name,
        formatComma(&cb5, mean),
        formatComma(&cb6, median),
        formatComma(&cb7, n50),
        formatComma(&cb8, @intCast(l50)),
        formatComma(&cb9, n90),
        formatComma(&cb10, @intCast(l90)),
        formatComma(&cb11, au),
        duplicates,
    }) catch {
        printErrorAndExit("error: write failed\n", .{});
    };

    // Tier 2: composition details (unless --index-only)
    if (comp) |c| {
        writer.print("\nComposition:\n", .{}) catch {};

        if (c.seq_type == .nucleotide) {
            // Nucleotide stats
            const a_count = c.counts['A'] + c.counts['a'];
            const c_count = c.counts['C'] + c.counts['c'];
            const g_count = c.counts['G'] + c.counts['g'];
            const t_count = c.counts['T'] + c.counts['t'];
            const n_count = c.counts['N'] + c.counts['n'];

            const acgt_total = a_count + c_count + g_count + t_count;
            const all_total = c.total_bases;

            // Other = total - A - C - G - T - N
            const other = if (all_total > a_count + c_count + g_count + t_count + n_count)
                all_total - a_count - c_count - g_count - t_count - n_count
            else
                0;

            const f_total: f64 = @floatFromInt(all_total);
            const f_acgt: f64 = @floatFromInt(acgt_total);

            writer.print("  A:   {d:.2}%\n", .{if (f_total > 0) @as(f64, @floatFromInt(a_count)) / f_total * 100.0 else 0.0}) catch {};
            writer.print("  C:   {d:.2}%\n", .{if (f_total > 0) @as(f64, @floatFromInt(c_count)) / f_total * 100.0 else 0.0}) catch {};
            writer.print("  G:   {d:.2}%\n", .{if (f_total > 0) @as(f64, @floatFromInt(g_count)) / f_total * 100.0 else 0.0}) catch {};
            writer.print("  T:   {d:.2}%\n", .{if (f_total > 0) @as(f64, @floatFromInt(t_count)) / f_total * 100.0 else 0.0}) catch {};
            writer.print("  N:   {d:.2}%\n", .{if (f_total > 0) @as(f64, @floatFromInt(n_count)) / f_total * 100.0 else 0.0}) catch {};

            if (other > 0) {
                writer.print("  Other: {d:.2}%\n", .{@as(f64, @floatFromInt(other)) / f_total * 100.0}) catch {};
            }

            // GC content (N excluded from denominator)
            const gc = @as(f64, @floatFromInt(g_count + c_count));
            const gc_pct = if (f_acgt > 0) gc / f_acgt * 100.0 else 0.0;
            const n_pct = if (f_total > 0) @as(f64, @floatFromInt(n_count)) / f_total * 100.0 else 0.0;

            if (n_pct > 1.0) {
                writer.print("  GC:  {d:.2}% (excl. N)\n", .{gc_pct}) catch {};
            } else {
                writer.print("  GC:  {d:.2}%\n", .{gc_pct}) catch {};
            }

            // GC skew
            const f_g: f64 = @floatFromInt(g_count);
            const f_c: f64 = @floatFromInt(c_count);
            const gc_sum = f_g + f_c;
            if (gc_sum > 0) {
                writer.print("  GC skew: {d:.3}\n", .{(f_g - f_c) / gc_sum}) catch {};
            }

            // N content
            writer.print("  N content: {s}\n", .{formatComma(&cb12, n_count)}) catch {};

            // Soft-masked
            const soft_pct = if (f_total > 0) @as(f64, @floatFromInt(c.lowercase_count)) / f_total * 100.0 else 0.0;
            writer.print("  Soft-masked: {d:.2}%\n", .{soft_pct}) catch {};
        } else {
            // Protein stats — top 3 most frequent amino acids
            const f_total: f64 = @floatFromInt(c.total_bases);

            // Collect amino acid counts (combine upper+lower)
            const AaCount = struct { code: u8, count: u64 };
            var aa_counts: [20]AaCount = undefined;
            for (amino_acid_names, 0..) |aa, i| {
                const lower = aa.code + 32;
                aa_counts[i] = .{
                    .code = aa.code,
                    .count = c.counts[aa.code] + c.counts[lower],
                };
            }

            // Sort by count descending
            std.mem.sort(AaCount, &aa_counts, {}, struct {
                fn cmp(_: void, a: AaCount, b: AaCount) bool {
                    return a.count > b.count;
                }
            }.cmp);

            // Print top 3
            for (0..@min(3, aa_counts.len)) |i| {
                const pct = if (f_total > 0) @as(f64, @floatFromInt(aa_counts[i].count)) / f_total * 100.0 else 0.0;
                writer.print("  {c}:   {d:.2}%  ({s})\n", .{
                    aa_counts[i].code,
                    pct,
                    getAminoAcidName(aa_counts[i].code),
                }) catch {};
            }
            writer.print("  (20 amino acids total)\n", .{}) catch {};

            // Lowercase fraction
            const soft_pct = if (f_total > 0) @as(f64, @floatFromInt(c.lowercase_count)) / f_total * 100.0 else 0.0;
            writer.print("  Lowercase: {d:.2}%\n", .{soft_pct}) catch {};
        }
    }

    stdout_fw.flush() catch {};
}

/// Get the name of a record, handling both .zfi and .fai sources.
fn getRecordName(idx: *const LoadedIndex, rec_idx: usize) []const u8 {
    const rec = idx.records[rec_idx];

    // For .zfi, names are slices into the FASTA mmap
    if (idx.source == .zfi and rec.name_len > 0) {
        return idx.fasta_data[rec.name_offset..][0..rec.name_len];
    }

    // For .fai, look up in name_map (iterate to find by index)
    var it = idx.name_map.iterator();
    while (it.next()) |entry| {
        if (entry.value_ptr.* == rec_idx) {
            return entry.key_ptr.*;
        }
    }
    return "?";
}

/// Scan all sequence regions in the FASTA for composition.
fn scanComposition(idx: *const LoadedIndex) CompositionStats {
    const fasta = idx.fasta_data;

    // Set madvise to sequential for full scan
    posix.madvise(
        @alignCast(@constCast(fasta.ptr)),
        fasta.len,
        posix.MADV.SEQUENTIAL,
    ) catch {};

    var counts: [256]u64 = .{0} ** 256;
    var lowercase_count: u64 = 0;
    var total_bases: u64 = 0;

    // Scan each sequence region
    for (idx.records) |rec| {
        if (rec.seq_len == 0) continue;

        const start: usize = @intCast(rec.seq_offset);

        if (tryScanFixedWidthRecord(fasta, rec, start, &counts, &total_bases, &lowercase_count)) {
            continue;
        }

        // Compute end of sequence region
        // end = seq_offset + (full_lines * line_bytes) + last_line_bytes
        const full_lines = rec.seq_len / rec.line_bases;
        const remainder = rec.seq_len % rec.line_bases;
        var end: usize = start;
        if (remainder > 0) {
            end = start + (full_lines * rec.line_bytes) + remainder;
            // Account for possible trailing newline
            if (end < fasta.len and (fasta[end] == '\n' or fasta[end] == '\r')) {
                end += 1;
                if (end < fasta.len and fasta[end - 1] == '\r' and fasta[end] == '\n') {
                    end += 1;
                }
            }
        } else {
            // seq_len is exact multiple of line_bases
            end = start + (full_lines * rec.line_bytes);
        }

        if (end > fasta.len) end = fasta.len;

        // Branchless byte counting
        const region = fasta[start..end];
        for (region) |byte| {
            if (byte != '\n' and byte != '\r') {
                counts[byte] += 1;
                total_bases += 1;
                if (byte >= 'a' and byte <= 'z') {
                    lowercase_count += 1;
                }
            }
        }
    }

    // Detect nucleotide vs protein
    const seq_type = detectType(&counts, total_bases);

    return CompositionStats{
        .counts = counts,
        .total_bases = total_bases,
        .seq_type = seq_type,
        .lowercase_count = lowercase_count,
    };
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
    for (data) |byte| {
        counts[byte] += 1;
        total_bases.* += 1;
        if (byte >= 'a' and byte <= 'z') {
            lowercase_count.* += 1;
        }
    }
}

/// Detect if sequences are nucleotide or protein based on first 100K bases.
pub fn detectType(counts: *const [256]u64, total: u64) SequenceType {
    if (total == 0) return .nucleotide;

    // Sum A, C, G, T, N, U (case-insensitive)
    const nuc_count = counts['A'] + counts['a'] +
        counts['C'] + counts['c'] +
        counts['G'] + counts['g'] +
        counts['T'] + counts['t'] +
        counts['N'] + counts['n'] +
        counts['U'] + counts['u'];

    // >90% nucleotide characters => nucleotide
    if (nuc_count * 10 > total * 9) {
        return .nucleotide;
    }
    return .protein;
}
