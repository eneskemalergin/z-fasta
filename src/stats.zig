//! Assembly and proteome statistics: index-only and full composition scan modes.
//!
//! Source duplicate extras are counted from the FASTA when scanning; index-only mode
//! reports `n/a` unless the index retained repeats (`--no-dedup`).

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

/// Extra repeated records by name: for each name with count k, add (k - 1).
pub fn countNameDuplicateExtras(map: *const std.StringHashMap(usize)) usize {
    var extras: usize = 0;
    var it = map.iterator();
    while (it.next()) |entry| {
        if (entry.value_ptr.* > 1) extras += entry.value_ptr.* - 1;
    }
    return extras;
}

/// Tally first-token FASTA header names (line-leading `>`). Names are slices into `fasta`.
pub fn tallyFastaHeaderNames(fasta: []const u8, map: *std.StringHashMap(usize)) !void {
    var at_line_start = true;
    var i: usize = 0;
    while (i < fasta.len) : (i += 1) {
        const c = fasta[i];
        if (at_line_start and c == '>') {
            const name_start = i + 1;
            var name_end = name_start;
            while (name_end < fasta.len and
                fasta[name_end] != ' ' and
                fasta[name_end] != '\t' and
                fasta[name_end] != '\n' and
                fasta[name_end] != '\r')
            {
                name_end += 1;
            }
            const name = fasta[name_start..name_end];
            const gop = try map.getOrPut(name);
            if (!gop.found_existing) gop.value_ptr.* = 0;
            gop.value_ptr.* += 1;
            at_line_start = false;
            continue;
        }
        at_line_start = (c == '\n');
    }
}

const DuplicateReport = union(enum) {
    /// Known source-level extras (`sum(k-1)`).
    count: usize,
    /// Index-only cannot prove absence of duplicates on a deduplicated index.
    unknown,
};

fn reportSourceDuplicates(
    idx: *LoadedIndex,
    io: std.Io,
    allocator: std.mem.Allocator,
    index_only: bool,
) DuplicateReport {
    var map = std.StringHashMap(usize).init(allocator);
    defer map.deinit();

    if (!index_only) {
        tallyFastaHeaderNames(idx.fasta_data, &map) catch {
            printErrorAndExit("error: out of memory\n", .{});
        };
        return .{ .count = countNameDuplicateExtras(&map) };
    }

    // Index-only: report only when the index itself retains repeated names
    // (--no-dedup). A fully unique index cannot distinguish "no source dups"
    // from "dedup dropped them", so never fabricate 0.
    for (0..idx.records.len) |ri| {
        const name = getRecordName(idx, io, ri);
        const gop = map.getOrPut(name) catch {
            printErrorAndExit("error: out of memory\n", .{});
        };
        if (!gop.found_existing) gop.value_ptr.* = 0;
        gop.value_ptr.* += 1;
    }
    var had_repeat = false;
    var it = map.iterator();
    while (it.next()) |entry| {
        if (entry.value_ptr.* > 1) {
            had_repeat = true;
            break;
        }
    }
    if (!had_repeat) return .unknown;
    return .{ .count = countNameDuplicateExtras(&map) };
}

/// Run the stats command.
pub fn runStats(io: std.Io, fasta_path: []const u8, index_only: bool) void {
    var idx = index_format.loadIndexWithMode(io, fasta_path, .stats_scan);
    defer idx.deinit(io);

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

    // Get names - handle both .zfi and .fai sources.
    const shortest_name = getRecordName(&idx, io, shortest_idx);
    const longest_name = getRecordName(&idx, io, longest_idx);

    const duplicates = reportSourceDuplicates(&idx, io, allocator, index_only);

    // Run composition scan early (if not --index-only) so we can include Type in the header
    const comp: ?CompositionStats = if (!index_only) scanComposition(&idx) else null;

    // Output
    var out_buf: [65536]u8 = undefined;
    var stdout_fw = std.Io.File.Writer.initStreaming(.stdout(), io, &out_buf);
    const writer = &stdout_fw.interface;

    var size_buf: [64]u8 = undefined;
    const size_str = formatSize(&size_buf, idx.fasta_size);

    var comma_buf: [64]u8 = undefined;

    const index_ext: []const u8 = switch (idx.source) {
        .zfi => ".zfi",
        .fai => ".fai",
    };

    const type_str: []const u8 = if (comp) |c| switch (c.seq_type) {
        .nucleotide => "Nucleotide",
        .protein => "Protein",
    } else "(run without --index-only for composition)";

    writer.print("File:           {s} ({s} on disk)\n", .{ fasta_path, size_str }) catch {
        printErrorAndExit("error: write failed\n", .{});
    };
    writer.print("Index:          {s}{s}\n", .{ fasta_path, index_ext }) catch {
        printErrorAndExit("error: write failed\n", .{});
    };
    writer.print("Type:           {s}\n", .{type_str}) catch {
        printErrorAndExit("error: write failed\n", .{});
    };
    writer.print("Sequences:      {s}\n", .{formatComma(&comma_buf, @intCast(num_seqs))}) catch {
        printErrorAndExit("error: write failed\n", .{});
    };
    writer.print("Total bases:    {s}\n", .{formatComma(&comma_buf, total_bases)}) catch {
        printErrorAndExit("error: write failed\n", .{});
    };
    writer.print("Shortest:       {s} ({s})\n", .{ formatComma(&comma_buf, records[shortest_idx].seq_len), shortest_name }) catch {
        printErrorAndExit("error: write failed\n", .{});
    };
    writer.print("Longest:        {s} ({s})\n", .{ formatComma(&comma_buf, records[longest_idx].seq_len), longest_name }) catch {
        printErrorAndExit("error: write failed\n", .{});
    };
    writer.print("Mean:           {s}\n", .{formatComma(&comma_buf, mean)}) catch {
        printErrorAndExit("error: write failed\n", .{});
    };
    writer.print("Median:         {s}\n", .{formatComma(&comma_buf, median)}) catch {
        printErrorAndExit("error: write failed\n", .{});
    };
    writer.print("N50:            {s}\n", .{formatComma(&comma_buf, n50)}) catch {
        printErrorAndExit("error: write failed\n", .{});
    };
    writer.print("L50:            {s}\n", .{formatComma(&comma_buf, @intCast(l50))}) catch {
        printErrorAndExit("error: write failed\n", .{});
    };
    writer.print("N90:            {s}\n", .{formatComma(&comma_buf, n90)}) catch {
        printErrorAndExit("error: write failed\n", .{});
    };
    writer.print("L90:            {s}\n", .{formatComma(&comma_buf, @intCast(l90))}) catch {
        printErrorAndExit("error: write failed\n", .{});
    };
    writer.print("AU:             {s}\n", .{formatComma(&comma_buf, au)}) catch {
        printErrorAndExit("error: write failed\n", .{});
    };
    switch (duplicates) {
        .count => |n| writer.print("Duplicates:     {d}\n", .{n}) catch {
            printErrorAndExit("error: write failed\n", .{});
        },
        .unknown => writer.print("Duplicates:     n/a (run without --index-only)\n", .{}) catch {
            printErrorAndExit("error: write failed\n", .{});
        },
    }

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
            writer.print("  N content: {s}\n", .{formatComma(&comma_buf, n_count)}) catch {};

            // Soft-masked
            const soft_pct = if (f_total > 0) @as(f64, @floatFromInt(c.lowercase_count)) / f_total * 100.0 else 0.0;
            writer.print("  Soft-masked: {d:.2}%\n", .{soft_pct}) catch {};
        } else {
            // Protein stats: top 3 most frequent amino acids
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
fn getRecordName(idx: *LoadedIndex, io: std.Io, rec_idx: usize) []const u8 {
    return idx.getRecordNameWithIo(io, rec_idx);
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

/// Inline composition scan for one record (hot path; kept in scanComposition for inlining).
fn scanCompositionRecordInline(
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
    const record_count = records.len;
    if (record_count == 0) {
        return CompositionStats{
            .counts = counts,
            .total_bases = total_bases,
            .seq_type = .nucleotide,
            .lowercase_count = lowercase_count,
        };
    }

    var release = CompositionReleaseCursor{};

    if (recordsInSeqOffsetOrder(records)) {
        var scan_end: usize = 0;
        for (records) |rec| {
            if (rec.seq_len == 0) continue;

            if (!rec.isUniformWidth()) {
                scanCompositionRecordInline(idx, rec, &counts, &total_bases, &lowercase_count, &release);
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
        const sorted_indices = std.heap.page_allocator.alloc(usize, record_count) catch {
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
            scanCompositionRecordInline(
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

    // Detect nucleotide vs protein
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
