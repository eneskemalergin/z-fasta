//! Stats unit and CLI tests: type detection and exact report behavior.
//!
//! Exercises complete stats output against fixture FASTAs.

const std = @import("std");
const builtin = @import("builtin");
const main = @import("main");
const detectType = main.stats.detectType;
const SequenceType = main.stats.SequenceType;

// ============================================================================
// detectType tests
// ============================================================================

test "detectType - all ACGT is nucleotide" {
    var counts: [256]u64 = .{0} ** 256;
    counts['A'] = 250;
    counts['C'] = 250;
    counts['G'] = 250;
    counts['T'] = 250;
    try std.testing.expectEqual(SequenceType.nucleotide, detectType(&counts, 1000));
}

test "detectType - ACGTN is nucleotide" {
    var counts: [256]u64 = .{0} ** 256;
    counts['A'] = 200;
    counts['C'] = 200;
    counts['G'] = 200;
    counts['T'] = 200;
    counts['N'] = 200;
    try std.testing.expectEqual(SequenceType.nucleotide, detectType(&counts, 1000));
}

test "detectType - full IUPAC ambiguity alphabet is nucleotide" {
    var counts: [256]u64 = .{0} ** 256;
    const letters = "ACGTURYSWKMBDHVNacgturyswkmbdhvnu";
    for (letters) |byte| counts[byte] += 1;
    try std.testing.expectEqual(SequenceType.nucleotide, detectType(&counts, letters.len));
}

test "detectType - mixed amino acids is protein" {
    var counts: [256]u64 = .{0} ** 256;
    counts['M'] = 100;
    counts['A'] = 100;
    counts['L'] = 100;
    counts['F'] = 100;
    counts['P'] = 100;
    counts['W'] = 100;
    counts['H'] = 100;
    counts['K'] = 100;
    counts['D'] = 100;
    counts['E'] = 100;
    try std.testing.expectEqual(SequenceType.protein, detectType(&counts, 1000));
}

test "detectType - below 90% threshold is protein" {
    var counts: [256]u64 = .{0} ** 256;
    counts['A'] = 200;
    counts['C'] = 200;
    counts['G'] = 200;
    // A+C+G = 600 out of 700, which is 85.7% => protein.
    // Use letters outside the IUPAC nucleotide alphabet for the remainder.
    counts['L'] = 50;
    counts['F'] = 50;
    try std.testing.expectEqual(SequenceType.protein, detectType(&counts, 700));
}

test "detectType - exactly 91% is nucleotide" {
    var counts: [256]u64 = .{0} ** 256;
    counts['A'] = 910;
    counts['L'] = 90;
    // ACGTN = 910/1000 = 91% > 90%
    try std.testing.expectEqual(SequenceType.nucleotide, detectType(&counts, 1000));
}

test "detectType - empty is nucleotide (default)" {
    var counts: [256]u64 = .{0} ** 256;
    try std.testing.expectEqual(SequenceType.nucleotide, detectType(&counts, 0));
}

test "detectType - lowercase nucleotides" {
    var counts: [256]u64 = .{0} ** 256;
    counts['a'] = 250;
    counts['c'] = 250;
    counts['g'] = 250;
    counts['t'] = 250;
    try std.testing.expectEqual(SequenceType.nucleotide, detectType(&counts, 1000));
}

test "detectType - large totals do not overflow the 90% threshold" {
    var counts: [256]u64 = .{0} ** 256;
    // total*9 and nuc*10 both overflow u64 at this scale; u128 keeps the compare correct.
    const total: u64 = 2_300_000_000_000_000_000;
    const nuc_hi: u64 = 2_070_000_000_000_000_001;
    counts['A'] = nuc_hi;
    counts['L'] = total - nuc_hi;
    try std.testing.expectEqual(SequenceType.nucleotide, detectType(&counts, total));

    counts = .{0} ** 256;
    const nuc_lo: u64 = 2_070_000_000_000_000_000;
    counts['A'] = nuc_lo;
    counts['L'] = total - nuc_lo;
    try std.testing.expectEqual(SequenceType.protein, detectType(&counts, total));
}

// ============================================================================
// Integration: stats via process spawn
// ============================================================================

const ZFASTA_BIN = if (builtin.os.tag == .windows) "zig-out\\bin\\z-fasta.exe" else "zig-out/bin/z-fasta";

fn runStatsAndCapture(allocator: std.mem.Allocator, fasta_path: []const u8) ![]u8 {
    var threaded = std.Io.Threaded.init(allocator, .{});
    defer threaded.deinit();
    const io = threaded.io();

    var proc = try std.process.spawn(io, .{
        .argv = &.{ ZFASTA_BIN, "stats", fasta_path },
        .stdout = .pipe,
    });
    defer proc.kill(io);

    var read_buf: [4096]u8 = undefined;
    var stdout_reader = proc.stdout.?.reader(io, &read_buf);
    const result = try stdout_reader.interface.allocRemaining(allocator, .limited(10 * 1024 * 1024));
    switch (try proc.wait(io)) {
        .exited => |code| if (code != 0) return error.ChildProcessFailed,
        else => return error.ChildProcessFailed,
    }
    return result;
}

test "stats renders the exact nucleotide report" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const output = try runStatsAndCapture(arena.allocator(), "tests/data/simple.fasta");

    const expected =
        \\File:
        \\  path: tests/data/simple.fasta
        \\  index: tests/data/simple.fasta.zfi
        \\  size_bytes: 82
        \\
        \\Lengths:
        \\  indexed_records: 2
        \\  total_symbols: 36
        \\  shortest_length: 12
        \\  shortest_name: seq2
        \\  longest_length: 24
        \\  longest_name: seq1
        \\  mean: 18
        \\  q1: 12
        \\  median: 18
        \\  q3: 24
        \\  range: 12
        \\
        \\Nx:
        \\  n50: 24
        \\  l50: 1
        \\  n90: 12
        \\  l90: 2
        \\  aun: 20.00
        \\
        \\Composition:
        \\  type: nucleotide_t
        \\  percent_denominator: total_symbols
        \\  a: 10 27.78%
        \\  c: 10 27.78%
        \\  g: 10 27.78%
        \\  t: 6 16.67%
        \\  u: 0 0.00%
        \\  n: 0 0.00%
        \\  iupac_ambiguous: 0 0.00%
        \\  invalid: 0 0.00%
        \\  gc: 55.56%
        \\  gc_skew: 0.000
        \\  lowercase: 0 0.00%
        \\
    ;

    try std.testing.expectEqualStrings(expected, output);
}

test "stats renders the exact complete protein field set" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const output = try runStatsAndCapture(arena.allocator(), "tests/data/proteome.fasta");

    const composition =
        \\Composition:
        \\  type: protein
        \\  percent_denominator: total_symbols
        \\  a_alanine: 8 11.27%
        \\  r_arginine: 2 2.82%
        \\  n_asparagine: 2 2.82%
        \\  d_aspartate: 3 4.23%
        \\  c_cysteine: 1 1.41%
        \\  e_glutamate: 4 5.63%
        \\  q_glutamine: 1 1.41%
        \\  g_glycine: 5 7.04%
        \\  h_histidine: 4 5.63%
        \\  i_isoleucine: 1 1.41%
        \\  l_leucine: 5 7.04%
        \\  k_lysine: 5 7.04%
        \\  m_methionine: 3 4.23%
        \\  f_phenylalanine: 5 7.04%
        \\  p_proline: 4 5.63%
        \\  s_serine: 4 5.63%
        \\  t_threonine: 5 7.04%
        \\  w_tryptophan: 2 2.82%
        \\  y_tyrosine: 3 4.23%
        \\  v_valine: 4 5.63%
        \\  b_asx: 0 0.00%
        \\  z_glx: 0 0.00%
        \\  j_xle: 0 0.00%
        \\  x_unknown: 0 0.00%
        \\  u_selenocysteine: 0 0.00%
        \\  o_pyrrolysine: 0 0.00%
        \\  stop: 0 0.00%
        \\  invalid: 0 0.00%
        \\  lowercase: 0 0.00%
        \\
    ;
    const start = std.mem.indexOf(u8, output, "Composition:\n") orelse return error.MissingComposition;

    try std.testing.expectEqualStrings(composition, output[start..]);
}
