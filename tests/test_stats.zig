//! Stats unit and CLI tests: type detection and exact report behavior.
//!
//! Exercises complete stats output against fixture FASTAs.

const std = @import("std");
const main = @import("main");
const utility = @import("utility.zig");
const stats = main.stats;
const io = std.testing.io;

const ZFASTA_BIN = utility.ZFASTA_BIN;
const expectCliFailure = utility.expectCliFailure;
const expectUnknownOptionRejected = utility.expectUnknownOptionRejected;

fn countSymbols(symbols: []const u8) [256]u64 {
    var counts: [256]u64 = .{0} ** 256;
    for (symbols) |symbol| counts[symbol] += 1;
    return counts;
}

fn runStatsAndCapture(allocator: std.mem.Allocator, fasta_path: []const u8) ![]u8 {
    var threaded = std.Io.Threaded.init(allocator, .{});
    defer threaded.deinit();

    const result = try std.process.run(allocator, threaded.io(), .{
        .argv = &.{ ZFASTA_BIN, "stats", fasta_path },
        .stdout_limit = .limited(10 * 1024 * 1024),
        .stderr_limit = .limited(64 * 1024),
    });
    errdefer allocator.free(result.stdout);
    defer allocator.free(result.stderr);

    switch (result.term) {
        .exited => |code| if (code != 0) return error.ChildProcessFailed,
        else => return error.ChildProcessFailed,
    }
    try std.testing.expectEqual(@as(usize, 0), result.stderr.len);
    return result.stdout;
}

// --- Type detection ---

test "detectType classifies representative alphabets" {
    const Case = struct {
        symbols: []const u8,
        expected: stats.SequenceType,
    };
    const cases = [_]Case{
        .{ .symbols = "", .expected = .nucleotide },
        .{ .symbols = "ACGTACGT", .expected = .nucleotide },
        .{ .symbols = "acgtacgt", .expected = .nucleotide },
        .{ .symbols = "ACGTURYSWKMBDHVNacgturyswkmbdhvn", .expected = .nucleotide },
        .{ .symbols = "EFLPQXZJO*eflpqxzjo*", .expected = .protein },
    };

    for (cases) |case| {
        const counts = countSymbols(case.symbols);

        try std.testing.expectEqual(case.expected, stats.detectType(&counts, case.symbols.len));
    }
}

test "detectType uses a strict overflow-safe 90 percent boundary" {
    const Case = struct {
        total: u64,
        nucleotide: u64,
        expected: stats.SequenceType,
    };
    const cases = [_]Case{
        .{ .total = 10, .nucleotide = 9, .expected = .protein },
        .{ .total = 10, .nucleotide = 10, .expected = .nucleotide },
        .{ .total = 11, .nucleotide = 10, .expected = .nucleotide },
        .{ .total = 1000, .nucleotide = 900, .expected = .protein },
        .{ .total = 1000, .nucleotide = 901, .expected = .nucleotide },
        .{
            .total = 2_300_000_000_000_000_000,
            .nucleotide = 2_070_000_000_000_000_000,
            .expected = .protein,
        },
        .{
            .total = 2_300_000_000_000_000_000,
            .nucleotide = 2_070_000_000_000_000_001,
            .expected = .nucleotide,
        },
    };

    for (cases) |case| {
        var counts: [256]u64 = .{0} ** 256;
        counts['A'] = case.nucleotide;
        counts['L'] = case.total - case.nucleotide;

        try std.testing.expectEqual(case.expected, stats.detectType(&counts, case.total));
    }
}

test "detectType matches its threshold across deterministic fuzzy inputs" {
    const nucleotide_codes = "ACGTURYSWKMBDHVNacgturyswkmbdhvn";
    const protein_only_codes = "EFLPQXZJO*eflpqxzjo*";
    var prng = std.Random.DefaultPrng.init(0x7a_fa_57_a7_5);
    const random = prng.random();

    for (0..4096) |i| {
        const total = random.intRangeAtMost(u64, 1, 1_000_000);
        const nucleotide = random.intRangeAtMost(u64, 0, total);
        var counts: [256]u64 = .{0} ** 256;
        counts[nucleotide_codes[i % nucleotide_codes.len]] = nucleotide;
        counts[protein_only_codes[i % protein_only_codes.len]] = total - nucleotide;
        const expected: stats.SequenceType = if (@as(u128, nucleotide) * 10 > @as(u128, total) * 9)
            .nucleotide
        else
            .protein;

        try std.testing.expectEqual(expected, stats.detectType(&counts, total));
    }
}

// --- CLI reports ---

test "[cli] - [stats]: rejects unknown options regardless of position" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const fasta = "tests/data/simple.fasta";

    try expectUnknownOptionRejected(allocator, &.{ ZFASTA_BIN, "stats", "--not-a-flag", fasta }, "--not-a-flag");
    try expectUnknownOptionRejected(allocator, &.{ ZFASTA_BIN, "stats", fasta, "--not-a-flag" }, "--not-a-flag");
    try expectUnknownOptionRejected(allocator, &.{ ZFASTA_BIN, "stats", "--index-only", fasta }, "--index-only");
}

test "[cli] - [stats]: rejects invalid FASTA path counts" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    try expectCliFailure(
        allocator,
        &.{ ZFASTA_BIN, "stats" },
        1,
        "error: usage: z-fasta stats <file.fasta>\n",
    );
    try expectCliFailure(
        allocator,
        &.{ ZFASTA_BIN, "stats", "tests/data/simple.fasta", "tests/data/single.fasta" },
        1,
        "error: stats accepts exactly one FASTA path\n",
    );
    try expectCliFailure(
        allocator,
        &.{ ZFASTA_BIN, "stats", "tests/data/definitely-missing-cli-failure.fasta" },
        1,
        "error: file not found: tests/data/definitely-missing-cli-failure.fasta\n",
    );
}

test "[cli] - [stats index selection]: invalid zfi blocks a valid fai" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta = try utility.uniqueArtifactPath(allocator, "stats-invalid-zfi", "fa");
    const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{fasta});
    const fai_path = try std.fmt.allocPrint(allocator, "{s}.fai", .{fasta});
    defer std.Io.Dir.cwd().deleteFile(io, fasta) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, zfi_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, fai_path) catch {};

    {
        const fasta_file = try std.Io.Dir.cwd().createFile(io, fasta, .{});
        defer fasta_file.close(io);
        try std.Io.File.writeStreamingAll(fasta_file, io, ">seq\nACGT\n");
    }
    {
        const fai_file = try std.Io.Dir.cwd().createFile(io, fai_path, .{});
        defer fai_file.close(io);
        try std.Io.File.writeStreamingAll(fai_file, io, "seq\t4\t5\t4\t5\n");
    }
    {
        const zfi_file = try std.Io.Dir.cwd().createFile(io, zfi_path, .{});
        defer zfi_file.close(io);
        try std.Io.File.writeStreamingAll(zfi_file, io, "not a zfi index");
    }

    const expected = try std.fmt.allocPrint(allocator, "error: corrupt index file for: {s}\n", .{fasta});
    try expectCliFailure(allocator, &.{ ZFASTA_BIN, "stats", fasta }, 1, expected);
}

test "stats renders the exact nucleotide report" {
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

    const output = try runStatsAndCapture(std.testing.allocator, "tests/data/simple.fasta");
    defer std.testing.allocator.free(output);

    try std.testing.expectEqualStrings(expected, output);
}

test "stats renders the exact protein report" {
    const expected =
        \\File:
        \\  path: tests/data/proteome.fasta
        \\  index: tests/data/proteome.fasta.zfi
        \\  size_bytes: 187
        \\
        \\Lengths:
        \\  indexed_records: 2
        \\  total_symbols: 71
        \\  shortest_length: 20
        \\  shortest_name: sp|Q98765|ANOT_MOUSE
        \\  longest_length: 51
        \\  longest_name: sp|P12345|PROT_HUMAN
        \\  mean: 35
        \\  q1: 20
        \\  median: 35
        \\  q3: 51
        \\  range: 31
        \\
        \\Nx:
        \\  n50: 51
        \\  l50: 1
        \\  n90: 20
        \\  l90: 2
        \\  aun: 42.27
        \\
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

    const output = try runStatsAndCapture(std.testing.allocator, "tests/data/proteome.fasta");
    defer std.testing.allocator.free(output);

    try std.testing.expectEqualStrings(expected, output);
}
