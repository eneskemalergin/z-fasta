//! Stats tests: type detection, exact reports, index parity, and CLI failures.
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
        .stdout_limit = .limited(64 * 1024),
        .stderr_limit = .limited(64 * 1024),
    });
    errdefer allocator.free(result.stdout);
    defer allocator.free(result.stderr);

    switch (result.term) {
        .exited => |code| try std.testing.expectEqual(@as(u8, 0), code),
        else => return error.ChildProcessFailed,
    }
    try std.testing.expectEqualStrings("", result.stderr);
    return result.stdout;
}

fn writeSingleLineFaiFixture(
    allocator: std.mem.Allocator,
    stem: []const u8,
    name: []const u8,
    sequence: []const u8,
) !struct { fasta_path: []u8, fai_path: []u8 } {
    const fasta_path = try utility.uniqueArtifactPath(allocator, stem, "fa");
    errdefer allocator.free(fasta_path);
    const fai_path = try std.fmt.allocPrint(allocator, "{s}.fai", .{fasta_path});
    errdefer allocator.free(fai_path);
    errdefer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    errdefer std.Io.Dir.cwd().deleteFile(io, fai_path) catch {};

    {
        const fasta_file = try std.Io.Dir.cwd().createFile(io, fasta_path, .{});
        defer fasta_file.close(io);
        try std.Io.File.writeStreamingAll(fasta_file, io, ">");
        try std.Io.File.writeStreamingAll(fasta_file, io, name);
        try std.Io.File.writeStreamingAll(fasta_file, io, "\n");
        try std.Io.File.writeStreamingAll(fasta_file, io, sequence);
        try std.Io.File.writeStreamingAll(fasta_file, io, "\n");
    }
    {
        const fai_file = try std.Io.Dir.cwd().createFile(io, fai_path, .{});
        defer fai_file.close(io);
        var buffer: [128]u8 = undefined;
        var writer = fai_file.writer(io, &buffer);
        try writer.interface.print("{s}\t{d}\t{d}\t{d}\t{d}\n", .{
            name,
            sequence.len,
            name.len + 2,
            sequence.len,
            sequence.len + 1,
        });
        try writer.interface.flush();
    }
    return .{ .fasta_path = fasta_path, .fai_path = fai_path };
}

fn expectComposition(expected: []const u8, report: []const u8) !void {
    const start = std.mem.indexOf(u8, report, "Composition:\n") orelse return error.TestExpectedEqual;
    try std.testing.expectEqualStrings(expected, report[start..]);
}

// --- Type detection ---

test "[unit] - [sequence type]: classifies representative alphabets" {
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

test "[edge] - [sequence type]: uses a strict overflow-safe 90 percent boundary" {
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

test "[property] - [sequence type]: matches the threshold across generated counts" {
    const nucleotide_codes = "ACGTURYSWKMBDHVNacgturyswkmbdhvn";
    const protein_only_codes = "EFLPQXZJO*eflpqxzjo*";
    var prng = std.Random.DefaultPrng.init(std.testing.random_seed);
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

        std.testing.expectEqual(expected, stats.detectType(&counts, total)) catch {
            std.debug.print(
                "sequence type mismatch: seed={d}, iteration={d}, total={d}, nucleotide={d}\n",
                .{ std.testing.random_seed, i, total, nucleotide },
            );
            return error.TestExpectedEqual;
        };
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

test "[cli] - [stats]: reports missing and corrupt indexes" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta = try utility.writeFastaArtifact(allocator, "stats-index-errors", ">seq\nACGT\n");
    defer std.Io.Dir.cwd().deleteFile(io, fasta) catch {};
    const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{fasta});
    defer std.Io.Dir.cwd().deleteFile(io, zfi_path) catch {};

    const missing = try std.fmt.allocPrint(
        allocator,
        "error: no index found for {s}. Run 'z-fasta index {s}' first.\n",
        .{ fasta, fasta },
    );
    try expectCliFailure(allocator, &.{ ZFASTA_BIN, "stats", fasta }, 1, missing);

    {
        const zfi_file = try std.Io.Dir.cwd().createFile(io, zfi_path, .{});
        defer zfi_file.close(io);
        try std.Io.File.writeStreamingAll(zfi_file, io, "not a zfi index");
    }

    const corrupt = try std.fmt.allocPrint(allocator, "error: corrupt index file for: {s}\n", .{fasta});
    try expectCliFailure(allocator, &.{ ZFASTA_BIN, "stats", fasta }, 1, corrupt);
}

test "[cli] - [stats report]: renders the exact nucleotide report" {
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

test "[integration] - [stats report]: ZFI and FAI render identical mixed nucleotide statistics" {
    const fasta_data = ">mixed\nacgtuACGTURYSWKMBDHVN?\n";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const paths = try writeSingleLineFaiFixture(
        allocator,
        "stats-index-parity",
        "mixed",
        "acgtuACGTURYSWKMBDHVN?",
    );
    const fasta_path = paths.fasta_path;
    const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{fasta_path});
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, zfi_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, paths.fai_path) catch {};
    {
        const fasta_file = try std.Io.Dir.cwd().openFile(io, fasta_path, .{});
        defer fasta_file.close(io);
        const fasta_stat = try fasta_file.stat(io);
        var index = try main.indexer.scanZfiData(fasta_data, true, allocator);
        defer index.deinit(allocator);
        try main.indexer.writeZfiIndexFile(
            io,
            zfi_path,
            &index,
            fasta_data.len,
            try main.index_format.timestampToNs(fasta_stat.mtime),
        );
    }

    const expected_composition =
        \\Composition:
        \\  type: nucleotide_mixed_tu
        \\  percent_denominator: total_symbols
        \\  a: 2 9.09%
        \\  c: 2 9.09%
        \\  g: 2 9.09%
        \\  t: 2 9.09%
        \\  u: 2 9.09%
        \\  n: 1 4.55%
        \\  iupac_ambiguous: 10 45.45%
        \\  invalid: 1 4.55%
        \\  gc: 40.00%
        \\  gc_skew: 0.000
        \\  lowercase: 5 22.73%
        \\
    ;

    const zfi_output = try runStatsAndCapture(allocator, fasta_path);
    defer allocator.free(zfi_output);
    try expectComposition(expected_composition, zfi_output);

    try std.Io.Dir.cwd().deleteFile(io, zfi_path);
    const fai_output = try runStatsAndCapture(allocator, fasta_path);
    defer allocator.free(fai_output);
    try expectComposition(expected_composition, fai_output);

    const zfi_index = try std.fmt.allocPrint(allocator, "  index: {s}.zfi\n", .{fasta_path});
    const fai_index = try std.fmt.allocPrint(allocator, "  index: {s}.fai\n", .{fasta_path});
    const zfi_index_start = std.mem.indexOf(u8, zfi_output, zfi_index) orelse return error.TestExpectedEqual;
    const fai_index_start = std.mem.indexOf(u8, fai_output, fai_index) orelse return error.TestExpectedEqual;
    try std.testing.expectEqualStrings(zfi_output[0..zfi_index_start], fai_output[0..fai_index_start]);
    try std.testing.expectEqualStrings(
        zfi_output[zfi_index_start + zfi_index.len ..],
        fai_output[fai_index_start + fai_index.len ..],
    );
}

test "[cli] - [stats report]: renders RNA and ambiguity-only nucleotide variants" {
    const Case = struct {
        stem: []const u8,
        name: []const u8,
        sequence: []const u8,
        composition: []const u8,
    };
    const cases = [_]Case{
        .{
            .stem = "stats-rna",
            .name = "rna",
            .sequence = "ACcuUU",
            .composition =
            \\Composition:
            \\  type: nucleotide_u
            \\  percent_denominator: total_symbols
            \\  a: 1 16.67%
            \\  c: 2 33.33%
            \\  g: 0 0.00%
            \\  t: 0 0.00%
            \\  u: 3 50.00%
            \\  n: 0 0.00%
            \\  iupac_ambiguous: 0 0.00%
            \\  invalid: 0 0.00%
            \\  gc: 33.33%
            \\  gc_skew: -1.000
            \\  lowercase: 2 33.33%
            \\
            ,
        },
        .{
            .stem = "stats-ambiguity",
            .name = "ambiguity",
            .sequence = "RYSWKMBDHVNN",
            .composition =
            \\Composition:
            \\  type: nucleotide
            \\  percent_denominator: total_symbols
            \\  a: 0 0.00%
            \\  c: 0 0.00%
            \\  g: 0 0.00%
            \\  t: 0 0.00%
            \\  u: 0 0.00%
            \\  n: 2 16.67%
            \\  iupac_ambiguous: 10 83.33%
            \\  invalid: 0 0.00%
            \\  gc: n/a
            \\  gc_skew: n/a
            \\  lowercase: 0 0.00%
            \\
            ,
        },
    };
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    for (cases) |case| {
        const paths = try writeSingleLineFaiFixture(
            allocator,
            case.stem,
            case.name,
            case.sequence,
        );
        defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
        defer std.Io.Dir.cwd().deleteFile(io, paths.fai_path) catch {};

        const output = try runStatsAndCapture(allocator, paths.fasta_path);
        defer allocator.free(output);
        try expectComposition(case.composition, output);
    }
}

test "[cli] - [stats report]: renders the exact protein report" {
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

test "[cli] - [stats report]: counts the complete protein alphabet and exceptional symbols" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const paths = try writeSingleLineFaiFixture(
        allocator,
        "stats-protein-categories",
        "protein",
        "ARNDCEQGHILKMFPSTWYVBZJXUOar*?",
    );
    defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, paths.fai_path) catch {};

    const expected =
        \\Composition:
        \\  type: protein
        \\  percent_denominator: total_symbols
        \\  a_alanine: 2 6.67%
        \\  r_arginine: 2 6.67%
        \\  n_asparagine: 1 3.33%
        \\  d_aspartate: 1 3.33%
        \\  c_cysteine: 1 3.33%
        \\  e_glutamate: 1 3.33%
        \\  q_glutamine: 1 3.33%
        \\  g_glycine: 1 3.33%
        \\  h_histidine: 1 3.33%
        \\  i_isoleucine: 1 3.33%
        \\  l_leucine: 1 3.33%
        \\  k_lysine: 1 3.33%
        \\  m_methionine: 1 3.33%
        \\  f_phenylalanine: 1 3.33%
        \\  p_proline: 1 3.33%
        \\  s_serine: 1 3.33%
        \\  t_threonine: 1 3.33%
        \\  w_tryptophan: 1 3.33%
        \\  y_tyrosine: 1 3.33%
        \\  v_valine: 1 3.33%
        \\  b_asx: 1 3.33%
        \\  z_glx: 1 3.33%
        \\  j_xle: 1 3.33%
        \\  x_unknown: 1 3.33%
        \\  u_selenocysteine: 1 3.33%
        \\  o_pyrrolysine: 1 3.33%
        \\  stop: 1 3.33%
        \\  invalid: 1 3.33%
        \\  lowercase: 2 6.67%
        \\
    ;

    const output = try runStatsAndCapture(allocator, paths.fasta_path);
    defer allocator.free(output);
    try expectComposition(expected, output);
}
