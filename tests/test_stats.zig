const std = @import("std");
const main = @import("main");
const formatComma = main.stats.formatComma;
const formatSize = main.stats.formatSize;
const detectType = main.stats.detectType;
const SequenceType = main.stats.SequenceType;

// ============================================================================
// formatComma tests
// ============================================================================

test "formatComma — zero" {
    var buf: [64]u8 = undefined;
    try std.testing.expectEqualStrings("0", formatComma(&buf, 0));
}

test "formatComma — single digit" {
    var buf: [64]u8 = undefined;
    try std.testing.expectEqualStrings("5", formatComma(&buf, 5));
}

test "formatComma — three digits" {
    var buf: [64]u8 = undefined;
    try std.testing.expectEqualStrings("999", formatComma(&buf, 999));
}

test "formatComma — four digits" {
    var buf: [64]u8 = undefined;
    try std.testing.expectEqualStrings("1,000", formatComma(&buf, 1000));
}

test "formatComma — six digits" {
    var buf: [64]u8 = undefined;
    try std.testing.expectEqualStrings("999,999", formatComma(&buf, 999_999));
}

test "formatComma — seven digits" {
    var buf: [64]u8 = undefined;
    try std.testing.expectEqualStrings("1,000,000", formatComma(&buf, 1_000_000));
}

test "formatComma — large number (billions)" {
    var buf: [64]u8 = undefined;
    try std.testing.expectEqualStrings("3,099,750,718", formatComma(&buf, 3_099_750_718));
}

test "formatComma — max u32" {
    var buf: [64]u8 = undefined;
    try std.testing.expectEqualStrings("4,294,967,295", formatComma(&buf, 4_294_967_295));
}

test "formatComma — large u64" {
    var buf: [64]u8 = undefined;
    try std.testing.expectEqualStrings("1,000,000,000,000", formatComma(&buf, 1_000_000_000_000));
}

test "formatComma — 12345" {
    var buf: [64]u8 = undefined;
    try std.testing.expectEqualStrings("12,345", formatComma(&buf, 12345));
}

test "formatComma — 100" {
    var buf: [64]u8 = undefined;
    try std.testing.expectEqualStrings("100", formatComma(&buf, 100));
}

// ============================================================================
// formatSize tests
// ============================================================================

test "formatSize — bytes" {
    var buf: [64]u8 = undefined;
    try std.testing.expectEqualStrings("512 B", formatSize(&buf, 512));
}

test "formatSize — kilobytes" {
    var buf: [64]u8 = undefined;
    const result = formatSize(&buf, 10 * 1024);
    try std.testing.expectEqualStrings("10.0 KB", result);
}

test "formatSize — megabytes" {
    var buf: [64]u8 = undefined;
    const result = formatSize(&buf, 66 * 1024 * 1024);
    try std.testing.expectEqualStrings("66.0 MB", result);
}

test "formatSize — gigabytes" {
    var buf: [64]u8 = undefined;
    const result = formatSize(&buf, 3 * 1024 * 1024 * 1024);
    try std.testing.expectEqualStrings("3.0 GB", result);
}

// ============================================================================
// detectType tests
// ============================================================================

test "detectType — all ACGT is nucleotide" {
    var counts: [256]u64 = .{0} ** 256;
    counts['A'] = 250;
    counts['C'] = 250;
    counts['G'] = 250;
    counts['T'] = 250;
    try std.testing.expectEqual(SequenceType.nucleotide, detectType(&counts, 1000));
}

test "detectType — ACGTN is nucleotide" {
    var counts: [256]u64 = .{0} ** 256;
    counts['A'] = 200;
    counts['C'] = 200;
    counts['G'] = 200;
    counts['T'] = 200;
    counts['N'] = 200;
    try std.testing.expectEqual(SequenceType.nucleotide, detectType(&counts, 1000));
}

test "detectType — mixed amino acids is protein" {
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

test "detectType — below 90% threshold is protein" {
    var counts: [256]u64 = .{0} ** 256;
    counts['A'] = 200;
    counts['C'] = 200;
    counts['G'] = 200;
    // A+C+G = 600 out of 700, which is 85.7% => protein
    counts['L'] = 50;
    counts['M'] = 50;
    try std.testing.expectEqual(SequenceType.protein, detectType(&counts, 700));
}

test "detectType — exactly 91% is nucleotide" {
    var counts: [256]u64 = .{0} ** 256;
    counts['A'] = 910;
    counts['L'] = 90;
    // ACGTN = 910/1000 = 91% > 90%
    try std.testing.expectEqual(SequenceType.nucleotide, detectType(&counts, 1000));
}

test "detectType — empty is nucleotide (default)" {
    var counts: [256]u64 = .{0} ** 256;
    try std.testing.expectEqual(SequenceType.nucleotide, detectType(&counts, 0));
}

test "detectType — lowercase nucleotides" {
    var counts: [256]u64 = .{0} ** 256;
    counts['a'] = 250;
    counts['c'] = 250;
    counts['g'] = 250;
    counts['t'] = 250;
    try std.testing.expectEqual(SequenceType.nucleotide, detectType(&counts, 1000));
}

// ============================================================================
// Integration: stats via process spawn
// ============================================================================

const ZFASTA_BIN = "zig-out/bin/z-fasta";

fn runStatsAndCapture(allocator: std.mem.Allocator, fasta_path: []const u8, index_only: bool) ![]u8 {
    var argv = std.ArrayList([]const u8).init(allocator);
    defer argv.deinit();
    try argv.append(ZFASTA_BIN);
    try argv.append("stats");
    if (index_only) try argv.append("--index-only");
    try argv.append(fasta_path);

    var proc = std.process.Child.init(argv.items, allocator);
    proc.stdout_behavior = .Pipe;
    proc.stderr_behavior = .Pipe;
    try proc.spawn();

    const result = try proc.stdout.?.readToEndAlloc(allocator, 10 * 1024 * 1024);
    _ = try proc.wait();
    return result;
}

test "stats — simple.fasta has 2 sequences" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const output = try runStatsAndCapture(arena.allocator(), "tests/data/simple.fasta", false);
    try std.testing.expect(std.mem.indexOf(u8, output, "Sequences:      2") != null);
}

test "stats — simple.fasta has 36 total bases" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const output = try runStatsAndCapture(arena.allocator(), "tests/data/simple.fasta", false);
    try std.testing.expect(std.mem.indexOf(u8, output, "Total bases:    36") != null);
}

test "stats — simple.fasta detected as Nucleotide" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const output = try runStatsAndCapture(arena.allocator(), "tests/data/simple.fasta", false);
    try std.testing.expect(std.mem.indexOf(u8, output, "Type:           Nucleotide") != null);
}

test "stats — proteome.fasta detected as Protein" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const output = try runStatsAndCapture(arena.allocator(), "tests/data/proteome.fasta", false);
    try std.testing.expect(std.mem.indexOf(u8, output, "Type:           Protein") != null);
}

test "stats — index-only does not show composition" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const output = try runStatsAndCapture(arena.allocator(), "tests/data/simple.fasta", true);
    try std.testing.expect(std.mem.indexOf(u8, output, "Composition:") == null);
}

test "stats — index-only shows placeholder type" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const output = try runStatsAndCapture(arena.allocator(), "tests/data/simple.fasta", true);
    try std.testing.expect(std.mem.indexOf(u8, output, "run without --index-only") != null);
}

test "stats — simple.fasta N50=24 L50=1" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const output = try runStatsAndCapture(arena.allocator(), "tests/data/simple.fasta", true);
    try std.testing.expect(std.mem.indexOf(u8, output, "N50:            24") != null);
    try std.testing.expect(std.mem.indexOf(u8, output, "L50:            1") != null);
}

test "stats — simple.fasta GC content" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const output = try runStatsAndCapture(arena.allocator(), "tests/data/simple.fasta", false);
    try std.testing.expect(std.mem.indexOf(u8, output, "GC:  55.56%") != null);
}

test "stats — mixed_widths.fasta has 3 sequences" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const output = try runStatsAndCapture(arena.allocator(), "tests/data/mixed_widths.fasta", true);
    try std.testing.expect(std.mem.indexOf(u8, output, "Sequences:      3") != null);
}

test "stats — proteome Composition shows amino acids" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const output = try runStatsAndCapture(arena.allocator(), "tests/data/proteome.fasta", false);
    try std.testing.expect(std.mem.indexOf(u8, output, "Alanine") != null);
    try std.testing.expect(std.mem.indexOf(u8, output, "20 amino acids total") != null);
}
