const std = @import("std");
const main = @import("main");
const parseRegion = main.getter.parseRegion;

// ============================================================================
// Region parsing tests
// ============================================================================

test "parseRegion — simple name" {
    const r = parseRegion("chr1");
    try std.testing.expectEqualStrings("chr1", r.name);
    try std.testing.expect(r.is_full);
}

test "parseRegion — name with range" {
    const r = parseRegion("chr1:100-200");
    try std.testing.expectEqualStrings("chr1", r.name);
    try std.testing.expectEqual(@as(u64, 100), r.start);
    try std.testing.expectEqual(@as(?u64, 200), r.end);
    try std.testing.expect(!r.is_full);
}

test "parseRegion — name with open end" {
    const r = parseRegion("chr1:100-");
    try std.testing.expectEqualStrings("chr1", r.name);
    try std.testing.expectEqual(@as(u64, 100), r.start);
    try std.testing.expectEqual(@as(?u64, null), r.end);
    try std.testing.expect(!r.is_full);
}

test "parseRegion — Ensembl colon name with range" {
    const r = parseRegion("chromosome:GRCh38:1:1:248956422:1:100-200");
    try std.testing.expectEqualStrings("chromosome:GRCh38:1:1:248956422:1", r.name);
    try std.testing.expectEqual(@as(u64, 100), r.start);
    try std.testing.expectEqual(@as(?u64, 200), r.end);
    try std.testing.expect(!r.is_full);
}

test "parseRegion — Ensembl colon name without range" {
    const r = parseRegion("chromosome:GRCh38:1:1:248956422:1");
    try std.testing.expectEqualStrings("chromosome:GRCh38:1:1:248956422:1", r.name);
    try std.testing.expect(r.is_full);
}

test "parseRegion — pipe-delimited protein name" {
    const r = parseRegion("sp|P12345|PROT_NAME:1-50");
    try std.testing.expectEqualStrings("sp|P12345|PROT_NAME", r.name);
    try std.testing.expectEqual(@as(u64, 1), r.start);
    try std.testing.expectEqual(@as(?u64, 50), r.end);
}

test "parseRegion — single base" {
    const r = parseRegion("chr1:1-1");
    try std.testing.expectEqualStrings("chr1", r.name);
    try std.testing.expectEqual(@as(u64, 1), r.start);
    try std.testing.expectEqual(@as(?u64, 1), r.end);
    try std.testing.expect(!r.is_full);
}

test "parseRegion — large coordinates" {
    const r = parseRegion("chr1:1000000-2000000");
    try std.testing.expectEqualStrings("chr1", r.name);
    try std.testing.expectEqual(@as(u64, 1_000_000), r.start);
    try std.testing.expectEqual(@as(?u64, 2_000_000), r.end);
}

test "parseRegion — name with underscore and range" {
    const r = parseRegion("KI270394.1:1-100");
    try std.testing.expectEqualStrings("KI270394.1", r.name);
    try std.testing.expectEqual(@as(u64, 1), r.start);
    try std.testing.expectEqual(@as(?u64, 100), r.end);
}

test "parseRegion — name with colon but no valid range" {
    // "1:1" looks like it could be parsed as name="1", start=1, end=null ...
    // but there's no dash, so parseRangeSuffix returns null.
    // The whole thing should be treated as a name.
    const r = parseRegion("1:1");
    try std.testing.expectEqualStrings("1:1", r.name);
    try std.testing.expect(r.is_full);
}

test "parseRegion — name only with dots" {
    const r = parseRegion("chr1.1.2.3");
    try std.testing.expectEqualStrings("chr1.1.2.3", r.name);
    try std.testing.expect(r.is_full);
}

// ============================================================================
// Index loading tests
// ============================================================================

test "loadIndex — .zfi file" {
    var idx = main.index_format.loadIndex("tests/data/simple.fasta");
    defer idx.deinit();

    try std.testing.expectEqual(@as(usize, 2), idx.records.len);
    try std.testing.expectEqual(@as(u64, 24), idx.records[0].seq_len);
    try std.testing.expectEqual(@as(u64, 12), idx.records[1].seq_len);

    const seq1_idx = idx.lookupName("seq1");
    try std.testing.expect(seq1_idx != null);
    try std.testing.expectEqual(@as(usize, 0), seq1_idx.?);

    const seq2_idx = idx.lookupName("seq2");
    try std.testing.expect(seq2_idx != null);
    try std.testing.expectEqual(@as(usize, 1), seq2_idx.?);

    try std.testing.expectEqual(@as(?usize, null), idx.lookupName("nonexistent"));
}

// ============================================================================
// Integration: get via process spawn
// ============================================================================

const ZFASTA_BIN = "zig-out/bin/z-fasta";

fn runGetAndCapture(allocator: std.mem.Allocator, args: []const []const u8) ![]u8 {
    var argv = std.ArrayList([]const u8).init(allocator);
    defer argv.deinit();
    try argv.append(ZFASTA_BIN);
    try argv.append("get");
    for (args) |a| try argv.append(a);

    var proc = std.process.Child.init(argv.items, allocator);
    proc.stdout_behavior = .Pipe;
    proc.stderr_behavior = .Pipe;
    try proc.spawn();

    var stdout_data = std.ArrayList(u8).init(allocator);
    var stderr_data = std.ArrayList(u8).init(allocator);

    // Read both stdout and stderr
    try collectOutput(proc.stdout.?, &stdout_data);
    try collectOutput(proc.stderr.?, &stderr_data);

    const term = try proc.wait();
    if (term.Exited != 0) {
        std.debug.print("z-fasta exited with code {d}: {s}\n", .{ term.Exited, stderr_data.items });
        return error.ProcessFailed;
    }

    return stdout_data.toOwnedSlice();
}

fn runSamtoolsAndCapture(allocator: std.mem.Allocator, args: []const []const u8) ![]u8 {
    var argv = std.ArrayList([]const u8).init(allocator);
    defer argv.deinit();
    try argv.append("samtools");
    try argv.append("faidx");
    for (args) |a| try argv.append(a);

    var proc = std.process.Child.init(argv.items, allocator);
    proc.stdout_behavior = .Pipe;
    proc.stderr_behavior = .Pipe;
    try proc.spawn();

    var stdout_data = std.ArrayList(u8).init(allocator);
    var stderr_data = std.ArrayList(u8).init(allocator);

    try collectOutput(proc.stdout.?, &stdout_data);
    try collectOutput(proc.stderr.?, &stderr_data);

    const term = try proc.wait();
    if (term.Exited != 0) {
        std.debug.print("samtools exited with code {d}: {s}\n", .{ term.Exited, stderr_data.items });
        return error.ProcessFailed;
    }

    return stdout_data.toOwnedSlice();
}

fn collectOutput(reader: std.fs.File, list: *std.ArrayList(u8)) !void {
    var buf: [4096]u8 = undefined;
    while (true) {
        const n = try reader.read(&buf);
        if (n == 0) break;
        try list.appendSlice(buf[0..n]);
    }
}

fn expectSameOutput(allocator: std.mem.Allocator, args: []const []const u8) !void {
    const zfasta = try runGetAndCapture(allocator, args);
    defer allocator.free(zfasta);
    const samtools = try runSamtoolsAndCapture(allocator, args);
    defer allocator.free(samtools);

    try std.testing.expectEqualStrings(samtools, zfasta);
}

test "get — full sequence matches samtools (seq1)" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    try expectSameOutput(arena.allocator(), &.{ "tests/data/simple.fasta", "seq1" });
}

test "get — full sequence matches samtools (seq2)" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    try expectSameOutput(arena.allocator(), &.{ "tests/data/simple.fasta", "seq2" });
}

test "get — region matches samtools (seq1:1-10)" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    try expectSameOutput(arena.allocator(), &.{ "tests/data/simple.fasta", "seq1:1-10" });
}

test "get — region matches samtools (seq1:5-18)" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    try expectSameOutput(arena.allocator(), &.{ "tests/data/simple.fasta", "seq1:5-18" });
}

test "get — single base matches samtools" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    try expectSameOutput(arena.allocator(), &.{ "tests/data/simple.fasta", "seq1:1-1" });
}

test "get — last base matches samtools" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    try expectSameOutput(arena.allocator(), &.{ "tests/data/simple.fasta", "seq1:24-24" });
}

test "get — full via region syntax matches samtools" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    try expectSameOutput(arena.allocator(), &.{ "tests/data/simple.fasta", "seq1:1-24" });
}

test "get — across line boundary matches samtools" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    try expectSameOutput(arena.allocator(), &.{ "tests/data/simple.fasta", "seq1:10-15" });
}

test "get — proteome pipe name matches samtools" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    try expectSameOutput(arena.allocator(), &.{ "tests/data/proteome.fasta", "sp|P12345|PROT_HUMAN" });
}

test "get — proteome region matches samtools" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    try expectSameOutput(arena.allocator(), &.{ "tests/data/proteome.fasta", "sp|P12345|PROT_HUMAN:1-10" });
}

test "get — single.fasta full matches samtools" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    try expectSameOutput(arena.allocator(), &.{ "tests/data/single.fasta", "single_sequence" });
}

test "get — mixed_widths mixed1 full matches samtools" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    try expectSameOutput(arena.allocator(), &.{ "tests/data/mixed_widths.fasta", "mixed1" });
}

test "get — mixed_widths mixed2 region matches samtools" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    try expectSameOutput(arena.allocator(), &.{ "tests/data/mixed_widths.fasta", "mixed2:100-500" });
}

test "get — mixed_widths mixed3 last base matches samtools" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    try expectSameOutput(arena.allocator(), &.{ "tests/data/mixed_widths.fasta", "mixed3:1658-1658" });
}

test "get — edge_cases lowercase matches samtools" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    try expectSameOutput(arena.allocator(), &.{ "tests/data/edge_cases.fasta", "lowercase" });
}

test "get — edge_cases single_line matches samtools" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    try expectSameOutput(arena.allocator(), &.{ "tests/data/edge_cases.fasta", "single_line" });
}

test "get — clamp END > seq_len header matches samtools" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    try expectSameOutput(arena.allocator(), &.{ "tests/data/simple.fasta", "seq1:1-1000" });
}
