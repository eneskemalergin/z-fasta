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
