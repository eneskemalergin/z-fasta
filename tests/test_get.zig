const std = @import("std");
const main = @import("main");
const parseRegion = main.getter.parseRegion;
const resolveRegion = main.getter.resolveRegion;
const io = std.Io.Threaded.global_single_threaded.io();

// ============================================================================
// Region parsing tests
// ============================================================================

test "parseRegion - simple name" {
    const r = parseRegion("chr1");
    try std.testing.expectEqualStrings("chr1", r.name);
    try std.testing.expect(r.is_full);
}

test "parseRegion - name with range" {
    const r = parseRegion("chr1:100-200");
    try std.testing.expectEqualStrings("chr1", r.name);
    try std.testing.expectEqual(@as(u64, 100), r.start);
    try std.testing.expectEqual(@as(?u64, 200), r.end);
    try std.testing.expect(!r.is_full);
}

test "parseRegion - name with open end" {
    const r = parseRegion("chr1:100-");
    try std.testing.expectEqualStrings("chr1", r.name);
    try std.testing.expectEqual(@as(u64, 100), r.start);
    try std.testing.expectEqual(@as(?u64, null), r.end);
    try std.testing.expect(!r.is_full);
}

test "parseRegion - Ensembl colon name with range" {
    const r = parseRegion("chromosome:GRCh38:1:1:248956422:1:100-200");
    try std.testing.expectEqualStrings("chromosome:GRCh38:1:1:248956422:1", r.name);
    try std.testing.expectEqual(@as(u64, 100), r.start);
    try std.testing.expectEqual(@as(?u64, 200), r.end);
    try std.testing.expect(!r.is_full);
}

test "parseRegion - Ensembl colon name without range" {
    const r = parseRegion("chromosome:GRCh38:1:1:248956422:1");
    try std.testing.expectEqualStrings("chromosome:GRCh38:1:1:248956422:1", r.name);
    try std.testing.expect(r.is_full);
}

test "parseRegion - pipe-delimited protein name" {
    const r = parseRegion("sp|P12345|PROT_NAME:1-50");
    try std.testing.expectEqualStrings("sp|P12345|PROT_NAME", r.name);
    try std.testing.expectEqual(@as(u64, 1), r.start);
    try std.testing.expectEqual(@as(?u64, 50), r.end);
}

test "parseRegion - single base" {
    const r = parseRegion("chr1:1-1");
    try std.testing.expectEqualStrings("chr1", r.name);
    try std.testing.expectEqual(@as(u64, 1), r.start);
    try std.testing.expectEqual(@as(?u64, 1), r.end);
    try std.testing.expect(!r.is_full);
}

test "parseRegion - large coordinates" {
    const r = parseRegion("chr1:1000000-2000000");
    try std.testing.expectEqualStrings("chr1", r.name);
    try std.testing.expectEqual(@as(u64, 1_000_000), r.start);
    try std.testing.expectEqual(@as(?u64, 2_000_000), r.end);
}

test "parseRegion - name with underscore and range" {
    const r = parseRegion("KI270394.1:1-100");
    try std.testing.expectEqualStrings("KI270394.1", r.name);
    try std.testing.expectEqual(@as(u64, 1), r.start);
    try std.testing.expectEqual(@as(?u64, 100), r.end);
}

test "parseRegion - name with colon but no valid range" {
    // "1:1" looks like it could be parsed as name="1", start=1, end=null ...
    // but there's no dash, so parseRangeSuffix returns null.
    // The whole thing should be treated as a name.
    const r = parseRegion("1:1");
    try std.testing.expectEqualStrings("1:1", r.name);
    try std.testing.expect(r.is_full);
}

test "parseRegion - name only with dots" {
    const r = parseRegion("chr1.1.2.3");
    try std.testing.expectEqualStrings("chr1.1.2.3", r.name);
    try std.testing.expect(r.is_full);
}

// ============================================================================
// Index loading tests
// ============================================================================

test "loadIndex - .zfi file" {
    var idx = main.index_format.loadIndex(io, "tests/data/simple.fasta");
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
// Multi-region resolution tests (v0.2.4)
// ============================================================================

test "resolveRegion - single region, full sequence" {
    var idx = main.index_format.loadIndex(io, "tests/data/simple.fasta");
    defer idx.deinit();

    const r = resolveRegion(&idx, "seq1", 0);
    try std.testing.expectEqualStrings("seq1", r.name);
    try std.testing.expect(r.is_full);
    try std.testing.expectEqual(@as(u64, 24), r.num_bases);
    try std.testing.expectEqual(@as(usize, 0), r.original_index);
}

test "resolveRegion - single region, sub-range" {
    var idx = main.index_format.loadIndex(io, "tests/data/simple.fasta");
    defer idx.deinit();

    const r = resolveRegion(&idx, "seq1:1-12", 0);
    try std.testing.expectEqualStrings("seq1", r.name);
    try std.testing.expect(!r.is_full);
    try std.testing.expectEqual(@as(u64, 12), r.num_bases);
    try std.testing.expectEqual(@as(u64, 1), r.start);
    try std.testing.expectEqual(@as(u64, 12), r.display_end);
}

test "resolveRegion - original_index preserved" {
    var idx = main.index_format.loadIndex(io, "tests/data/simple.fasta");
    defer idx.deinit();

    const r0 = resolveRegion(&idx, "seq1", 0);
    const r1 = resolveRegion(&idx, "seq2", 1);
    const r2 = resolveRegion(&idx, "seq1:1-5", 2);

    try std.testing.expectEqual(@as(usize, 0), r0.original_index);
    try std.testing.expectEqual(@as(usize, 1), r1.original_index);
    try std.testing.expectEqual(@as(usize, 2), r2.original_index);
}

test "resolveRegion - end clamped silently" {
    var idx = main.index_format.loadIndex(io, "tests/data/simple.fasta");
    defer idx.deinit();

    // seq1 has 24 bases; request end=9999 should clamp to 24
    const r = resolveRegion(&idx, "seq1:1-9999", 0);
    try std.testing.expectEqual(@as(u64, 24), r.num_bases);
    // display_end should be the user-supplied value (before clamping)
    try std.testing.expectEqual(@as(u64, 9999), r.display_end);
}

test "resolveRegion - byte offset for first base" {
    var idx = main.index_format.loadIndex(io, "tests/data/simple.fasta");
    defer idx.deinit();

    // seq1 starts immediately after ">seq1 test sequence\n"
    // Verify start_byte is the offset of the first base character
    const r = resolveRegion(&idx, "seq1:1-1", 0);
    try std.testing.expectEqual(@as(u64, 1), r.num_bases);
    // The byte at start_byte in fasta_data should be 'A'
    try std.testing.expectEqual(@as(u8, 'A'), idx.fasta_data[r.start_byte]);
}

test "resolveRegion - duplicate region allowed, same start_byte" {
    var idx = main.index_format.loadIndex(io, "tests/data/simple.fasta");
    defer idx.deinit();

    const r0 = resolveRegion(&idx, "seq1:1-5", 0);
    const r1 = resolveRegion(&idx, "seq1:1-5", 1);

    // Same region twice should resolve to identical start_byte and num_bases
    try std.testing.expectEqual(r0.start_byte, r1.start_byte);
    try std.testing.expectEqual(r0.num_bases, r1.num_bases);
    try std.testing.expectEqual(@as(usize, 0), r0.original_index);
    try std.testing.expectEqual(@as(usize, 1), r1.original_index);
}

// ============================================================================
// Edge case tests for resolveRegion (v0.2.4 bug hunt)
// ============================================================================

test "resolveRegion - open-ended region (NAME:START-) uses seq_len as end" {
    var idx = main.index_format.loadIndex(io, "tests/data/simple.fasta");
    defer idx.deinit();

    // seq1 has 24 bases; NAME:13- should return bases 13..24 = 12 bases
    const r = resolveRegion(&idx, "seq1:13-", 0);
    try std.testing.expect(!r.is_full);
    try std.testing.expectEqual(@as(u64, 12), r.num_bases);
    // display_end should be seq_len (24), not null
    try std.testing.expectEqual(@as(u64, 24), r.display_end);
}

test "resolveRegion - single-base region" {
    var idx = main.index_format.loadIndex(io, "tests/data/simple.fasta");
    defer idx.deinit();

    const r = resolveRegion(&idx, "seq1:1-1", 0);
    try std.testing.expectEqual(@as(u64, 1), r.num_bases);
    // First base of seq1 should be 'A'
    try std.testing.expectEqual(@as(u8, 'A'), idx.fasta_data[r.start_byte]);
}

test "resolveRegion - last-base region" {
    var idx = main.index_format.loadIndex(io, "tests/data/simple.fasta");
    defer idx.deinit();

    // seq1: ACGTACGTACGTACGTACGTACGT (24 bases), last base = 'T'
    const r = resolveRegion(&idx, "seq1:24-24", 0);
    try std.testing.expectEqual(@as(u64, 1), r.num_bases);
    // Last base of seq1 should be 'T'
    const byte = idx.fasta_data[r.start_byte];
    try std.testing.expectEqual(@as(u8, 'T'), byte);
}

test "resolveRegion - cross-line region (starts line 1, ends line 2)" {
    var idx = main.index_format.loadIndex(io, "tests/data/simple.fasta");
    defer idx.deinit();

    // seq1 wraps at 12 bases per line; region 10-15 crosses the line boundary
    const r = resolveRegion(&idx, "seq1:10-15", 0);
    try std.testing.expectEqual(@as(u64, 6), r.num_bases);
}

test "resolveRegion - full sequence start_byte points to first base character" {
    var idx = main.index_format.loadIndex(io, "tests/data/simple.fasta");
    defer idx.deinit();

    const r_full = resolveRegion(&idx, "seq1", 0);
    const r_range = resolveRegion(&idx, "seq1:1-24", 1);
    // Both forms should land on the same start byte
    try std.testing.expectEqual(r_full.start_byte, r_range.start_byte);
    try std.testing.expectEqual(r_full.num_bases, r_range.num_bases);
}

test "resolveRegion - display_end before clamp, num_bases after clamp" {
    var idx = main.index_format.loadIndex(io, "tests/data/simple.fasta");
    defer idx.deinit();

    const r = resolveRegion(&idx, "seq1:1-9999", 0);
    // num_bases clamped to seq_len
    try std.testing.expectEqual(@as(u64, 24), r.num_bases);
    // display_end preserves user value
    try std.testing.expectEqual(@as(u64, 9999), r.display_end);
}

test "resolveRegion - proteome pipe-delimited name" {
    var idx = main.index_format.loadIndex(io, "tests/data/proteome.fasta");
    defer idx.deinit();

    const r = resolveRegion(&idx, "sp|P12345|PROT_HUMAN:1-10", 0);
    try std.testing.expectEqualStrings("sp|P12345|PROT_HUMAN", r.name);
    try std.testing.expectEqual(@as(u64, 10), r.num_bases);
}

test "resolveRegion - long header name (200-char sequence name)" {
    var idx = main.index_format.loadIndex(io, "tests/data/edge_cases.fasta");
    defer idx.deinit();

    const long_name = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    const r = resolveRegion(&idx, long_name, 0);
    try std.testing.expectEqualStrings(long_name, r.name);
    try std.testing.expect(r.is_full);
    try std.testing.expectEqual(@as(u64, 8), r.num_bases);
}

test "resolveRegion - lowercase bases preserved (byte offset still correct)" {
    var idx = main.index_format.loadIndex(io, "tests/data/edge_cases.fasta");
    defer idx.deinit();

    // 'lowercase' in edge_cases.fasta: acgtACGTacgt (12 bases)
    const r = resolveRegion(&idx, "lowercase:1-1", 0);
    try std.testing.expectEqual(@as(u64, 1), r.num_bases);
    // First base should be 'a' (lowercase)
    try std.testing.expectEqual(@as(u8, 'a'), idx.fasta_data[r.start_byte]);
}

test "resolveRegion - ordering: seq2 has higher file offset than seq1" {
    var idx = main.index_format.loadIndex(io, "tests/data/simple.fasta");
    defer idx.deinit();

    const r1 = resolveRegion(&idx, "seq1:1-1", 0);
    const r2 = resolveRegion(&idx, "seq2:1-1", 1);
    // seq2 appears after seq1 in the file, so its start_byte must be greater
    try std.testing.expect(r2.start_byte > r1.start_byte);
}

test "resolveRegion - reversed CLI order: seq2 before seq1 in args" {
    var idx = main.index_format.loadIndex(io, "tests/data/simple.fasta");
    defer idx.deinit();

    // Even if caller passes seq2 first, original_index tracks CLI position
    const r2 = resolveRegion(&idx, "seq2:1-1", 0); // position 0 in args
    const r1 = resolveRegion(&idx, "seq1:1-1", 1); // position 1 in args
    try std.testing.expectEqual(@as(usize, 0), r2.original_index);
    try std.testing.expectEqual(@as(usize, 1), r1.original_index);
    // But seq1 still has lower file offset
    try std.testing.expect(r1.start_byte < r2.start_byte);
}

test "resolveRegion - nonstandard characters in sequence (stars/dashes)" {
    var idx = main.index_format.loadIndex(io, "tests/data/edge_cases.fasta");
    defer idx.deinit();

    // 'nonstandard' has ACG*-NACGT (10 chars)
    const r = resolveRegion(&idx, "nonstandard:1-10", 0);
    try std.testing.expectEqual(@as(u64, 10), r.num_bases);
    try std.testing.expectEqual(@as(u8, 'A'), idx.fasta_data[r.start_byte]);
}
