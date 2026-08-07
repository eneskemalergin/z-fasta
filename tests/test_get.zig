//! GET unit and CLI tests: region parsing, extraction, and failure-path subprocess checks.
//!
//! Includes indexed extraction and exact-exit CLI failure contracts.

const std = @import("std");
const builtin = @import("builtin");
const main = @import("main");
const parseRegion = main.getter.parseRegion;
const resolveRegion = main.getter.resolveRegion;
const io = std.testing.io;

const ZFASTA_BIN = if (builtin.os.tag == .windows) "zig-out\\bin\\z-fasta.exe" else "zig-out/bin/z-fasta";

// ============================================================================
// Region parsing tests
// ============================================================================

test "parseRegion - simple name" {
    const r = parseRegion("chr1");
    try std.testing.expectEqualStrings("chr1", r.name);
    try std.testing.expect(r.is_full);
}

test "parseRegion - empty identifier" {
    const r = parseRegion("");
    try std.testing.expectEqualStrings("", r.name);
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
    defer idx.deinit(io);

    var first_byte: [1]u8 = undefined;
    try std.testing.expectEqual(
        @as(usize, 1),
        try std.Io.File.readPositionalAll(idx.fasta_file, io, &first_byte, 0),
    );
    try std.testing.expectEqual(@as(u8, '>'), first_byte[0]);

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
// Multi-region resolution tests
// ============================================================================

test "resolveRegion - single region, full sequence" {
    var idx = main.index_format.loadIndex(io, "tests/data/simple.fasta");
    defer idx.deinit(io);

    const r = resolveRegion(&idx, "seq1", 0);
    try std.testing.expectEqualStrings("seq1", r.name);
    try std.testing.expect(r.is_full);
    try std.testing.expectEqual(@as(u64, 24), r.num_bases);
    try std.testing.expectEqual(@as(usize, 0), r.original_index);
}

test "resolveRegion - single region, sub-range" {
    var idx = main.index_format.loadIndex(io, "tests/data/simple.fasta");
    defer idx.deinit(io);

    const r = resolveRegion(&idx, "seq1:1-12", 0);
    try std.testing.expectEqualStrings("seq1", r.name);
    try std.testing.expect(!r.is_full);
    try std.testing.expectEqual(@as(u64, 12), r.num_bases);
    try std.testing.expectEqual(@as(u64, 1), r.start);
    try std.testing.expectEqual(@as(u64, 12), r.display_end);
}

test "resolveRegion - original_index preserved" {
    var idx = main.index_format.loadIndex(io, "tests/data/simple.fasta");
    defer idx.deinit(io);

    const r0 = resolveRegion(&idx, "seq1", 0);
    const r1 = resolveRegion(&idx, "seq2", 1);
    const r2 = resolveRegion(&idx, "seq1:1-5", 2);

    try std.testing.expectEqual(@as(usize, 0), r0.original_index);
    try std.testing.expectEqual(@as(usize, 1), r1.original_index);
    try std.testing.expectEqual(@as(usize, 2), r2.original_index);
}

test "resolveRegion - end clamped silently" {
    var idx = main.index_format.loadIndex(io, "tests/data/simple.fasta");
    defer idx.deinit(io);

    // seq1 has 24 bases; request end=9999 should clamp to 24
    const r = resolveRegion(&idx, "seq1:1-9999", 0);
    try std.testing.expectEqual(@as(u64, 24), r.num_bases);
    // display_end should be the user-supplied value (before clamping)
    try std.testing.expectEqual(@as(u64, 9999), r.display_end);
}

test "resolveRegion - byte offset for first base" {
    var idx = main.index_format.loadIndex(io, "tests/data/simple.fasta");
    defer idx.deinit(io);

    // seq1 starts immediately after ">seq1 test sequence\n"
    // Verify start_byte is the offset of the first base character
    const r = resolveRegion(&idx, "seq1:1-1", 0);
    try std.testing.expectEqual(@as(u64, 1), r.num_bases);
    // The byte at start_byte in fasta_data should be 'A'
    try std.testing.expectEqual(@as(u8, 'A'), idx.fasta_data[r.start_byte]);
}

test "resolveRegion - duplicate region allowed, same start_byte" {
    var idx = main.index_format.loadIndex(io, "tests/data/simple.fasta");
    defer idx.deinit(io);

    const r0 = resolveRegion(&idx, "seq1:1-5", 0);
    const r1 = resolveRegion(&idx, "seq1:1-5", 1);

    // Same region twice should resolve to identical start_byte and num_bases
    try std.testing.expectEqual(r0.start_byte, r1.start_byte);
    try std.testing.expectEqual(r0.num_bases, r1.num_bases);
    try std.testing.expectEqual(@as(usize, 0), r0.original_index);
    try std.testing.expectEqual(@as(usize, 1), r1.original_index);
}

// ============================================================================
// resolveRegion edge cases
// ============================================================================

test "resolveRegion - open-ended region (NAME:START-) uses seq_len as end" {
    var idx = main.index_format.loadIndex(io, "tests/data/simple.fasta");
    defer idx.deinit(io);

    // seq1 has 24 bases; NAME:13- should return bases 13..24 = 12 bases
    const r = resolveRegion(&idx, "seq1:13-", 0);
    try std.testing.expect(!r.is_full);
    try std.testing.expectEqual(@as(u64, 12), r.num_bases);
    // display_end should be seq_len (24), not null
    try std.testing.expectEqual(@as(u64, 24), r.display_end);
}

test "resolveRegion - single-base region" {
    var idx = main.index_format.loadIndex(io, "tests/data/simple.fasta");
    defer idx.deinit(io);

    const r = resolveRegion(&idx, "seq1:1-1", 0);
    try std.testing.expectEqual(@as(u64, 1), r.num_bases);
    // First base of seq1 should be 'A'
    try std.testing.expectEqual(@as(u8, 'A'), idx.fasta_data[r.start_byte]);
}

test "resolveRegion - last-base region" {
    var idx = main.index_format.loadIndex(io, "tests/data/simple.fasta");
    defer idx.deinit(io);

    // seq1: ACGTACGTACGTACGTACGTACGT (24 bases), last base = 'T'
    const r = resolveRegion(&idx, "seq1:24-24", 0);
    try std.testing.expectEqual(@as(u64, 1), r.num_bases);
    // Last base of seq1 should be 'T'
    const byte = idx.fasta_data[r.start_byte];
    try std.testing.expectEqual(@as(u8, 'T'), byte);
}

test "resolveRegion - cross-line region (starts line 1, ends line 2)" {
    var idx = main.index_format.loadIndex(io, "tests/data/simple.fasta");
    defer idx.deinit(io);

    // seq1 wraps at 12 bases per line; region 10-15 crosses the line boundary
    const r = resolveRegion(&idx, "seq1:10-15", 0);
    try std.testing.expectEqual(@as(u64, 6), r.num_bases);
}

test "resolveRegion - full sequence start_byte points to first base character" {
    var idx = main.index_format.loadIndex(io, "tests/data/simple.fasta");
    defer idx.deinit(io);

    const r_full = resolveRegion(&idx, "seq1", 0);
    const r_range = resolveRegion(&idx, "seq1:1-24", 1);
    // Both forms should land on the same start byte
    try std.testing.expectEqual(r_full.start_byte, r_range.start_byte);
    try std.testing.expectEqual(r_full.num_bases, r_range.num_bases);
}

test "resolveRegion - display_end before clamp, num_bases after clamp" {
    var idx = main.index_format.loadIndex(io, "tests/data/simple.fasta");
    defer idx.deinit(io);

    const r = resolveRegion(&idx, "seq1:1-9999", 0);
    // num_bases clamped to seq_len
    try std.testing.expectEqual(@as(u64, 24), r.num_bases);
    // display_end preserves user value
    try std.testing.expectEqual(@as(u64, 9999), r.display_end);
}

test "resolveRegion - proteome pipe-delimited name" {
    var idx = main.index_format.loadIndex(io, "tests/data/proteome.fasta");
    defer idx.deinit(io);

    const r = resolveRegion(&idx, "sp|P12345|PROT_HUMAN:1-10", 0);
    try std.testing.expectEqualStrings("sp|P12345|PROT_HUMAN", r.name);
    try std.testing.expectEqual(@as(u64, 10), r.num_bases);
}

test "resolveRegion - long header name (200-char sequence name)" {
    var idx = main.index_format.loadIndex(io, "tests/data/edge_cases.fasta");
    defer idx.deinit(io);

    const long_name = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    const r = resolveRegion(&idx, long_name, 0);
    try std.testing.expectEqualStrings(long_name, r.name);
    try std.testing.expect(r.is_full);
    try std.testing.expectEqual(@as(u64, 8), r.num_bases);
}

test "resolveRegion - lowercase bases preserved (byte offset still correct)" {
    var idx = main.index_format.loadIndex(io, "tests/data/edge_cases.fasta");
    defer idx.deinit(io);

    // 'lowercase' in edge_cases.fasta: acgtACGTacgt (12 bases)
    const r = resolveRegion(&idx, "lowercase:1-1", 0);
    try std.testing.expectEqual(@as(u64, 1), r.num_bases);
    // First base should be 'a' (lowercase)
    try std.testing.expectEqual(@as(u8, 'a'), idx.fasta_data[r.start_byte]);
}

test "resolveRegion - ordering: seq2 has higher file offset than seq1" {
    var idx = main.index_format.loadIndex(io, "tests/data/simple.fasta");
    defer idx.deinit(io);

    const r1 = resolveRegion(&idx, "seq1:1-1", 0);
    const r2 = resolveRegion(&idx, "seq2:1-1", 1);
    // seq2 appears after seq1 in the file, so its start_byte must be greater
    try std.testing.expect(r2.start_byte > r1.start_byte);
}

test "resolveRegion - reversed CLI order: seq2 before seq1 in args" {
    var idx = main.index_format.loadIndex(io, "tests/data/simple.fasta");
    defer idx.deinit(io);

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
    defer idx.deinit(io);

    // 'nonstandard' has ACG*-NACGT (10 chars)
    const r = resolveRegion(&idx, "nonstandard:1-10", 0);
    try std.testing.expectEqual(@as(u64, 10), r.num_bases);
    try std.testing.expectEqual(@as(u8, 'A'), idx.fasta_data[r.start_byte]);
}

// ============================================================================
// Index-backed get fixtures
// ============================================================================

fn uniqueArtifactPath(allocator: std.mem.Allocator, stem: []const u8, ext: []const u8) ![]u8 {
    try std.Io.Dir.cwd().createDirPath(io, "zig-cache/test-artifacts");
    const now = std.Io.Clock.Timestamp.now(io, .awake);
    const nanos: u64 = @intCast(now.raw.toNanoseconds());
    return std.fmt.allocPrint(allocator, "zig-cache/test-artifacts/{s}-{d}.{s}", .{
        stem,
        nanos,
        ext,
    });
}

fn writeFastaArtifact(allocator: std.mem.Allocator, stem: []const u8, data: []const u8) ![]const u8 {
    const path = try uniqueArtifactPath(allocator, stem, "fa");
    const file = try std.Io.Dir.cwd().createFile(io, path, .{});
    defer file.close(io);
    try std.Io.File.writeStreamingAll(file, io, data);
    return path;
}

fn writeZfi(allocator: std.mem.Allocator, fasta_path: []const u8, data: []const u8, enable_dedup: bool) !void {
    var arena = std.heap.ArenaAllocator.init(allocator);
    defer arena.deinit();

    var index = try main.indexer.scanZfiData(data, enable_dedup, arena.allocator());
    defer index.deinit(arena.allocator());

    const fasta_file = try std.Io.Dir.cwd().openFile(io, fasta_path, .{});
    defer fasta_file.close(io);
    const mtime_ns = main.index_format.timestampToNs((try fasta_file.stat(io)).mtime);

    var zfi_path_buf: [4096]u8 = undefined;
    const zfi_path = try std.fmt.bufPrint(&zfi_path_buf, "{s}.zfi", .{fasta_path});
    try main.indexer.writeZfiIndexFile(io, zfi_path, &index, data.len, mtime_ns);
}

fn captureExtractRegion(allocator: std.mem.Allocator, fasta_path: []const u8, region: []const u8) ![]u8 {
    var idx = try main.index_format.loadIndexChecked(io, fasta_path);
    defer idx.deinit(io);

    var out = std.Io.Writer.Allocating.init(allocator);
    main.getter.extractRegion(&idx, region, &out.writer);
    return out.toOwnedSlice();
}

test "get on zfi output for simple.fasta seq1:1-10" {
    const data = @embedFile("data/simple.fasta");
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const path = try writeFastaArtifact(allocator, "get-sam", data);
    const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{path});
    defer std.Io.Dir.cwd().deleteFile(io, path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, zfi_path) catch {};

    try writeZfi(allocator, path, data, true);

    const got = try captureExtractRegion(allocator, path, "seq1:1-10");
    const expected =
        \\>seq1:1-10
        \\ACGTACGTAC
        \\
    ;
    try std.testing.expectEqualStrings(expected, got);
}

test "get on messy mixed_widths after indexing" {
    const data =
        \\>mixed_widths internal line widths vary
        \\AAAACCCCGGGG
        \\TTTTAA
        \\AACCCCGGGGTT
        \\TT
        \\
    ;
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const path = try writeFastaArtifact(allocator, "get-messy", data);
    const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{path});
    defer std.Io.Dir.cwd().deleteFile(io, path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, zfi_path) catch {};

    try writeZfi(allocator, path, data, true);

    const region = "mixed_widths:3-24";
    const output = try captureExtractRegion(allocator, path, region);

    const expected =
        \\>mixed_widths:3-24
        \\AACCCCGGGGTTTTAAAACCCC
        \\
    ;
    try std.testing.expectEqualStrings(expected, output);
}

test "CLI get selects the first empty identifier through zfi and fai" {
    const data = ">\nAAAA\n>\nCCCCCC\n";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_path = try writeFastaArtifact(allocator, "get-empty-name", data);
    const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{fasta_path});
    const fai_path = try std.fmt.allocPrint(allocator, "{s}.fai", .{fasta_path});
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, zfi_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, fai_path) catch {};

    try writeZfi(allocator, fasta_path, data, false);
    try expectCliSuccess(allocator, &.{ ZFASTA_BIN, "get", fasta_path, "" }, ">\nAAAA\n");

    try std.Io.Dir.cwd().deleteFile(io, zfi_path);
    {
        const fai_file = try std.Io.Dir.cwd().createFile(io, fai_path, .{});
        defer fai_file.close(io);
        try std.Io.File.writeStreamingAll(fai_file, io, "\t4\t2\t4\t5\n\t6\t9\t6\t7\n");
    }
    try expectCliSuccess(allocator, &.{ ZFASTA_BIN, "get", fasta_path, "" }, ">\nAAAA\n");
}

/// Hard CLI failures: exact exit code, exact stderr, empty stdout (no partial success).
fn expectCliFailure(
    allocator: std.mem.Allocator,
    argv: []const []const u8,
    exit_code: u8,
    expected_stderr: []const u8,
) !void {
    var threaded = std.Io.Threaded.init(allocator, .{});
    defer threaded.deinit();
    const spawn_io = threaded.io();

    const result = try std.process.run(allocator, spawn_io, .{
        .argv = argv,
        .stdout_limit = .limited(64 * 1024),
        .stderr_limit = .limited(64 * 1024),
    });
    defer allocator.free(result.stdout);
    defer allocator.free(result.stderr);

    switch (result.term) {
        .exited => |code| try std.testing.expectEqual(exit_code, code),
        else => return error.ChildProcessFailed,
    }
    try std.testing.expectEqual(@as(usize, 0), result.stdout.len);
    try std.testing.expectEqualStrings(expected_stderr, result.stderr);
}

fn expectCliSuccess(allocator: std.mem.Allocator, argv: []const []const u8, expected_stdout: []const u8) !void {
    var threaded = std.Io.Threaded.init(allocator, .{});
    defer threaded.deinit();

    const result = try std.process.run(allocator, threaded.io(), .{
        .argv = argv,
        .stdout_limit = .limited(128 * 1024),
        .stderr_limit = .limited(64 * 1024),
    });
    defer allocator.free(result.stdout);
    defer allocator.free(result.stderr);

    switch (result.term) {
        .exited => |code| try std.testing.expectEqual(@as(u8, 0), code),
        else => return error.ChildProcessFailed,
    }
    try std.testing.expectEqualStrings(expected_stdout, result.stdout);
    try std.testing.expectEqual(@as(usize, 0), result.stderr.len);
}

/// Usage dumps: exit 1, empty stdout, stderr matches `z-fasta --help` stdout exactly.
fn expectCliUsageFailure(allocator: std.mem.Allocator, argv: []const []const u8) !void {
    var threaded = std.Io.Threaded.init(allocator, .{});
    defer threaded.deinit();
    const spawn_io = threaded.io();

    const help = try std.process.run(allocator, spawn_io, .{
        .argv = &.{ ZFASTA_BIN, "--help" },
        .stdout_limit = .limited(64 * 1024),
        .stderr_limit = .limited(64 * 1024),
    });
    defer allocator.free(help.stdout);
    defer allocator.free(help.stderr);
    switch (help.term) {
        .exited => |code| try std.testing.expectEqual(@as(u8, 0), code),
        else => return error.ChildProcessFailed,
    }
    try std.testing.expectEqual(@as(usize, 0), help.stderr.len);
    try std.testing.expect(help.stdout.len > 0);

    const result = try std.process.run(allocator, spawn_io, .{
        .argv = argv,
        .stdout_limit = .limited(64 * 1024),
        .stderr_limit = .limited(64 * 1024),
    });
    defer allocator.free(result.stdout);
    defer allocator.free(result.stderr);

    switch (result.term) {
        .exited => |code| try std.testing.expectEqual(@as(u8, 1), code),
        else => return error.ChildProcessFailed,
    }
    try std.testing.expectEqual(@as(usize, 0), result.stdout.len);
    try std.testing.expectEqualStrings(help.stdout, result.stderr);
}

fn expectUnknownOptionRejected(allocator: std.mem.Allocator, argv: []const []const u8, unknown: []const u8) !void {
    const expected = try std.fmt.allocPrint(allocator, "error: unknown option: {s}\n", .{unknown});
    try expectCliFailure(allocator, argv, 1, expected);
}

test "get help describes the request source contract" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    var threaded = std.Io.Threaded.init(allocator, .{});
    defer threaded.deinit();

    const result = try std.process.run(allocator, threaded.io(), .{
        .argv = &.{ ZFASTA_BIN, "get", "--help" },
        .stdout_limit = .limited(64 * 1024),
        .stderr_limit = .limited(64 * 1024),
    });
    defer allocator.free(result.stdout);
    defer allocator.free(result.stderr);

    switch (result.term) {
        .exited => |code| try std.testing.expectEqual(@as(u8, 0), code),
        else => return error.ChildProcessFailed,
    }
    try std.testing.expectEqual(@as(usize, 0), result.stderr.len);
    try std.testing.expect(std.mem.indexOf(u8, result.stdout, "--chunk-size") == null);
    try std.testing.expect(std.mem.indexOf(u8, result.stdout, "Choose exactly one request source") != null);
    try std.testing.expect(std.mem.indexOf(u8, result.stdout, "Names and BED input stream.") != null);
}

const get_usage_stderr =
    \\error: usage: z-fasta get <file.fasta> (<region>... | --names file.txt|- | --bed file.bed|-)
    \\
;

test "index rejects unknown options before and after FASTA path" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const fasta = "tests/data/simple.fasta";
    const unknown = "--not-a-flag";

    try expectUnknownOptionRejected(allocator, &.{ ZFASTA_BIN, "index", unknown, fasta }, unknown);
    try expectUnknownOptionRejected(allocator, &.{ ZFASTA_BIN, "index", fasta, unknown }, unknown);
    try expectUnknownOptionRejected(allocator, &.{ ZFASTA_BIN, "index", "-z", fasta }, "-z");
    try expectUnknownOptionRejected(allocator, &.{ ZFASTA_BIN, "index", "--low-mem", fasta }, "--low-mem");
}

test "get rejects unknown options before, between, and after positionals" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const fasta = "tests/data/simple.fasta";
    const unknown = "--not-a-flag";

    try expectUnknownOptionRejected(allocator, &.{ ZFASTA_BIN, "get", unknown, fasta, "seq1" }, unknown);
    try expectUnknownOptionRejected(allocator, &.{ ZFASTA_BIN, "get", fasta, unknown, "seq1" }, unknown);
    try expectUnknownOptionRejected(allocator, &.{ ZFASTA_BIN, "get", fasta, "seq1", unknown }, unknown);
    try expectUnknownOptionRejected(allocator, &.{ ZFASTA_BIN, "get", fasta, "--chunk-size", "1", "seq1" }, "--chunk-size");
    try expectUnknownOptionRejected(allocator, &.{ ZFASTA_BIN, "get", fasta, "--bed", "tests/data/simple.bed", "--chunk-size=-1" }, "--chunk-size=-1");
}

test "get rejects repeated request source flags" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const fasta = "tests/data/simple.fasta";

    try expectCliFailure(
        allocator,
        &.{ ZFASTA_BIN, "get", fasta, "--bed", "first.bed", "--bed", "second.bed" },
        1,
        "error: --bed may be specified only once\n",
    );
    try expectCliFailure(
        allocator,
        &.{ ZFASTA_BIN, "get", fasta, "--names", "first.txt", "--names", "second.txt" },
        1,
        "error: --names may be specified only once\n",
    );
}

test "get rejects mixed request source families" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const fasta = "tests/data/simple.fasta";
    const expected = "error: choose one request source: positional regions, --names, or --bed\n";

    try expectCliFailure(allocator, &.{ ZFASTA_BIN, "get", fasta, "seq1", "--names", "ids.txt" }, 1, expected);
    try expectCliFailure(allocator, &.{ ZFASTA_BIN, "get", fasta, "seq1", "--bed", "regions.bed" }, 1, expected);
    try expectCliFailure(allocator, &.{ ZFASTA_BIN, "get", fasta, "--names", "ids.txt", "--bed", "regions.bed" }, 1, expected);
}

test "get rejects strand handling outside BED input" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const fasta = "tests/data/simple.fasta";
    const expected = "error: --strand-aware requires --bed\n";

    try expectCliFailure(allocator, &.{ ZFASTA_BIN, "get", fasta, "seq1", "--strand-aware" }, 1, expected);
    try expectCliFailure(allocator, &.{ ZFASTA_BIN, "get", fasta, "--names", "ids.txt", "--honor-strand" }, 1, expected);
}

test "names accepts the index name-length boundary in both formats" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const name = try allocator.alloc(u8, main.indexer.max_index_name_len);
    @memset(name, 'A');
    var fasta = std.ArrayList(u8).empty;
    try fasta.ensureTotalCapacity(allocator, name.len + 4);
    try fasta.append(allocator, '>');
    try fasta.appendSlice(allocator, name);
    try fasta.appendSlice(allocator, "\nA\n");
    const names_path = try uniqueArtifactPath(allocator, "names-max", "txt");
    defer std.Io.Dir.cwd().deleteFile(io, names_path) catch {};
    {
        const names_file = try std.Io.Dir.cwd().createFile(io, names_path, .{});
        defer names_file.close(io);
        try std.Io.File.writeStreamingAll(names_file, io, name);
        try std.Io.File.writeStreamingAll(names_file, io, "\r\n");
    }
    var expected = std.ArrayList(u8).empty;
    try expected.ensureTotalCapacity(allocator, name.len + 5);
    try expected.append(allocator, '>');
    try expected.appendSlice(allocator, name);
    try expected.appendSlice(allocator, "\nA\n");

    const zfi_fasta = try writeFastaArtifact(allocator, "names-max-zfi", fasta.items);
    const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{zfi_fasta});
    defer std.Io.Dir.cwd().deleteFile(io, zfi_fasta) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, zfi_path) catch {};
    try writeZfi(allocator, zfi_fasta, fasta.items, true);

    const fai_fasta = try writeFastaArtifact(allocator, "names-max-fai", fasta.items);
    const fai_path = try std.fmt.allocPrint(allocator, "{s}.fai", .{fai_fasta});
    defer std.Io.Dir.cwd().deleteFile(io, fai_fasta) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, fai_path) catch {};
    {
        const fai_file = try std.Io.Dir.cwd().createFile(io, fai_path, .{});
        defer fai_file.close(io);
        var buf: [128]u8 = undefined;
        const fields = try std.fmt.bufPrint(&buf, "\t1\t{d}\t1\t2\n", .{name.len + 2});
        try std.Io.File.writeStreamingAll(fai_file, io, name);
        try std.Io.File.writeStreamingAll(fai_file, io, fields);
    }

    try expectCliSuccess(allocator, &.{ ZFASTA_BIN, "get", zfi_fasta, "--names", names_path }, expected.items);
    try expectCliSuccess(allocator, &.{ ZFASTA_BIN, "get", fai_fasta, "--names", names_path }, expected.items);
}

test "names rejects a request name above the index limit" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const names_path = try uniqueArtifactPath(allocator, "names-too-long", "txt");
    defer std.Io.Dir.cwd().deleteFile(io, names_path) catch {};
    const names_file = try std.Io.Dir.cwd().createFile(io, names_path, .{});
    const name = try allocator.alloc(u8, main.indexer.max_index_name_len + 1);
    @memset(name, 'A');
    try std.Io.File.writeStreamingAll(names_file, io, name);
    names_file.close(io);

    try expectCliFailure(
        allocator,
        &.{ ZFASTA_BIN, "get", "tests/data/simple.fasta", "--names", names_path },
        1,
        "error: name at line 1 exceeds 65535 bytes\n",
    );
}

test "stats rejects unknown options before and after FASTA path" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const fasta = "tests/data/simple.fasta";
    const unknown = "--not-a-flag";

    try expectUnknownOptionRejected(allocator, &.{ ZFASTA_BIN, "stats", unknown, fasta }, unknown);
    try expectUnknownOptionRejected(allocator, &.{ ZFASTA_BIN, "stats", fasta, unknown }, unknown);
    try expectUnknownOptionRejected(allocator, &.{ ZFASTA_BIN, "stats", "--index-only", fasta }, "--index-only");
}

test "validate rejects unknown options before and after FASTA path" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const fasta = "tests/data/simple.fasta";
    const unknown = "--not-a-flag";

    try expectUnknownOptionRejected(allocator, &.{ ZFASTA_BIN, "validate", unknown, fasta }, unknown);
    try expectUnknownOptionRejected(allocator, &.{ ZFASTA_BIN, "validate", fasta, unknown }, unknown);
}

test "CLI failures: missing command and missing subcommand args print usage on stderr" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    try expectCliUsageFailure(allocator, &.{ZFASTA_BIN});
    try expectCliUsageFailure(allocator, &.{ ZFASTA_BIN, "index" });
    try expectCliUsageFailure(allocator, &.{ ZFASTA_BIN, "not-a-command" });
}

test "CLI failures: missing FASTA path is exit 1 with exact stderr" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const missing = "tests/data/definitely-missing-cli-failure.fasta";
    const expected = "error: file not found: tests/data/definitely-missing-cli-failure.fasta\n";

    try expectCliFailure(allocator, &.{ ZFASTA_BIN, "index", missing }, 1, expected);
    try expectCliFailure(allocator, &.{ ZFASTA_BIN, "stats", missing }, 1, expected);
    try expectCliFailure(allocator, &.{ ZFASTA_BIN, "validate", missing }, 1, expected);
    try expectCliFailure(allocator, &.{ ZFASTA_BIN, "get", missing, "seq1" }, 1, expected);
}

test "get validates FASTA and sidecar before opening a request source" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    try expectCliFailure(
        allocator,
        &.{ ZFASTA_BIN, "get", "tests/data/definitely-missing-source-order.fasta", "--names", "also-missing.txt" },
        1,
        "error: file not found: tests/data/definitely-missing-source-order.fasta\n",
    );
    try expectCliFailure(
        allocator,
        &.{ ZFASTA_BIN, "get", "tests/data/simple.fasta", "--names", "tests/data/also-missing.txt" },
        1,
        "error: file not found: tests/data/also-missing.txt\n",
    );
}

test "CLI failures: invalid present zfi blocks valid fai for get and stats" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta = try writeFastaArtifact(allocator, "invalid-zfi-cli", ">seq\nACGT\n");
    const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{fasta});
    const fai_path = try std.fmt.allocPrint(allocator, "{s}.fai", .{fasta});
    defer std.Io.Dir.cwd().deleteFile(io, fasta) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, zfi_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, fai_path) catch {};

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
    try expectCliFailure(allocator, &.{ ZFASTA_BIN, "get", fasta, "seq" }, 1, expected);
    try expectCliFailure(allocator, &.{ ZFASTA_BIN, "get", fasta, "--names", "missing.txt" }, 1, expected);
    try expectCliFailure(allocator, &.{ ZFASTA_BIN, "stats", fasta }, 1, expected);
}

test "streamed names failures suppress summary and expose only valid output prefixes" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const early_path = try uniqueArtifactPath(allocator, "names-early-failure", "txt");
    const late_path = try uniqueArtifactPath(allocator, "names-late-failure", "txt");
    defer std.Io.Dir.cwd().deleteFile(io, early_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, late_path) catch {};
    {
        const file = try std.Io.Dir.cwd().createFile(io, early_path, .{});
        defer file.close(io);
        try std.Io.File.writeStreamingAll(file, io, "missing\nseq1\n");
    }
    var late_names = std.ArrayList(u8).empty;
    try late_names.ensureTotalCapacity(allocator, 5 * (65536 + 1) + 3);
    for (0..65536) |_| try late_names.appendSlice(allocator, "seq1\n");
    try late_names.appendSlice(allocator, "missing\n");
    {
        const file = try std.Io.Dir.cwd().createFile(io, late_path, .{});
        defer file.close(io);
        try std.Io.File.writeStreamingAll(file, io, late_names.items);
    }

    try expectCliFailure(
        allocator,
        &.{ ZFASTA_BIN, "get", "tests/data/simple.fasta", "--names", early_path, "--summary" },
        1,
        "error: sequence not found: missing\n",
    );

    var threaded = std.Io.Threaded.init(allocator, .{});
    defer threaded.deinit();
    const result = try std.process.run(allocator, threaded.io(), .{
        .argv = &.{ ZFASTA_BIN, "get", "tests/data/simple.fasta", "--names", late_path, "--summary" },
        .stdout_limit = .limited(4 * 1024 * 1024),
        .stderr_limit = .limited(64 * 1024),
    });
    defer allocator.free(result.stdout);
    defer allocator.free(result.stderr);

    switch (result.term) {
        .exited => |code| try std.testing.expectEqual(@as(u8, 1), code),
        else => return error.ChildProcessFailed,
    }
    try std.testing.expect(result.stdout.len > 0);
    try std.testing.expectEqualStrings("error: sequence not found: missing\n", result.stderr);
    const one_output = ">seq1\nACGTACGTACGTACGTACGTACGT\n";
    for (result.stdout, 0..) |byte, i| {
        try std.testing.expectEqual(one_output[i % one_output.len], byte);
    }
}

test "CLI failures: get usage, conflicts, and missing sequence" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta = try writeFastaArtifact(allocator, "get-cli-errors", @embedFile("data/simple.fasta"));
    const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{fasta});
    defer std.Io.Dir.cwd().deleteFile(io, fasta) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, zfi_path) catch {};

    {
        var threaded = std.Io.Threaded.init(allocator, .{});
        defer threaded.deinit();
        const spawn_io = threaded.io();
        const indexed = try std.process.run(allocator, spawn_io, .{
            .argv = &.{ ZFASTA_BIN, "index", fasta },
            .stdout_limit = .limited(64 * 1024),
            .stderr_limit = .limited(64 * 1024),
        });
        defer allocator.free(indexed.stdout);
        defer allocator.free(indexed.stderr);
        switch (indexed.term) {
            .exited => |code| try std.testing.expectEqual(@as(u8, 0), code),
            else => return error.ChildProcessFailed,
        }
    }

    try expectCliFailure(allocator, &.{ ZFASTA_BIN, "get", fasta }, 1, get_usage_stderr);
    try expectCliFailure(
        allocator,
        &.{ ZFASTA_BIN, "get", fasta, "missing" },
        1,
        "error: sequence not found: missing\n",
    );
    try expectCliFailure(
        allocator,
        &.{ ZFASTA_BIN, "get", fasta, "seq1", "--annotate-rc" },
        1,
        "error: --annotate-rc requires a transform or --strand-aware BED input\n",
    );
    try expectCliFailure(
        allocator,
        &.{ ZFASTA_BIN, "get", fasta, "--rc", "--complement-only", "seq1" },
        1,
        "error: --rc, --complement-only, and --reverse-only are mutually exclusive\n",
    );
}

test "BED annotations describe final composed orientation" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const bed_path = try uniqueArtifactPath(allocator, "bed-annotation", "bed");
    defer std.Io.Dir.cwd().deleteFile(io, bed_path) catch {};
    const bed_file = try std.Io.Dir.cwd().createFile(io, bed_path, .{});
    try std.Io.File.writeStreamingAll(
        bed_file,
        io,
        "seq1\t0\t5\tplus\t0\t+\nseq1\t0\t5\tminus\t0\t-\n",
    );
    bed_file.close(io);

    try expectCliSuccess(
        allocator,
        &.{ ZFASTA_BIN, "get", "tests/data/simple.fasta", "--bed", bed_path, "--strand-aware", "--annotate-rc" },
        ">seq1:1-5\nACGTA\n>seq1:1-5 (reverse complement)\nTACGT\n",
    );
    try expectCliSuccess(
        allocator,
        &.{ ZFASTA_BIN, "get", "tests/data/simple.fasta", "--bed", bed_path, "--strand-aware", "--rc", "--annotate-rc" },
        ">seq1:1-5 (reverse complement)\nTACGT\n>seq1:1-5\nACGTA\n",
    );
}

test "CLI failures: stats and validate usage errors" {
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
        &.{ ZFASTA_BIN, "validate", "--summary", "tests/data/simple.fasta" },
        1,
        "error: validate --summary requires --json\n",
    );
    try expectCliFailure(
        allocator,
        &.{ ZFASTA_BIN, "validate" },
        1,
        "error: usage: z-fasta validate [options] <file.fasta>\n",
    );
}

test "CLI validate warnings: exit 2, empty stderr, exact WARNING stdout" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const path = try writeFastaArtifact(allocator, "cli-validate-warn", ">empty_rec\n");
    defer std.Io.Dir.cwd().deleteFile(io, path) catch {};

    var threaded = std.Io.Threaded.init(allocator, .{});
    defer threaded.deinit();
    const spawn_io = threaded.io();

    const result = try std.process.run(allocator, spawn_io, .{
        .argv = &.{ ZFASTA_BIN, "validate", path },
        .stdout_limit = .limited(64 * 1024),
        .stderr_limit = .limited(64 * 1024),
    });
    defer allocator.free(result.stdout);
    defer allocator.free(result.stderr);

    switch (result.term) {
        .exited => |code| try std.testing.expectEqual(@as(u8, 2), code),
        else => return error.ChildProcessFailed,
    }
    try std.testing.expectEqual(@as(usize, 0), result.stderr.len);
    try std.testing.expectEqualStrings("WARNING: line 1: empty sequence 'empty_rec'\n", result.stdout);
}
