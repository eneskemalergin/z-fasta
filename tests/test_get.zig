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
// Low-mem index get
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

fn writeZfi(allocator: std.mem.Allocator, fasta_path: []const u8, data: []const u8) !void {
    var arena = std.heap.ArenaAllocator.init(allocator);
    defer arena.deinit();

    var index = try main.indexer.scanZfiData(data, true, arena.allocator());
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

    try writeZfi(allocator, path, data);

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

    try writeZfi(allocator, path, data);

    const region = "mixed_widths:3-24";
    const output = try captureExtractRegion(allocator, path, region);

    const expected =
        \\>mixed_widths:3-24
        \\AACCCCGGGGTTTTAAAACCCC
        \\
    ;
    try std.testing.expectEqualStrings(expected, output);
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

const get_usage_stderr =
    \\error: usage: z-fasta get <file.fasta> [--bed file.bed|-] [--names file.txt] [--strand-aware] [--summary] [--rc|--complement-only|--reverse-only] [--annotate-rc] [--chunk-size N|-1] <region> [region ...]
    \\
;

test "--names rejects file over max_input_file_bytes with clear error" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const names_path = try uniqueArtifactPath(allocator, "names-oversize", "txt");
    defer std.Io.Dir.cwd().deleteFile(io, names_path) catch {};

    const names_file = try std.Io.Dir.cwd().createFile(io, names_path, .{});
    try names_file.setLength(io, main.getter.max_input_file_bytes + 1);
    names_file.close(io);

    const expected = try std.fmt.allocPrint(
        allocator,
        "error: names file exceeds {d} MiB limit: {s} (--chunk-size does not stream --names)\n",
        .{ main.getter.max_input_file_mib, names_path },
    );
    try expectCliFailure(
        allocator,
        &.{ ZFASTA_BIN, "get", "tests/data/simple.fasta", "--names", names_path },
        1,
        expected,
    );
}

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
}

test "stats rejects unknown options before and after FASTA path" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const fasta = "tests/data/simple.fasta";
    const unknown = "--not-a-flag";

    try expectUnknownOptionRejected(allocator, &.{ ZFASTA_BIN, "stats", unknown, fasta }, unknown);
    try expectUnknownOptionRejected(allocator, &.{ ZFASTA_BIN, "stats", fasta, unknown }, unknown);
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
        "error: --annotate-rc requires --rc, --complement-only, or --reverse-only\n",
    );
    try expectCliFailure(
        allocator,
        &.{ ZFASTA_BIN, "get", fasta, "--rc", "--complement-only", "seq1" },
        1,
        "error: --rc, --complement-only, and --reverse-only are mutually exclusive\n",
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
        "error: usage: z-fasta stats [--index-only] <file.fasta>\n",
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
