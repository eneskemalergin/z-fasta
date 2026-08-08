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

// --- Region parsing tests ---

test "parseRegion recognizes names and coordinate suffixes" {
    const Case = struct {
        input: []const u8,
        name: []const u8,
        start: u64 = 1,
        end: ?u64 = null,
        is_full: bool = true,
    };
    const cases = [_]Case{
        .{ .input = "chr1", .name = "chr1" },
        .{ .input = "", .name = "" },
        .{ .input = "chr1:100-200", .name = "chr1", .start = 100, .end = 200, .is_full = false },
        .{ .input = "chr1:100-", .name = "chr1", .start = 100, .is_full = false },
        .{
            .input = "chromosome:GRCh38:1:1:248956422:1:100-200",
            .name = "chromosome:GRCh38:1:1:248956422:1",
            .start = 100,
            .end = 200,
            .is_full = false,
        },
        .{
            .input = "chromosome:GRCh38:1:1:248956422:1",
            .name = "chromosome:GRCh38:1:1:248956422:1",
        },
        .{ .input = "sp|P12345|PROT_NAME:1-50", .name = "sp|P12345|PROT_NAME", .end = 50, .is_full = false },
        .{ .input = "chr1:1-1", .name = "chr1", .end = 1, .is_full = false },
        .{ .input = "chr1:1000000-2000000", .name = "chr1", .start = 1_000_000, .end = 2_000_000, .is_full = false },
        .{ .input = "KI270394.1:1-100", .name = "KI270394.1", .end = 100, .is_full = false },
        .{ .input = "1:1", .name = "1:1" },
        .{ .input = "chr1.1.2.3", .name = "chr1.1.2.3" },
    };

    for (cases) |case| {
        const region = parseRegion(case.input);

        try std.testing.expectEqualStrings(case.name, region.name);
        try std.testing.expectEqual(case.start, region.start);
        try std.testing.expectEqual(case.end, region.end);
        try std.testing.expectEqual(case.is_full, region.is_full);
    }
}

// --- Index loading tests ---

test "loadIndex - .zfi file" {
    var idx = main.index_format.loadIndex(std.testing.allocator, io, "tests/data/simple.fasta");
    defer idx.deinit();

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

// --- Multi-region resolution tests ---

test "resolveRegion - single region, full sequence" {
    var idx = main.index_format.loadIndex(std.testing.allocator, io, "tests/data/simple.fasta");
    defer idx.deinit();

    const r = resolveRegion(&idx, "seq1");

    try std.testing.expectEqualStrings("seq1", r.name);
    try std.testing.expect(r.is_full);
    try std.testing.expectEqual(@as(u64, 24), r.num_bases);
}

test "resolveRegion - single region, sub-range" {
    var idx = main.index_format.loadIndex(std.testing.allocator, io, "tests/data/simple.fasta");
    defer idx.deinit();

    const r = resolveRegion(&idx, "seq1:1-12");

    try std.testing.expectEqualStrings("seq1", r.name);
    try std.testing.expect(!r.is_full);
    try std.testing.expectEqual(@as(u64, 12), r.num_bases);
    try std.testing.expectEqual(@as(u64, 1), r.start);
    try std.testing.expectEqual(@as(u64, 12), r.display_end);
}

test "resolveRegion - end clamped silently" {
    var idx = main.index_format.loadIndex(std.testing.allocator, io, "tests/data/simple.fasta");
    defer idx.deinit();

    const r = resolveRegion(&idx, "seq1:1-9999");

    try std.testing.expectEqual(@as(u64, 24), r.num_bases);
    try std.testing.expectEqual(@as(u64, 9999), r.display_end);
}

test "resolveRegion - first base resolves one symbol" {
    var idx = main.index_format.loadIndex(std.testing.allocator, io, "tests/data/simple.fasta");
    defer idx.deinit();

    const r = resolveRegion(&idx, "seq1:1-1");

    try std.testing.expectEqual(@as(u64, 1), r.num_bases);
}

// --- resolveRegion edge cases ---

test "resolveRegion - open-ended region (NAME:START-) uses seq_len as end" {
    var idx = main.index_format.loadIndex(std.testing.allocator, io, "tests/data/simple.fasta");
    defer idx.deinit();

    const r = resolveRegion(&idx, "seq1:13-");

    try std.testing.expect(!r.is_full);
    try std.testing.expectEqual(@as(u64, 12), r.num_bases);
    try std.testing.expectEqual(@as(u64, 24), r.display_end);
}

test "resolveRegion - last-base region" {
    var idx = main.index_format.loadIndex(std.testing.allocator, io, "tests/data/simple.fasta");
    defer idx.deinit();

    const r = resolveRegion(&idx, "seq1:24-24");

    try std.testing.expectEqual(@as(u64, 1), r.num_bases);
}

test "resolveRegion - full and explicit ranges preserve geometry" {
    var idx = main.index_format.loadIndex(std.testing.allocator, io, "tests/data/simple.fasta");
    defer idx.deinit();

    const r_full = resolveRegion(&idx, "seq1");
    const r_range = resolveRegion(&idx, "seq1:1-24");

    try std.testing.expectEqual(r_full.seq_offset, r_range.seq_offset);
    try std.testing.expectEqual(r_full.num_bases, r_range.num_bases);
}

test "resolveRegion - proteome pipe-delimited name" {
    var idx = main.index_format.loadIndex(std.testing.allocator, io, "tests/data/proteome.fasta");
    defer idx.deinit();

    const r = resolveRegion(&idx, "sp|P12345|PROT_HUMAN:1-10");

    try std.testing.expectEqualStrings("sp|P12345|PROT_HUMAN", r.name);
    try std.testing.expectEqual(@as(u64, 10), r.num_bases);
}

test "resolveRegion - long header name (200-char sequence name)" {
    var idx = main.index_format.loadIndex(std.testing.allocator, io, "tests/data/edge_cases.fasta");
    defer idx.deinit();

    const long_name = "A" ** 200;
    const r = resolveRegion(&idx, long_name);

    try std.testing.expectEqualStrings(long_name, r.name);
    try std.testing.expect(r.is_full);
    try std.testing.expectEqual(@as(u64, 8), r.num_bases);
}

// --- Index-backed get fixtures ---

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
    var idx = try main.index_format.loadIndexChecked(allocator, io, fasta_path);
    defer idx.deinit();

    var out = std.Io.Writer.Allocating.init(allocator);
    errdefer out.deinit();
    main.getter.extractRegion(&idx, region, &out.writer);
    return out.toOwnedSlice();
}

test "get extracts a region across uniform FASTA lines" {
    const data = @embedFile("data/simple.fasta");
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const path = try writeFastaArtifact(allocator, "get-sam", data);
    const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{path});
    defer std.Io.Dir.cwd().deleteFile(io, path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, zfi_path) catch {};

    try writeZfi(allocator, path, data, true);

    const got = try captureExtractRegion(allocator, path, "seq1:10-15");
    const expected =
        \\>seq1:10-15
        \\CGTACG
        \\
    ;
    try std.testing.expectEqualStrings(expected, got);
}

test "get preserves lowercase and nonstandard sequence bytes" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const lowercase = try captureExtractRegion(allocator, "tests/data/edge_cases.fasta", "lowercase");
    const nonstandard = try captureExtractRegion(allocator, "tests/data/edge_cases.fasta", "nonstandard");

    try std.testing.expectEqualStrings(">lowercase\nacgtACGTacgt\n", lowercase);
    try std.testing.expectEqualStrings(">nonstandard\nACG*-NACGT\n", nonstandard);
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

// Hard CLI failures require exact exit code and stderr with no partial stdout.
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

// Usage failures must reproduce the top-level help on stderr exactly.
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
    defer allocator.free(expected);
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

const GET_USAGE_STDERR =
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

test "index rejects multiple FASTA paths" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();

    try expectCliFailure(
        arena.allocator(),
        &.{ ZFASTA_BIN, "index", "tests/data/simple.fasta", "tests/data/single.fasta" },
        1,
        "error: index accepts exactly one FASTA path\n",
    );
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

    const name = try allocator.alloc(u8, main.indexer.MAX_INDEX_NAME_LEN);
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
    const name = try allocator.alloc(u8, main.indexer.MAX_INDEX_NAME_LEN + 1);
    @memset(name, 'A');
    {
        const names_file = try std.Io.Dir.cwd().createFile(io, names_path, .{});
        defer names_file.close(io);
        try std.Io.File.writeStreamingAll(names_file, io, name);
    }

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

test "stats rejects multiple FASTA paths" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();

    try expectCliFailure(
        arena.allocator(),
        &.{ ZFASTA_BIN, "stats", "tests/data/simple.fasta", "tests/data/single.fasta" },
        1,
        "error: stats accepts exactly one FASTA path\n",
    );
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

    try expectCliFailure(allocator, &.{ ZFASTA_BIN, "get", fasta }, 1, GET_USAGE_STDERR);
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
    {
        const bed_file = try std.Io.Dir.cwd().createFile(io, bed_path, .{});
        defer bed_file.close(io);
        try std.Io.File.writeStreamingAll(
            bed_file,
            io,
            "seq1\t0\t5\tplus\t0\t+\nseq1\t0\t5\tminus\t0\t-\n",
        );
    }

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
        &.{ ZFASTA_BIN, "validate", "--fix-format-only", "tests/data/simple.fasta" },
        1,
        "error: validate --fix-format-only requires --fix\n",
    );
    try expectCliFailure(
        allocator,
        &.{ ZFASTA_BIN, "validate", "-o", "zig-cache/test-artifacts/unused.fasta", "tests/data/simple.fasta" },
        1,
        "error: validate -o requires --fix\n",
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
