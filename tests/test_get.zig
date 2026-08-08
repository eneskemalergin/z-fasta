//! GET integration and CLI tests: region resolution, extraction, and failure paths.
//!
//! Includes indexed extraction and exact-exit CLI failure contracts.

const std = @import("std");
const main = @import("main");
const utility = @import("utility.zig");
const resolveRegion = main.getter.resolveRegion;
const io = std.testing.io;

const ZFASTA_BIN = utility.ZFASTA_BIN;
const expectCliFailure = utility.expectCliFailure;
const expectCliSuccess = utility.expectCliSuccess;
const expectUnknownOptionRejected = utility.expectUnknownOptionRejected;
const uniqueArtifactPath = utility.uniqueArtifactPath;
const writeFastaArtifact = utility.writeFastaArtifact;
const writeZfi = utility.writeZfi;
const captureExtractRegion = utility.captureExtractRegion;

test "[integration] - [region resolution]: preserves indexed coordinate geometry" {
    var idx = main.index_format.loadIndex(std.testing.allocator, io, "tests/data/simple.fasta");
    defer idx.deinit();

    const Case = struct {
        request: []const u8,
        start: u64,
        display_end: u64,
        num_bases: u64,
        is_full: bool = false,
    };
    const cases = [_]Case{
        .{ .request = "seq1", .start = 1, .display_end = 24, .num_bases = 24, .is_full = true },
        .{ .request = "seq1:1-12", .start = 1, .display_end = 12, .num_bases = 12 },
        .{ .request = "seq1:1-9999", .start = 1, .display_end = 9999, .num_bases = 24 },
        .{ .request = "seq1:1-1", .start = 1, .display_end = 1, .num_bases = 1 },
        .{ .request = "seq1:13-", .start = 13, .display_end = 24, .num_bases = 12 },
        .{ .request = "seq1:24-24", .start = 24, .display_end = 24, .num_bases = 1 },
        .{ .request = "seq1:1-24", .start = 1, .display_end = 24, .num_bases = 24 },
    };

    for (cases) |case| {
        const resolved = resolveRegion(&idx, case.request);

        try std.testing.expectEqualStrings("seq1", resolved.name);
        try std.testing.expectEqual(case.start, resolved.start);
        try std.testing.expectEqual(case.display_end, resolved.display_end);
        try std.testing.expectEqual(case.num_bases, resolved.num_bases);
        try std.testing.expectEqual(case.is_full, resolved.is_full);
        try std.testing.expectEqual(idx.records[0].seq_offset, resolved.seq_offset);
    }
}

test "[integration] - [region resolution]: accepts a pipe-delimited protein identifier" {
    var idx = main.index_format.loadIndex(std.testing.allocator, io, "tests/data/proteome.fasta");
    defer idx.deinit();

    const r = resolveRegion(&idx, "sp|P12345|PROT_HUMAN:1-10");

    try std.testing.expectEqualStrings("sp|P12345|PROT_HUMAN", r.name);
    try std.testing.expectEqual(@as(u64, 10), r.num_bases);
}

test "[integration] - [region resolution]: accepts a 200-byte identifier" {
    var idx = main.index_format.loadIndex(std.testing.allocator, io, "tests/data/edge_cases.fasta");
    defer idx.deinit();

    const long_name = "A" ** 200;
    const r = resolveRegion(&idx, long_name);

    try std.testing.expectEqualStrings(long_name, r.name);
    try std.testing.expect(r.is_full);
    try std.testing.expectEqual(@as(u64, 8), r.num_bases);
}

// --- Index-backed get fixtures ---

test "[property] - [get extraction]: zfi and fai match across uniform FASTA lines" {
    const data = @embedFile("data/simple.fasta");
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const path = try writeFastaArtifact(allocator, "get-sam", data);
    defer std.Io.Dir.cwd().deleteFile(io, path) catch {};
    const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{path});
    const fai_path = try std.fmt.allocPrint(allocator, "{s}.fai", .{path});
    defer std.Io.Dir.cwd().deleteFile(io, zfi_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, fai_path) catch {};

    try writeZfi(allocator, path, data, true);

    const zfi_output = try captureExtractRegion(allocator, path, "seq1:10-15");
    const expected =
        \\>seq1:10-15
        \\CGTACG
        \\
    ;
    try std.testing.expectEqualStrings(expected, zfi_output);

    try std.Io.Dir.cwd().deleteFile(io, zfi_path);
    {
        const fai_file = try std.Io.Dir.cwd().createFile(io, fai_path, .{});
        defer fai_file.close(io);
        try std.Io.File.writeStreamingAll(
            fai_file,
            io,
            "seq1\t24\t20\t12\t13\nseq2\t12\t69\t12\t13\n",
        );
    }
    const fai_output = try captureExtractRegion(allocator, path, "seq1:10-15");

    try std.testing.expectEqualStrings(expected, fai_output);
    try std.testing.expectEqualStrings(zfi_output, fai_output);
}

test "[integration] - [get extraction]: preserves lowercase and nonstandard bytes" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const lowercase = try captureExtractRegion(allocator, "tests/data/edge_cases.fasta", "lowercase");
    const nonstandard = try captureExtractRegion(allocator, "tests/data/edge_cases.fasta", "nonstandard");

    try std.testing.expectEqualStrings(">lowercase\nacgtACGTacgt\n", lowercase);
    try std.testing.expectEqualStrings(">nonstandard\nACG*-NACGT\n", nonstandard);
}

test "[integration] - [get extraction]: extracts across non-uniform FASTA lines" {
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
    defer std.Io.Dir.cwd().deleteFile(io, path) catch {};
    const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{path});
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

test "[cli] - [get]: selects the first empty identifier through zfi and fai" {
    const data = ">\nAAAA\n>\nCCCCCC\n";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_path = try writeFastaArtifact(allocator, "get-empty-name", data);
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{fasta_path});
    const fai_path = try std.fmt.allocPrint(allocator, "{s}.fai", .{fasta_path});
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

test "[cli] - [get help]: describes the request source contract" {
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

test "[cli] - [get]: rejects unknown options regardless of position" {
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

test "[cli] - [get]: rejects repeated request source flags" {
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

test "[cli] - [get]: rejects mixed request source families" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const fasta = "tests/data/simple.fasta";
    const expected = "error: choose one request source: positional regions, --names, or --bed\n";

    try expectCliFailure(allocator, &.{ ZFASTA_BIN, "get", fasta, "seq1", "--names", "ids.txt" }, 1, expected);
    try expectCliFailure(allocator, &.{ ZFASTA_BIN, "get", fasta, "seq1", "--bed", "regions.bed" }, 1, expected);
    try expectCliFailure(allocator, &.{ ZFASTA_BIN, "get", fasta, "--names", "ids.txt", "--bed", "regions.bed" }, 1, expected);
}

test "[cli] - [get]: rejects strand handling outside BED input" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const fasta = "tests/data/simple.fasta";
    const expected = "error: --strand-aware requires --bed\n";

    try expectCliFailure(allocator, &.{ ZFASTA_BIN, "get", fasta, "seq1", "--strand-aware" }, 1, expected);
    try expectCliFailure(allocator, &.{ ZFASTA_BIN, "get", fasta, "--names", "ids.txt", "--honor-strand" }, 1, expected);
}

test "[cli] - [get names]: accepts the index name-length boundary in both formats" {
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
    defer std.Io.Dir.cwd().deleteFile(io, zfi_fasta) catch {};
    const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{zfi_fasta});
    defer std.Io.Dir.cwd().deleteFile(io, zfi_path) catch {};
    try writeZfi(allocator, zfi_fasta, fasta.items, true);

    const fai_fasta = try writeFastaArtifact(allocator, "names-max-fai", fasta.items);
    defer std.Io.Dir.cwd().deleteFile(io, fai_fasta) catch {};
    const fai_path = try std.fmt.allocPrint(allocator, "{s}.fai", .{fai_fasta});
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

test "[cli] - [get names]: rejects a request name above the index limit" {
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

test "[cli] - [get]: missing FASTA path returns an exact error" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const missing = "tests/data/definitely-missing-cli-failure.fasta";
    const expected = "error: file not found: tests/data/definitely-missing-cli-failure.fasta\n";

    try expectCliFailure(allocator, &.{ ZFASTA_BIN, "get", missing, "seq1" }, 1, expected);
}

test "[cli] - [get]: validates FASTA and sidecar before request input" {
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

test "[cli] - [get index selection]: invalid zfi blocks a valid fai" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta = try writeFastaArtifact(allocator, "invalid-zfi-cli", ">seq\nACGT\n");
    defer std.Io.Dir.cwd().deleteFile(io, fasta) catch {};
    const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{fasta});
    const fai_path = try std.fmt.allocPrint(allocator, "{s}.fai", .{fasta});
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
}

test "[cli] - [get names]: failures preserve only valid output prefixes" {
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
    try std.testing.expect(result.stdout.len <= one_output.len * 65536);
    for (result.stdout, 0..) |byte, i| {
        try std.testing.expectEqual(one_output[i % one_output.len], byte);
    }
}

test "[cli] - [get]: rejects invalid requests, missing names, and transform conflicts" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta = try writeFastaArtifact(allocator, "get-cli-errors", @embedFile("data/simple.fasta"));
    defer std.Io.Dir.cwd().deleteFile(io, fasta) catch {};
    const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{fasta});
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
        &.{ ZFASTA_BIN, "get", fasta, "seq1:0-1" },
        1,
        "error: start position must be >= 1\n",
    );
    try expectCliFailure(
        allocator,
        &.{ ZFASTA_BIN, "get", fasta, "seq1:25-25" },
        1,
        "error: start position 25 exceeds sequence length 24\n",
    );
    try expectCliFailure(
        allocator,
        &.{ ZFASTA_BIN, "get", fasta, "seq1:10-9" },
        1,
        "error: end position must be >= start position\n",
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

test "[cli] - [get BED]: annotations describe the composed orientation" {
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
