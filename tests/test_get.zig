//! GET integration and CLI tests: region resolution, extraction, and failure paths.
//!
//! Includes indexed extraction and exact-exit CLI failure contracts.

const std = @import("std");
const main = @import("main");
const utility = @import("utility.zig");

const parseRegion = main.getter.parseRegion;
fn resolveRegion(idx: *const main.index_format.LoadedIndex, input: []const u8) !main.getter.ResolvedRegion {
    var failure: main.index_format.diagnostics.Failure = undefined;
    return main.getter.resolveRegion(idx, input, &failure);
}
const io = std.testing.io;

const ZFASTA_BIN = utility.ZFASTA_BIN;
const expectCliFailure = utility.expectCliFailure;
const expectCliSuccess = utility.expectCliSuccess;
const expectUnknownOptionRejected = utility.expectUnknownOptionRejected;
const uniqueArtifactPath = utility.uniqueArtifactPath;
const writeFastaArtifact = utility.writeFastaArtifact;
const writeZfi = utility.writeZfi;
const captureExtractRegion = utility.captureExtractRegion;

const GET_USAGE_STDERR = "error [usage.missing_argument] get: A required argument is missing. argument=request_source\nhint: Correct the command arguments; use --help for accepted options.\n";

fn checkGetOperationAllocations(allocator: std.mem.Allocator, source: main.getter.RequestSource) !void {
    const diagnostics = main.index_format.diagnostics;
    var failure: diagnostics.Failure = undefined;
    _ = main.getter.runGetWithOptions(allocator, io, "tests/data/simple.fasta", .{ .source = source }, &failure) catch |err| {
        if (err == error.OutOfMemory) {
            try std.testing.expectEqual(diagnostics.Code.@"resource.out_of_memory", failure.primary().?.code());
            return err;
        }
        try std.testing.expectEqual(error.InvalidRequest, err);
        try std.testing.expectEqual(diagnostics.Code.@"request.sequence_not_found", failure.primary().?.code());
        return;
    };
    return error.UnexpectedSuccess;
}

test "[unit] - [region parser]: distinguishes literal names and rightmost coordinate suffixes" {
    const Case = struct {
        input: []const u8,
        name: []const u8,
        start: u64 = 1,
        end: ?u64 = null,
        is_full: bool = true,
    };
    const max = std.math.maxInt(u64);
    const cases = [_]Case{
        .{ .input = "seq1", .name = "seq1" },
        .{ .input = "seq1:3-8", .name = "seq1", .start = 3, .end = 8, .is_full = false },
        .{ .input = "seq1:3-", .name = "seq1", .start = 3, .is_full = false },
        .{ .input = ":1-2", .name = "", .end = 2, .is_full = false },
        .{
            .input = "chromosome:GRCh38:1:1:248956422:1",
            .name = "chromosome:GRCh38:1:1:248956422:1",
        },
        .{
            .input = "chromosome:GRCh38:1:1:248956422:1:100-200",
            .name = "chromosome:GRCh38:1:1:248956422:1",
            .start = 100,
            .end = 200,
            .is_full = false,
        },
        .{ .input = "name:1-2:bad", .name = "name:1-2:bad" },
        .{ .input = "name:18446744073709551616-1", .name = "name:18446744073709551616-1" },
        .{ .input = "name:1-18446744073709551616", .name = "name:1-18446744073709551616" },
        .{
            .input = "name:18446744073709551615-18446744073709551615",
            .name = "name",
            .start = max,
            .end = max,
            .is_full = false,
        },
    };

    for (cases) |case| {
        const region = parseRegion(case.input);
        try std.testing.expectEqualStrings(case.name, region.name);
        try std.testing.expectEqual(case.start, region.start);
        try std.testing.expectEqual(case.end, region.end);
        try std.testing.expectEqual(case.is_full, region.is_full);
    }
}

test "[unit] - [get orientation]: composes independent reverse and complement operations" {
    const reverse_complement = main.getter.Orientation.reverseComplement();
    const reverse = main.getter.Orientation.reverseOnly();
    const complement = main.getter.Orientation.complementOnly();

    try std.testing.expect(reverse_complement.compose(reverse_complement).isIdentity());
    try std.testing.expectEqual(complement, reverse_complement.compose(reverse));
    try std.testing.expectEqual(reverse, reverse_complement.compose(complement));
}

test "[integration] - [region resolution]: preserves indexed coordinate geometry" {
    var idx = try main.index_format.loadIndexCheckedLegacy(std.testing.allocator, io, "tests/data/simple.fasta");
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
        const resolved = try resolveRegion(&idx, case.request);

        try std.testing.expectEqualStrings("seq1", resolved.name);
        try std.testing.expectEqual(case.start, resolved.start);
        try std.testing.expectEqual(case.display_end, resolved.display_end);
        try std.testing.expectEqual(case.num_bases, resolved.num_bases);
        try std.testing.expectEqual(case.is_full, resolved.is_full);
        try std.testing.expectEqual(idx.records[0].seq_offset, resolved.seq_offset);
    }
}

test "[integration] - [region resolution]: accepts a pipe-delimited protein identifier" {
    var idx = try main.index_format.loadIndexCheckedLegacy(std.testing.allocator, io, "tests/data/proteome.fasta");
    defer idx.deinit();

    const r = try resolveRegion(&idx, "sp|P12345|PROT_HUMAN:1-10");

    try std.testing.expectEqualStrings("sp|P12345|PROT_HUMAN", r.name);
    try std.testing.expectEqual(@as(u64, 10), r.num_bases);
}

test "[integration] - [region resolution]: accepts a 200-byte identifier" {
    var idx = try main.index_format.loadIndexCheckedLegacy(std.testing.allocator, io, "tests/data/edge_cases.fasta");
    defer idx.deinit();

    const long_name = "A" ** 200;
    const r = try resolveRegion(&idx, long_name);

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
        "error [usage.repeated_option] get: This option cannot be repeated. option=bed\nhint: Correct the command arguments; use --help for accepted options.\n",
    );
    try expectCliFailure(
        allocator,
        &.{ ZFASTA_BIN, "get", fasta, "--names", "first.txt", "--names", "second.txt" },
        1,
        "error [usage.repeated_option] get: This option cannot be repeated. option=names\nhint: Correct the command arguments; use --help for accepted options.\n",
    );
}

test "[cli] - [get]: rejects mixed request source families" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const fasta = "tests/data/simple.fasta";
    try utility.expectCliFailureCode(allocator, &.{ ZFASTA_BIN, "get", fasta, "seq1", "--names", "ids.txt" }, .@"usage.conflicting_options");
    try utility.expectCliFailureCode(allocator, &.{ ZFASTA_BIN, "get", fasta, "seq1", "--bed", "regions.bed" }, .@"usage.conflicting_options");
    try utility.expectCliFailureCode(allocator, &.{ ZFASTA_BIN, "get", fasta, "--names", "ids.txt", "--bed", "regions.bed" }, .@"usage.conflicting_options");
}

test "[cli] - [get]: rejects strand handling outside BED input" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const fasta = "tests/data/simple.fasta";
    const expected = "error [usage.option_requires] get: An option prerequisite is missing. option=strand_aware required=bed\nhint: Correct the command arguments; use --help for accepted options.\n";

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

    const expected_delivery = try std.fmt.allocPrint(allocator, "error [request.identifier_too_long] get: \"zig-cache/test-artifacts\" ({d} bytes):1 column=1 offset=0 request=1: Request identifier exceeds its byte limit. identifier_bytes=65536 limit_bytes=65535\nhint: Correct the reported request or select the intended index.\n", .{names_path.len});
    defer allocator.free(expected_delivery);
    try expectCliFailure(
        allocator,
        &.{ ZFASTA_BIN, "get", "tests/data/simple.fasta", "--names", names_path },
        1,
        expected_delivery,
    );
}

test "[cli] - [get]: missing FASTA path returns an exact error" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const missing = "tests/data/definitely-missing-cli-failure.fasta";
    const expected = "error [io.file_not_found] get: \"tests/data/definitely-mi\" (47 bytes): Required file was not found.; cause=\"FileNotFound\"\nhint: Check the named input/output and retained I/O cause.\n";

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
        "error [io.file_not_found] get: \"tests/data/definitely-mi\" (48 bytes): Required file was not found.; cause=\"FileNotFound\"\nhint: Check the named input/output and retained I/O cause.\n",
    );
    try expectCliFailure(
        allocator,
        &.{ ZFASTA_BIN, "get", "tests/data/simple.fasta", "--names", "tests/data/also-missing.txt" },
        1,
        "error [io.file_not_found] get: \"tests/data/also-missing.\" (27 bytes): Required file was not found.; cause=\"FileNotFound\"\nhint: Check the named input/output and retained I/O cause.\n",
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

    try utility.expectCliFailureCode(allocator, &.{ ZFASTA_BIN, "get", fasta, "seq" }, .@"index.truncated");
    try utility.expectCliFailureCode(allocator, &.{ ZFASTA_BIN, "get", fasta, "--names", "missing.txt" }, .@"index.truncated");
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

    try utility.expectCliFailureCode(allocator, &.{ ZFASTA_BIN, "get", "tests/data/simple.fasta", "--names", early_path, "--summary" }, .@"request.sequence_not_found");

    var threaded = std.Io.Threaded.init(allocator, .{});
    defer threaded.deinit();
    const result = try std.process.run(allocator, threaded.io(), .{
        .argv = &.{ ZFASTA_BIN, "get", "--error-format=json", "tests/data/simple.fasta", "--names", late_path, "--summary" },
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
    try utility.expectDiagnostic(allocator, result.stderr, .@"request.sequence_not_found");
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
        "error [request.sequence_not_found] get: \"positional requests\" request=1: Sequence was not found in the selected lookup. lookup_scope=complete_catalog index_source=zfi identifier=\"missing\"\nhint: Correct the reported request or select the intended index.\n",
    );
    try expectCliFailure(
        allocator,
        &.{ ZFASTA_BIN, "get", fasta, "seq1:0-1" },
        1,
        "error [request.start_before_one] get: \"positional requests\" request=1: Inclusive request start must be at least one. start=0 coordinate_system=one_based_inclusive\nhint: Correct the reported request or select the intended index.\n",
    );
    try expectCliFailure(
        allocator,
        &.{ ZFASTA_BIN, "get", fasta, "seq1:25-25" },
        1,
        "error [request.start_after_sequence] get: \"positional requests\" request=1: Request start exceeds the indexed sequence length. start=25 sequence_length=24\nhint: Correct the reported request or select the intended index.\n",
    );
    try expectCliFailure(
        allocator,
        &.{ ZFASTA_BIN, "get", fasta, "seq1:10-9" },
        1,
        "error [request.end_before_start] get: \"positional requests\" request=1: Request end is before its start. start=10 end=9\nhint: Correct the reported request or select the intended index.\n",
    );
    try expectCliFailure(
        allocator,
        &.{ ZFASTA_BIN, "get", fasta, "seq1", "--annotate-rc" },
        1,
        "error [usage.option_requires] get: An option prerequisite is missing. option=annotate_rc required=transform_or_strand\nhint: Correct the command arguments; use --help for accepted options.\n",
    );
    try expectCliFailure(
        allocator,
        &.{ ZFASTA_BIN, "get", fasta, "--rc", "--complement-only", "seq1" },
        1,
        "error [usage.conflicting_options] get: Selected options conflict. arguments=rc,complement_only\nhint: Correct the command arguments; use --help for accepted options.\n",
    );
    const empty_names = try writeFastaArtifact(allocator, "get-empty-names", "");
    defer std.Io.Dir.cwd().deleteFile(io, empty_names) catch {};
    try utility.expectCliFailureCode(allocator, &.{ ZFASTA_BIN, "get", fasta, "--names", empty_names }, .@"request.no_regions");
}

test "[cli] - [get BED]: converts coordinates and skips metadata lines" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const bed_path = try uniqueArtifactPath(allocator, "bed-parser-success", "bed");
    defer std.Io.Dir.cwd().deleteFile(io, bed_path) catch {};
    {
        const bed_file = try std.Io.Dir.cwd().createFile(io, bed_path, .{});
        defer bed_file.close(io);
        try std.Io.File.writeStreamingAll(
            bed_file,
            io,
            "track name=ignored\n# ignored\nbrowser position chr1:1-10\nseq1\t0\t1\nseq1\t5\t9\r\n",
        );
    }

    try expectCliSuccess(
        allocator,
        &.{ ZFASTA_BIN, "get", "tests/data/simple.fasta", "--bed", bed_path },
        ">seq1:1-1\nA\n>seq1:6-9\nCGTA\n",
    );
}

test "[cli] - [get BED]: reports parser boundary errors" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const bed_path = try uniqueArtifactPath(allocator, "bed-parser-failures", "bed");
    defer std.Io.Dir.cwd().deleteFile(io, bed_path) catch {};
    const Case = struct {
        line: []const u8,
        code: main.index_format.diagnostics.Code,
    };
    const cases = [_]Case{
        .{ .line = "\t0\t1\n", .code = .@"request.bed.field_missing" },
        .{ .line = "seq1\n", .code = .@"request.bed.field_missing" },
        .{ .line = "seq1\t0\n", .code = .@"request.bed.field_missing" },
        .{ .line = "seq1\t-0\t10\n", .code = .@"request.bed.invalid_coordinate" },
        .{ .line = "seq1\t+1\t10\n", .code = .@"request.bed.invalid_coordinate" },
        .{ .line = "seq1\t1_0\t20\n", .code = .@"request.bed.invalid_coordinate" },
        .{ .line = "seq1\t18446744073709551616\t20\n", .code = .@"request.bed.coordinate_out_of_range" },
        .{ .line = "seq1\t0\t+10\n", .code = .@"request.bed.invalid_coordinate" },
        .{ .line = "seq1\t0\t1_0\n", .code = .@"request.bed.invalid_coordinate" },
        .{ .line = "seq1\t0\t18446744073709551616\n", .code = .@"request.bed.coordinate_out_of_range" },
        .{ .line = "seq1\t10\t10\n", .code = .@"request.bed.invalid_interval" },
        .{ .line = "seq1\t10\t9\n", .code = .@"request.bed.invalid_interval" },
        .{
            .line = "seq1\t18446744073709551614\t18446744073709551615\n",
            .code = .@"request.start_after_sequence",
        },
    };

    for (cases) |case| {
        {
            const bed_file = try std.Io.Dir.cwd().createFile(io, bed_path, .{ .truncate = true });
            defer bed_file.close(io);
            try std.Io.File.writeStreamingAll(bed_file, io, case.line);
        }
        try utility.expectCliFailureCode(
            allocator,
            &.{ ZFASTA_BIN, "get", "tests/data/simple.fasta", "--bed", bed_path },
            case.code,
        );
    }

    {
        const bed_file = try std.Io.Dir.cwd().createFile(io, bed_path, .{ .truncate = true });
        defer bed_file.close(io);
        try std.Io.File.writeStreamingAll(bed_file, io, "seq1\t0\t1\tname\t0\t?\n");
    }
    try expectCliSuccess(
        allocator,
        &.{ ZFASTA_BIN, "get", "tests/data/simple.fasta", "--bed", bed_path },
        ">seq1:1-1\nA\n",
    );
    const expected_delivery = try std.fmt.allocPrint(allocator, "error [request.bed.invalid_strand] get: \"zig-cache/test-artifacts\" ({d} bytes):1 column=17 offset=16 request=1: Selected BED strand is invalid. strand=\"?\" accepted=plus_minus_dot_or_absent\nhint: Correct the reported request or select the intended index.\n", .{bed_path.len});
    defer allocator.free(expected_delivery);
    try expectCliFailure(
        allocator,
        &.{ ZFASTA_BIN, "get", "tests/data/simple.fasta", "--bed", bed_path, "--strand-aware" },
        1,
        expected_delivery,
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

test "[failure] - [BED parser]: retains field locations and distinguishes decimal overflow" {
    const diagnostics = main.index_format.diagnostics;
    const cases = [_]struct {
        row: []const u8,
        code: diagnostics.Code,
        column: u64,
        legacy: main.bed_parser.ParseError,
        value: ?[]const u8 = null,
    }{
        .{ .row = "\t0\t1", .code = .@"request.bed.field_missing", .column = 1, .legacy = error.MissingChrom },
        .{ .row = "seq", .code = .@"request.bed.field_missing", .column = 4, .legacy = error.MissingStart },
        .{ .row = "seq\t\t1", .code = .@"request.bed.field_missing", .column = 5, .legacy = error.MissingStart },
        .{ .row = "seq\t0", .code = .@"request.bed.field_missing", .column = 6, .legacy = error.MissingEnd },
        .{ .row = "seq\t0\t", .code = .@"request.bed.field_missing", .column = 7, .legacy = error.MissingEnd },
        .{ .row = "seq\tx\t1", .code = .@"request.bed.invalid_coordinate", .column = 5, .legacy = error.InvalidStart, .value = "x" },
        .{ .row = "seq\t0\t-1", .code = .@"request.bed.invalid_coordinate", .column = 7, .legacy = error.InvalidEnd, .value = "-1" },
        .{ .row = "seq\t18446744073709551616\t1", .code = .@"request.bed.coordinate_out_of_range", .column = 5, .legacy = error.InvalidStart, .value = "18446744073709551616" },
        .{ .row = "seq\t0\t18446744073709551616", .code = .@"request.bed.coordinate_out_of_range", .column = 7, .legacy = error.InvalidEnd, .value = "18446744073709551616" },
        .{ .row = "seq\t18446744073709551616000x\t1", .code = .@"request.bed.invalid_coordinate", .column = 5, .legacy = error.InvalidStart, .value = "18446744073709551616000x" },
        .{ .row = "seq\t2\t2", .code = .@"request.bed.invalid_interval", .column = 5, .legacy = error.EmptyInterval },
        .{ .row = "seq\t2\t1", .code = .@"request.bed.invalid_interval", .column = 5, .legacy = error.EmptyInterval },
        .{ .row = "seq\tx", .code = .@"request.bed.field_missing", .column = 6, .legacy = error.MissingEnd },
    };
    for (cases) |case| {
        var failure: diagnostics.Failure = undefined;
        failure.reset(.get);
        try std.testing.expectError(error.InvalidRequest, main.bed_parser.parseBedLineChecked(case.row, false, .{
            .source = .bed,
            .location = .{ .offset = 100, .line = 12 },
            .request = 8,
        }, &failure));
        const primary = failure.primary().?;
        try std.testing.expectEqual(case.code, primary.code());
        try std.testing.expectEqual(@as(?u64, case.column), primary.at.location.column);
        try std.testing.expectEqual(@as(?u64, 100 + case.column - 1), primary.at.location.offset);
        try std.testing.expectEqual(@as(?u64, 12), primary.at.location.line);
        try std.testing.expectEqual(@as(?u64, 8), primary.at.request);
        try std.testing.expectEqual(@as(?diagnostics.SourceId, .bed), primary.at.source);
        try std.testing.expectEqual(diagnostics.Stage.parse, primary.stage);
        try std.testing.expect(primary.cause == null);
        switch (primary.details) {
            .@"request.bed.field_missing" => |detail| {
                const expected: diagnostics.BedField = switch (case.legacy) {
                    error.MissingChrom => .chrom,
                    error.MissingStart => .start,
                    else => .end,
                };
                try std.testing.expectEqual(expected, detail.field);
                try std.testing.expectEqual(@as(u8, switch (expected) {
                    .chrom => 1,
                    .start => 2,
                    .end => 3,
                }), detail.column_number);
            },
            .@"request.bed.invalid_coordinate" => |detail| {
                try std.testing.expectEqual(diagnostics.CoordinateSystem.zero_based_half_open, detail.coordinate_system);
                try std.testing.expectEqual(if (case.legacy == error.InvalidStart) diagnostics.BedCoordinate.start else .end, detail.field);
                try std.testing.expectEqualStrings(case.value.?, failure.preview(detail.value.?));
            },
            .@"request.bed.coordinate_out_of_range" => |detail| {
                try std.testing.expectEqual(std.math.maxInt(u64), detail.max_value);
                try std.testing.expectEqual(diagnostics.CoordinateSystem.zero_based_half_open, detail.coordinate_system);
                try std.testing.expectEqual(if (case.legacy == error.InvalidStart) diagnostics.BedCoordinate.start else .end, detail.field);
                try std.testing.expectEqualStrings(case.value.?, failure.preview(detail.value.?));
            },
            .@"request.bed.invalid_interval" => |detail| {
                try std.testing.expectEqual(@as(u64, 2), detail.start);
                try std.testing.expectEqual(diagnostics.CoordinateSystem.zero_based_half_open, detail.coordinate_system);
                try std.testing.expectEqual(@as(u64, if (std.mem.endsWith(u8, case.row, "2")) 2 else 1), detail.end);
            },
            else => return error.UnexpectedDiagnostic,
        }
        try std.testing.expectError(case.legacy, main.bed_parser.parseBedLine(case.row));
    }
}

test "[edge] - [BED parser]: preserves skips, leading zeros, maximum coordinates and optional strand" {
    const diagnostics = main.index_format.diagnostics;
    for ([_][]const u8{ "", "\r", "# comment", "track", "track name=x", "browser", "browser position chr1" }) |row| {
        var failure: diagnostics.Failure = undefined;
        failure.reset(.get);
        const parsed = try main.bed_parser.parseBedLineChecked(row, true, .{}, &failure);
        try std.testing.expect(parsed == .skip);
        try std.testing.expect(failure.primary() == null);
    }
    const rows = [_][]const u8{
        "seq\t0000\t0001\tname\t0\t-\r",
        "seq\t18446744073709551614\t18446744073709551615",
        "seq\t0\t1\tname\t0\t?",
        "seq\t0\t1\tname\t0\t",
    };
    for (rows, 0..) |row, i| {
        var failure: diagnostics.Failure = undefined;
        failure.reset(.get);
        const parsed = try main.bed_parser.parseBedLineChecked(row, false, .{}, &failure);
        const legacy = try main.bed_parser.parseBedLine(row);
        try std.testing.expectEqualDeep(legacy, parsed);
        try std.testing.expect(failure.primary() == null);
        try std.testing.expectEqual(@as(u64, if (i == 1) std.math.maxInt(u64) else 1), parsed.region.end_1based);
    }
    var failure: diagnostics.Failure = undefined;
    failure.reset(.get);
    try std.testing.expectError(error.InvalidRequest, main.bed_parser.parseBedLineChecked(rows[2], true, .{
        .source = .bed,
        .location = .{ .line = 9, .offset = 40 },
    }, &failure));
    const primary = failure.primary().?;
    try std.testing.expectEqual(diagnostics.Code.@"request.bed.invalid_strand", primary.code());
    try std.testing.expectEqualStrings("?", failure.preview(primary.details.@"request.bed.invalid_strand".strand.?));
    try std.testing.expectEqual(@as(?u64, 16), primary.at.location.column);
    try std.testing.expectEqual(@as(?u64, 55), primary.at.location.offset);
    try std.testing.expectEqual(@as(?u64, 9), primary.at.location.line);
}

test "[failure] - [BED parser]: owns bounded previews after the input is reused" {
    const diagnostics = main.index_format.diagnostics;
    var row = [_]u8{'9'} ** 128;
    @memcpy(row[0..6], "seq\t0\t");
    var failure: diagnostics.Failure = undefined;
    failure.reset(.get);
    try std.testing.expectError(error.InvalidRequest, main.bed_parser.parseBedLineChecked(&row, false, .{ .source = .bed }, &failure));
    @memset(&row, 0);
    const primary = failure.primary().?;
    const preview = primary.details.@"request.bed.coordinate_out_of_range".value.?;
    try std.testing.expectEqual(@as(u64, 122), preview.original_bytes);
    try std.testing.expect(preview.truncated);
    try std.testing.expectEqualStrings("99999999999999999999999999999999", failure.preview(preview));
    try std.testing.expect(primary.at.location.offset == null);
    try std.testing.expectEqual(@as(?u64, 7), primary.at.location.column);
}

test "[failure] - [GET diagnostics]: returns exact request facts without exiting" {
    const diagnostics = main.index_format.diagnostics;
    var idx = try main.index_format.loadIndexCheckedLegacy(std.testing.allocator, io, "tests/data/simple.fasta");
    defer idx.deinit();
    const cases = [_]struct { request: []const u8, code: diagnostics.Code }{
        .{ .request = "missing", .code = .@"request.sequence_not_found" },
        .{ .request = "seq1:0-1", .code = .@"request.start_before_one" },
        .{ .request = "seq1:25-25", .code = .@"request.start_after_sequence" },
        .{ .request = "seq1:10-9", .code = .@"request.end_before_start" },
    };
    var failure: diagnostics.Failure = undefined;
    for (cases) |case| {
        try std.testing.expectError(error.InvalidRequest, main.getter.resolveRegion(&idx, case.request, &failure));
        const primary = failure.primary().?;
        try std.testing.expectEqual(case.code, primary.code());
        try std.testing.expectEqual(@as(?u64, 1), primary.at.request);
        try std.testing.expect(primary.cause == null);
        try std.testing.expectEqual(diagnostics.Stage.resolve, primary.stage);
        try std.testing.expectEqual(diagnostics.SourceId.request_snapshot, primary.at.source.?);
        try std.testing.expect(!failure.work_complete and !failure.complete);
        switch (primary.details) {
            .@"request.sequence_not_found" => |detail| {
                try std.testing.expectEqualStrings("missing", failure.preview(detail.identifier.?));
                try std.testing.expectEqual(diagnostics.LookupScope.retained_matches, detail.lookup_scope);
                try std.testing.expectEqual(diagnostics.SourceId.zfi, detail.index_source);
            },
            .@"request.start_before_one" => |detail| {
                try std.testing.expectEqual(@as(u64, 0), detail.start);
                try std.testing.expectEqual(diagnostics.CoordinateSystem.one_based_inclusive, detail.coordinate_system);
            },
            .@"request.start_after_sequence" => |detail| {
                try std.testing.expectEqual(@as(u64, 25), detail.start);
                try std.testing.expectEqual(@as(u64, 24), detail.sequence_length);
            },
            .@"request.end_before_start" => |detail| {
                try std.testing.expectEqual(@as(u64, 10), detail.start);
                try std.testing.expectEqual(@as(u64, 9), detail.end);
            },
            else => return error.UnexpectedDiagnostic,
        }
    }
    const good = try main.getter.resolveRegion(&idx, "seq1:1-999", &failure);
    try std.testing.expectEqual(@as(u64, 24), good.num_bases);
    try std.testing.expect(failure.primary() == null);
    try std.testing.expect(failure.complete);
}

test "[failure] - [GET diagnostics]: retains FASTA read causes and exact short reads" {
    const diagnostics = main.index_format.diagnostics;
    const Fault = struct {
        fn read(_: ?*anyopaque, _: std.Io.File, _: []const []u8, _: u64) std.Io.File.ReadPositionalError!usize {
            return error.InputOutput;
        }
        fn short(_: ?*anyopaque, _: std.Io.File, _: []const []u8, _: u64) std.Io.File.ReadPositionalError!usize {
            return 0;
        }
    };
    var idx = try main.index_format.loadIndexCheckedLegacy(std.testing.allocator, io, "tests/data/simple.fasta");
    defer idx.deinit();
    for (0..2) |i| {
        var vtable = io.vtable.*;
        vtable.fileReadPositional = if (i == 0) Fault.read else Fault.short;
        idx.io = .{ .userdata = io.userdata, .vtable = &vtable };
        defer idx.io = io;
        var output = std.Io.Writer.Allocating.init(std.testing.allocator);
        defer output.deinit();
        var failure: diagnostics.Failure = undefined;
        try std.testing.expectError(error.Io, main.getter.extractRegion(&idx, "seq1:1-2", &output.writer, &failure));
        const primary = failure.primary().?;
        try std.testing.expectEqual(@as(?u64, idx.records[0].seq_offset), primary.at.location.offset);
        try std.testing.expectEqual(@as(?diagnostics.SourceId, .fasta), primary.at.source);
        try std.testing.expectEqual(@as(?u64, 1), primary.related.?.request);
        if (i == 0) {
            try std.testing.expectEqual(diagnostics.Code.@"io.read_failed", primary.code());
            try std.testing.expectEqual(error.InputOutput, primary.cause.?.zig);
            try std.testing.expectEqual(@as(?u64, 2), primary.details.@"io.read_failed".requested_bytes);
        } else {
            const detail = primary.details.@"io.unexpected_eof".exact;
            try std.testing.expectEqual(@as(u64, 2), detail.requested_bytes);
            try std.testing.expectEqual(@as(u64, 0), detail.received_bytes);
            try std.testing.expect(primary.cause == null);
        }
        try std.testing.expectEqualStrings(">seq1:1-2\n", output.written());
        try std.testing.expect(failure.outputs.file.may_be_partial);
    }
}

test "[failure] - [GET diagnostics]: returns partial writer failure with no final flush" {
    const diagnostics = main.index_format.diagnostics;
    var idx = try main.index_format.loadIndexCheckedLegacy(std.testing.allocator, io, "tests/data/simple.fasta");
    defer idx.deinit();
    var buffer: [1]u8 = undefined;
    var writer = std.Io.Writer.fixed(&buffer);
    var failure: diagnostics.Failure = undefined;
    try std.testing.expectError(error.Io, main.getter.extractRegion(&idx, "seq1", &writer, &failure));
    try std.testing.expectEqual(diagnostics.Code.@"io.write_failed", failure.primary().?.code());
    try std.testing.expectEqual(error.WriteFailed, failure.primary().?.cause.?.zig);
    try std.testing.expectEqual(@as(?diagnostics.SourceId, .output), failure.primary().?.at.source);
    try std.testing.expectEqualStrings(">", writer.buffered());
    try std.testing.expectEqual(.not_attempted, failure.outputs.file.final_flush);
    try std.testing.expect(failure.outputs.file.may_be_partial);
    try std.testing.expect(!failure.complete);
}

test "[failure] - [GET diagnostics]: captures the actual protein sample and source priority" {
    const diagnostics = main.index_format.diagnostics;
    var failure: diagnostics.Failure = undefined;
    try std.testing.expectError(error.UnsupportedTransform, main.getter.runGetWithOptions(std.testing.allocator, io, "tests/data/proteome.fasta", .{
        .source = .{ .positional = &.{"sp|P12345|PROT_HUMAN:1-10"} },
        .orientation = main.getter.Orientation.complementOnly(),
    }, &failure));
    const detail = failure.primary().?.details.@"transform.complement_on_protein";
    try std.testing.expectEqual(@as(u32, 10), detail.sample.symbols);
    try std.testing.expectEqual(@as(?u64, 1), detail.sample.request);
    try std.testing.expectEqual(diagnostics.Interval{ .start = 1, .end = 10 }, detail.sample.interval.?);
    try std.testing.expect(detail.transform.complement and !detail.transform.reverse);
    try std.testing.expectEqualStrings("sp|P12345|PROT_HUMAN", failure.preview(detail.identifier.?));
    try std.testing.expectEqual(.not_started, failure.outputs.stdout.state);
    try std.testing.expectError(error.Io, main.getter.runGetWithOptions(std.testing.allocator, io, "tests/data/absent-native-get.fasta", .{
        .source = .{ .names = "absent-native-requests.txt" },
    }, &failure));
    try std.testing.expectEqual(diagnostics.Code.@"io.file_not_found", failure.primary().?.code());
    try std.testing.expectEqual(@as(?diagnostics.SourceId, .fasta), failure.primary().?.at.source);
}

test "[cli] - [GET abort]: discards pending stdout after a later malformed BED row" {
    const allocator = std.testing.allocator;
    const path = try uniqueArtifactPath(allocator, "get-abort-bed", "bed");
    defer allocator.free(path);
    defer std.Io.Dir.cwd().deleteFile(io, path) catch {};
    var rows = std.Io.Writer.Allocating.init(allocator);
    defer rows.deinit();
    for (0..8192) |_| try rows.writer.writeAll("seq1\t0\t1\n");
    try rows.writer.writeAll("seq1\n");
    const file = try std.Io.Dir.cwd().createFile(io, path, .{});
    try std.Io.File.writeStreamingAll(file, io, rows.written());
    file.close(io);
    var expected = std.Io.Writer.Allocating.init(allocator);
    defer expected.deinit();
    for (0..8192) |_| try expected.writer.writeAll(">seq1:1-1\nA\n");
    const expected_delivery = try std.fmt.allocPrint(allocator, "error [request.bed.field_missing] get: \"zig-cache/test-artifacts\" ({d} bytes):8193 column=5 offset=73732 request=8193: A required BED field is missing or empty. field=start column_number=2\nhint: Correct the reported request or select the intended index.\n", .{path.len});
    defer allocator.free(expected_delivery);
    try utility.expectCliResult(allocator, &.{ ZFASTA_BIN, "get", "tests/data/simple.fasta", "--bed", path, "--summary" }, 1, expected.written()[0..65537], expected_delivery);
}

test "[failure] - [GET operation]: releases allocation failures in positional and cached file modes" {
    const allocator = std.testing.allocator;
    const names = try uniqueArtifactPath(allocator, "get-oom-names", "txt");
    defer allocator.free(names);
    defer std.Io.Dir.cwd().deleteFile(io, names) catch {};
    const bed = try uniqueArtifactPath(allocator, "get-oom-bed", "bed");
    defer allocator.free(bed);
    defer std.Io.Dir.cwd().deleteFile(io, bed) catch {};
    for ([_]struct { path: []const u8, bytes: []const u8 }{
        .{ .path = names, .bytes = "missing\n" }, .{ .path = bed, .bytes = "missing\t0\t1\n" },
    }) |fixture| {
        const file = try std.Io.Dir.cwd().createFile(io, fixture.path, .{});
        defer file.close(io);
        try std.Io.File.writeStreamingAll(file, io, fixture.bytes);
    }
    for ([_]main.getter.RequestSource{ .{ .positional = &.{"missing"} }, .{ .names = names }, .{ .bed = bed } }) |source|
        try std.testing.checkAllAllocationFailures(allocator, checkGetOperationAllocations, .{source});
}

test "[failure] - [GET stdout]: distinguishes writes from final flush and retains the backend cause" {
    const diagnostics = main.index_format.diagnostics;
    const Fault = struct {
        fn operate(userdata: ?*anyopaque, operation: std.Io.Operation) std.Io.Cancelable!std.Io.Operation.Result {
            return switch (operation) {
                .file_write_streaming => .{ .file_write_streaming = error.BrokenPipe },
                else => std.testing.io.vtable.operate(userdata, operation),
            };
        }
    };
    var vtable = io.vtable.*;
    vtable.operate = Fault.operate;
    const fault_io = std.Io{ .userdata = io.userdata, .vtable = &vtable };
    var failure: diagnostics.Failure = undefined;
    try std.testing.expectError(error.Io, main.getter.runGetWithOptions(std.testing.allocator, fault_io, "tests/data/simple.fasta", .{
        .source = .{ .positional = &.{"seq1"} },
    }, &failure));
    try std.testing.expectEqual(diagnostics.Code.@"io.flush_failed", failure.primary().?.code());
    try std.testing.expectEqual(error.BrokenPipe, failure.primary().?.cause.?.zig);
    try std.testing.expect(failure.work_complete and !failure.complete);
    try std.testing.expectEqual(.failed, failure.outputs.stdout.final_flush);
    try std.testing.expectEqual(@as(?u64, 1), failure.outputs.stdout.records);

    const path = try uniqueArtifactPath(std.testing.allocator, "get-write-fault", "txt");
    defer std.testing.allocator.free(path);
    defer std.Io.Dir.cwd().deleteFile(io, path) catch {};
    var rows = std.Io.Writer.Allocating.init(std.testing.allocator);
    defer rows.deinit();
    for (0..3000) |_| try rows.writer.writeAll("seq1\n");
    const file = try std.Io.Dir.cwd().createFile(io, path, .{});
    try std.Io.File.writeStreamingAll(file, io, rows.written());
    file.close(io);
    try std.testing.expectError(error.Io, main.getter.runGetWithOptions(std.testing.allocator, fault_io, "tests/data/simple.fasta", .{
        .source = .{ .names = path },
    }, &failure));
    try std.testing.expectEqual(diagnostics.Code.@"io.write_failed", failure.primary().?.code());
    try std.testing.expectEqual(error.BrokenPipe, failure.primary().?.cause.?.zig);
    try std.testing.expect(!failure.work_complete and !failure.complete);
    try std.testing.expectEqual(.not_attempted, failure.outputs.stdout.final_flush);
    try std.testing.expect(failure.outputs.stdout.records.? > 0);
}

test "[cli] - [GET missing identifiers]: reports original lengths and bounded prefixes" {
    const allocator = std.testing.allocator;
    const name = "absent" ** 10;
    const path = try uniqueArtifactPath(allocator, "get-long-miss", "txt");
    defer allocator.free(path);
    defer std.Io.Dir.cwd().deleteFile(io, path) catch {};
    const file = try std.Io.Dir.cwd().createFile(io, path, .{});
    try std.Io.File.writeStreamingAll(file, io, name ++ "\n");
    file.close(io);
    var threaded = std.Io.Threaded.init(allocator, .{});
    defer threaded.deinit();
    for ([_]struct { argv: []const []const u8, source: []const u8 }{
        .{ .argv = &.{ ZFASTA_BIN, "get", "tests/data/simple.fasta", name, "--error-format=json" }, .source = "request_snapshot" },
        .{ .argv = &.{ ZFASTA_BIN, "get", "tests/data/simple.fasta", "--names", path, "--error-format=json" }, .source = "names" },
    }) |case| {
        const result = try std.process.run(allocator, threaded.io(), .{ .argv = case.argv });
        defer allocator.free(result.stdout);
        defer allocator.free(result.stderr);
        try std.testing.expectEqual(@as(u8, 1), result.term.exited);
        try std.testing.expectEqualStrings("", result.stdout);
        try utility.expectDiagnostic(allocator, result.stderr, .@"request.sequence_not_found");
        const parsed = try std.json.parseFromSlice(std.json.Value, allocator, result.stderr, .{});
        defer parsed.deinit();
        const primary = parsed.value.object.get("failure").?.object.get("primary").?.object;
        try std.testing.expectEqualStrings(case.source, primary.get("source").?.string);
        try std.testing.expectEqualStrings("1", primary.get("request").?.string);
        const identifier = primary.get("data").?.object.get("identifier").?.object;
        try std.testing.expectEqualStrings("utf8", identifier.get("encoding").?.string);
        try std.testing.expectEqualStrings(name[0..32], identifier.get("value").?.string);
        try std.testing.expectEqualStrings("60", identifier.get("original_bytes").?.string);
        try std.testing.expect(identifier.get("truncated").?.bool);
    }
}

test "[cli] - [get machine summary]: frames metrics without changing FASTA output" {
    const allocator = std.testing.allocator;
    var threaded = std.Io.Threaded.init(allocator, .{});
    defer threaded.deinit();
    const result = try std.process.run(allocator, threaded.io(), .{ .argv = &.{ ZFASTA_BIN, "get", "--error-format=json", "tests/data/simple.fasta", "seq1:1-4", "--summary" } });
    defer allocator.free(result.stdout);
    defer allocator.free(result.stderr);
    try std.testing.expectEqual(@as(u8, 0), result.term.exited);
    try std.testing.expectEqualStrings(">seq1:1-4\nACGT\n", result.stdout);
    const parsed = try std.json.parseFromSlice(std.json.Value, allocator, result.stderr, .{});
    defer parsed.deinit();
    const object = parsed.value.object;
    try std.testing.expectEqualStrings("summary", object.get("kind").?.string);
    try std.testing.expect(object.get("accepted").?.bool);
    try std.testing.expect(object.get("operation").?.object.get("complete").?.bool);
    try std.testing.expect(object.get("failure").? == .null);
    const metrics = object.get("metrics").?.object;
    try std.testing.expectEqualStrings("1", metrics.get("region_count").?.string);
    try std.testing.expectEqualStrings("4", metrics.get("total_bases").?.string);
    _ = try std.fmt.parseInt(u64, metrics.get("elapsed_ns").?.string, 10);
}

test "[cli] - [get machine failure]: preserves opaque identifier bytes" {
    const allocator = std.testing.allocator;
    const names = try writeFastaArtifact(allocator, "get-machine-opaque", "\xff\n");
    defer allocator.free(names);
    defer std.Io.Dir.cwd().deleteFile(io, names) catch {};
    var threaded = std.Io.Threaded.init(allocator, .{});
    defer threaded.deinit();
    const result = try std.process.run(allocator, threaded.io(), .{ .argv = &.{ ZFASTA_BIN, "get", "tests/data/simple.fasta", "--names", names, "--error-format=json" } });
    defer allocator.free(result.stdout);
    defer allocator.free(result.stderr);
    try std.testing.expectEqual(@as(u8, 1), result.term.exited);
    try std.testing.expectEqualStrings("", result.stdout);
    try utility.expectDiagnostic(allocator, result.stderr, .@"request.sequence_not_found");
    const parsed = try std.json.parseFromSlice(std.json.Value, allocator, result.stderr, .{});
    defer parsed.deinit();
    const primary = parsed.value.object.get("failure").?.object.get("primary").?.object;
    try std.testing.expectEqualStrings("names", primary.get("source").?.string);
    try std.testing.expectEqualStrings("1", primary.get("location").?.object.get("line").?.string);
    try std.testing.expectEqualStrings("1", primary.get("request").?.string);
    const identifier = primary.get("data").?.object.get("identifier").?.object;
    try std.testing.expectEqualStrings("base64", identifier.get("encoding").?.string);
    try std.testing.expectEqualStrings("/w==", identifier.get("value").?.string);
    try std.testing.expect(!identifier.get("truncated").?.bool);
}

test "[failure] - [GET diagnostics]: classifies an overlong FASTA path" {
    const diagnostics = main.index_format.diagnostics;
    const path: [256]u8 = @splat('x');
    var failure: diagnostics.Failure = undefined;
    try std.testing.expectError(error.Io, main.getter.runGetWithOptions(std.testing.allocator, io, &path, .{ .source = .{ .positional = &.{"seq1"} } }, &failure));
    const primary = failure.primary().?;
    try std.testing.expectEqual(diagnostics.Code.@"io.invalid_path", primary.code());
    try std.testing.expectEqual(diagnostics.Stage.open, primary.stage);
    try std.testing.expectEqual(diagnostics.SourceId.fasta, primary.at.source);
    try std.testing.expectEqual(error.NameTooLong, primary.cause.?.zig);
    const details = primary.details.@"io.invalid_path";
    try std.testing.expectEqual(diagnostics.PathReason.too_long, details.reason);
    try std.testing.expect(details.limit_bytes == null);
    try std.testing.expectEqual(diagnostics.Operation.get, failure.operation);
    try std.testing.expectEqualStrings(&path, failure.sourceLabel(.fasta).?);
    try std.testing.expectEqual(diagnostics.SourceKind.path, failure.sources.get(.fasta).?.kind);
    try std.testing.expect(!failure.work_complete and !failure.complete);
    try std.testing.expectEqual(.not_started, failure.outputs.stdout.state);
    try std.testing.expectEqual(.not_attempted, failure.outputs.stdout.final_flush);
    try std.testing.expectEqual(@as(usize, 0), failure.secondaryItems().len);
}
