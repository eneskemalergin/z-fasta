//! Shared test mechanics for command suites.

const std = @import("std");
const builtin = @import("builtin");
const main = @import("main");
const diagnostics = main.index_format.diagnostics;

const test_options = @import("test_options");
pub const ZFASTA_BIN = if (std.fs.path.isAbsolute(test_options.binary_path)) test_options.binary_path else test_options.build_root ++ std.fs.path.sep_str ++ test_options.binary_path;

/// Loads with automatic selection, for tests that need a loaded index; a failure fails the test.
pub fn loadIndex(allocator: std.mem.Allocator, io: std.Io, fasta_path: []const u8, mode: main.index_format.LoadMode) diagnostics.Error!main.index_format.LoadedIndex {
    var failure: diagnostics.Failure = undefined;
    failure.reset(if (mode == .stats_scan) .stats else .get);
    return main.index_format.loadIndex(allocator, io, fasta_path, .automatic, mode, &failure);
}

pub fn expectCliResult(
    allocator: std.mem.Allocator,
    argv: []const []const u8,
    exit_code: u8,
    expected_stdout: []const u8,
    expected_stderr: []const u8,
) !void {
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
        .exited => |code| try std.testing.expectEqual(exit_code, code),
        else => return error.ChildProcessFailed,
    }
    try std.testing.expectEqualStrings(expected_stdout, result.stdout);
    try std.testing.expectEqualStrings(expected_stderr, result.stderr);
}

/// Runs `argv` with `stdin` as its standard input.
pub fn expectCliStdinResult(
    allocator: std.mem.Allocator,
    argv: []const []const u8,
    stdin: []const u8,
    exit_code: u8,
    expected_stdout: []const u8,
    expected_stderr: []const u8,
) !void {
    var threaded = std.Io.Threaded.init(allocator, .{});
    defer threaded.deinit();
    const io = threaded.io();
    var child = try std.process.spawn(io, .{ .argv = argv, .stdin = .pipe, .stdout = .pipe, .stderr = .pipe });
    try child.stdin.?.writeStreamingAll(io, stdin);
    child.stdin.?.close(io);
    child.stdin = null;
    var stdout_buffer: [4096]u8 = undefined;
    var stdout_reader = child.stdout.?.readerStreaming(io, &stdout_buffer);
    const stdout = try stdout_reader.interface.allocRemaining(allocator, .limited(64 * 1024));
    defer allocator.free(stdout);
    var stderr_buffer: [4096]u8 = undefined;
    var stderr_reader = child.stderr.?.readerStreaming(io, &stderr_buffer);
    const stderr = try stderr_reader.interface.allocRemaining(allocator, .limited(64 * 1024));
    defer allocator.free(stderr);
    try std.testing.expectEqual(std.process.Child.Term{ .exited = exit_code }, try child.wait(io));
    try std.testing.expectEqualStrings(expected_stdout, stdout);
    try std.testing.expectEqualStrings(expected_stderr, stderr);
}

pub fn expectCliFailure(
    allocator: std.mem.Allocator,
    argv: []const []const u8,
    exit_code: u8,
    expected_stderr: []const u8,
) !void {
    try expectCliResult(allocator, argv, exit_code, "", expected_stderr);
}

pub fn expectDiagnostic(allocator: std.mem.Allocator, bytes: []const u8, code: main.index_format.diagnostics.Code) !void {
    const parsed = try std.json.parseFromSlice(std.json.Value, allocator, bytes, .{});
    defer parsed.deinit();
    const object = parsed.value.object;
    try std.testing.expectEqualStrings("0.4.0-1", object.get("catalog_revision").?.string);
    try std.testing.expectEqualStrings("failure", object.get("kind").?.string);
    try std.testing.expect(!object.get("accepted").?.bool);
    try std.testing.expect(!object.get("operation").?.object.get("complete").?.bool);
    const primary = object.get("failure").?.object.get("primary").?.object;
    try std.testing.expectEqualStrings(@tagName(code), primary.get("code").?.string);
    try std.testing.expectEqualStrings("error", primary.get("severity").?.string);
    try std.testing.expectEqualStrings(@tagName(main.index_format.diagnostics.metadata(code).category), primary.get("category").?.string);
}

pub fn expectCliCodeResult(allocator: std.mem.Allocator, argv: []const []const u8, exit_code: u8, stdout: []const u8, code: main.index_format.diagnostics.Code) !void {
    const args = try allocator.alloc([]const u8, argv.len + 1);
    defer allocator.free(args);
    @memcpy(args[0..argv.len], argv);
    args[argv.len] = "--error-format=json";
    var threaded = std.Io.Threaded.init(allocator, .{});
    defer threaded.deinit();
    const result = try std.process.run(allocator, threaded.io(), .{ .argv = args, .stdout_limit = .limited(128 * 1024), .stderr_limit = .limited(64 * 1024) });
    defer allocator.free(result.stdout);
    defer allocator.free(result.stderr);
    switch (result.term) {
        .exited => |status| try std.testing.expectEqual(exit_code, status),
        else => return error.ChildProcessFailed,
    }
    try std.testing.expectEqualStrings(stdout, result.stdout);
    try expectDiagnostic(allocator, result.stderr, code);
}

pub fn expectCliFailureCode(allocator: std.mem.Allocator, argv: []const []const u8, code: main.index_format.diagnostics.Code) !void {
    try expectCliCodeResult(allocator, argv, 1, "", code);
}

pub fn expectCliSuccess(
    allocator: std.mem.Allocator,
    argv: []const []const u8,
    expected_stdout: []const u8,
) !void {
    try expectCliResult(allocator, argv, 0, expected_stdout, "");
}

pub fn expectUnknownOptionRejected(
    allocator: std.mem.Allocator,
    argv: []const []const u8,
    unknown: []const u8,
) !void {
    const expected = try std.fmt.allocPrint(allocator, "error [usage.unknown_option] {s}: Unknown option. option=\"{s}\"\nhint: Correct the command arguments; use --help for accepted options.\n", .{ argv[1], unknown });
    defer allocator.free(expected);
    try expectCliFailure(allocator, argv, 1, expected);
}

pub fn uniqueArtifactPath(allocator: std.mem.Allocator, stem: []const u8, ext: []const u8) ![]u8 {
    const io = std.testing.io;
    try std.Io.Dir.cwd().createDirPath(io, "zig-cache/test-artifacts");
    var nonce_bytes: [16]u8 = undefined;
    io.random(&nonce_bytes);
    const nonce = std.mem.readInt(u128, &nonce_bytes, .little);
    return std.fmt.allocPrint(allocator, "zig-cache/test-artifacts/{s}-{x}.{s}", .{ stem, nonce, ext });
}

pub fn writeFastaArtifact(allocator: std.mem.Allocator, stem: []const u8, data: []const u8) ![]u8 {
    const io = std.testing.io;
    const path = try uniqueArtifactPath(allocator, stem, "fa");
    errdefer allocator.free(path);
    errdefer std.Io.Dir.cwd().deleteFile(io, path) catch {};
    const file = try std.Io.Dir.cwd().createFile(io, path, .{});
    defer file.close(io);
    try std.Io.File.writeStreamingAll(file, io, data);
    return path;
}

pub fn writeZfi(
    allocator: std.mem.Allocator,
    fasta_path: []const u8,
    data: []const u8,
) !void {
    const io = std.testing.io;
    var index = try main.indexer.scanZfiData(allocator, data);
    defer index.deinit(allocator);

    var path_buffer: [4096]u8 = undefined;
    const zfi_path = try std.fmt.bufPrint(&path_buffer, "{s}.zfi", .{fasta_path});
    try main.indexer.writeZfiIndexFile(io, zfi_path, &index, data.len, try statMtimeNs(fasta_path));
}

pub fn captureExtractRegion(
    allocator: std.mem.Allocator,
    fasta_path: []const u8,
    region: []const u8,
) ![]u8 {
    var index = try loadIndex(allocator, std.testing.io, fasta_path, .lookup_full_map);
    defer index.deinit();

    var output = std.Io.Writer.Allocating.init(allocator);
    errdefer output.deinit();
    var failure: main.index_format.diagnostics.Failure = undefined;
    try main.getter.extractRegion(&index, region, &output.writer, &failure);
    return output.toOwnedSlice();
}

pub fn statMtimeNs(path: []const u8) !u64 {
    const io = std.testing.io;
    const file = try std.Io.Dir.cwd().openFile(io, path, .{});
    defer file.close(io);
    return main.index_format.timestampToNs((try file.stat(io)).mtime);
}

pub fn readRequiredFile(allocator: std.mem.Allocator, path: []const u8) ![]u8 {
    const io = std.testing.io;
    const file = std.Io.Dir.cwd().openFile(io, path, .{}) catch |err| switch (err) {
        error.FileNotFound => {
            std.debug.print("required test fixture missing: {s}\n", .{path});
            return error.RequiredFixtureMissing;
        },
        else => |e| return e,
    };
    defer file.close(io);
    const stat = try file.stat(io);
    const size = std.math.cast(usize, stat.size) orelse return error.FileTooBig;
    const data = try allocator.alloc(u8, size);
    errdefer allocator.free(data);
    const bytes_read = try file.readPositionalAll(io, data, 0);
    if (bytes_read != data.len) return error.UnexpectedEndOfFile;
    return data;
}
