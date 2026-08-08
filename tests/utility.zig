//! Shared test mechanics for command suites.

const std = @import("std");
const builtin = @import("builtin");
const main = @import("main");

pub const ZFASTA_BIN = if (builtin.os.tag == .windows) "zig-out\\bin\\z-fasta.exe" else "zig-out/bin/z-fasta";

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

pub fn expectCliFailure(
    allocator: std.mem.Allocator,
    argv: []const []const u8,
    exit_code: u8,
    expected_stderr: []const u8,
) !void {
    try expectCliResult(allocator, argv, exit_code, "", expected_stderr);
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
    const expected = try std.fmt.allocPrint(allocator, "error: unknown option: {s}\n", .{unknown});
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
    enable_dedup: bool,
) !void {
    const io = std.testing.io;
    var index = try main.indexer.scanZfiData(data, enable_dedup, allocator);
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
    var index = try main.index_format.loadIndexChecked(allocator, std.testing.io, fasta_path);
    defer index.deinit();

    var output = std.Io.Writer.Allocating.init(allocator);
    errdefer output.deinit();
    main.getter.extractRegion(&index, region, &output.writer);
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
