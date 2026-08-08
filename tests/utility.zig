//! Shared process assertions and artifact paths for command test suites.

const std = @import("std");
const builtin = @import("builtin");

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
    const now = std.Io.Clock.Timestamp.now(io, .awake);
    const nanos: u64 = @intCast(now.raw.toNanoseconds());
    return std.fmt.allocPrint(allocator, "zig-cache/test-artifacts/{s}-{d}.{s}", .{ stem, nanos, ext });
}
