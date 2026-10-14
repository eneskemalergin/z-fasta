//! Read-only file mapping, sidecar opening, and memory advice, with the platform differences kept here.

const std = @import("std");
const builtin = @import("builtin");

pub const MappedBytes = []align(std.heap.page_size_min) const u8;

/// Windows needs `.populate = true` (Zig 0.16 fails with `STATUS_INVALID_PARAMETER` otherwise).
pub fn mapFileReadOnly(io: std.Io, file: std.Io.File, len: usize) std.Io.File.MemoryMap.CreateError!std.Io.File.MemoryMap {
    return file.createMemoryMap(io, .{
        .len = len,
        .protection = .{ .read = true, .write = false },
        .populate = builtin.os.tag == .windows,
    });
}

/// No-op on Windows.
pub fn adviseSequential(memory: MappedBytes) void {
    if (comptime builtin.os.tag == .windows) return;
    if (memory.len == 0) return;
    std.posix.madvise(
        @constCast(memory.ptr),
        memory.len,
        std.posix.MADV.SEQUENTIAL,
    ) catch {};
}

pub const OpenSidecarError = std.posix.OpenError || std.posix.UnexpectedError;

/// Opens an index path derived from the FASTA. A FIFO there opens without waiting for a writer;
/// the returned descriptor blocks as `openFile` would.
pub fn openSidecar(io: std.Io, path: []const u8) OpenSidecarError!std.Io.File {
    const handle = try std.posix.openat(std.posix.AT.FDCWD, path, .{ .ACCMODE = .RDONLY, .NONBLOCK = true, .CLOEXEC = true }, 0);
    const file = std.Io.File{ .handle = handle, .flags = .{ .nonblocking = false } };
    errdefer file.close(io);
    const status = std.posix.system.fcntl(handle, std.posix.F.GETFL, @as(usize, 0));
    switch (std.posix.errno(status)) {
        .SUCCESS => {},
        else => |err| return std.posix.unexpectedErrno(err),
    }
    const blocking = @as(usize, @intCast(status)) & ~@as(usize, 1 << @bitOffsetOf(std.posix.O, "NONBLOCK"));
    switch (std.posix.errno(std.posix.system.fcntl(handle, std.posix.F.SETFL, blocking))) {
        .SUCCESS => {},
        else => |err| return std.posix.unexpectedErrno(err),
    }
    return file;
}

test "[integration] - [file mapping]: exposes file contents" {
    const contents = "mapped bytes";
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();

    const file = try tmp.dir.createFile(std.testing.io, "input", .{ .read = true });
    defer file.close(std.testing.io);
    try std.Io.File.writeStreamingAll(file, std.testing.io, contents);

    var map = try mapFileReadOnly(std.testing.io, file, contents.len);
    defer map.destroy(std.testing.io);
    adviseSequential(map.memory);

    try std.testing.expectEqualStrings(contents, map.memory);
}
