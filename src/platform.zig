//! Read-only file mapping and sidecar opening: the POSIX calls the commands share.

const std = @import("std");
const builtin = @import("builtin");

comptime {
    if (builtin.os.tag == .windows) @compileError("z-fasta supports Linux and macOS; on Windows, run the Linux build through WSL");
}

pub const MappedBytes = []align(std.heap.page_size_min) const u8;

/// Pages load on first access: the loaders read only the parts they check.
pub fn mapFileReadOnly(io: std.Io, file: std.Io.File, len: usize) std.Io.File.MemoryMap.CreateError!std.Io.File.MemoryMap {
    return file.createMemoryMap(io, .{
        .len = len,
        .protection = .{ .read = true, .write = false },
        .populate = false,
    });
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

    try std.testing.expectEqualStrings(contents, map.memory);
}
