//! Read-only file mapping and memory advice, with the platform differences kept here.

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
