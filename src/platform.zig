//! Cross-platform policy for read-only file mapping and memory advice.

const std = @import("std");
const builtin = @import("builtin");

/// Read-only mapped bytes aligned for platform mapping APIs. The corresponding
/// `MemoryMap` owns the lifetime.
pub const MappedBytes = []align(std.heap.page_size_min) const u8;

/// Maps exactly `len` nonzero bytes read-only without prefaulting pages on POSIX.
///
/// Zig 0.16's Windows `NtCreateSection` path requires `.populate = true`; false
/// produces `STATUS_INVALID_PARAMETER`. Call `MemoryMap.destroy` when finished.
pub fn mapFileReadOnly(io: std.Io, file: std.Io.File, len: usize) error{MmapFailed}!std.Io.File.MemoryMap {
    return file.createMemoryMap(io, .{
        .len = len,
        .protection = .{ .read = true, .write = false },
        .populate = builtin.os.tag == .windows,
    }) catch error.MmapFailed;
}

/// Hints that a mapped region will be read sequentially. No-op on Windows.
pub fn adviseSequential(memory: MappedBytes) void {
    if (comptime builtin.os.tag == .windows) return;
    if (memory.len == 0) return;
    std.posix.madvise(
        @constCast(memory.ptr),
        memory.len,
        std.posix.MADV.SEQUENTIAL,
    ) catch {};
}

test "[integration] - [read-only mapping]: exposes file contents" {
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
