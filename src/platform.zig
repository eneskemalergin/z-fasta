//! Portable file mapping and memory advice.
//!
//! Core modules must not call `std.posix.mmap` / `munmap` / `madvise` directly.
//! Mapping goes through `std.Io.File.MemoryMap` (native Windows path included).
//! Advice is a best-effort optimization: real on Linux/macOS, no-op on Windows.

const std = @import("std");
const builtin = @import("builtin");

pub const mapped_align = std.heap.page_size_min;
pub const MappedBytes = []align(mapped_align) const u8;

pub const Advice = enum {
    sequential,
    random,
    dont_need,
};

/// Owns a whole-file mapping.
///
/// On POSIX, `.populate = false` keeps pages lazy. On Windows, Zig 0.16's
/// `NtCreateSection` path passes `{ .COMMIT = populate }` as allocation
/// attributes; `populate = false` yields attributes 0 and
/// `STATUS_INVALID_PARAMETER` (0xc00000f4), so Windows must use `populate = true`.
pub const FileView = struct {
    map: std.Io.File.MemoryMap,

    pub fn mapFile(io: std.Io, file: std.Io.File, len: usize) MapError!FileView {
        const map = file.createMemoryMap(io, .{
            .len = len,
            .protection = .{ .read = true, .write = false },
            .populate = builtin.os.tag == .windows,
        }) catch return error.MmapFailed;
        return .{ .map = map };
    }

    pub fn bytes(self: *const FileView) MappedBytes {
        return self.map.memory;
    }

    pub fn destroy(self: *FileView, io: std.Io) void {
        self.map.destroy(io);
    }
};

pub const MapError = error{MmapFailed};

/// Optional access hint for a mapped region. No-op on Windows.
pub fn advise(memory: []const u8, kind: Advice) void {
    if (comptime builtin.os.tag == .windows) return;
    if (memory.len == 0) return;
    const len = std.mem.alignBackward(usize, memory.len, mapped_align);
    if (len == 0) return;
    const madv: u32 = switch (kind) {
        .sequential => std.posix.MADV.SEQUENTIAL,
        .random => std.posix.MADV.RANDOM,
        .dont_need => std.posix.MADV.DONTNEED,
    };
    std.posix.madvise(@alignCast(@constCast(memory.ptr)), len, madv) catch {};
}

/// Release (or hint release of) pages covering `[start, end_exclusive)`.
pub fn releaseSpan(memory: []const u8, start: usize, end_exclusive: usize) void {
    if (comptime builtin.os.tag == .windows) return;
    if (start >= end_exclusive or memory.len == 0) return;
    const drop_start = std.mem.alignBackward(usize, start, mapped_align);
    const drop_end = @min(memory.len, std.mem.alignForward(usize, end_exclusive, mapped_align));
    if (drop_end <= drop_start) return;
    std.posix.madvise(
        @alignCast(@constCast(memory.ptr + drop_start)),
        drop_end - drop_start,
        std.posix.MADV.DONTNEED,
    ) catch {};
}
