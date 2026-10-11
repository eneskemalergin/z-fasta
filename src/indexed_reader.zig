//! Open-once indexed FASTA reading: record metadata and strict interval fetches.
//!
//! Intervals are zero-based and half-open, and are never clipped. A reader owns its
//! source descriptor, index storage and one read buffer until `deinit`.

const std = @import("std");

const diagnostics = @import("diagnostics");
const getter = @import("getter.zig");
const index_format = @import("index_format.zig");

/// Zero-based half-open symbol interval `[start, end)` within one record.
pub const Interval = struct {
    start: u64,
    end: u64,
};

/// One retained index entry. `name` borrows reader storage until `deinit`.
pub const RecordInfo = struct {
    ordinal: usize,
    name: []const u8,
    length: u64,
};

/// Sidecar choice: automatic discovery keeps a present `.zfi` authoritative and
/// falls back to `.fai` only when `.zfi` is absent; an explicit path never falls back.
pub const IndexSelection = index_format.IndexSelection;

/// Reads intervals from one FASTA through its `.zfi` or `.fai` sidecar. Not safe for
/// concurrent use; open one reader per thread.
pub const IndexedReader = struct {
    index: index_format.LoadedIndex,
    read_buffer: []u8,

    pub const Options = struct {
        index: IndexSelection = .automatic,
    };

    /// Opens `path` and its selected sidecar and checks index structure and source
    /// identity. Failure storage resets once. `path` and explicit sidecar paths are
    /// borrowed only during the call; `io` must outlive the reader.
    pub fn open(allocator: std.mem.Allocator, io: std.Io, path: []const u8, options: Options, failure: *diagnostics.Failure) diagnostics.Error!IndexedReader {
        failure.reset(.get);
        var index = try index_format.loadIndexSelected(allocator, io, path, options.index, .lookup_full_map, failure);
        errdefer index.deinit();
        try rejectCompressed(index, failure);
        const read_buffer = index.arena.allocator().alloc(u8, READ_BUFFER_BYTES) catch return failure.fail(.{
            .details = .{ .@"resource.out_of_memory" = .{ .allocation_role = .reader_buffer, .requested_bytes = READ_BUFFER_BYTES } },
            .stage = .allocate,
            .at = .{ .source = .fasta },
            .cause = .{ .zig = error.OutOfMemory },
        });
        failure.succeed();
        return .{ .index = index, .read_buffer = read_buffer };
    }

    /// Closes the source and releases all reader storage.
    pub fn deinit(self: *IndexedReader) void {
        self.index.deinit();
        self.* = undefined;
    }

    /// Number of retained index entries, including repeated names.
    pub fn recordCount(self: *const IndexedReader) usize {
        return self.index.records.len;
    }

    /// Returns the entry at `ordinal`, or null when it is out of range.
    pub fn record(self: *const IndexedReader, ordinal: usize) ?RecordInfo {
        if (ordinal >= self.index.records.len) return null;
        return .{
            .ordinal = ordinal,
            .name = self.index.recordName(ordinal) orelse "",
            .length = self.index.records[ordinal].seq_len,
        };
    }

    /// Returns the first entry whose name equals `name` byte for byte.
    pub fn lookup(self: *const IndexedReader, name: []const u8) ?RecordInfo {
        return self.record(self.index.lookupName(name) orelse return null);
    }

    /// Copies the symbols of `interval` in record `ordinal` to the front of `dst`
    /// and returns that prefix, preserving case and omitting line breaks. Failure
    /// storage resets once. Ordinal, interval and capacity are checked before `dst`
    /// changes; a later read failure can leave a changed prefix. `dst` must not
    /// alias reader storage.
    pub fn fetchInto(self: *IndexedReader, ordinal: usize, interval: Interval, dst: []u8, failure: *diagnostics.Failure) diagnostics.Error![]u8 {
        return getter.fetchSymbols(&self.index, ordinal, interval.start, interval.end, dst, self.read_buffer, failure);
    }
};

const READ_BUFFER_BYTES = 256 * 1024;

fn rejectCompressed(index: index_format.LoadedIndex, failure: *diagnostics.Failure) diagnostics.Error!void {
    var magic: [2]u8 = undefined;
    const got = index.fasta_file.readPositionalAll(index.io, &magic, 0) catch |err| return failure.fail(.{
        .details = .{ .@"io.read_failed" = .{ .requested_bytes = magic.len } },
        .stage = .read,
        .at = .{ .source = .fasta, .location = .{ .offset = 0 } },
        .cause = .{ .zig = err },
    });
    if (got == magic.len and magic[0] == 0x1f and magic[1] == 0x8b) return failure.fail(.{
        .details = .{ .@"input.unsupported_format" = .{ .detected = .gzip } },
        .stage = .read,
        .at = .{ .source = .fasta, .location = .{ .offset = 0 } },
    });
}
