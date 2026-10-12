//! The `z-fasta` package. Only its public declarations are supported. Calls fill a caller-owned
//! `diagnostics.Failure` and never print or exit.
//! `IndexedReader` opens a FASTA and its sidecar once and fetches many intervals from it.
//! Intervals are zero-based half-open and never clipped. A reader is not safe for concurrent use,
//! and its `io` must outlive it.

const std = @import("std");

pub const diagnostics = @import("diagnostics");
const getter = @import("getter.zig");
const index_format = @import("index_format.zig");

pub const Interval = struct {
    start: u64,
    end: u64,
};

/// `name` borrows reader storage until `deinit`.
pub const RecordInfo = struct {
    ordinal: usize,
    name: []const u8,
    length: u64,
};

/// Automatic prefers a present `.zfi`; an explicit path never falls back.
pub const IndexSelection = index_format.IndexSelection;

pub const IndexedReader = struct {
    index: index_format.LoadedIndex,
    read_buffer: []u8,

    pub const Options = struct {
        index: IndexSelection = .automatic,
    };

    pub fn open(allocator: std.mem.Allocator, io: std.Io, path: []const u8, options: Options, failure: *diagnostics.Failure) diagnostics.Error!IndexedReader {
        failure.reset(.get);
        var index = try index_format.loadIndex(allocator, io, path, options.index, .lookup_full_map, failure);
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

    pub fn deinit(self: *IndexedReader) void {
        self.index.deinit();
        self.* = undefined;
    }

    /// Counts repeated names separately.
    pub fn recordCount(self: *const IndexedReader) usize {
        return self.index.records.len;
    }

    pub fn record(self: *const IndexedReader, ordinal: usize) ?RecordInfo {
        if (ordinal >= self.index.records.len) return null;
        return .{
            .ordinal = ordinal,
            .name = self.index.recordName(ordinal) orelse "",
            .length = self.index.records[ordinal].seq_len,
        };
    }

    /// First match wins.
    pub fn lookup(self: *const IndexedReader, name: []const u8) ?RecordInfo {
        return self.record(self.index.lookupName(name) orelse return null);
    }

    /// Checks everything before writing `dst`; a later read failure can leave a prefix.
    pub fn fetchInto(self: *IndexedReader, ordinal: usize, interval: Interval, dst: []u8, failure: *diagnostics.Failure) diagnostics.Error![]u8 {
        return getter.fetchSymbols(&self.index, ordinal, interval.start, interval.end, dst, self.read_buffer, failure);
    }
};

const READ_BUFFER_BYTES = 256 * 1024;

fn rejectCompressed(index: index_format.LoadedIndex, failure: *diagnostics.Failure) diagnostics.Error!void {
    var magic: [2]u8 = undefined;
    const got = index.fasta_file.readPositionalAll(index.io, &magic, 0) catch |err| return failure.fail(.{
        .details = diagnostics.readFailure(magic.len, null),
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
