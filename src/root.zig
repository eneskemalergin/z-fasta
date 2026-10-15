//! The `z-fasta` package. Only its public declarations are supported. Calls fill a caller-owned
//! `diagnostics.Failure` and never print or exit.
//! `IndexedReader` opens a FASTA and its sidecar once and fetches many intervals from it; a
//! `.zfi` is checked as its records are used, so opening costs the same for any index size.
//! Intervals are zero-based half-open and never clipped. A reader is not safe for concurrent use,
//! and its `io` must outlive it.
//! `RecordReader` streams every record of a FASTA from any `std.Io.Reader`, by `index`'s rules.
//! `validation` checks a FASTA stream in one pass, and can write a repaired copy in the same pass.

const std = @import("std");

pub const diagnostics = @import("diagnostics");
const fasta = @import("fasta.zig");
const getter = @import("getter.zig");
const index_format = @import("index_format.zig");
const validator = @import("validator.zig");

pub const RecordReader = fasta.RecordReader;

pub const validation = struct {
    pub const Rules = validator.Rules;
    pub const Report = validator.Report;
    pub const NameSource = validator.NameSource;
    pub const READ_BUFFER_BYTES = validator.READ_BUFFER_BYTES;
    pub const validate = validator.validate;
    pub const repair = validator.repair;
};

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
        var index = try index_format.loadIndex(allocator, io, path, options.index, .on_access, failure);
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

    pub fn recordCount(self: *const IndexedReader) usize {
        return self.index.records.len();
    }

    pub fn record(self: *IndexedReader, ordinal: usize, failure: *diagnostics.Failure) diagnostics.Error!?RecordInfo {
        failure.reset(.get);
        if (ordinal >= self.index.records.len()) {
            failure.succeed();
            return null;
        }
        const rec = try self.index.checkedRecord(ordinal, failure);
        failure.succeed();
        return .{ .ordinal = ordinal, .name = self.index.recordName(ordinal) orelse "", .length = rec.seq_len };
    }

    /// First match wins.
    pub fn lookup(self: *IndexedReader, name: []const u8, failure: *diagnostics.Failure) diagnostics.Error!?RecordInfo {
        failure.reset(.get);
        const ordinal = try self.index.findName(name, failure) orelse {
            failure.succeed();
            return null;
        };
        return self.record(ordinal, failure);
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
