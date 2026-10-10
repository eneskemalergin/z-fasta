//! FASTA indexing through one bounded reader for `.zfi` and `.fai` output.
//!
//! The parser skips empty lines for geometry, sets `seq_offset` to the first base-bearing
//! line, invents a separator when the final line has no newline, and rejects names longer
//! than `MAX_INDEX_NAME_LEN`. Proven dense blocks update the same `ChunkParseState`; the
//! first anomaly resumes the general line machine at the exact unconsumed byte.
//!
//! Duplicate-name filtering is first-wins and collision-safe: lookup may use a hash, but
//! identity is always full-string equality, including source-backed name references.

const std = @import("std");
const builtin = @import("builtin");
const index_format = @import("index_format.zig");
const diagnostics = @import("diagnostics");

/// Size in bytes of the scanner's sequence read buffer.
pub const INDEX_READ_BUFFER_SIZE = 256 * 1024;

/// Maximum identifier length in bytes for both index formats.
pub const MAX_INDEX_NAME_LEN = std.math.maxInt(u16);

/// Reader policy; identifier limits may be stricter than `MAX_INDEX_NAME_LEN`.
pub const ScanOptions = struct {
    enable_dedup: bool,
    max_name_len: ?usize = null,
    collect_side_tables: bool = false,
    require_initial_header: bool = false,
    /// Borrowed diagnostic destination, initialized by the owning operation.
    failure: ?*diagnostics.Failure = null,
};

/// Selects CLI output format and first-occurrence duplicate filtering.
pub const Options = struct {
    emit_fai: bool = false,
    enable_dedup: bool = true,
};

/// Number of retained nonempty records after publication and required cleanup.
pub const Result = struct { record_count: u32 };

/// First-wins name set whose decisions use full-string equality despite hash collisions.
/// FAI owns stable keys in a child arena; `.zfi` shares names with its output blob.
pub const NameDedup = NameDedupWith(std.hash_map.StringContext);

/// One parsed FASTA record passed to reader emit callbacks.
/// Name and side-table slices are borrowed only for the callback's duration.
pub const FastaRecordEmit = struct {
    record: IndexRecord,
    name: []const u8,
    /// Original FASTA position; borrowed only for the callback duration.
    source_reference: ?*const diagnostics.Reference = null,
    nonuniform_reason: diagnostics.NonuniformReason = .nonuniform_layout,
    uses_uniform_formula: bool,
    /// Incremental side-table bytes for `.zfi` (line count plus `SideTableLine` rows).
    side_table: []const u8 = &.{},
    /// When true, `record` already has blob-relative `name_offset`/`name_len` and `NAME_IN_ZFI_FLAG`
    /// (`.zfi` blob-backed `NameDedup`); skip `embedZfiName`.
    name_embedded: bool = false,
};

/// Owned records, side tables, and names; release with `deinit` using the scan allocator.
pub const ZfiIndex = struct {
    records: std.ArrayList(IndexRecord),
    side_tables: std.ArrayList(u8),
    name_blob: std.ArrayList(u8),

    /// Releases all catalog storage with the allocator used to create it.
    pub fn deinit(self: *ZfiIndex, allocator: std.mem.Allocator) void {
        self.records.deinit(allocator);
        self.side_tables.deinit(allocator);
        self.name_blob.deinit(allocator);
    }
};

/// Collects an owned `.zfi` catalog; `read_buf` must be nonempty and remains caller-owned.
pub fn scanZfiReader(
    allocator: std.mem.Allocator,
    reader: *std.Io.Reader,
    read_buf: []u8,
    enable_dedup: bool,
) !ZfiIndex {
    return scanZfiReaderWithOptions(
        allocator,
        reader,
        read_buf,
        .{ .enable_dedup = enable_dedup },
    );
}

/// Collects an owned `.zfi` catalog from borrowed FASTA bytes.
pub fn scanZfiData(
    allocator: std.mem.Allocator,
    data: []const u8,
    enable_dedup: bool,
) !ZfiIndex {
    var r = std.Io.Reader.fixed(data);
    var read_buf: [INDEX_READ_BUFFER_SIZE]u8 = undefined;
    return scanZfiReader(allocator, &r, &read_buf, enable_dedup);
}

/// Writes a collected catalog directly to `path`; the CLI owns atomic publication.
pub fn writeZfiIndexFile(
    io: std.Io,
    path: []const u8,
    index: *const ZfiIndex,
    source_size: u64,
    source_mtime_ns: u64,
) !void {
    _ = try checkedZfiRecordCount(index.records.items.len, .{});
    const file = try std.Io.Dir.cwd().createFile(io, path, .{});
    defer file.close(io);

    var file_buf: [INDEX_OUTPUT_BUFFER_SIZE]u8 = undefined;
    var file_fw = file.writer(io, &file_buf);
    try writeZfiIndex(&file_fw.interface, index, source_size, source_mtime_ns);
    try file_fw.flush();
}

/// Caller owns the returned production-layout bytes.
pub fn zfiIndexToBytes(
    allocator: std.mem.Allocator,
    index: *const ZfiIndex,
    source_size: u64,
    source_mtime_ns: u64,
) ![]u8 {
    return serializeZfiIndex(allocator, index, source_size, source_mtime_ns, .{});
}

/// Serializes an owned index without resetting the caller's diagnostic destination.
/// The returned bytes belong to `allocator`; rejects counts above the ZFI u32 limit.
/// Allocation failures retain OOM facts without allocating diagnostic storage.
pub fn zfiIndexToBytesChecked(allocator: std.mem.Allocator, index: *const ZfiIndex, source_size: u64, source_mtime_ns: u64, failure: *diagnostics.Failure) diagnostics.Error![]u8 {
    return serializeZfiIndex(allocator, index, source_size, source_mtime_ns, .{ .failure = failure }) catch return diagnostics.metadata(failure.primary().?.code()).returned_error;
}

/// Scans records through a nonempty caller-owned buffer and returns the retained count.
/// Callback slices expire on return; an external name set is borrowed and may be updated.
/// Returns `TooManySequences` before emitting more than 4294967295 retained records.
/// Rejects empty buffers and buffers larger than 4294967295 bytes before reading.
pub fn scanFastaReader(
    allocator: std.mem.Allocator,
    reader: *std.Io.Reader,
    read_buf: []u8,
    options: ScanOptions,
    external_dedup: ?*NameDedup,
    ctx: anytype,
    comptime emitRecord: fn (@TypeOf(ctx), FastaRecordEmit) anyerror!void,
) !u32 {
    try checkReadBufferSize(read_buf.len, options.failure);
    return scanFastaReaderWithDedup(allocator, reader, read_buf, options, external_dedup, ctx, emitRecord);
}

/// Writes `.fai` rows to `writer`; errors may leave a prefix in this low-level API.
pub fn scanFaiData(
    allocator: std.mem.Allocator,
    data: []const u8,
    writer: anytype,
    enable_dedup: bool,
) !u32 {
    var r = std.Io.Reader.fixed(data);
    var read_buf: [INDEX_READ_BUFFER_SIZE]u8 = undefined;
    return scanFaiReader(allocator, &r, &read_buf, writer, .{ .enable_dedup = enable_dedup });
}

/// Indexes a file and returns after publication and required temporary-file removal.
/// The caller owns `failure` and `path` through diagnostic handling. This entry resets
/// `failure` once. ZFI replaces its sidecar atomically; FAI stdout starts after a valid scan.
pub fn runIndex(
    io: std.Io,
    environ: std.process.Environ,
    path: []const u8,
    options: Options,
    failure: *diagnostics.Failure,
) diagnostics.Error!Result {
    failure.reset(.index);
    failure.sources.set(.fasta, .{ .label = .{ .borrowed = path }, .kind = .path });
    failure.sources.set(.zfi, .{ .label = .{ .borrowed = path }, .kind = .path, .suffix = .zfi });
    failure.sources.set(.stdout, .{ .label = .{ .borrowed = "stdout" }, .kind = .stream });
    const diag = Diagnostics{ .failure = failure, .format = if (options.emit_fai) .fai else .zfi };
    const file = std.Io.Dir.cwd().openFile(io, path, .{}) catch |err| return failure.fail(.{
        .details = switch (err) {
            error.FileNotFound => .{ .@"io.file_not_found" = {} },
            error.AccessDenied => .{ .@"io.access_denied" = {} },
            error.NameTooLong => .{ .@"io.invalid_path" = .{ .reason = .too_long } },
            else => .{ .@"io.open_failed" = {} },
        },
        .stage = .open,
        .at = .{ .source = .fasta },
        .cause = .{ .zig = err },
    });
    defer file.close(io);
    const stat = file.stat(io) catch |err| return failure.fail(.{ .details = .{ .@"io.stat_failed" = {} }, .stage = .stat, .at = .{ .source = .fasta }, .cause = .{ .zig = err } });
    if (stat.size == 0) return failure.fail(.{ .details = .{ .@"fasta.no_records" = {} }, .stage = .stat, .at = .{ .source = .fasta } });
    const count = (if (stat.size <= SMALL_INDEX_MAX_BYTES)
        runIndexSmall(io, environ, path, options, file, stat, diag)
    else
        runIndexLarge(io, environ, path, options, file, stat, diag)) catch |err| {
        if (failure.primary()) |primary| return diagnostics.metadata(primary.code()).returned_error;
        return failure.fail(.{ .details = .{ .@"internal.consistency_failed" = .{ .component = .indexer, .check = .codec_state } }, .stage = .scan, .at = .{ .source = .fasta }, .cause = .{ .zig = err } });
    };
    failure.work_complete = true;
    try failure.finishDeferred();
    failure.succeed();
    return .{ .record_count = count };
}

const IndexRecord = index_format.IndexRecord;
const ZfiHeader = index_format.ZfiHeader;
const ZFI_MAGIC = index_format.ZFI_MAGIC;

const SIMD_CHUNK_SIZE = 32;
const SimdVec = @Vector(SIMD_CHUNK_SIZE, u8);
const STRIDE_VECTOR_LEN = std.simd.suggestVectorLength(u8) orelse 1;
const STRIDE_VALIDATION_BLOCK_LINES = 256;

const FILE_IO_BUF_SIZE = 8 * 1024;
const INDEX_OUTPUT_BUFFER_SIZE = 64 * 1024;
const INDEX_PATH_BUFFER_SIZE = std.Io.Dir.max_path_bytes;
const SMALL_INDEX_MAX_BYTES = 32 * 1024 * 1024;
const ZFI_INITIAL_RECORD_CAPACITY = 4096;
const DEDUP_INITIAL_CAPACITY: u32 = 16384;
const FAI_SPOOL_NAME_BUFFER_SIZE = 64;
const FAI_U64_DECIMAL_DIGITS = 20;
const FAI_SUFFIX_BUFFER_SIZE = 4 * (1 + FAI_U64_DECIMAL_DIGITS) + 1;
const FAI_SPOOL_CREATE_ATTEMPTS = 16;

const IndexReference = struct {
    source: ?diagnostics.SourceId = null,
    location: struct { offset: ?u64 = null } = .{},
    record: ?u64 = null,

    fn reference(self: @This()) diagnostics.Reference {
        return .{ .source = self.source, .location = .{ .offset = self.location.offset }, .record = self.record };
    }
};

const Diagnostics = struct {
    failure: ?*diagnostics.Failure = null,
    format: diagnostics.Format = .zfi,
    deferred: bool = false,
    at: IndexReference = .{ .source = .fasta },

    fn forRecord(self: @This(), info: FastaRecordEmit) @This() {
        var diag = self;
        if (info.source_reference) |reference| diag.at = .{ .source = reference.source, .location = .{ .offset = reference.location.offset }, .record = reference.record };
        return diag;
    }

    fn failReadAt(self: @This(), source: diagnostics.SourceId, offset: u64, requested: u64, received: ?u64, cause: ?anyerror, legacy: anyerror) anyerror {
        @branchHint(.cold);
        var diag = self;
        diag.at = .{ .source = source, .location = .{ .offset = offset } };
        const details: diagnostics.Details = if (received) |actual|
            .{ .@"io.unexpected_eof" = .{ .exact = .{ .requested_bytes = requested, .received_bytes = actual } } }
        else
            .{ .@"io.read_failed" = .{ .requested_bytes = requested } };
        return diag.fail(details, .read, source, cause, legacy);
    }

    fn fail(self: @This(), details: diagnostics.Details, stage: diagnostics.Stage, source: diagnostics.SourceId, cause: ?anyerror, legacy: anyerror) anyerror {
        @branchHint(.cold);
        const failure = self.failure orelse return legacy;
        var at = self.at;
        if (at.source != source) at = .{};
        at.source = source;
        const diagnostic = diagnostics.Diagnostic{ .details = details, .stage = stage, .at = at.reference(), .cause = if (cause) |err| .{ .zig = err } else null };
        if (self.deferred) {
            failure.deferFailure(diagnostic);
            return diagnostics.metadata(diagnostic.code()).returned_error;
        }
        return failure.fail(diagnostic);
    }

    fn oom(self: @This(), role: diagnostics.AllocationRole, bytes: ?u64) error{OutOfMemory} {
        @branchHint(.cold);
        if (self.failure) |failure| {
            switch (failure.fail(.{ .details = .{ .@"resource.out_of_memory" = .{ .allocation_role = role, .requested_bytes = bytes } }, .stage = .allocate, .at = self.at.reference(), .cause = .{ .zig = error.OutOfMemory } })) {
                else => {},
            }
        }
        return error.OutOfMemory;
    }

    fn capacity(self: @This(), role: diagnostics.CapacityRole, limit: u128, required: u128, legacy: anyerror) anyerror {
        return self.fail(.{ .@"resource.capacity_limit_exceeded" = .{ .capacity_role = role, .limit = limit, .required = .{ .exact = required } } }, .allocate, self.at.source orelse .fasta, null, legacy);
    }

    fn consistency(self: @This(), check: diagnostics.Check, legacy: anyerror) anyerror {
        return self.fail(.{ .@"internal.consistency_failed" = .{ .component = .indexer, .check = check } }, .scan, self.at.source orelse .fasta, null, legacy);
    }

    fn lineLimit(self: @This(), field: diagnostics.IndexField, observed: u64, width_role: diagnostics.LineWidthRole) anyerror {
        return self.fail(.{ .@"index.sequence_line_too_long" = .{ .field = field, .observed = .{ .exact = observed }, .limit = std.math.maxInt(u32), .width_role = width_role, .format = self.format } }, .scan, .fasta, null, error.SequenceLineTooLong);
    }

    fn recordLimit(self: @This(), count: u64) anyerror {
        return self.fail(.{ .@"index.record_count_limit_exceeded" = .{ .record_count = .{ .exact = count }, .max_records = std.math.maxInt(u32), .format = self.format, .limit_role = .scanner_count } }, .scan, .fasta, null, error.TooManySequences);
    }

    fn markNonUniform(self: @This(), record: *IndexRecord, offset: u64, role: diagnostics.OffsetRole) !void {
        record.markNonUniform(offset) catch return self.fail(.{ .@"index.side_table_offset_limit_exceeded" = .{ .offset = offset, .max_offset = index_format.MAX_SIDE_TABLE_OFFSET, .offset_role = role } }, .scan, .zfi, null, error.SideTableOffsetTooLarge);
    }
};

const NameBlobRef = struct { offset: u32, len: u32 };

fn NameDedupWith(comptime Context: type) type {
    return struct {
        diag: Diagnostics = .{},
        map: std.HashMap([]const u8, void, Context, std.hash_map.default_max_load_percentage),
        key_arena: ?std.heap.ArenaAllocator = null,
        name_blob: ?*std.ArrayList(u8) = null,
        ref_map: ?BlobRefMap = null,
        last_offset: u32 = 0,
        last_len: u16 = 0,

        const BlobHashContext = struct {
            blob: *std.ArrayList(u8),
            pub fn hash(self: @This(), key: NameBlobRef) u64 {
                return std.hash_map.hashString(self.blob.items[key.offset..][0..key.len]);
            }
            pub fn eql(self: @This(), a: NameBlobRef, b: NameBlobRef) bool {
                if (a.len != b.len) return false;
                return std.mem.eql(
                    u8,
                    self.blob.items[a.offset..][0..a.len],
                    self.blob.items[b.offset..][0..b.len],
                );
            }
        };

        const BlobRefMap = std.HashMap(
            NameBlobRef,
            void,
            BlobHashContext,
            std.hash_map.default_max_load_percentage,
        );

        const SliceLookup = struct {
            blob: *std.ArrayList(u8),
            pub fn hash(_: @This(), name: []const u8) u64 {
                return std.hash_map.hashString(name);
            }
            pub fn eql(self: @This(), name: []const u8, ref: NameBlobRef) bool {
                if (name.len != ref.len) return false;
                return std.mem.eql(u8, name, self.blob.items[ref.offset..][0..ref.len]);
            }
        };

        const Self = @This();

        /// Owns a name set using `allocator` until `deinit`.
        pub fn init(allocator: std.mem.Allocator) Self {
            return .{
                .map = std.HashMap([]const u8, void, Context, std.hash_map.default_max_load_percentage).init(allocator),
                .key_arena = std.heap.ArenaAllocator.init(allocator),
            };
        }

        fn observeAt(self: *Self, name: []const u8, offset: u64) !bool {
            _ = offset;
            return self.observe(name);
        }

        fn initOwningBlob(allocator: std.mem.Allocator, blob: *std.ArrayList(u8)) Self {
            return .{
                .map = std.HashMap([]const u8, void, Context, std.hash_map.default_max_load_percentage).init(allocator),
                .name_blob = blob,
                .ref_map = BlobRefMap.initContext(allocator, .{ .blob = blob }),
            };
        }

        fn bindsToBlob(self: *const Self) bool {
            return self.name_blob != null;
        }

        fn ensureCapacity(self: *Self, min_cap: u32) !void {
            if (self.ref_map) |*rm| {
                rm.ensureTotalCapacity(min_cap) catch return self.diag.oom(.dedup_table, null);
                return;
            }
            self.map.ensureTotalCapacity(min_cap) catch return self.diag.oom(.dedup_table, null);
        }

        /// Releases the map and its owned names.
        pub fn deinit(self: *Self) void {
            if (self.ref_map) |*rm| {
                rm.deinit();
                self.ref_map = null;
            }
            self.map.deinit();
            if (self.key_arena) |*arena| arena.deinit();
            self.key_arena = null;
            self.name_blob = null;
        }

        /// Returns true for an existing exact name; otherwise retains a stable copy.
        pub fn observe(self: *Self, name: []const u8) !bool {
            if (self.name_blob) |blob| {
                var rm = &(self.ref_map orelse unreachable);
                const lookup = SliceLookup{ .blob = blob };
                const gop = rm.getOrPutAdapted(name, lookup) catch return self.diag.oom(.dedup_table, null);
                if (gop.found_existing) {
                    self.last_offset = gop.key_ptr.offset;
                    self.last_len = @intCast(gop.key_ptr.len);
                    return true;
                }
                const off = std.math.cast(u32, blob.items.len) orelse {
                    _ = rm.removeByPtr(gop.key_ptr);
                    return self.diag.capacity(.name_offset, std.math.maxInt(u32), blob.items.len, error.OutOfMemory);
                };
                blob.appendSlice(rm.allocator, name) catch {
                    _ = rm.removeByPtr(gop.key_ptr);
                    return self.diag.oom(.names, name.len);
                };
                gop.key_ptr.* = .{ .offset = off, .len = @intCast(name.len) };
                self.last_offset = off;
                self.last_len = @intCast(name.len);
                return false;
            }
            const arena = &(self.key_arena orelse unreachable);
            const gop = self.map.getOrPut(name) catch return self.diag.oom(.dedup_table, null);
            if (gop.found_existing) return true;
            const owned = arena.allocator().dupe(u8, name) catch {
                _ = self.map.remove(name);
                return self.diag.oom(.names, name.len);
            };
            gop.key_ptr.* = owned;
            return false;
        }
    };
}

const LineMetrics = struct {
    seq_len: u64 = 0,
    line_bases: u32 = 0,
    line_bytes: u32 = 0,
    is_uniform_width: bool = true,
    line_count: u64 = 0,
    first_content_dense: bool = false,
};

const SequenceLengthInfo = struct {
    seq_len: u64,
    uses_uniform_formula: bool,
};

fn countBases(data: []const u8) u64 {
    var count: u64 = 0;
    var pos: usize = 0;
    while (pos + SIMD_CHUNK_SIZE <= data.len) {
        const chunk: SimdVec = data[pos..][0..SIMD_CHUNK_SIZE].*;
        const space_char: SimdVec = @splat(' ');
        count += std.simd.countTrues(chunk > space_char);
        pos += SIMD_CHUNK_SIZE;
    }
    while (pos < data.len) {
        if (data[pos] > ' ') count += 1;
        pos += 1;
    }
    return count;
}

fn hasLineSeparatorAt(data: []const u8, pos: usize, sep_len: u32) bool {
    return switch (sep_len) {
        1 => pos < data.len and data[pos] == '\n',
        2 => pos + 1 < data.len and data[pos] == '\r' and data[pos + 1] == '\n',
        else => false,
    };
}

fn lineSliceHasWhitespace(data: []const u8, start: usize, end: usize) bool {
    var pos = start;
    while (pos + SIMD_CHUNK_SIZE <= end) {
        const chunk: SimdVec = data[pos..][0..SIMD_CHUNK_SIZE].*;
        const space_char: SimdVec = @splat(' ');
        const mask = chunk <= space_char;
        if (@reduce(.Or, mask)) return true;
        pos += SIMD_CHUNK_SIZE;
    }
    while (pos < end) : (pos += 1) {
        if (data[pos] <= ' ') return true;
    }
    return false;
}

fn countBasesInLineSlice(data: []const u8, start: usize, end: usize) u32 {
    if (start >= end) return 0;
    if (!lineSliceHasWhitespace(data, start, end)) {
        return @intCast(end - start);
    }
    return @intCast(countBases(data[start..end]));
}

fn firstLineIsDense(line_bases: u32, line_bytes: u32) bool {
    if (line_bases == 0 or line_bytes <= line_bases) return false;
    const sep_len = line_bytes - line_bases;
    return (sep_len == 1 or sep_len == 2) and line_bases + sep_len == line_bytes;
}

fn countWhitespaceDirect(comptime vector_len: usize, data: []const u8) usize {
    const Vec = @Vector(vector_len, u8);
    var count: usize = 0;
    var pos: usize = 0;
    while (pos + vector_len <= data.len) : (pos += vector_len) {
        const chunk: Vec = data[pos..][0..vector_len].*;
        count += std.simd.countTrues(chunk <= @as(Vec, @splat(' ')));
    }
    while (pos < data.len) : (pos += 1) {
        count += @intFromBool(data[pos] <= ' ');
    }
    return count;
}

fn countWhitespace(comptime vector_len: usize, data: []const u8) usize {
    if (data.len < 2048) return countWhitespaceDirect(vector_len, data);
    const Vec = @Vector(vector_len, u8);
    const Wide = @Vector(vector_len, u16);
    const group = 128;
    var count: usize = 0;
    var pos: usize = 0;
    while (data.len - pos >= vector_len * group) {
        var lanes: Vec = @splat(0);
        for (0..group) |_| {
            const bytes: Vec = data[pos..][0..vector_len].*;
            lanes += @select(u8, bytes <= @as(Vec, @splat(' ')), @as(Vec, @splat(1)), @as(Vec, @splat(0)));
            pos += vector_len;
        }
        count += @reduce(.Add, @as(Wide, lanes));
    }
    var lanes: Vec = @splat(0);
    while (data.len - pos >= vector_len) : (pos += vector_len) {
        const bytes: Vec = data[pos..][0..vector_len].*;
        lanes += @select(u8, bytes <= @as(Vec, @splat(' ')), @as(Vec, @splat(1)), @as(Vec, @splat(0)));
    }
    count += @reduce(.Add, @as(Wide, lanes));
    while (pos < data.len) : (pos += 1) {
        count += @intFromBool(data[pos] <= ' ');
    }
    return count;
}

fn validatedStrideBlock(
    comptime vector_len: usize,
    data: []const u8,
    start: usize,
    line_count: usize,
    line_bases: usize,
    line_bytes: usize,
) usize {
    const sep_len = line_bytes - line_bases;
    var line_index: usize = 0;
    while (line_index < line_count) : (line_index += 1) {
        const line_start = start + line_index * line_bytes;
        if (data[line_start] == '>') break;
        if (!hasLineSeparatorAt(data, line_start + line_bases, @intCast(sep_len))) {
            break;
        }
    }
    const candidate_lines = line_index;
    if (candidate_lines == 0) return 0;
    const block_end = start + candidate_lines * line_bytes;
    if (countWhitespace(vector_len, data[start..block_end]) == candidate_lines * sep_len) {
        return candidate_lines;
    }

    line_index = 0;
    while (line_index < candidate_lines) : (line_index += 1) {
        const line_start = start + line_index * line_bytes;
        if (countWhitespace(vector_len, data[line_start .. line_start + line_bases]) != 0) {
            return line_index;
        }
    }
    return 0;
}

fn validatedStrideRun(
    comptime vector_len: usize,
    data: []const u8,
    start: usize,
    line_bases_u32: u32,
    line_bytes_u32: u32,
) usize {
    const line_bases: usize = line_bases_u32;
    const line_bytes: usize = line_bytes_u32;
    if (line_bases == 0 or line_bytes <= line_bases) return 0;
    const sep_len = line_bytes - line_bases;
    if (sep_len != 1 and sep_len != 2) return 0;

    var cursor = start;
    var accepted_lines: usize = 0;
    while (data.len - cursor >= line_bytes) {
        const available_lines = (data.len - cursor) / line_bytes;
        const block_lines = @min(available_lines, STRIDE_VALIDATION_BLOCK_LINES);
        const valid_lines = validatedStrideBlock(
            vector_len,
            data,
            cursor,
            block_lines,
            line_bases,
            line_bytes,
        );
        accepted_lines += valid_lines;
        if (valid_lines != block_lines) return accepted_lines;
        cursor += block_lines * line_bytes;
    }
    return accepted_lines;
}

fn findNameDelimiter(data: []const u8, start: usize) usize {
    const width = 16;
    const Vec = @Vector(width, u8);
    var position = start;
    while (data.len - position >= width) : (position += width) {
        const bytes: Vec = data[position..][0..width].*;
        const separators = (bytes == @as(Vec, @splat(' '))) |
            (bytes == @as(Vec, @splat('\t'))) |
            (bytes == @as(Vec, @splat('\r'))) |
            (bytes == @as(Vec, @splat('\n')));
        if (std.simd.firstTrue(separators)) |offset| return position + offset;
    }
    while (position < data.len) : (position += 1) {
        switch (data[position]) {
            ' ', '\t', '\r', '\n' => return position,
            else => {},
        }
    }
    return position;
}

fn findNextNewline(data: []const u8, start: usize, end: usize) usize {
    const relative = std.mem.findScalar(u8, data[start..end], '\n') orelse return end;
    return start + relative;
}

const LineMetricsBuilder = struct {
    metrics: LineMetrics = .{},
    nonuniform_reason: diagnostics.NonuniformReason = .nonuniform_layout,
    first_actual_bytes: u32 = 0,
    first_content_len: u32 = 0,
    pending_bases: u32 = 0,
    pending_actual_bytes: u32 = 0,
    pending_content_len: u32 = 0,
    have_pending: bool = false,
    blank_after_bases: bool = false,

    fn ingestLine(self: *LineMetricsBuilder, bases: u32, actual_bytes: u32, has_lf: bool, content_len: u32) void {
        if (bases == 0) {
            if (self.metrics.line_count > 0) {
                self.blank_after_bases = true;
            }
            return;
        }
        if (self.blank_after_bases) {
            self.metrics.is_uniform_width = false;
            self.nonuniform_reason = .interior_blank_line;
            self.blank_after_bases = false;
        }
        self.metrics.seq_len += bases;
        self.metrics.line_count += 1;
        if (bases != content_len) {
            if (self.metrics.is_uniform_width) self.nonuniform_reason = .sequence_whitespace;
            self.metrics.is_uniform_width = false;
        }

        if (self.metrics.line_count == 1) {
            self.first_actual_bytes = actual_bytes;
            self.first_content_len = content_len;
            self.metrics.line_bases = bases;
            self.metrics.line_bytes = if (has_lf) actual_bytes else actual_bytes + 1;
            const sep_len = if (has_lf) actual_bytes -| content_len else @as(u32, 1);
            self.metrics.first_content_dense = bases == content_len and (sep_len == 1 or sep_len == 2);
        } else if (self.have_pending) {
            if (self.pending_bases != self.metrics.line_bases or
                self.pending_actual_bytes != self.first_actual_bytes)
            {
                if (self.metrics.is_uniform_width) self.nonuniform_reason = if (self.pending_bases < self.metrics.line_bases) .interior_short_line else if (self.pending_bases > self.metrics.line_bases) .wider_line else .mixed_line_endings;
                self.metrics.is_uniform_width = false;
            }
        }

        self.pending_bases = bases;
        self.pending_actual_bytes = actual_bytes;
        self.pending_content_len = content_len;
        self.have_pending = true;
    }

    fn finish(self: *LineMetricsBuilder) LineMetrics {
        self.blank_after_bases = false;
        if (self.metrics.line_count > 1 and self.have_pending) {
            if (self.pending_bases > self.metrics.line_bases or
                self.pending_content_len > self.first_content_len)
            {
                if (self.metrics.is_uniform_width) self.nonuniform_reason = .wider_line;
                self.metrics.is_uniform_width = false;
            }
        }
        return self.metrics;
    }
};

fn sequenceLengthInfoFromMetrics(metrics: LineMetrics) SequenceLengthInfo {
    if (metrics.seq_len == 0) return .{ .seq_len = 0, .uses_uniform_formula = false };
    if (!metrics.is_uniform_width or metrics.line_bases == 0 or !metrics.first_content_dense) {
        return .{ .seq_len = metrics.seq_len, .uses_uniform_formula = false };
    }
    return .{
        .seq_len = metrics.seq_len,
        .uses_uniform_formula = firstLineIsDense(metrics.line_bases, metrics.line_bytes),
    };
}

fn finalizeSideTableOffsets(index: *ZfiIndex, diag: Diagnostics) !void {
    const side_table_base = @sizeOf(index_format.ZfiHeader) + index.records.items.len * @sizeOf(IndexRecord);
    for (index.records.items) |*rec| {
        if (!rec.isUniformWidth()) {
            const local_offset = rec.sideTableOffset();
            try diag.markNonUniform(rec, side_table_base + local_offset, .wire);
        }
    }
}

fn embedZfiName(
    allocator: std.mem.Allocator,
    rec: *IndexRecord,
    name: []const u8,
    name_blob: *std.ArrayList(u8),
    diag: Diagnostics,
) !void {
    rec.name_offset = name_blob.items.len;
    rec.name_len = @intCast(name.len);
    name_blob.appendSlice(allocator, name) catch return diag.oom(.names, name.len);
    rec.flags_and_offset[0] |= index_format.NAME_IN_ZFI_FLAG;
}

fn scanZfiReaderWithOptions(
    allocator: std.mem.Allocator,
    reader: *std.Io.Reader,
    read_buf: []u8,
    options: ScanOptions,
) !ZfiIndex {
    var index = ZfiIndex{
        .records = .empty,
        .side_tables = .empty,
        .name_blob = .empty,
    };
    errdefer index.deinit(allocator);
    const diag = Diagnostics{ .failure = options.failure };
    try checkReadBufferSize(read_buf.len, options.failure);
    index.records.ensureTotalCapacity(allocator, ZFI_INITIAL_RECORD_CAPACITY) catch return diag.oom(.records, ZFI_INITIAL_RECORD_CAPACITY * @sizeOf(IndexRecord));

    const Ctx = struct {
        index: *ZfiIndex,
        allocator: std.mem.Allocator,
        diag: Diagnostics,

        fn emit(ctx: *@This(), emit_info: FastaRecordEmit) !void {
            var rec = emit_info.record;
            if (!emit_info.uses_uniform_formula) {
                const local_offset = ctx.index.side_tables.items.len;
                ctx.index.side_tables.appendSlice(ctx.allocator, emit_info.side_table) catch return ctx.diag.forRecord(emit_info).oom(.side_tables, emit_info.side_table.len);
                try ctx.diag.forRecord(emit_info).markNonUniform(&rec, local_offset, .local);
            }
            if (!emit_info.name_embedded) {
                try embedZfiName(ctx.allocator, &rec, emit_info.name, &ctx.index.name_blob, ctx.diag.forRecord(emit_info));
            }
            ctx.index.records.append(ctx.allocator, rec) catch return ctx.diag.forRecord(emit_info).oom(.records, @sizeOf(IndexRecord));
        }
    };

    var ctx = Ctx{ .index = &index, .allocator = allocator, .diag = diag };
    var blob_dedup: ?NameDedup = null;
    defer if (blob_dedup) |*seen| seen.deinit();
    if (options.enable_dedup) {
        blob_dedup = NameDedup.initOwningBlob(allocator, &index.name_blob);
        blob_dedup.?.diag = diag;
        try blob_dedup.?.ensureCapacity(DEDUP_INITIAL_CAPACITY);
    }
    const external_dedup: ?*NameDedup = if (blob_dedup) |*s| s else null;
    _ = try scanFastaReader(
        allocator,
        reader,
        read_buf,
        .{
            .enable_dedup = options.enable_dedup,
            .collect_side_tables = true,
            .require_initial_header = options.require_initial_header,
            .failure = options.failure,
        },
        external_dedup,
        &ctx,
        Ctx.emit,
    );
    try finalizeSideTableOffsets(&index, diag);
    return index;
}

fn writeZfiIndex(
    writer: *std.Io.Writer,
    index: *const ZfiIndex,
    source_size: u64,
    source_mtime_ns: u64,
) !void {
    const header = ZfiHeader{
        .magic = ZFI_MAGIC,
        .record_count = @intCast(index.records.items.len),
        .source_size = source_size,
    };
    try writer.writeAll(std.mem.asBytes(&header));
    try writer.writeAll(std.mem.sliceAsBytes(index.records.items));
    try writer.writeAll(index.side_tables.items);
    try writer.writeAll(index.name_blob.items);
    try writer.writeAll(&index_format.encodeZfiSourceId(source_mtime_ns));
    try writer.writeAll(&index_format.encodeZfiNameFooter(index.name_blob.items.len));
}

fn serializeZfiIndex(allocator: std.mem.Allocator, index: *const ZfiIndex, source_size: u64, source_mtime_ns: u64, diag: Diagnostics) ![]u8 {
    _ = try checkedZfiRecordCount(index.records.items.len, diag);
    var aw: std.Io.Writer.Allocating = .init(allocator);
    errdefer aw.deinit();
    writeZfiIndex(&aw.writer, index, source_size, source_mtime_ns) catch return diag.oom(.serialized_index, null);
    return aw.toOwnedSlice() catch return diag.oom(.serialized_index, null);
}

fn checkedZfiRecordCount(count: usize, diag: Diagnostics) !u32 {
    return std.math.cast(u32, count) orelse return diag.fail(.{ .@"index.record_count_limit_exceeded" = .{ .record_count = .{ .exact = count }, .max_records = std.math.maxInt(u32), .format = .zfi, .limit_role = .zfi_header } }, .write, .zfi, null, error.TooManySequences);
}

const ChunkParseState = struct {
    diag: Diagnostics = .{},
    allocator: std.mem.Allocator,
    collect_side_tables: bool = false,
    active: bool = false,
    in_header: bool = false,
    parsing_name: bool = false,
    at_line_start: bool = true,
    header_start_offset: u64 = 0,
    source_record: u64 = 0,
    source_reference: diagnostics.Reference = .{ .source = .fasta },
    name: std.ArrayList(u8) = .empty,
    seq_offset: u64 = 0,
    seq_offset_set: bool = false,
    line_builder: LineMetricsBuilder = .{},
    in_pending_line: bool = false,
    pending_line_bases: u32 = 0,
    pending_line_bytes: u32 = 0,
    pending_pre_lf_bytes: u32 = 0,
    pending_pre_lf_ends_with_cr: bool = false,
    pending_line_file_offset: u64 = 0,
    side_table: std.ArrayList(u8) = .empty,
    side_base_start: u64 = 0,
    side_table_active: bool = false,
    deferred_side_line: ?struct {
        bases: u32,
        actual_bytes: u32,
        line_file_offset: u64,
    } = null,
    deferred_short_tail: ?struct {
        bases: u32,
        actual_bytes: u32,
        line_file_offset: u64,
    } = null,

    fn recordDiagnostics(self: *const ChunkParseState) Diagnostics {
        var diag = self.diag;
        diag.at.location.offset = self.header_start_offset;
        diag.at.record = if (self.source_record == 0) null else self.source_record;
        return diag;
    }

    fn ensureSideTableHeader(self: *ChunkParseState) !void {
        if (self.side_table.items.len >= @sizeOf(u64)) return;
        var line_count: u64 = 0;
        self.side_table.appendSlice(self.allocator, std.mem.asBytes(&line_count)) catch return self.recordDiagnostics().oom(.side_tables, @sizeOf(u64));
    }

    fn activateSideTable(self: *ChunkParseState, lines_to_materialize: u64) !void {
        if (self.side_table_active) return;
        self.side_table_active = true;
        try self.ensureSideTableHeader();
        if (lines_to_materialize == 0) {
            self.deferred_side_line = null;
            return;
        }
        const first = self.deferred_side_line orelse return self.recordDiagnostics().consistency(.missing_side_table, error.MissingSideTable);
        var i: u64 = 0;
        while (i < lines_to_materialize) : (i += 1) {
            const off = first.line_file_offset + i * @as(u64, first.actual_bytes);
            try self.appendSideTableLine(first.bases, first.actual_bytes, off);
        }
        self.deferred_side_line = null;
    }

    fn startHeader(self: *ChunkParseState, header_start_offset: u64) !void {
        self.active = true;
        self.in_header = true;
        self.parsing_name = true;
        self.at_line_start = false;
        self.header_start_offset = header_start_offset;
        self.source_record += 1;
        self.source_reference.location.offset = header_start_offset;
        self.source_reference.record = self.source_record;
        self.name.clearRetainingCapacity();
        self.seq_offset = 0;
        self.seq_offset_set = false;
        self.line_builder = .{};
        self.in_pending_line = false;
        self.pending_line_bases = 0;
        self.pending_line_bytes = 0;
        self.pending_pre_lf_bytes = 0;
        self.pending_pre_lf_ends_with_cr = false;
        self.pending_line_file_offset = 0;
        self.side_base_start = 0;
        self.side_table.clearRetainingCapacity();
        self.side_table_active = false;
        self.deferred_side_line = null;
        self.deferred_short_tail = null;
    }

    fn appendSideTableLine(
        self: *ChunkParseState,
        bases: u32,
        actual_bytes: u32,
        line_file_offset: u64,
    ) !void {
        const entry = index_format.SideTableLine{
            .base_start = self.side_base_start,
            .byte_offset = line_file_offset,
            .line_bytes = actual_bytes,
            .line_bases = bases,
        };
        self.side_table.appendSlice(self.allocator, std.mem.asBytes(&entry)) catch return self.recordDiagnostics().oom(.side_tables, @sizeOf(index_format.SideTableLine));
        self.side_base_start += bases;
    }

    fn flushDeferredShortTailThen(
        self: *ChunkParseState,
        metrics_line_count: u64,
        bases: u32,
        actual_bytes: u32,
        line_file_offset: u64,
    ) !void {
        const short = self.deferred_short_tail orelse return self.recordDiagnostics().consistency(.missing_side_table, error.MissingSideTable);
        self.deferred_short_tail = null;
        if (metrics_line_count < 2) return self.recordDiagnostics().consistency(.missing_side_table, error.MissingSideTable);
        try self.activateSideTable(metrics_line_count - 2);
        try self.appendSideTableLine(short.bases, short.actual_bytes, short.line_file_offset);
        try self.appendSideTableLine(bases, actual_bytes, line_file_offset);
    }

    fn recordSequenceLine(
        self: *ChunkParseState,
        bases: u32,
        actual_bytes: u32,
        has_lf: bool,
        content_len: u32,
        line_file_offset: u64,
    ) !void {
        self.line_builder.ingestLine(bases, actual_bytes, has_lf, content_len);

        if (!self.collect_side_tables) return;

        if (bases == 0) return;

        const metrics = self.line_builder.metrics;
        const matches_stride = metrics.line_bases > 0 and
            bases == metrics.line_bases and
            actual_bytes == self.line_builder.first_actual_bytes;
        const dense = metrics.first_content_dense and firstLineIsDense(metrics.line_bases, metrics.line_bytes);
        const can_use_formula = metrics.is_uniform_width and matches_stride and dense;

        if (can_use_formula) {
            if (metrics.line_count == 1) {
                self.deferred_side_line = .{
                    .bases = bases,
                    .actual_bytes = actual_bytes,
                    .line_file_offset = line_file_offset,
                };
            }
            return;
        }

        if (self.deferred_short_tail != null) {
            try self.flushDeferredShortTailThen(metrics.line_count, bases, actual_bytes, line_file_offset);
            return;
        }

        const short_ok_if_final = metrics.is_uniform_width and
            metrics.line_bases > 0 and
            dense and
            bases <= metrics.line_bases and
            content_len <= self.line_builder.first_content_len and
            actual_bytes <= self.line_builder.first_actual_bytes + 1;
        if (short_ok_if_final) {
            self.deferred_short_tail = .{
                .bases = bases,
                .actual_bytes = actual_bytes,
                .line_file_offset = line_file_offset,
            };
            return;
        }

        try self.activateSideTable(metrics.line_count - 1);
        try self.appendSideTableLine(bases, actual_bytes, line_file_offset);
    }

    fn noteSeqOffset(self: *ChunkParseState, line_file_offset: u64) void {
        if (self.seq_offset_set) return;
        self.seq_offset = line_file_offset;
        self.seq_offset_set = true;
    }

    fn strideEligible(self: *const ChunkParseState) bool {
        if (!self.active or self.in_header or self.in_pending_line) return false;
        if (self.side_table_active or self.deferred_short_tail != null) return false;
        if (self.line_builder.blank_after_bases) return false;
        const m = self.line_builder.metrics;
        if (m.line_count == 0 or !m.is_uniform_width) return false;
        if (!m.first_content_dense or !firstLineIsDense(m.line_bases, m.line_bytes)) return false;
        if (self.line_builder.first_actual_bytes != m.line_bytes) return false;
        if (self.line_builder.have_pending and
            (self.line_builder.pending_bases != m.line_bases or
                self.line_builder.pending_actual_bytes != self.line_builder.first_actual_bytes))
        {
            return false;
        }
        return true;
    }

    fn ingestStrideRun(self: *ChunkParseState, n_lines: u64) void {
        std.debug.assert(n_lines >= 1);
        std.debug.assert(self.line_builder.metrics.line_count >= 1);
        const bases = self.line_builder.metrics.line_bases;
        const line_bytes = self.line_builder.metrics.line_bytes;
        self.line_builder.metrics.seq_len += n_lines * @as(u64, bases);
        self.line_builder.metrics.line_count += n_lines;
        self.line_builder.pending_bases = bases;
        self.line_builder.pending_actual_bytes = line_bytes;
        self.line_builder.pending_content_len = self.line_builder.first_content_len;
        self.line_builder.have_pending = true;
        self.at_line_start = true;
    }

    fn ingestSequenceLine(
        self: *ChunkParseState,
        data: []const u8,
        line_start: usize,
        line_end: usize,
        has_lf: bool,
        line_file_offset: u64,
    ) !void {
        var content_end = line_end;
        if (content_end > line_start and data[content_end - 1] == '\r') {
            content_end -= 1;
        }

        const actual_bytes: u32 = @intCast((if (has_lf) line_end + 1 else line_end) - line_start);
        const content_len: u32 = @intCast(content_end - line_start);
        const bases = countBasesInLineSlice(data, line_start, content_end);
        if (bases > 0) self.noteSeqOffset(line_file_offset);

        try self.recordSequenceLine(bases, actual_bytes, has_lf, content_len, line_file_offset);
        self.at_line_start = has_lf;
        self.in_pending_line = false;
        self.pending_line_bases = 0;
        self.pending_line_bytes = 0;
        self.pending_pre_lf_bytes = 0;
        self.pending_pre_lf_ends_with_cr = false;
    }

    fn appendPendingFragment(self: *ChunkParseState, fragment: []const u8, fragment_file_offset: u64) !void {
        if (!self.in_pending_line) {
            self.in_pending_line = true;
            self.pending_line_bases = 0;
            self.pending_line_bytes = 0;
            self.pending_pre_lf_bytes = 0;
            self.pending_pre_lf_ends_with_cr = false;
            self.pending_line_file_offset = fragment_file_offset;
        }
        const fragment_bases = countBasesInLineSlice(fragment, 0, fragment.len);
        var diag = self.recordDiagnostics();
        diag.at.location.offset = self.pending_line_file_offset;
        const next_bases = std.math.add(u32, self.pending_line_bases, fragment_bases) catch
            return diag.lineLimit(.line_bases, @as(u64, self.pending_line_bases) + fragment_bases, .index_geometry);
        const next_bytes = std.math.add(u32, self.pending_line_bytes, @intCast(fragment.len)) catch
            return diag.lineLimit(.line_bytes, @as(u64, self.pending_line_bytes) + fragment.len, .index_geometry);
        self.pending_line_bases = next_bases;
        self.pending_line_bytes = next_bytes;
        self.pending_pre_lf_bytes = next_bytes;
        if (fragment.len > 0) {
            self.pending_pre_lf_ends_with_cr = fragment[fragment.len - 1] == '\r';
        }
    }

    fn commitPendingLine(self: *ChunkParseState, has_lf: bool) !void {
        if (!self.in_pending_line and self.pending_line_bytes == 0) return;
        var diag = self.recordDiagnostics();
        diag.at.location.offset = self.pending_line_file_offset;
        if (!has_lf and self.pending_line_bytes == std.math.maxInt(u32)) {
            return diag.lineLimit(.line_bytes, @as(u64, self.pending_line_bytes) + 1, .synthetic_separator);
        }
        if (self.pending_line_bases > 0) self.noteSeqOffset(self.pending_line_file_offset);
        var content_len = if (has_lf) self.pending_pre_lf_bytes else self.pending_line_bytes;
        if (content_len > 0 and self.pending_pre_lf_ends_with_cr) {
            content_len -= 1;
        }
        try self.recordSequenceLine(
            self.pending_line_bases,
            self.pending_line_bytes,
            has_lf,
            content_len,
            self.pending_line_file_offset,
        );
        self.at_line_start = has_lf;
        self.in_pending_line = false;
        self.pending_line_bases = 0;
        self.pending_line_bytes = 0;
        self.pending_pre_lf_bytes = 0;
        self.pending_pre_lf_ends_with_cr = false;
    }

    fn finalizeHeader(self: *ChunkParseState) void {
        self.in_header = false;
        self.at_line_start = true;
    }

    fn emitRecordIfReady(
        self: *ChunkParseState,
        ctx: anytype,
        seen_names: anytype,
        comptime emitRecord: fn (@TypeOf(ctx), FastaRecordEmit) anyerror!void,
        record_count: *u32,
    ) !void {
        try self.commitPendingLine(false);
        const metrics = self.line_builder.finish();
        const seq_info = sequenceLengthInfoFromMetrics(metrics);
        if (seq_info.seq_len == 0) {
            self.deferred_short_tail = null;
            return;
        }

        var name_embedded = false;
        var name_offset: u64 = self.header_start_offset + 1;
        var out_name_len: u16 = @intCast(self.name.items.len);
        var name_pad: [6]u8 = .{0} ** 6;
        if (seen_names) |seen| {
            if (try seen.observeAt(self.name.items, self.header_start_offset + 1)) {
                self.deferred_short_tail = null;
                return;
            }
            if (comptime @hasField(@TypeOf(seen.*), "last_offset")) {
                if (seen.bindsToBlob()) {
                    name_offset = seen.last_offset;
                    out_name_len = seen.last_len;
                    name_pad[0] = index_format.NAME_IN_ZFI_FLAG;
                    name_embedded = true;
                }
            }
        }

        if (record_count.* == std.math.maxInt(u32)) {
            @branchHint(.cold);
            return self.recordDiagnostics().recordLimit(@as(u64, record_count.*) + 1);
        }

        var side_table_slice: []const u8 = &.{};
        if (seq_info.uses_uniform_formula) {
            self.deferred_short_tail = null;
        } else if (self.collect_side_tables) {
            if (self.deferred_short_tail) |short| {
                if (!self.side_table_active) {
                    if (metrics.line_count < 1) return self.recordDiagnostics().consistency(.missing_side_table, error.MissingSideTable);
                    try self.activateSideTable(metrics.line_count - 1);
                    try self.appendSideTableLine(short.bases, short.actual_bytes, short.line_file_offset);
                }
                self.deferred_short_tail = null;
            } else if (!self.side_table_active) {
                try self.activateSideTable(metrics.line_count);
            }
            if (!self.side_table_active) {
                return self.recordDiagnostics().consistency(.missing_side_table, error.MissingSideTable);
            }
            @memcpy(self.side_table.items[0..@sizeOf(u64)], std.mem.asBytes(&metrics.line_count));
            side_table_slice = self.side_table.items;
        }
        try emitRecord(ctx, .{
            .record = .{
                .name_offset = name_offset,
                .name_len = out_name_len,
                .seq_offset = self.seq_offset,
                .seq_len = seq_info.seq_len,
                .line_bases = metrics.line_bases,
                .line_bytes = metrics.line_bytes,
                .flags_and_offset = name_pad,
            },
            .name = self.name.items,
            .source_reference = &self.source_reference,
            .nonuniform_reason = if (!metrics.first_content_dense) .sequence_whitespace else self.line_builder.nonuniform_reason,
            .uses_uniform_formula = seq_info.uses_uniform_formula,
            .side_table = side_table_slice,
            .name_embedded = name_embedded,
        });
        record_count.* += 1;
    }
};

fn processChunkBytes(
    state: *ChunkParseState,
    data: []const u8,
    file_offset: u64,
    max_name_len: ?usize,
    seen_names: anytype,
    ctx: anytype,
    comptime emitRecord: fn (@TypeOf(ctx), FastaRecordEmit) anyerror!void,
    record_count: *u32,
) !void {
    var i: usize = 0;

    if (state.in_pending_line) {
        const nl = findNextNewline(data, 0, data.len);
        if (nl >= data.len) {
            try state.appendPendingFragment(data, file_offset);
            return;
        }

        try state.appendPendingFragment(data[0..nl], file_offset);
        state.pending_line_bytes = std.math.add(u32, state.pending_line_bytes, 1) catch {
            var diag = state.diag;
            diag.at.location.offset = state.pending_line_file_offset;
            return diag.lineLimit(.line_bytes, @as(u64, state.pending_line_bytes) + 1, .index_geometry);
        };
        try state.commitPendingLine(true);
        i = nl + 1;
    }

    while (i < data.len) {
        const byte_offset = file_offset + @as(u64, @intCast(i));
        const byte = data[i];

        if (state.at_line_start and byte == '>') {
            if (state.active) {
                try state.emitRecordIfReady(ctx, seen_names, emitRecord, record_count);
            }
            try state.startHeader(byte_offset);
            i += 1;
            continue;
        }

        if (state.in_header) {
            if (state.parsing_name) {
                const name_end = findNameDelimiter(data, i);

                const fragment = data[i..name_end];
                if (state.name.items.len > MAX_INDEX_NAME_LEN or
                    fragment.len > MAX_INDEX_NAME_LEN - state.name.items.len)
                {
                    return state.recordDiagnostics().fail(.{ .@"index.identifier_too_long" = .{ .observed_bytes = .{ .at_least = state.name.items.len + fragment.len }, .limit_bytes = MAX_INDEX_NAME_LEN, .format = state.diag.format } }, .scan, .fasta, null, error.HeaderTooLong);
                }
                if (max_name_len) |limit| {
                    if (state.name.items.len > limit or fragment.len > limit - state.name.items.len) {
                        return state.recordDiagnostics().fail(.{ .@"index.identifier_too_long" = .{ .observed_bytes = .{ .at_least = state.name.items.len + fragment.len }, .limit_bytes = limit, .format = state.diag.format } }, .scan, .fasta, null, error.HeaderTooLong);
                    }
                }
                state.name.appendSlice(state.allocator, fragment) catch return state.diag.oom(.names, fragment.len);
                i = name_end;
                if (i >= data.len) return;
                if (data[i] == '\n') {
                    state.finalizeHeader();
                    i += 1;
                    continue;
                }
                state.parsing_name = false;
            }

            const header_end = findNextNewline(data, i, data.len);
            if (header_end >= data.len) return;
            state.finalizeHeader();
            i = header_end + 1;
            continue;
        }

        if (!state.active) {
            state.at_line_start = byte == '\n';
            i += 1;
            continue;
        }

        if (state.at_line_start and state.strideEligible()) {
            const line_bases = state.line_builder.metrics.line_bases;
            const line_bytes = state.line_builder.metrics.line_bytes;
            const run_lines = validatedStrideRun(
                STRIDE_VECTOR_LEN,
                data,
                i,
                line_bases,
                line_bytes,
            );
            if (run_lines > 0) {
                state.ingestStrideRun(run_lines);
                i += run_lines * @as(usize, line_bytes);
            }
            if (i >= data.len) return;
            if (data[i] == '>') continue;
            if (data.len - i < line_bytes) {
                if (findNextNewline(data, i, data.len) >= data.len) {
                    try state.appendPendingFragment(
                        data[i..],
                        file_offset + @as(u64, @intCast(i)),
                    );
                    return;
                }
            }
        }

        const line_start = i;
        const nl = findNextNewline(data, i, data.len);
        if (nl >= data.len) {
            try state.appendPendingFragment(
                data[line_start..],
                file_offset + @as(u64, @intCast(line_start)),
            );
            return;
        }

        try state.ingestSequenceLine(
            data,
            line_start,
            nl,
            true,
            file_offset + @as(u64, @intCast(line_start)),
        );
        i = nl + 1;
    }
}

fn scanFastaReaderWithDedup(
    allocator: std.mem.Allocator,
    reader: *std.Io.Reader,
    read_buf: []u8,
    options: ScanOptions,
    external_dedup: anytype,
    ctx: anytype,
    comptime emitRecord: fn (@TypeOf(ctx), FastaRecordEmit) anyerror!void,
) !u32 {
    var owned_original: ?NameDedup = null;
    defer if (owned_original) |*seen| seen.deinit();
    const Dedup = if (@TypeOf(external_dedup) == ?*SourceAdapter) SourceAdapter else if (@TypeOf(external_dedup) == ?*PackedNames) PackedNames else NameDedup;
    const external: ?*Dedup = external_dedup;
    const dedup_ptr: ?*Dedup = blk: {
        if (external) |pointer| break :blk pointer;
        if (comptime Dedup == NameDedup) {
            if (options.enable_dedup) {
                owned_original = NameDedup.init(allocator);
                owned_original.?.diag = .{ .failure = options.failure };
                try owned_original.?.ensureCapacity(DEDUP_INITIAL_CAPACITY);
                break :blk &owned_original.?;
            }
        } else {
            std.debug.assert(!options.enable_dedup);
        }
        break :blk null;
    };

    var state = ChunkParseState{
        .allocator = allocator,
        .collect_side_tables = options.collect_side_tables,
        .diag = .{ .failure = options.failure, .format = if (options.collect_side_tables) .zfi else .fai },
    };
    defer state.name.deinit(allocator);
    defer if (options.collect_side_tables) state.side_table.deinit(allocator);
    var record_count: u32 = 0;
    var file_offset: u64 = 0;

    while (true) {
        var read_diag = state.diag;
        read_diag.at = .{ .source = .fasta, .location = .{ .offset = file_offset } };
        const n = reader.readSliceShort(read_buf) catch |err| return read_diag.fail(.{ .@"io.read_failed" = .{ .requested_bytes = read_buf.len } }, .read, .fasta, err, error.SourceReadFailed);
        if (n == 0) break;
        if (file_offset == 0 and options.require_initial_header and read_buf[0] != '>') {
            return read_diag.fail(.{ .@"index.initial_header_required" = .{ .byte = read_buf[0] } }, .scan, .fasta, null, error.NotFasta);
        }
        const next_offset = std.math.add(u64, file_offset, n) catch return read_diag.fail(.{ .@"resource.count_overflow" = .{ .counter = .source_bytes, .max_value = std.math.maxInt(u64), .current = file_offset, .increment = n } }, .scan, .fasta, null, error.SourceTooLong);
        if (comptime Dedup == SourceAdapter) {
            if (dedup_ptr) |seen| seen.source_names.setChunk(read_buf[0..n], file_offset);
        }
        try processChunkBytes(
            &state,
            read_buf[0..n],
            file_offset,
            options.max_name_len,
            dedup_ptr,
            ctx,
            emitRecord,
            &record_count,
        );
        file_offset = next_offset;
    }

    if (state.active) {
        try state.emitRecordIfReady(ctx, dedup_ptr, emitRecord, &record_count);
    }

    return record_count;
}

fn checkReadBufferSize(len: usize, failure: ?*diagnostics.Failure) !void {
    if (len == 0 or len > std.math.maxInt(u32)) {
        @branchHint(.cold);
        const diag = Diagnostics{ .failure = failure };
        if (len == 0) return diag.fail(.{ .@"usage.invalid_option_value" = .{ .option = .read_buffer, .expected = .nonempty_buffer } }, .scan, .fasta, null, error.InvalidArgument);
        return diag.capacity(.reader_buffer, std.math.maxInt(u32), len, error.CapacityExceeded);
    }
}

const FaiEmitBuffer = struct {
    diag: Diagnostics = .{ .format = .fai },
    writer: *std.Io.Writer,
    buf: [INDEX_OUTPUT_BUFFER_SIZE]u8 = undefined,
    len: usize = 0,

    fn flush(self: *FaiEmitBuffer) !void {
        if (self.len == 0) return;
        self.writer.writeAll(self.buf[0..self.len]) catch |err| return self.diag.fail(.{ .@"io.write_failed" = {} }, .write, .fai_spool, err, error.FaiSpoolWriteFailed);
        self.len = 0;
    }

    fn emitRecord(self: *FaiEmitBuffer, emit_info: FastaRecordEmit) !void {
        if (!emit_info.uses_uniform_formula) {
            const scan_diag = self.diag.forRecord(emit_info);
            return scan_diag.fail(.{ .@"index.nonuniform_fai" = .{ .reason = emit_info.nonuniform_reason, .observed = .{ .line_bases = emit_info.record.line_bases, .line_bytes = emit_info.record.line_bytes } } }, .scan, .fasta, null, error.NonUniformFai);
        }
        const rec = emit_info.record;
        var suffix_buf: [FAI_SUFFIX_BUFFER_SIZE]u8 = undefined;
        const suffix = formatFaiSuffix(&suffix_buf, rec);
        const total = emit_info.name.len + suffix.len;
        if (total > self.buf.len) {
            try self.flush();
            self.writer.writeAll(emit_info.name) catch |err| return self.diag.fail(.{ .@"io.write_failed" = {} }, .write, .fai_spool, err, error.FaiSpoolWriteFailed);
            self.writer.writeAll(suffix) catch |err| return self.diag.fail(.{ .@"io.write_failed" = {} }, .write, .fai_spool, err, error.FaiSpoolWriteFailed);
            return;
        }
        if (self.len + total > self.buf.len) try self.flush();
        @memcpy(self.buf[self.len..][0..emit_info.name.len], emit_info.name);
        self.len += emit_info.name.len;
        @memcpy(self.buf[self.len..][0..suffix.len], suffix);
        self.len += suffix.len;
    }
};

fn appendUnsignedDecimal(buffer: []u8, value: u64) usize {
    var reversed: [FAI_U64_DECIMAL_DIGITS]u8 = undefined;
    var remaining = value;
    var digit_count: usize = 0;
    while (true) {
        reversed[digit_count] = @intCast('0' + remaining % 10);
        digit_count += 1;
        remaining /= 10;
        if (remaining == 0) break;
    }

    for (0..digit_count) |index| {
        buffer[index] = reversed[digit_count - index - 1];
    }
    return digit_count;
}

fn formatFaiSuffix(buffer: *[FAI_SUFFIX_BUFFER_SIZE]u8, rec: IndexRecord) []const u8 {
    const values = [_]u64{ rec.seq_len, rec.seq_offset, rec.line_bases, rec.line_bytes };
    var len: usize = 0;
    for (values) |value| {
        buffer[len] = '\t';
        len += 1;
        len += appendUnsignedDecimal(buffer[len..], value);
    }
    buffer[len] = '\n';
    return buffer[0 .. len + 1];
}

const FaiSpool = struct {
    diag: Diagnostics = .{},
    role: diagnostics.SourceId = .fai_spool,
    dir: std.Io.Dir,
    close_dir: bool,
    file: std.Io.File,
    name: [FAI_SPOOL_NAME_BUFFER_SIZE]u8,
    name_len: usize,

    fn deinit(self: *FaiSpool, io: std.Io) void {
        self.file.close(io);
        self.dir.deleteFile(io, self.name[0..self.name_len]) catch |err| {
            if (err != error.FileNotFound) if (self.diag.failure) |failure| {
                switch (failure.fail(.{ .details = .{ .@"io.cleanup_failed" = .{ .resource_role = self.role } }, .stage = .cleanup, .at = .{ .source = self.role }, .cause = .{ .zig = err } })) {
                    else => {},
                }
            };
        };
        if (self.close_dir) self.dir.close(io);
    }
};

fn registerTemporaryPath(failure: *diagnostics.Failure, role: diagnostics.SourceId, directory: []const u8, name: []const u8) void {
    if (failure.sources.get(role) == null) {
        var prefix: [64]u8 = undefined;
        var retained: usize = 0;
        for ([_][]const u8{ directory, "/", name }) |part| {
            const n = @min(part.len, prefix.len - retained);
            @memcpy(prefix[retained..][0..n], part[0..n]);
            retained += n;
        }
        failure.sources.set(role, .{ .label = .{ .retained = failure.retainPrefix(prefix[0..retained], .{ .exact = directory.len + 1 + name.len }, 64) }, .kind = .path });
    }
}

fn tryCreateFaiSpool(io: std.Io, dir: std.Io.Dir, close_dir: bool) ?FaiSpool {
    var attempts: diagnostics.TempAttempts = .{};
    return tryCreateSpool(io, dir, close_dir, .{}, .fai_spool, .current_directory, ".", &attempts);
}

fn tryCreateSpool(io: std.Io, dir: std.Io.Dir, close_dir: bool, diag: Diagnostics, role: diagnostics.SourceId, place: diagnostics.TemporaryPlace, directory: []const u8, attempts: *diagnostics.TempAttempts) ?FaiSpool {
    for (0..FAI_SPOOL_CREATE_ATTEMPTS) |_| {
        var nonce: u128 = undefined;
        std.Io.randomSecure(io, std.mem.asBytes(&nonce)) catch |err| {
            attempts.append(.{ .place = place, .stage = .create, .reason = .randomness, .cause = .{ .zig = err } });
            return null;
        };

        var name_buf: [FAI_SPOOL_NAME_BUFFER_SIZE]u8 = undefined;
        const name = std.fmt.bufPrint(&name_buf, "z-fasta-{x}.fai.tmp", .{nonce}) catch unreachable;
        const file = dir.createFile(io, name, .{ .read = true, .exclusive = true }) catch |err| switch (err) {
            error.PathAlreadyExists => continue,
            else => {
                attempts.append(.{ .place = place, .stage = .create, .reason = .creation, .cause = .{ .zig = err } });
                return null;
            },
        };

        var spool = FaiSpool{
            .diag = diag,
            .role = role,
            .dir = dir,
            .close_dir = close_dir,
            .file = file,
            .name = undefined,
            .name_len = name.len,
        };
        @memcpy(spool.name[0..name.len], name);
        if (diag.failure) |failure| {
            registerTemporaryPath(failure, role, directory, name);
        }
        return spool;
    }
    attempts.append(.{ .place = place, .stage = .create, .reason = .collisions, .cause = .{ .zig = error.PathAlreadyExists } });
    return null;
}

fn tryOpenSpoolDir(io: std.Io, path: []const u8, diag: Diagnostics, role: diagnostics.SourceId, place: diagnostics.TemporaryPlace, attempts: *diagnostics.TempAttempts) ?FaiSpool {
    const dir = std.Io.Dir.cwd().openDir(io, path, .{}) catch |err| {
        attempts.append(.{ .place = place, .stage = .open, .reason = .directory_open, .cause = .{ .zig = err } });
        return null;
    };
    return tryCreateSpool(io, dir, true, diag, role, place, path, attempts) orelse {
        dir.close(io);
        return null;
    };
}

fn openFaiSpool(allocator: std.mem.Allocator, io: std.Io, environ: std.process.Environ, diag: Diagnostics, role: diagnostics.SourceId) !FaiSpool {
    var attempts: diagnostics.TempAttempts = .{};
    const env_names = [_][]const u8{ "TMPDIR", "TEMP", "TMP" };
    if (comptime @hasDecl(std.process.Environ.Block, "view")) {
        for (env_names) |env_name| {
            const path = std.process.Environ.getPosix(environ, env_name) orelse continue;
            if (path.len == 0) continue;
            if (tryOpenSpoolDir(io, path, diag, role, .environment, &attempts)) |spool| return spool;
        }
    } else {
        var env_map = environ.createMap(allocator) catch return diag.oom(.arena, null);
        defer env_map.deinit();
        for (env_names) |env_name| {
            const path = env_map.get(env_name) orelse continue;
            if (path.len == 0) continue;
            if (tryOpenSpoolDir(io, path, diag, role, .environment, &attempts)) |spool| return spool;
        }
    }
    if (builtin.os.tag != .windows) {
        if (tryOpenSpoolDir(io, "/tmp", diag, role, .system, &attempts)) |spool| return spool;
    }
    if (tryCreateSpool(io, .cwd(), false, diag, role, .current_directory, ".", &attempts)) |spool| return spool;
    return diag.fail(.{ .@"io.temp_unavailable" = .{ .purpose = switch (role) {
        .name_spool => .names,
        .side_spool => .side_tables,
        .catalog_spool => .catalog,
        else => .fai,
    }, .attempts = attempts } }, .create, role, null, error.NoUsableFaiSpool);
}

fn cleanupAtomic(file: *std.Io.File.Atomic, io: std.Io, failure: ?*diagnostics.Failure) void {
    if (file.file_open) {
        file.file.close(io);
        file.file_open = false;
    }
    if (file.file_exists) {
        var removed = true;
        const basename = std.fmt.hex(file.file_basename_hex);
        file.dir.deleteFile(io, &basename) catch |err| {
            removed = err == error.FileNotFound;
            if (err != error.FileNotFound) if (failure) |destination| {
                registerTemporaryPath(destination, .staged_output, std.fs.path.dirname(destination.sourceLabel(.zfi) orelse "") orelse ".", &basename);
                switch (destination.fail(.{ .details = .{ .@"io.cleanup_failed" = .{ .resource_role = .staged_output } }, .stage = .cleanup, .at = .{ .source = .staged_output }, .cause = .{ .zig = err } })) {
                    else => {},
                }
            };
        };
        if (failure) |destination| {
            if (destination.outputs.file.state == .staged) destination.outputs.file.state = if (removed) .aborted else .unknown;
        }
        file.file_exists = false;
    }
    file.deinit(io);
}

fn checkSourceUnchanged(io: std.Io, dir: std.Io.Dir, path: []const u8, before: std.Io.File.Stat, diag: Diagnostics) !void {
    const after = dir.statFile(io, path, .{}) catch |err| {
        if (err == error.FileNotFound) return diag.fail(.{ .@"input.changed" = .{ .changed_fields = .{ .path = true }, .before = .{ .bytes = before.size, .mtime_ns = before.mtime.nanoseconds }, .after = null, .disappeared = true } }, .compare_source, .fasta, err, error.SourceChanged);
        return diag.fail(.{ .@"io.stat_failed" = {} }, .compare_source, .fasta, err, error.SourceStatFailed);
    };
    if (!sourceMetadataMatches(before, after)) return diag.fail(.{ .@"input.changed" = .{ .changed_fields = .{ .size = before.size != after.size, .mtime = before.mtime.nanoseconds != after.mtime.nanoseconds }, .before = .{ .bytes = before.size, .mtime_ns = before.mtime.nanoseconds }, .after = .{ .bytes = after.size, .mtime_ns = after.mtime.nanoseconds } } }, .compare_source, .fasta, null, error.SourceChanged);
}

fn sourceMetadataMatches(before: std.Io.File.Stat, after: std.Io.File.Stat) bool {
    return before.size == after.size and before.mtime.nanoseconds == after.mtime.nanoseconds;
}

fn scanFaiReader(
    allocator: std.mem.Allocator,
    reader: *std.Io.Reader,
    read_buf: []u8,
    writer: anytype,
    options: ScanOptions,
) !u32 {
    return scanFaiReaderWithDedup(allocator, reader, read_buf, writer, options, @as(?*NameDedup, null));
}

fn scanFaiReaderWithDedup(
    allocator: std.mem.Allocator,
    reader: *std.Io.Reader,
    read_buf: []u8,
    writer: anytype,
    options: ScanOptions,
    dedup: anytype,
) !u32 {
    try checkReadBufferSize(read_buf.len, options.failure);
    var fai_buf = FaiEmitBuffer{ .writer = writer, .diag = .{ .failure = options.failure, .format = .fai } };
    const Ctx = struct {
        buffer: *FaiEmitBuffer,

        fn emit(ctx: *@This(), emit_info: FastaRecordEmit) !void {
            try ctx.buffer.emitRecord(emit_info);
        }
    };

    var replay = FaiNameReplay{ .buffer = &fai_buf, .file = undefined, .io = undefined };
    if (comptime @TypeOf(dedup) == ?*SourceAdapter) {
        if (dedup) |seen| {
            replay.file = seen.source_names.replay_file.?;
            replay.io = seen.source_names.io;
            seen.source_names.replay_context = &replay;
            seen.source_names.replay_fn = FaiNameReplay.replay;
        }
    }
    var ctx = Ctx{ .buffer = &fai_buf };
    const count = try scanFastaReaderWithDedup(allocator, reader, read_buf, options, dedup, &ctx, Ctx.emit);
    try fai_buf.flush();
    return count;
}

fn sourceTimestamp(stat: std.Io.File.Stat, diag: Diagnostics) !u64 {
    return index_format.timestampToNs(stat.mtime) catch return diag.fail(.{ .@"index.source_timestamp_out_of_range" = .{ .timestamp_ns = stat.mtime.nanoseconds, .min_ns = 0, .max_ns = std.math.maxInt(u64) } }, .stat, .fasta, null, error.UnsupportedTimestamp);
}

fn noIndexableRecords(options: Options, diag: Diagnostics) anyerror {
    return diag.fail(.{ .@"index.no_indexable_records" = .{ .naming_policy = if (options.enable_dedup) .first_nonempty else .all_nonempty } }, .scan, .fasta, null, error.NoValidSequences);
}

fn zfiOutputPath(buffer: []u8, path: []const u8, diag: Diagnostics) ![]const u8 {
    return std.fmt.bufPrint(buffer, "{s}.zfi", .{path}) catch return diag.fail(.{ .@"io.invalid_path" = .{ .reason = .too_long, .limit_bytes = buffer.len } }, .create, .zfi, null, error.OutputPathTooLong);
}

fn retainScanCause(failure: ?*diagnostics.Failure, reader: *std.Io.File.Reader, writer: ?*std.Io.File.Writer) void {
    const destination = failure orelse return;
    if (!destination.has_primary) return;
    const primary = &destination.primary_storage;
    if (primary.stage == .read and primary.at.source == .fasta) {
        if (reader.err) |err| primary.cause = .{ .zig = err };
    } else if ((primary.stage == .write or primary.stage == .flush) and primary.at.source == .fai_spool) {
        if (writer) |output| if (output.err) |err| {
            primary.cause = .{ .zig = err };
        };
    }
}

fn indexFai(
    allocator: std.mem.Allocator,
    io: std.Io,
    environ: std.process.Environ,
    path: []const u8,
    options: Options,
    stat: std.Io.File.Stat,
    reader: *std.Io.File.Reader,
    read_buf: []u8,
    dedup: anytype,
    diag: Diagnostics,
) !u32 {
    var spool = try openFaiSpool(allocator, io, environ, diag, .fai_spool);
    defer spool.deinit(io);
    if (comptime @TypeOf(dedup) == ?*SourceAdapter) {
        if (dedup) |seen| seen.source_names.replay_file = spool.file;
    }
    var tmp_io_buf: [INDEX_OUTPUT_BUFFER_SIZE]u8 = undefined;
    var tmp_fw = spool.file.writer(io, &tmp_io_buf);
    const count = scanFaiReaderWithDedup(allocator, &reader.interface, read_buf, &tmp_fw.interface, .{ .enable_dedup = options.enable_dedup, .require_initial_header = true, .failure = diag.failure }, dedup) catch |err| {
        retainScanCause(diag.failure, reader, &tmp_fw);
        return err;
    };
    tmp_fw.interface.flush() catch |err| return diag.fail(.{ .@"io.flush_failed" = {} }, .flush, .fai_spool, if (tmp_fw.err) |cause| cause else err, error.FaiSpoolWriteFailed);
    if (count == 0) return noIndexableRecords(options, diag);
    try checkSourceUnchanged(io, .cwd(), path, stat, diag);
    const tmp_size = (spool.file.stat(io) catch |err| return diag.fail(.{ .@"io.stat_failed" = {} }, .stat, .fai_spool, err, error.FaiSpoolReadFailed)).size;
    var out_buf: [INDEX_OUTPUT_BUFFER_SIZE]u8 = undefined;
    var stdout_fw = std.Io.File.Writer.initStreaming(.stdout(), io, &out_buf);
    var offset: u64 = 0;
    var copy_buf: [INDEX_OUTPUT_BUFFER_SIZE]u8 = undefined;
    while (offset < tmp_size) {
        const want: usize = @intCast(@min(copy_buf.len, tmp_size - offset));
        const n = spool.file.readPositionalAll(io, copy_buf[0..want], offset) catch |err| return diag.failReadAt(.fai_spool, offset, want, null, err, error.FaiSpoolReadFailed);
        if (n == 0) return diag.failReadAt(.fai_spool, offset, want, n, null, error.FaiSpoolReadFailed);
        if (diag.failure) |failure| {
            failure.outputs.stdout.state = .streaming;
            failure.outputs.stdout.may_be_partial = true;
        }
        stdout_fw.interface.writeAll(copy_buf[0..n]) catch |err| return diag.fail(.{ .@"io.write_failed" = {} }, .write, .stdout, if (stdout_fw.err) |cause| cause else err, error.StdoutReplayFailed);
        offset += n;
        if (diag.failure) |failure| failure.outputs.stdout.logical_bytes = offset;
    }
    stdout_fw.flush() catch |err| {
        if (diag.failure) |failure| failure.outputs.stdout.final_flush = .failed;
        return diag.fail(.{ .@"io.flush_failed" = {} }, .flush, .stdout, if (stdout_fw.err) |cause| cause else err, error.StdoutFlushFailed);
    };
    if (diag.failure) |failure| {
        failure.outputs.stdout.final_flush = .succeeded;
        failure.outputs.stdout.records = count;
        failure.outputs.stdout.may_be_partial = false;
        failure.work_complete = true;
    }
    return count;
}

fn runIndexSmall(io: std.Io, environ: std.process.Environ, path: []const u8, options: Options, file: std.Io.File, stat: std.Io.File.Stat, diag: Diagnostics) !u32 {
    const allocator = std.heap.page_allocator;
    var io_buf: [FILE_IO_BUF_SIZE]u8 = undefined;
    var read_buf: [INDEX_READ_BUFFER_SIZE]u8 = undefined;
    var file_reader = file.reader(io, &io_buf);
    if (options.emit_fai) {
        var seen = PackedNames.init(allocator);
        seen.diag = diag;
        defer seen.deinit();
        if (options.enable_dedup) try seen.ensureCapacity(DEDUP_INITIAL_CAPACITY);
        const dedup: ?*PackedNames = if (options.enable_dedup) &seen else null;
        return indexFai(allocator, io, environ, path, options, stat, &file_reader, &read_buf, dedup, diag);
    }
    var path_buf: [INDEX_PATH_BUFFER_SIZE]u8 = undefined;
    const zfi_path = try zfiOutputPath(&path_buf, path, diag);
    var index = scanZfiReaderWithOptions(allocator, &file_reader.interface, &read_buf, .{ .enable_dedup = options.enable_dedup, .require_initial_header = true, .failure = diag.failure }) catch |err| {
        retainScanCause(diag.failure, &file_reader, null);
        return err;
    };
    defer index.deinit(allocator);
    const count: u32 = @intCast(index.records.items.len);
    if (count == 0) return noIndexableRecords(options, diag);
    const mtime = try sourceTimestamp(stat, diag);
    var atomic = std.Io.Dir.cwd().createFileAtomic(io, zfi_path, .{ .replace = true }) catch |err| return diag.fail(.{ .@"io.create_failed" = {} }, .create, .zfi, err, error.ZfiWriteFailed);
    defer cleanupAtomic(&atomic, io, diag.failure);
    if (diag.failure) |failure| failure.outputs.file.state = .staged;
    var out_buf: [INDEX_OUTPUT_BUFFER_SIZE]u8 = undefined;
    var output = atomic.file.writer(io, &out_buf);
    writeZfiIndex(&output.interface, &index, stat.size, mtime) catch |err| return diag.fail(.{ .@"io.write_failed" = {} }, .write, .zfi, if (output.err) |cause| cause else err, error.ZfiWriteFailed);
    output.flush() catch |err| {
        if (diag.failure) |failure| failure.outputs.file.final_flush = .failed;
        return diag.fail(.{ .@"io.flush_failed" = {} }, .flush, .zfi, if (output.err) |cause| cause else err, error.ZfiWriteFailed);
    };
    if (diag.failure) |failure| failure.outputs.file.final_flush = .succeeded;
    try checkSourceUnchanged(io, .cwd(), path, stat, diag);
    atomic.replace(io) catch |err| {
        if (diag.failure) |failure| failure.outputs.file.state = .unknown;
        return diag.fail(.{ .@"io.publish_failed" = {} }, .publish, .zfi, err, error.ZfiFinalizeFailed);
    };
    if (diag.failure) |failure| {
        failure.outputs.file.state = .published;
        failure.outputs.file.records = count;
        failure.work_complete = true;
    }
    return count;
}

fn runIndexLarge(io: std.Io, environ: std.process.Environ, path: []const u8, options: Options, file: std.Io.File, stat: std.Io.File.Stat, diag: Diagnostics) !u32 {
    const allocator = std.heap.page_allocator;
    var names = SourceDedup.init(allocator, io, file, stat.size);
    names.diag = diag;
    defer names.deinit();
    var adapter = SourceAdapter{ .source_names = &names };
    const dedup: ?*SourceAdapter = if (options.enable_dedup) &adapter else null;
    var io_buf: [FILE_IO_BUF_SIZE]u8 = undefined;
    var read_buf: [INDEX_READ_BUFFER_SIZE]u8 = undefined;
    var reader = file.reader(io, &io_buf);
    if (options.emit_fai) return indexFai(allocator, io, environ, path, options, stat, &reader, &read_buf, dedup, diag);
    var path_buf: [INDEX_PATH_BUFFER_SIZE]u8 = undefined;
    const zfi_path = try zfiOutputPath(&path_buf, path, diag);
    var deferred = diag;
    deferred.deferred = true;
    deferred.at = .{ .source = .fasta };
    const storage = CatalogBytes{ .allocator = allocator, .io = io, .environ = environ, .diag = deferred, .role = .name_spool };
    var catalog = StreamedCatalog{ .io = io, .path = zfi_path, .diag = deferred, .names = .{ .raw = storage, .compressed = .{ .encoded = storage } }, .sides = storage };
    catalog.sides.role = .side_spool;
    defer catalog.deinit();
    names.replay_context = &catalog;
    names.replay_fn = StreamedCatalog.replayNames;
    const count = scanFastaReaderWithDedup(allocator, &reader.interface, &read_buf, .{ .enable_dedup = options.enable_dedup, .require_initial_header = true, .collect_side_tables = true, .failure = diag.failure }, dedup, &catalog, StreamedCatalog.emit) catch |err| {
        retainScanCause(diag.failure, &reader, null);
        return err;
    };
    if (count == 0) return noIndexableRecords(options, diag);
    const mtime = try sourceTimestamp(stat, diag);
    if (diag.failure) |failure| try failure.finishDeferred();
    catalog.diag.deferred = false;
    catalog.names.raw.diag.deferred = false;
    catalog.names.compressed.encoded.diag.deferred = false;
    catalog.sides.diag.deferred = false;
    try catalog.finish(count, stat.size, mtime);
    try checkSourceUnchanged(io, .cwd(), path, stat, diag);
    catalog.atomic_file.?.replace(io) catch |err| {
        if (diag.failure) |failure| failure.outputs.file.state = .unknown;
        return diag.fail(.{ .@"io.publish_failed" = {} }, .publish, .zfi, err, error.ZfiFinalizeFailed);
    };
    if (diag.failure) |failure| {
        failure.outputs.file.state = .published;
        failure.outputs.file.records = count;
        failure.work_complete = true;
    }
    return count;
}

const PackedNames = struct {
    diag: Diagnostics = .{},
    allocator: std.mem.Allocator,
    blocks: std.ArrayList([]u8) = .empty,
    position: usize = 0,
    map: Map,
    const BLOCK_BYTES = 128 * 1024;
    const Ref = struct { offset: u32, len: u32 };
    const Context = struct {
        owner: *PackedNames,
        pub fn hash(self: @This(), ref: Ref) u64 {
            return std.hash_map.hashString(self.owner.bytes(ref));
        }
        pub fn eql(self: @This(), a: Ref, b: Ref) bool {
            return a.len == b.len and std.mem.eql(u8, self.owner.bytes(a), self.owner.bytes(b));
        }
    };
    const Lookup = struct {
        owner: *PackedNames,
        pub fn hash(_: @This(), name: []const u8) u64 {
            return std.hash_map.hashString(name);
        }
        pub fn eql(self: @This(), name: []const u8, ref: Ref) bool {
            return name.len == ref.len and std.mem.eql(u8, name, self.owner.bytes(ref));
        }
    };
    const Map = std.HashMap(Ref, void, Context, std.hash_map.default_max_load_percentage);

    fn init(allocator: std.mem.Allocator) @This() {
        return .{ .allocator = allocator, .map = Map.initContext(allocator, .{ .owner = undefined }) };
    }

    fn ensureCapacity(self: *@This(), capacity: u32) !void {
        self.map.ctx.owner = self;
        self.map.ensureTotalCapacity(capacity) catch return self.diag.oom(.dedup_table, null);
    }

    fn deinit(self: *@This()) void {
        self.map.deinit();
        for (self.blocks.items) |block| self.allocator.free(block);
        self.blocks.deinit(self.allocator);
    }

    fn bytes(self: *@This(), ref: Ref) []const u8 {
        if (ref.len == 0) return &.{};
        return self.blocks.items[ref.offset / BLOCK_BYTES][ref.offset % BLOCK_BYTES ..][0..ref.len];
    }

    fn observeAt(self: *@This(), name: []const u8, source_offset: u64) !bool {
        _ = source_offset;
        self.map.ctx.owner = self;
        const entry = self.map.getOrPutAdapted(name, Lookup{ .owner = self }) catch return self.diag.oom(.dedup_table, null);
        if (entry.found_existing) return true;
        errdefer self.map.removeByPtr(entry.key_ptr);
        const remaining = BLOCK_BYTES - self.position % BLOCK_BYTES;
        if (name.len > remaining) {
            self.position = std.math.add(usize, self.position, remaining) catch return self.diag.capacity(.name_offset, std.math.maxInt(usize), @as(u128, self.position) + remaining, error.OutOfMemory);
        }
        const offset = std.math.cast(u32, self.position) orelse return self.diag.capacity(.name_offset, std.math.maxInt(u32), self.position, error.OutOfMemory);
        const next_position = std.math.add(usize, self.position, name.len) catch return self.diag.capacity(.name_offset, std.math.maxInt(usize), @as(u128, self.position) + name.len, error.OutOfMemory);
        if (name.len != 0 and self.position / BLOCK_BYTES == self.blocks.items.len) {
            const block = self.allocator.alloc(u8, BLOCK_BYTES) catch return self.diag.oom(.names, BLOCK_BYTES);
            self.blocks.append(self.allocator, block) catch {
                self.allocator.free(block);
                return self.diag.oom(.names, null);
            };
        }
        const reference = Ref{ .offset = offset, .len = @intCast(name.len) };
        if (name.len != 0) @memcpy(self.blocks.items[self.position / BLOCK_BYTES][self.position % BLOCK_BYTES ..][0..name.len], name);
        self.position = next_position;
        entry.key_ptr.* = reference;
        return false;
    }
};

const CACHE_BYTES = 512 * 1024;
const CACHE_BLOCK_BYTES: usize = 128 * 1024;

const SourceDedup = struct {
    diag: Diagnostics = .{},
    allocator: std.mem.Allocator,
    io: std.Io,
    file: std.Io.File,
    tables: union(enum) { narrow: Tables(u32), wide: Tables(u64) },
    chunk: []const u8 = &.{},
    chunk_offset: u64 = 0,
    cache: std.ArrayList(u8) = .empty,
    scratch: std.ArrayList(u8) = .empty,
    cache_full: bool = false,
    allocation_failed: bool = false,
    source_match: bool = false,
    source_name_hits: u32 = 0,
    hydrated: bool = false,
    hot_blocks: std.ArrayList([]u8) = .empty,
    hot_len: usize = 0,
    replay_file: ?std.Io.File = null,
    replay_context: ?*anyopaque = null,
    replay_fn: ?*const fn (*anyopaque, *SourceDedup) anyerror!void = null,
    read_failed: bool = false,

    fn init(allocator: std.mem.Allocator, io: std.Io, file: std.Io.File, size: u64) @This() {
        return .{
            .allocator = allocator,
            .io = io,
            .file = file,
            .tables = if (size <= std.math.maxInt(u31) - CACHE_BYTES) .{ .narrow = Tables(u32).init(allocator) } else .{ .wide = Tables(u64).init(allocator) },
        };
    }

    fn deinit(self: *@This()) void {
        switch (self.tables) {
            .narrow => |*set| set.deinit(),
            .wide => |*set| set.deinit(),
        }
        self.cache.deinit(self.allocator);
        self.scratch.deinit(self.allocator);
        for (self.hot_blocks.items) |block| self.allocator.free(block);
        self.hot_blocks.deinit(self.allocator);
    }

    fn setChunk(self: *@This(), data: []const u8, offset: u64) void {
        self.chunk = data;
        self.chunk_offset = offset;
    }

    fn observe(self: *@This(), name: []const u8, offset: u64) !bool {
        return switch (self.tables) {
            .narrow => |*set| if (offset >= Tables(u32).CACHE_FLAG)
                self.diag.capacity(.name_cache_offset, Tables(u32).CACHE_FLAG - 1, offset, error.SourceChanged)
            else
                set.observe(self, name, @intCast(offset)),
            .wide => |*set| if (offset >= Tables(u64).CACHE_FLAG)
                self.diag.capacity(.name_cache_offset, Tables(u64).CACHE_FLAG - 1, offset, error.SourceChanged)
            else
                set.observe(self, name, offset),
        };
    }

    const CacheEncoding = struct { header: [7]u8, size: usize, middle: []const u8 };

    fn encodeCacheName(name: []const u8, prefix: usize, suffix: usize) CacheEncoding {
        const middle = name[prefix .. name.len - suffix];
        const short = prefix <= 255 and suffix <= 255 and middle.len <= 255;
        const header_size: usize = if (short) 4 else 7;
        var header: [7]u8 = undefined;
        if (name.len <= 253 and name.len + 1 <= header_size + middle.len) {
            header[0] = @intCast(name.len);
            return .{ .header = header, .size = 1, .middle = name };
        }
        if (short) {
            header[0] = 254;
            header[1] = @intCast(prefix);
            header[2] = @intCast(suffix);
            header[3] = @intCast(middle.len);
        } else {
            header[0] = 255;
            std.mem.writeInt(u16, header[1..3], @intCast(prefix), .little);
            std.mem.writeInt(u16, header[3..5], @intCast(suffix), .little);
            std.mem.writeInt(u16, header[5..7], @intCast(middle.len), .little);
        }
        return .{ .header = header, .size = header_size, .middle = middle };
    }

    fn cacheAnchor(self: *@This(), base: usize) []const u8 {
        return switch (self.cache.items[base]) {
            0...253 => |len| self.cache.items[base + 1 ..][0..len],
            254 => self.cache.items[base + 4 ..][0..self.cache.items[base + 3]],
            255 => self.cache.items[base + 7 ..][0..std.mem.readInt(u16, self.cache.items[base + 5 ..][0..2], .little)],
        };
    }

    fn remember(self: *@This(), name: []const u8, offset: u64, flag: u64) !u64 {
        if (self.cache_full or self.cache.items.len >= CACHE_BYTES) return offset;
        var prefix: usize = 0;
        var suffix: usize = 0;
        if (name.len >= 16 and self.cache.items.len % CACHE_BLOCK_BYTES != 0) {
            const base = (self.cache.items.len - 1) & ~(CACHE_BLOCK_BYTES - 1);
            const anchor = self.cacheAnchor(base);
            prefix = commonPrefix(name, anchor);
            suffix = commonSuffix(name, anchor, prefix);
        }
        var encoded = encodeCacheName(name, prefix, suffix);
        var start = self.cache.items.len;
        const remaining = CACHE_BLOCK_BYTES - (start % CACHE_BLOCK_BYTES);
        if (encoded.size + encoded.middle.len > remaining) {
            start += remaining;
            encoded = encodeCacheName(name, 0, 0);
        }
        if (start + encoded.size + encoded.middle.len > CACHE_BYTES) {
            self.cache_full = true;
            return offset;
        }
        if (self.cache.capacity == 0) self.cache.ensureTotalCapacity(self.allocator, CACHE_BYTES) catch return self.diag.oom(.name_cache, CACHE_BYTES);
        self.cache.items.len = start;
        self.cache.appendSliceAssumeCapacity(encoded.header[0..encoded.size]);
        self.cache.appendSliceAssumeCapacity(encoded.middle);
        return flag | start;
    }

    fn sameCached(self: *@This(), name: []const u8, start: usize) bool {
        const tag = self.cache.items[start];
        if (tag <= 253) {
            return name.len == tag and std.mem.eql(u8, name, self.cache.items[start + 1 ..][0..tag]);
        }
        const prefix: usize = if (tag == 254) self.cache.items[start + 1] else std.mem.readInt(u16, self.cache.items[start + 1 ..][0..2], .little);
        const suffix: usize = if (tag == 254) self.cache.items[start + 2] else std.mem.readInt(u16, self.cache.items[start + 3 ..][0..2], .little);
        const middle_len: usize = if (tag == 254) self.cache.items[start + 3] else std.mem.readInt(u16, self.cache.items[start + 5 ..][0..2], .little);
        const header_len: usize = if (tag == 254) 4 else 7;
        if (name.len != prefix + suffix + middle_len) return false;
        const anchor = self.cacheAnchor(start & ~(CACHE_BLOCK_BYTES - 1));
        return std.mem.eql(u8, name[0..prefix], anchor[0..prefix]) and
            std.mem.eql(u8, name[prefix..][0..middle_len], self.cache.items[start + header_len ..][0..middle_len]) and
            std.mem.eql(u8, name[name.len - suffix ..], anchor[anchor.len - suffix ..]);
    }

    fn retainAt(self: *@This(), ordinal: u32, name: []const u8) !void {
        switch (self.tables) {
            .narrow => |*set| try set.retainAt(self, ordinal, name),
            .wide => |*set| try set.retainAt(self, ordinal, name),
        }
    }

    fn uniqueCount(self: *const @This()) u32 {
        return switch (self.tables) {
            .narrow => |*set| set.count(),
            .wide => |*set| set.count(),
        };
    }

    fn appendHot(self: *@This(), input: []const u8) !void {
        var remaining = input;
        while (remaining.len != 0) {
            const block_index = self.hot_len / CACHE_BLOCK_BYTES;
            const within = self.hot_len % CACHE_BLOCK_BYTES;
            if (block_index == self.hot_blocks.items.len) {
                const block = self.allocator.alloc(u8, CACHE_BLOCK_BYTES) catch return self.diag.oom(.retained_names, CACHE_BLOCK_BYTES);
                errdefer self.allocator.free(block);
                self.hot_blocks.append(self.allocator, block) catch return self.diag.oom(.retained_names, null);
            }
            const n = @min(remaining.len, CACHE_BLOCK_BYTES - within);
            @memcpy(self.hot_blocks.items[block_index][within..][0..n], remaining[0..n]);
            self.hot_len += n;
            remaining = remaining[n..];
        }
    }

    fn rememberHot(self: *@This(), name: []const u8, flag: u64) !u64 {
        if (self.hot_len % CACHE_BLOCK_BYTES == CACHE_BLOCK_BYTES - 1) try self.appendHot(&.{0});
        const start = self.hot_len;
        if (@as(u128, CACHE_BYTES) + start >= flag) return self.diag.capacity(.name_cache_offset, flag - 1, @as(u128, CACHE_BYTES) + start, error.OutOfMemory);
        var header: [2]u8 = undefined;
        std.mem.writeInt(u16, &header, @intCast(name.len), .little);
        try self.appendHot(&header);
        try self.appendHot(name);
        return flag | (CACHE_BYTES + start);
    }

    fn sameHot(self: *@This(), name: []const u8, start: usize) bool {
        const header_block = self.hot_blocks.items[start / CACHE_BLOCK_BYTES];
        const len = std.mem.readInt(u16, header_block[start % CACHE_BLOCK_BYTES ..][0..2], .little);
        if (len != name.len) return false;
        if (name.len == 0) return true;
        const bytes_start = start + 2;
        const within = bytes_start % CACHE_BLOCK_BYTES;
        const n = @min(name.len, CACHE_BLOCK_BYTES - within);
        const block_index = bytes_start / CACHE_BLOCK_BYTES;
        if (!std.mem.eql(u8, name[0..n], self.hot_blocks.items[block_index][within..][0..n])) return false;
        return n == name.len or std.mem.eql(u8, name[n..], self.hot_blocks.items[block_index + 1][0 .. name.len - n]);
    }

    fn sameName(self: *@This(), name: []const u8, offset: u64, flag: u64) bool {
        if (offset & flag != 0) {
            const start: usize = @intCast(offset & ~flag);
            return if (start < CACHE_BYTES) self.sameCached(name, start) else self.sameHot(name, start - CACHE_BYTES);
        }
        const needed = name.len + 1;
        if (offset >= self.chunk_offset and offset - self.chunk_offset <= self.chunk.len) {
            const start: usize = @intCast(offset - self.chunk_offset);
            if (needed <= self.chunk.len - start) return sameBytes(name, self.chunk[start..][0..needed]);
        }
        self.scratch.resize(self.allocator, needed) catch {
            switch (self.diag.oom(.retained_names, needed)) {
                else => {},
            }
            self.allocation_failed = true;
            return false;
        };
        const read = self.file.readPositionalAll(self.io, self.scratch.items[0..needed], offset) catch |err| {
            var diag = self.diag;
            diag.at.location.offset = offset;
            switch (diag.fail(.{ .@"io.read_failed" = .{ .requested_bytes = needed } }, .read, .fasta, err, error.SourceReadFailed)) {
                else => {},
            }
            self.read_failed = true;
            return false;
        };
        if (read != needed) {
            for (self.scratch.items[0..read]) |byte| {
                switch (byte) {
                    ' ', '\t', '\r', '\n' => return false,
                    else => {},
                }
            }
            var diag = self.diag;
            diag.at.location.offset = offset;
            switch (diag.fail(.{ .@"io.unexpected_eof" = .{ .exact = .{ .requested_bytes = needed, .received_bytes = read } } }, .read, .fasta, null, error.SourceReadFailed)) {
                else => {},
            }
            self.read_failed = true;
            return false;
        }
        const matches = sameBytes(name, self.scratch.items[0..needed]);
        if (matches) self.source_match = true;
        return matches;
    }

    fn sameBytes(name: []const u8, bytes: []const u8) bool {
        const terminator = bytes[name.len];
        return (terminator == ' ' or terminator == '\t' or terminator == '\r' or terminator == '\n') and
            std.mem.eql(u8, name, bytes[0..name.len]);
    }
};

fn commonPrefix(a: []const u8, b: []const u8) usize {
    const Vec = @Vector(16, u8);
    const limit = @min(a.len, b.len);
    var position: usize = 0;
    while (limit - position >= 16) : (position += 16) {
        const av: Vec = a[position..][0..16].*;
        const bv: Vec = b[position..][0..16].*;
        if (std.simd.firstTrue(av != bv)) |different| return position + different;
    }
    while (position < limit and a[position] == b[position]) : (position += 1) {}
    return position;
}

fn commonSuffix(a: []const u8, b: []const u8, prefix: usize) usize {
    const Vec = @Vector(16, u8);
    const limit = @min(a.len, b.len) - prefix;
    var matched: usize = 0;
    while (limit - matched >= 16) : (matched += 16) {
        const av: Vec = a[a.len - matched - 16 ..][0..16].*;
        const bv: Vec = b[b.len - matched - 16 ..][0..16].*;
        if (std.simd.lastTrue(av != bv)) |different| return matched + 15 - different;
    }
    while (matched < limit and a[a.len - matched - 1] == b[b.len - matched - 1]) : (matched += 1) {}
    return matched;
}

fn Tables(comptime Offset: type) type {
    return struct {
        const CACHE_FLAG: u64 = @as(u64, 1) << (@bitSizeOf(Offset) - 1);
        const Key = struct { fingerprint: u32, offset: Offset };
        const BLOCK_ENTRIES = 4096;
        const Vec = @Vector(16, u8);
        allocator: std.mem.Allocator,
        buckets: []u32 = &.{},
        metadata: []u8 = &.{},
        blocks: std.ArrayList([]Key) = .empty,
        len: u32 = 0,

        fn init(allocator: std.mem.Allocator) @This() {
            return .{ .allocator = allocator };
        }

        fn deinit(self: *@This()) void {
            self.allocator.free(self.buckets);
            self.allocator.free(self.metadata);
            for (self.blocks.items) |block| self.allocator.free(block);
            self.blocks.deinit(self.allocator);
        }

        fn count(self: *const @This()) u32 {
            return self.len;
        }

        fn keyAt(self: *@This(), index: u32) *Key {
            return &self.blocks.items[index / BLOCK_ENTRIES][index % BLOCK_ENTRIES];
        }

        fn tag(fingerprint: u32) u8 {
            return 0x80 | @as(u8, @truncate(fingerprint >> 25));
        }

        fn emptySlot(metadata: []const u8, fingerprint: u32) usize {
            const mask = metadata.len - 1;
            var group = spread(fingerprint) & mask & ~@as(usize, 15);
            while (true) : (group = (group + 16) & mask) {
                const bytes: Vec = metadata[group..][0..16].*;
                if (std.simd.firstTrue(bytes == @as(Vec, @splat(0)))) |index| return group + index;
            }
        }

        fn grow(self: *@This(), diag: Diagnostics) !void {
            var required: u128 = if (self.buckets.len == 0) 512 else @as(u128, self.buckets.len) * 2;
            while ((@as(u128, self.len) + 1) * 5 > required * 4) required *= 2;
            const capacity = std.math.cast(usize, required) orelse return diag.capacity(.name_map, std.math.maxInt(usize), required, error.OutOfMemory);
            if (required * @sizeOf(u32) > std.math.maxInt(usize)) return diag.capacity(.name_map, std.math.maxInt(usize), required * @sizeOf(u32), error.OutOfMemory);
            self.allocator.free(self.buckets);
            self.allocator.free(self.metadata);
            self.buckets = &.{};
            self.metadata = &.{};
            const metadata = self.allocator.alloc(u8, capacity) catch return diag.oom(.dedup_table, capacity);
            errdefer self.allocator.free(metadata);
            @memset(metadata, 0);
            const buckets = self.allocator.alloc(u32, capacity) catch return diag.oom(.dedup_table, capacity * @sizeOf(u32));
            for (0..self.len) |entry_index| {
                const key = self.keyAt(@intCast(entry_index));
                const position = emptySlot(metadata, key.fingerprint);
                buckets[position] = @intCast(entry_index + 1);
                metadata[position] = tag(key.fingerprint);
            }
            self.buckets = buckets;
            self.metadata = metadata;
        }

        fn observe(self: *@This(), source: *SourceDedup, name: []const u8, offset: Offset) !bool {
            return self.observeWithFingerprint(source, name, offset, @truncate(std.hash_map.hashString(name)));
        }

        fn observeWithFingerprint(self: *@This(), source: *SourceDedup, name: []const u8, offset: Offset, fingerprint: u32) !bool {
            if (self.buckets.len == 0) try self.grow(source.diag);
            const mask = self.buckets.len - 1;
            var group = spread(fingerprint) & mask & ~@as(usize, 15);
            const wanted: Vec = @splat(tag(fingerprint));
            source.source_match = false;
            var position: usize = undefined;
            while (true) : (group = (group + 16) & mask) {
                const bytes: Vec = self.metadata[group..][0..16].*;
                var matches: u16 = @bitCast(bytes == wanted);
                while (matches != 0) {
                    const candidate = group + @ctz(matches);
                    matches &= matches - 1;
                    const key = self.keyAt(self.buckets[candidate] - 1);
                    if (key.fingerprint != fingerprint) continue;
                    if (source.sameName(name, key.offset, CACHE_FLAG)) {
                        if (source.source_match) {
                            key.offset = @intCast(try source.rememberHot(name, CACHE_FLAG));
                            source.source_name_hits += 1;
                            if (!source.hydrated and source.source_name_hits >= 64) {
                                if (source.replay_fn) |replay| {
                                    source.hydrated = true;
                                    try replay(source.replay_context.?, source);
                                }
                            }
                        }
                        return true;
                    }
                    if (source.read_failed) return error.SourceReadFailed;
                    if (source.allocation_failed) return error.OutOfMemory;
                }
                if (std.simd.firstTrue(bytes == @as(Vec, @splat(0)))) |empty| {
                    position = group + empty;
                    break;
                }
            }
            if (self.len == std.math.maxInt(u32)) {
                @branchHint(.cold);
                return source.diag.recordLimit(@as(u64, self.len) + 1);
            }
            if ((@as(u64, self.len) + 1) * 5 > @as(u64, self.buckets.len) * 4) {
                try self.grow(source.diag);
                position = emptySlot(self.metadata, fingerprint);
            }
            if (self.len == self.blocks.items.len * BLOCK_ENTRIES) {
                const block = self.allocator.alloc(Key, BLOCK_ENTRIES) catch return source.diag.oom(.dedup_table, BLOCK_ENTRIES * @sizeOf(Key));
                errdefer self.allocator.free(block);
                self.blocks.append(self.allocator, block) catch return source.diag.oom(.dedup_table, null);
            }
            const reference = try source.remember(name, offset, CACHE_FLAG);
            self.keyAt(self.len).* = .{ .fingerprint = fingerprint, .offset = @intCast(reference) };
            self.buckets[position] = self.len + 1;
            self.metadata[position] = tag(fingerprint);
            self.len += 1;
            return false;
        }

        fn retainAt(self: *@This(), source: *SourceDedup, ordinal: u32, name: []const u8) !void {
            if (ordinal >= self.len) return source.diag.consistency(.record_ordinal, error.SourceReadFailed);
            const key = self.keyAt(ordinal);
            if (key.offset & CACHE_FLAG != 0) return;
            key.offset = @intCast(try source.rememberHot(name, CACHE_FLAG));
        }

        fn spread(fingerprint: u32) usize {
            return @truncate(@as(u64, fingerprint) *% 0x9e3779b97f4a7c15);
        }
    };
}

const CatalogBytes = struct {
    diag: Diagnostics = .{},
    allocator: std.mem.Allocator,
    io: std.Io,
    environ: std.process.Environ,
    role: diagnostics.SourceId = .catalog_spool,
    bytes: std.ArrayList(u8) = .empty,
    spool: ?FaiSpool = null,
    buffer: [INDEX_OUTPUT_BUFFER_SIZE]u8 = undefined,
    buffered: usize = 0,
    len: u64 = 0,
    const MEMORY_BYTES = 1024 * 1024;

    fn deinit(self: *@This()) void {
        self.bytes.deinit(self.allocator);
        if (self.spool) |*spool| spool.deinit(self.io);
    }

    fn append(self: *@This(), data: []const u8) !void {
        if (self.spool == null and self.len + data.len <= MEMORY_BYTES) {
            self.bytes.appendSlice(self.allocator, data) catch return self.diag.oom(if (self.role == .side_spool) .side_tables else .names, data.len);
            self.len += data.len;
            return;
        }
        if (self.spool == null) {
            self.spool = try openFaiSpool(self.allocator, self.io, self.environ, self.diag, self.role);
            self.spool.?.file.writePositionalAll(self.io, self.bytes.items, 0) catch |err| return self.diag.fail(.{ .@"io.write_failed" = {} }, .write, self.role, err, error.WriteFailed);
            self.bytes.deinit(self.allocator);
            self.bytes = .empty;
        }
        var remaining = data;
        while (remaining.len > 0) {
            const n = @min(self.buffer.len - self.buffered, remaining.len);
            @memcpy(self.buffer[self.buffered..][0..n], remaining[0..n]);
            self.buffered += n;
            self.len += n;
            remaining = remaining[n..];
            if (self.buffered == self.buffer.len) try self.flush();
        }
    }

    fn flush(self: *@This()) !void {
        if (self.buffered == 0) return;
        self.spool.?.file.writePositionalAll(self.io, self.buffer[0..self.buffered], self.len - self.buffered) catch |err| return self.diag.fail(.{ .@"io.write_failed" = {} }, .write, self.role, err, error.WriteFailed);
        self.buffered = 0;
    }

    fn take(self: *@This(), reader: *std.Io.Reader, file_reader: ?*std.Io.File.Reader, count: usize) ![]const u8 {
        return reader.take(count) catch |err| {
            if (self.spool == null) return self.diag.consistency(.codec_state, err);
            const cause: anyerror = if (file_reader) |input| if (input.err) |failure| failure else err else err;
            if (err == error.EndOfStream) return self.diag.fail(.{ .@"io.unexpected_eof" = .{ .exact = .{ .requested_bytes = count, .received_bytes = reader.bufferedLen() } } }, .read, self.role, null, err);
            return self.diag.fail(.{ .@"io.read_failed" = .{ .requested_bytes = count } }, .read, self.role, cause, err);
        };
    }

    fn copyTo(self: *@This(), writer: *std.Io.Writer) !void {
        try self.flush();
        if (self.spool == null) {
            writer.writeAll(self.bytes.items) catch |err| return self.diag.fail(.{ .@"io.write_failed" = {} }, .write, .zfi, err, error.WriteFailed);
            return;
        }
        var offset: u64 = 0;
        var buffer: [INDEX_OUTPUT_BUFFER_SIZE]u8 = undefined;
        while (offset < self.len) {
            const wanted: usize = @intCast(@min(buffer.len, self.len - offset));
            const n = self.spool.?.file.readPositionalAll(self.io, buffer[0..wanted], offset) catch |err| return self.diag.failReadAt(self.role, offset, wanted, null, err, error.ReadFailed);
            if (n != wanted) return self.diag.failReadAt(self.role, offset, wanted, n, null, error.ReadFailed);
            writer.writeAll(buffer[0..n]) catch |err| return self.diag.fail(.{ .@"io.write_failed" = {} }, .write, .zfi, err, error.WriteFailed);
            offset += n;
        }
    }
};

const CompressedNames = struct {
    encoded: CatalogBytes,
    anchor: [MAX_INDEX_NAME_LEN]u8 = undefined,
    anchor_len: usize = 0,
    block_bytes: usize = 0,
    count: u32 = 0,

    fn deinit(self: *@This()) void {
        self.encoded.deinit();
    }

    fn append(self: *@This(), name: []const u8) !void {
        var prefix: usize = 0;
        var suffix: usize = 0;
        if (self.count != 0) {
            prefix = commonPrefix(name, self.anchor[0..self.anchor_len]);
            suffix = commonSuffix(name, self.anchor[0..self.anchor_len], prefix);
        }
        var middle = name[prefix .. name.len - suffix];
        if (self.count == 0 or self.block_bytes + 6 + middle.len > 128 * 1024) {
            prefix = MAX_INDEX_NAME_LEN;
            suffix = MAX_INDEX_NAME_LEN;
            middle = name;
            @memcpy(self.anchor[0..name.len], name);
            self.anchor_len = name.len;
            self.block_bytes = 0;
        }
        var header: [6]u8 = undefined;
        std.mem.writeInt(u16, header[0..2], @intCast(prefix), .little);
        std.mem.writeInt(u16, header[2..4], @intCast(suffix), .little);
        std.mem.writeInt(u16, header[4..6], @intCast(middle.len), .little);
        try self.encoded.append(&header);
        try self.encoded.append(middle);
        self.block_bytes += 6 + middle.len;
        self.count += 1;
    }

    fn retain(self: *@This(), names: *SourceDedup) !void {
        try self.encoded.flush();
        var read_buffer: [INDEX_OUTPUT_BUFFER_SIZE]u8 = undefined;
        var memory_reader = std.Io.Reader.fixed(self.encoded.bytes.items);
        var file_reader: std.Io.File.Reader = undefined;
        const reader = if (self.encoded.spool) |spool| blk: {
            file_reader = spool.file.reader(self.encoded.io, &read_buffer);
            break :blk &file_reader.interface;
        } else &memory_reader;
        var anchor_len: usize = 0;
        for (0..self.count) |ordinal| {
            const header = (try self.encoded.take(reader, if (self.encoded.spool != null) &file_reader else null, 6))[0..6].*;
            const prefix = std.mem.readInt(u16, header[0..2], .little);
            const suffix = std.mem.readInt(u16, header[2..4], .little);
            const middle_len = std.mem.readInt(u16, header[4..6], .little);
            const middle = try self.encoded.take(reader, if (self.encoded.spool != null) &file_reader else null, middle_len);
            if (prefix == MAX_INDEX_NAME_LEN and suffix == MAX_INDEX_NAME_LEN) {
                @memcpy(self.anchor[0..middle.len], middle);
                anchor_len = middle.len;
                try names.retainAt(@intCast(ordinal), middle);
            } else {
                const full_len = @as(usize, prefix) + middle_len + suffix;
                names.scratch.resize(names.allocator, full_len) catch return names.diag.oom(.retained_names, full_len);
                @memcpy(names.scratch.items[0..prefix], self.anchor[0..prefix]);
                @memcpy(names.scratch.items[prefix..][0..middle_len], middle);
                @memcpy(names.scratch.items[@as(usize, prefix) + middle_len ..][0..suffix], self.anchor[anchor_len - suffix .. anchor_len]);
                try names.retainAt(@intCast(ordinal), names.scratch.items[0..full_len]);
            }
        }
    }

    fn write(self: *@This(), writer: *std.Io.Writer) !void {
        try self.encoded.flush();
        var read_buffer: [INDEX_OUTPUT_BUFFER_SIZE]u8 = undefined;
        var memory_reader = std.Io.Reader.fixed(self.encoded.bytes.items);
        var file_reader: std.Io.File.Reader = undefined;
        const reader = if (self.encoded.spool) |spool| blk: {
            file_reader = spool.file.reader(self.encoded.io, &read_buffer);
            break :blk &file_reader.interface;
        } else &memory_reader;
        var anchor_len: usize = 0;
        for (0..self.count) |_| {
            const header = (try self.encoded.take(reader, if (self.encoded.spool != null) &file_reader else null, 6))[0..6].*;
            const prefix = std.mem.readInt(u16, header[0..2], .little);
            const suffix = std.mem.readInt(u16, header[2..4], .little);
            const middle_len = std.mem.readInt(u16, header[4..6], .little);
            const middle = try self.encoded.take(reader, if (self.encoded.spool != null) &file_reader else null, middle_len);
            if (prefix == MAX_INDEX_NAME_LEN and suffix == MAX_INDEX_NAME_LEN) {
                @memcpy(self.anchor[0..middle.len], middle);
                anchor_len = middle.len;
                writer.writeAll(middle) catch |err| return self.encoded.diag.fail(.{ .@"io.write_failed" = {} }, .write, .zfi, err, error.WriteFailed);
            } else {
                writer.writeAll(self.anchor[0..prefix]) catch |err| return self.encoded.diag.fail(.{ .@"io.write_failed" = {} }, .write, .zfi, err, error.WriteFailed);
                writer.writeAll(middle) catch |err| return self.encoded.diag.fail(.{ .@"io.write_failed" = {} }, .write, .zfi, err, error.WriteFailed);
                writer.writeAll(self.anchor[anchor_len - suffix .. anchor_len]) catch |err| return self.encoded.diag.fail(.{ .@"io.write_failed" = {} }, .write, .zfi, err, error.WriteFailed);
            }
        }
    }
};

const StreamedCatalog = struct {
    diag: Diagnostics = .{},
    io: std.Io,
    dir: std.Io.Dir = .cwd(),
    path: []const u8,
    names: StagedNames,
    sides: CatalogBytes,
    atomic_file: ?std.Io.File.Atomic = null,
    records: [INDEX_OUTPUT_BUFFER_SIZE]u8 = undefined,
    buffered: usize = 0,
    record_bytes: u64 = 0,
    write_failed: bool = false,

    fn deinit(self: *@This()) void {
        self.names.deinit();
        self.sides.deinit();
        if (self.atomic_file) |*file| cleanupAtomic(file, self.io, self.diag.failure);
    }

    fn flushRecords(self: *@This()) !void {
        if (self.buffered == 0) return;
        if (self.atomic_file == null) {
            self.atomic_file = self.dir.createFileAtomic(self.io, self.path, .{ .replace = true }) catch |err| return self.diag.fail(.{ .@"io.create_failed" = {} }, .create, .zfi, err, error.ZfiWriteFailed);
            if (self.diag.failure) |failure| failure.outputs.file.state = .staged;
        }
        self.atomic_file.?.file.writePositionalAll(self.io, self.records[0..self.buffered], @sizeOf(ZfiHeader) + self.record_bytes - self.buffered) catch |err| return self.diag.fail(.{ .@"io.write_failed" = {} }, .write, .zfi, err, error.ZfiWriteFailed);
        self.buffered = 0;
    }

    fn openRecordsReader(self: *@This()) !std.Io.File {
        const atomic = &self.atomic_file.?;
        std.debug.assert(atomic.file_exists);
        const basename = std.fmt.hex(atomic.file_basename_hex);
        return atomic.dir.openFile(self.io, &basename, .{ .mode = .read_only }) catch |err| return self.diag.fail(.{ .@"io.open_failed" = {} }, .open, .catalog_spool, err, error.ZfiWriteFailed);
    }

    fn replayNames(context: *anyopaque, names: *SourceDedup) !void {
        const self: *@This() = @ptrCast(@alignCast(context));
        if (self.write_failed) return;
        self.flushRecords() catch {
            self.write_failed = true;
            return;
        };
        const read_file = self.openRecordsReader() catch {
            self.write_failed = true;
            return;
        };
        defer read_file.close(self.io);
        self.names.retain(names, read_file) catch |err| {
            if (err == error.OutOfMemory) return err;
            self.write_failed = true;
        };
    }

    fn emit(self: *@This(), info: FastaRecordEmit) !void {
        if (self.write_failed) return;
        self.emitInner(info) catch |err| {
            if (err == error.OutOfMemory) return err;
            self.write_failed = true;
        };
    }

    fn emitInner(self: *@This(), info: FastaRecordEmit) !void {
        var record = info.record;
        record.name_offset = self.names.raw_len;
        record.name_len = @intCast(info.name.len);
        record.flags_and_offset[0] |= index_format.NAME_IN_ZFI_FLAG;
        try self.names.append(info.name);
        if (!info.uses_uniform_formula) {
            try self.diag.markNonUniform(&record, self.sides.len, .local);
            try self.sides.append(info.side_table);
        }
        if (self.buffered + @sizeOf(IndexRecord) > self.records.len) try self.flushRecords();
        @memcpy(self.records[self.buffered..][0..@sizeOf(IndexRecord)], std.mem.asBytes(&record));
        self.buffered += @sizeOf(IndexRecord);
        self.record_bytes += @sizeOf(IndexRecord);
    }

    fn finish(self: *@This(), count: u32, size: u64, mtime: u64) !void {
        if (self.write_failed) return error.ZfiWriteFailed;
        try self.flushRecords();
        const file = self.atomic_file.?.file;
        const header = ZfiHeader{ .magic = ZFI_MAGIC, .record_count = count, .source_size = size };
        file.writePositionalAll(self.io, std.mem.asBytes(&header), 0) catch |err| return self.diag.fail(.{ .@"io.write_failed" = {} }, .write, .zfi, err, error.ZfiWriteFailed);
        if (self.sides.len != 0) {
            const read_file = try self.openRecordsReader();
            defer read_file.close(self.io);
            var records: [1024]IndexRecord = undefined;
            const bytes = std.mem.sliceAsBytes(&records);
            const base = @sizeOf(ZfiHeader) + self.record_bytes;
            var offset: u64 = 0;
            while (offset < self.record_bytes) {
                const wanted: usize = @intCast(@min(bytes.len, self.record_bytes - offset));
                const actual = read_file.readPositionalAll(self.io, bytes[0..wanted], @sizeOf(ZfiHeader) + offset) catch |err| return self.diag.failReadAt(.catalog_spool, @sizeOf(ZfiHeader) + offset, wanted, null, err, error.ReadFailed);
                if (actual != wanted) return self.diag.failReadAt(.catalog_spool, @sizeOf(ZfiHeader) + offset, wanted, actual, null, error.ReadFailed);
                for (records[0 .. wanted / @sizeOf(IndexRecord)]) |*record| {
                    if (!record.isUniformWidth()) try self.diag.markNonUniform(record, base + record.sideTableOffset(), .wire);
                }
                file.writePositionalAll(self.io, bytes[0..wanted], @sizeOf(ZfiHeader) + offset) catch |err| return self.diag.fail(.{ .@"io.write_failed" = {} }, .write, .zfi, err, error.ZfiWriteFailed);
                offset += wanted;
            }
        }
        var buffer: [INDEX_OUTPUT_BUFFER_SIZE]u8 = undefined;
        var output = file.writer(self.io, &buffer);
        output.seekTo(@sizeOf(ZfiHeader) + self.record_bytes) catch |err| return self.diag.fail(.{ .@"io.seek_failed" = .{ .requested_offset = @sizeOf(ZfiHeader) + self.record_bytes } }, .seek, .zfi, if (output.seek_err) |cause| cause else err, error.ZfiWriteFailed);
        self.sides.copyTo(&output.interface) catch |err| {
            self.retainWriteCause(output.err);
            return err;
        };
        self.names.write(&output.interface) catch |err| {
            self.retainWriteCause(output.err);
            return err;
        };
        output.interface.writeAll(&index_format.encodeZfiSourceId(mtime)) catch |err| return self.diag.fail(.{ .@"io.write_failed" = {} }, .write, .zfi, if (output.err) |cause| cause else err, error.ZfiWriteFailed);
        output.interface.writeAll(&index_format.encodeZfiNameFooter(self.names.raw_len)) catch |err| return self.diag.fail(.{ .@"io.write_failed" = {} }, .write, .zfi, if (output.err) |cause| cause else err, error.ZfiWriteFailed);
        output.flush() catch |err| {
            if (self.diag.failure) |failure| failure.outputs.file.final_flush = .failed;
            return self.diag.fail(.{ .@"io.flush_failed" = {} }, .flush, .zfi, if (output.err) |cause| cause else err, error.ZfiWriteFailed);
        };
        if (self.diag.failure) |failure| failure.outputs.file.final_flush = .succeeded;
    }

    fn retainWriteCause(self: *@This(), cause: ?std.Io.File.Writer.Error) void {
        const failure = self.diag.failure orelse return;
        if (failure.has_primary and failure.primary_storage.stage == .write and failure.primary_storage.at.source == .zfi) {
            if (cause) |err| failure.primary_storage.cause = .{ .zig = err };
        }
    }
};

const SourceAdapter = struct {
    source_names: *SourceDedup,

    fn observeAt(self: *@This(), name: []const u8, offset: u64) !bool {
        return self.source_names.observe(name, offset);
    }
};

const FaiNameReplay = struct {
    buffer: *FaiEmitBuffer,
    file: std.Io.File,
    io: std.Io,

    fn replay(context: *anyopaque, names: *SourceDedup) !void {
        const self: *@This() = @ptrCast(@alignCast(context));
        try self.buffer.flush();
        self.buffer.writer.flush() catch |err| return self.buffer.diag.fail(.{ .@"io.flush_failed" = {} }, .flush, .fai_spool, err, error.FaiSpoolWriteFailed);
        var read_buffer: [INDEX_OUTPUT_BUFFER_SIZE]u8 = undefined;
        var input = self.file.reader(self.io, &read_buffer);
        for (0..names.uniqueCount()) |ordinal| {
            const name = input.interface.takeDelimiterInclusive('\t') catch |err| return self.readFailure(&input, err, '\t');
            try names.retainAt(@intCast(ordinal), name[0 .. name.len - 1]);
            _ = input.interface.takeDelimiterInclusive('\n') catch |err| return self.readFailure(&input, err, '\n');
        }
    }

    fn readFailure(self: *@This(), input: *std.Io.File.Reader, err: anyerror, delimiter: u8) anyerror {
        if (err == error.EndOfStream) return self.buffer.diag.fail(.{ .@"io.unexpected_eof" = .{ .delimiter = .{ .delimiter = delimiter } } }, .read, .fai_spool, null, error.FaiSpoolReadFailed);
        return self.buffer.diag.fail(.{ .@"io.read_failed" = .{} }, .read, .fai_spool, if (input.err) |cause| cause else err, error.FaiSpoolReadFailed);
    }
};

const StagedNames = struct {
    raw: CatalogBytes,
    compressed: CompressedNames,
    mode: enum { sampling, raw, compressed } = .sampling,
    sample_lengths: [8]u16 = undefined,
    sample_count: usize = 0,
    estimated_bytes: usize = 0,
    raw_len: u64 = 0,
    count: u32 = 0,

    fn deinit(self: *@This()) void {
        self.raw.deinit();
        self.compressed.deinit();
    }

    fn append(self: *@This(), name: []const u8) !void {
        switch (self.mode) {
            .raw => try self.raw.append(name),
            .compressed => try self.compressed.append(name),
            .sampling => {
                if (self.sample_count == 0) {
                    @memcpy(self.compressed.anchor[0..name.len], name);
                    self.compressed.anchor_len = name.len;
                    self.estimated_bytes = 6 + name.len;
                } else {
                    const anchor = self.compressed.anchor[0..self.compressed.anchor_len];
                    const prefix = commonPrefix(name, anchor);
                    const suffix = commonSuffix(name, anchor, prefix);
                    self.estimated_bytes += 6 + name.len - prefix - suffix;
                }
                self.sample_lengths[self.sample_count] = @intCast(name.len);
                self.sample_count += 1;
                try self.raw.append(name);
                if (self.sample_count == self.sample_lengths.len) {
                    if (self.estimated_bytes * 4 < self.raw.len * 3) {
                        var offset: usize = 0;
                        for (self.sample_lengths) |len| {
                            try self.compressed.append(self.raw.bytes.items[offset..][0..len]);
                            offset += len;
                        }
                        self.raw.deinit();
                        self.raw.bytes = .empty;
                        self.raw.spool = null;
                        self.raw.buffered = 0;
                        self.mode = .compressed;
                    } else {
                        self.mode = .raw;
                    }
                }
            },
        }
        self.raw_len += name.len;
        self.count += 1;
    }

    fn write(self: *@This(), writer: *std.Io.Writer) !void {
        if (self.mode == .compressed) try self.compressed.write(writer) else try self.raw.copyTo(writer);
    }

    fn retain(self: *@This(), names: *SourceDedup, records_file: std.Io.File) !void {
        if (self.mode == .compressed) return self.compressed.retain(names);
        try self.raw.flush();
        var name_buffer: [INDEX_OUTPUT_BUFFER_SIZE]u8 = undefined;
        var record_buffer: [INDEX_OUTPUT_BUFFER_SIZE]u8 = undefined;
        var memory_reader = std.Io.Reader.fixed(self.raw.bytes.items);
        var file_reader: std.Io.File.Reader = undefined;
        const reader = if (self.raw.spool) |spool| blk: {
            file_reader = spool.file.reader(self.raw.io, &name_buffer);
            break :blk &file_reader.interface;
        } else &memory_reader;
        var records = records_file.reader(self.raw.io, &record_buffer);
        records.seekTo(@sizeOf(ZfiHeader)) catch |err| return self.raw.diag.fail(.{ .@"io.seek_failed" = .{ .requested_offset = @sizeOf(ZfiHeader) } }, .seek, .catalog_spool, if (records.seek_err) |cause| cause else err, error.ReadFailed);
        for (0..self.count) |ordinal| {
            const record = records.interface.takeArray(@sizeOf(IndexRecord)) catch |err| {
                if (err == error.EndOfStream) return self.raw.diag.fail(.{ .@"io.unexpected_eof" = .{ .exact = .{ .requested_bytes = @sizeOf(IndexRecord), .received_bytes = records.interface.bufferedLen() } } }, .read, .catalog_spool, null, err);
                return self.raw.diag.fail(.{ .@"io.read_failed" = .{ .requested_bytes = @sizeOf(IndexRecord) } }, .read, .catalog_spool, if (records.err) |cause| cause else err, err);
            };
            const start = @offsetOf(IndexRecord, "name_len");
            const len = std.mem.readInt(u16, record[start..][0..2], .little);
            try names.retainAt(@intCast(ordinal), try self.raw.take(reader, if (self.raw.spool != null) &file_reader else null, len));
        }
    }
};

fn checkPackedNameAllocations(allocator: std.mem.Allocator) !void {
    var failure: diagnostics.Failure = undefined;
    failure.reset(.index);
    const diag = Diagnostics{ .failure = &failure };
    errdefer |err| {
        if (err == error.OutOfMemory) {
            std.testing.expectEqual(diagnostics.Code.@"resource.out_of_memory", failure.primary().?.code()) catch @panic("missing allocation diagnostic");
            std.testing.expectEqual(diagnostics.Stage.allocate, failure.primary().?.stage) catch @panic("wrong allocation stage");
        }
    }

    var names = PackedNames.init(allocator);
    names.diag = diag;
    defer names.deinit();
    try names.ensureCapacity(16);
    var name = [_]u8{'A'} ** MAX_INDEX_NAME_LEN;
    for (0..3) |i| {
        name[0] = @intCast('A' + i);
        try std.testing.expect(!try names.observeAt(&name, 0));
    }
    for (0..300) |i| {
        var buffer: [32]u8 = undefined;
        const short = try std.fmt.bufPrint(&buffer, "name-{d}", .{i});
        try std.testing.expect(!try names.observeAt(short, 0));
    }
    for (0..3) |i| {
        name[0] = @intCast('A' + i);
        try std.testing.expect(try names.observeAt(&name, 0));
    }
}

fn checkSourceNameAllocations(allocator: std.mem.Allocator, data: []const u8) !void {
    var failure: diagnostics.Failure = undefined;
    failure.reset(.index);
    const diag = Diagnostics{ .failure = &failure };
    errdefer |err| {
        if (err == error.OutOfMemory) {
            std.testing.expectEqual(diagnostics.Code.@"resource.out_of_memory", failure.primary().?.code()) catch @panic("missing allocation diagnostic");
            std.testing.expectEqual(diagnostics.Stage.allocate, failure.primary().?.stage) catch @panic("wrong allocation stage");
        }
    }

    var names = SourceDedup.init(allocator, std.testing.io, undefined, data.len);
    names.diag = diag;
    defer names.deinit();
    names.setChunk(data, 0);
    for ([_]bool{ false, true }) |duplicate| {
        var lines = std.mem.tokenizeScalar(u8, data, '\n');
        while (lines.next()) |name| {
            const offset = @intFromPtr(name.ptr) - @intFromPtr(data.ptr);
            try std.testing.expectEqual(duplicate, try names.observe(name, offset));
        }
    }
    try std.testing.expectEqual(@as(u32, 5000), names.uniqueCount());
}

fn checkStagedNameReplay(allocator: std.mem.Allocator, data: []const u8, compressed: bool) !void {
    var failure: diagnostics.Failure = undefined;
    failure.reset(.index);
    const diag = Diagnostics{ .failure = &failure };
    errdefer |err| {
        if (err == error.OutOfMemory) {
            std.testing.expectEqual(diagnostics.Code.@"resource.out_of_memory", failure.primary().?.code()) catch @panic("missing allocation diagnostic");
            std.testing.expectEqual(diagnostics.Stage.allocate, failure.primary().?.stage) catch @panic("wrong allocation stage");
        }
    }

    const io = std.testing.io;
    var temporary = std.testing.tmpDir(.{});
    defer temporary.cleanup();
    const file = try temporary.dir.createFile(io, "source.fa", .{ .read = true });
    defer file.close(io);
    try file.writePositionalAll(io, data, 0);
    var names = SourceDedup.init(allocator, io, file, data.len);
    names.diag = diag;
    defer names.deinit();
    names.cache_full = true;
    var adapter = SourceAdapter{ .source_names = &names };
    const storage = CatalogBytes{ .allocator = allocator, .io = io, .environ = .empty, .diag = diag, .role = .name_spool };
    var catalog = StreamedCatalog{ .io = io, .dir = temporary.dir, .path = "out.zfi", .diag = diag, .names = .{ .raw = storage, .compressed = .{ .encoded = storage } }, .sides = storage };
    catalog.sides.role = .side_spool;
    defer catalog.deinit();
    names.replay_context = &catalog;
    names.replay_fn = StreamedCatalog.replayNames;
    var reader = std.Io.Reader.fixed(data);
    var read_buffer: [31]u8 = undefined;
    const count = try scanFastaReaderWithDedup(
        allocator,
        &reader,
        &read_buffer,
        .{ .enable_dedup = true, .collect_side_tables = true, .failure = &failure },
        @as(?*SourceAdapter, &adapter),
        &catalog,
        StreamedCatalog.emit,
    );
    try std.testing.expectEqual(@as(u32, 128), count);
    try std.testing.expect(names.hydrated);
    try std.testing.expectEqual(compressed, catalog.names.mode == .compressed);

    names.setChunk(&.{}, 0);
    try file.setLength(io, 0);
    var lines = std.mem.splitScalar(u8, data, '\n');
    while (lines.next()) |line| {
        if (line.len != 0 and line[0] == '>') {
            try std.testing.expect(try names.observe(line[1..], 0));
        }
    }
}

fn checkCollectedIndexAllocations(allocator: std.mem.Allocator) !void {
    var failure: diagnostics.Failure = undefined;
    failure.reset(.index);
    errdefer |err| {
        if (err == error.OutOfMemory) {
            std.testing.expectEqual(diagnostics.Code.@"resource.out_of_memory", failure.primary().?.code()) catch @panic("missing collected-index allocation diagnostic");
            std.testing.expectEqual(error.OutOfMemory, failure.primary().?.cause.?.zig) catch @panic("wrong allocation cause");
        }
    }
    const data = ">one\nACGT\nAC\nTTTT\n>two\nA C\n>one\nTT\n";
    var reader = std.Io.Reader.fixed(data);
    var buffer: [7]u8 = undefined;
    var index = try scanZfiReaderWithOptions(allocator, &reader, &buffer, .{ .enable_dedup = true, .failure = &failure });
    defer index.deinit(allocator);
    const bytes = try zfiIndexToBytesChecked(allocator, &index, data.len, 0, &failure);
    defer allocator.free(bytes);
    try std.testing.expectEqual(@as(usize, 2), index.records.items.len);
    try std.testing.expect(failure.primary() == null);
}

test "[property] - [stride scanner]: matches scalar structural boundaries" {
    const line_bases: u32 = 4;
    const line_count = 520;
    const target_line = 256;
    const cases = [_][]const u8{
        "AAAA\n",
        "AAAA\r\n",
    };
    const mutations = [_]u8{ ' ', '\t', '\r', '\n', 0, '>' };

    var storage: [line_count * 6]u8 = undefined;
    for (cases) |line| {
        const line_bytes: u32 = @intCast(line.len);
        const block = storage[0 .. line_count * line.len];
        for (0..line_count) |line_index| {
            @memcpy(block[line_index * line.len ..][0..line.len], line);
        }

        for (0..line.len + 1) |line_offset| {
            const mutation_offset = target_line * line.len + line_offset;
            const original = block[mutation_offset];
            for (mutations) |mutation| {
                block[mutation_offset] = mutation;

                const native = validatedStrideRun(
                    STRIDE_VECTOR_LEN,
                    block,
                    0,
                    line_bases,
                    line_bytes,
                );
                const scalar = reference: {
                    for (0..line_count) |line_index| {
                        const row = block[line_index * line.len ..][0..line.len];
                        if (row[0] == '>') break :reference line_index;
                        for (row[0..line_bases]) |byte| {
                            if (byte <= ' ') break :reference line_index;
                        }
                        if (!std.mem.eql(u8, row[line_bases..], line[line_bases..])) break :reference line_index;
                    }
                    break :reference line_count;
                };

                try std.testing.expectEqual(scalar, native);
                const at_next_line = line_offset == line.len;
                const invalid = if (line_offset < line_bases)
                    mutation <= ' ' or (line_offset == 0 and mutation == '>')
                else if (!at_next_line)
                    mutation != line[line_offset]
                else
                    mutation <= ' ' or mutation == '>';
                const stop_line = target_line + @as(usize, @intFromBool(at_next_line));
                try std.testing.expectEqual(if (invalid) stop_line else line_count, native);
            }
            block[mutation_offset] = original;
        }
    }
}

test "[regression] - [stride scanner]: checks whitespace before a later structural stop" {
    const cases = [_][]const u8{
        "AAAA\nAA A\n>b x\n",
        "AAAA\nAA A\nAAAAx",
        "AAAA\r\nAA A\r\n>b x\r\n",
        "AAAA\r\nAA A\r\nAAAAx\n",
    };
    for (cases) |data| {
        const line_bytes: u32 = if (data[4] == '\r') 6 else 5;
        try std.testing.expectEqual(@as(usize, 1), validatedStrideRun(1, data, 0, 4, line_bytes));
        try std.testing.expectEqual(@as(usize, 1), validatedStrideRun(STRIDE_VECTOR_LEN, data, 0, 4, line_bytes));
    }
}

test "[property] - [FAI formatter]: matches standard formatting at integer boundaries" {
    const cases = [_]IndexRecord{
        .{
            .name_offset = 0,
            .name_len = 0,
            .seq_len = 0,
            .seq_offset = 0,
            .line_bases = 0,
            .line_bytes = 0,
        },
        .{
            .name_offset = 0,
            .name_len = 0,
            .seq_len = std.math.maxInt(u64),
            .seq_offset = std.math.maxInt(u64),
            .line_bases = std.math.maxInt(u32),
            .line_bytes = std.math.maxInt(u32),
        },
        .{
            .name_offset = 0,
            .name_len = 0,
            .seq_len = 0,
            .seq_offset = std.math.maxInt(u64),
            .line_bases = 0,
            .line_bytes = std.math.maxInt(u32),
        },
        .{
            .name_offset = 0,
            .name_len = 0,
            .seq_len = std.math.maxInt(u64),
            .seq_offset = 0,
            .line_bases = std.math.maxInt(u32),
            .line_bytes = 0,
        },
    };

    for (cases) |rec| {
        var actual_buf: [FAI_SUFFIX_BUFFER_SIZE]u8 = undefined;
        const actual = formatFaiSuffix(&actual_buf, rec);
        var expected_buf: [128]u8 = undefined;
        const expected = try std.fmt.bufPrint(
            &expected_buf,
            "\t{d}\t{d}\t{d}\t{d}\n",
            .{ rec.seq_len, rec.seq_offset, rec.line_bases, rec.line_bytes },
        );

        try std.testing.expectEqualStrings(expected, actual);
    }
}

test "[unit] - [FASTA scanner]: gates initial-header validation by option" {
    const data = "not-fasta\n";
    var read_buf: [4]u8 = undefined;
    var output = std.Io.Writer.Allocating.init(std.testing.allocator);
    defer output.deinit();

    var unchecked_reader = std.Io.Reader.fixed(data);
    try std.testing.expectEqual(
        @as(u32, 0),
        try scanFaiReader(
            std.testing.allocator,
            &unchecked_reader,
            &read_buf,
            &output.writer,
            .{ .enable_dedup = true },
        ),
    );

    var checked_reader = std.Io.Reader.fixed(data);
    try std.testing.expectError(
        error.NotFasta,
        scanFaiReader(
            std.testing.allocator,
            &checked_reader,
            &read_buf,
            &output.writer,
            .{ .enable_dedup = true, .require_initial_header = true },
        ),
    );
}

test "[failure] - [FASTA scanner]: rejects sequence lines above u32 geometry" {
    var state = ChunkParseState{
        .allocator = std.testing.allocator,
        .in_pending_line = true,
        .pending_line_bases = std.math.maxInt(u32),
    };
    try std.testing.expectError(error.SequenceLineTooLong, state.appendPendingFragment("A", 0));

    state.pending_line_bases = 0;
    state.pending_line_bytes = std.math.maxInt(u32);
    const Ctx = struct {
        fn emit(_: *@This(), _: FastaRecordEmit) !void {
            unreachable;
        }
    };
    var ctx = Ctx{};
    var record_count: u32 = 0;
    try std.testing.expectError(
        error.SequenceLineTooLong,
        processChunkBytes(&state, "\n", 0, null, null, &ctx, Ctx.emit, &record_count),
    );

    state.pending_line_bases = 1;
    try std.testing.expectError(error.SequenceLineTooLong, state.commitPendingLine(false));
}

test "[edge] - [index count]: rejects excess records before emitting" {
    var state = ChunkParseState{ .allocator = std.testing.allocator };
    defer state.name.deinit(std.testing.allocator);
    defer state.side_table.deinit(std.testing.allocator);
    try state.name.appendSlice(std.testing.allocator, "limit");
    try state.recordSequenceLine(1, 2, true, 1, 7);
    const Ctx = struct {
        count: usize = 0,
        fn emit(self: *@This(), _: FastaRecordEmit) !void {
            self.count += 1;
        }
    };
    var ctx = Ctx{};
    var count: u32 = std.math.maxInt(u32) - 1;
    try state.emitRecordIfReady(&ctx, null, Ctx.emit, &count);
    try std.testing.expectEqual(std.math.maxInt(u32), count);
    try std.testing.expectError(error.TooManySequences, state.emitRecordIfReady(&ctx, null, Ctx.emit, &count));
    try std.testing.expectEqual(@as(usize, 1), ctx.count);
    try std.testing.expectEqual(std.math.maxInt(u32), count);

    var seen = NameDedup.init(std.testing.allocator);
    defer seen.deinit();
    try std.testing.expect(!try seen.observe("limit"));
    try state.emitRecordIfReady(&ctx, @as(?*NameDedup, &seen), Ctx.emit, &count);
    state.line_builder = .{};
    try state.emitRecordIfReady(&ctx, null, Ctx.emit, &count);
    try std.testing.expectEqual(@as(usize, 1), ctx.count);
    try std.testing.expectEqual(std.math.maxInt(u32), count);
}

test "[edge] - [index count]: full source tables retain duplicate lookups" {
    inline for (.{ u32, u64 }) |Offset| {
        var failing = std.testing.FailingAllocator.init(std.testing.allocator, .{});
        var source = SourceDedup.init(failing.allocator(), std.testing.io, undefined, if (Offset == u32) 0 else std.math.maxInt(u32));
        defer source.deinit();
        source.cache_full = true;
        source.setChunk("alpha\nbeta\n", 0);
        const set = if (Offset == u32) &source.tables.narrow else &source.tables.wide;
        try std.testing.expect(!try set.observeWithFingerprint(&source, "alpha", 0, 0));
        set.len = std.math.maxInt(u32);
        failing.fail_index = failing.alloc_index;
        try std.testing.expect(try set.observeWithFingerprint(&source, "alpha", 0, 0));
        try std.testing.expectError(error.TooManySequences, set.observeWithFingerprint(&source, "beta", 6, 0));
        try std.testing.expect(!failing.has_induced_failure);
        try std.testing.expectEqual(std.math.maxInt(u32), set.len);
    }
}

test "[property] - [name deduplication]: remains exact under hash collisions" {
    const CollisionContext = struct {
        pub fn hash(_: @This(), _: []const u8) u64 {
            return 0;
        }

        pub fn eql(_: @This(), a: []const u8, b: []const u8) bool {
            return std.mem.eql(u8, a, b);
        }
    };
    const CollidingDedup = NameDedupWith(CollisionContext);
    const names = [_][]const u8{ "alpha", "beta", "alpha", "gamma", "beta" };
    var seen = CollidingDedup.init(std.testing.allocator);
    defer seen.deinit();

    var kept: usize = 0;
    for (names) |name| {
        if (!(try seen.observe(name))) kept += 1;
    }

    try std.testing.expectEqual(@as(usize, 3), kept);
    try std.testing.expectEqual(@as(usize, 3), seen.map.count());
}

test "[failure] - [FASTA scanner]: maps reader failures to SourceReadFailed" {
    var reader: std.Io.Reader = .failing;
    var read_buf: [4]u8 = undefined;
    var output = std.Io.Writer.Allocating.init(std.testing.allocator);
    defer output.deinit();

    try std.testing.expectError(
        error.SourceReadFailed,
        scanFaiReader(
            std.testing.allocator,
            &reader,
            &read_buf,
            &output.writer,
            .{ .enable_dedup = true },
        ),
    );
}

test "[failure] - [FAI scanner]: maps writer failures to FaiSpoolWriteFailed" {
    var reader = std.Io.Reader.fixed(">seq\nA\n");
    var read_buf: [4]u8 = undefined;
    var writer: std.Io.Writer = .failing;

    try std.testing.expectError(
        error.FaiSpoolWriteFailed,
        scanFaiReader(
            std.testing.allocator,
            &reader,
            &read_buf,
            &writer,
            .{ .enable_dedup = true },
        ),
    );
}

test "[integration] - [FAI spool]: uses exclusive files and cleans them" {
    var tmp = std.testing.tmpDir(.{ .iterate = true });
    defer tmp.cleanup();

    {
        var first = tryCreateFaiSpool(std.testing.io, tmp.dir, false) orelse
            return error.TestUnexpectedResult;
        defer first.deinit(std.testing.io);
        var second = tryCreateFaiSpool(std.testing.io, tmp.dir, false) orelse
            return error.TestUnexpectedResult;
        defer second.deinit(std.testing.io);

        try std.testing.expect(!std.mem.eql(
            u8,
            first.name[0..first.name_len],
            second.name[0..second.name_len],
        ));
    }

    var entries = tmp.dir.iterate();
    try std.testing.expectEqual(null, try entries.next(std.testing.io));
}

test "[integration] - [source identity]: detects path replacement" {
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();

    const original = try tmp.dir.createFile(std.testing.io, "source.fa", .{});
    try std.Io.File.writeStreamingAll(original, std.testing.io, ">seq\nAAAA\n");
    const before = try original.stat(std.testing.io);
    original.close(std.testing.io);

    try tmp.dir.rename("source.fa", tmp.dir, "old.fa", std.testing.io);
    const replacement = try tmp.dir.createFile(std.testing.io, "source.fa", .{});
    try std.Io.File.writeStreamingAll(replacement, std.testing.io, ">seq\nA\n");
    replacement.close(std.testing.io);

    var failure: diagnostics.Failure = undefined;
    failure.reset(.index);
    try std.testing.expectError(error.InputChanged, checkSourceUnchanged(std.testing.io, tmp.dir, "source.fa", before, .{ .failure = &failure }));
    const details = failure.primary().?.details.@"input.changed";
    try std.testing.expect(details.changed_fields.size or details.changed_fields.mtime);
    try std.testing.expectEqual(before.size, details.before.bytes.?);
    try std.testing.expectEqual(before.mtime.nanoseconds, details.before.mtime_ns.?);
    try std.testing.expect(details.after != null);
    try std.testing.expect(!details.disappeared);
}

test "[property] - [source names]: complete names resolve fingerprint collisions" {
    var temporary = std.testing.tmpDir(.{});
    defer temporary.cleanup();
    const io = std.testing.io;
    const file = try temporary.dir.createFile(io, "names.fa", .{ .read = true });
    defer file.close(io);
    const long = "a-much-longer-different-name";
    const data = long ++ "\nalpha\nalphabet\nbeta\nalpha\n";
    try file.writePositionalAll(io, data, 0);
    var source = SourceDedup.init(std.testing.allocator, io, file, data.len);
    defer source.deinit();
    source.cache_full = true;
    const set = &source.tables.narrow;
    const start = long.len + 1;
    try std.testing.expect(!try set.observeWithFingerprint(&source, "alpha", start, 0));
    try std.testing.expect(!try set.observeWithFingerprint(&source, "alphabet", start + 6, 0));
    try std.testing.expect(!try set.observeWithFingerprint(&source, "beta", start + 15, 0));
    try std.testing.expect(try set.observeWithFingerprint(&source, "alpha", start + 20, 0));
    try std.testing.expect(!try set.observeWithFingerprint(&source, long, 0, 0));
}

test "[unit] - [source names]: exact comparisons use the current input chunk" {
    var temporary = std.testing.tmpDir(.{});
    defer temporary.cleanup();
    const io = std.testing.io;
    const file = try temporary.dir.createFile(io, "names.fa", .{ .read = true });
    defer file.close(io);
    const data = "alpha\nalphabet\nbeta\nalpha\n";
    var source = SourceDedup.init(std.testing.allocator, io, file, data.len);
    defer source.deinit();
    source.cache_full = true;
    source.setChunk(data, 0);
    const set = &source.tables.narrow;
    try std.testing.expect(!try set.observeWithFingerprint(&source, "alpha", 0, 0));
    try std.testing.expect(!try set.observeWithFingerprint(&source, "alphabet", 6, 0));
    try std.testing.expect(try set.observeWithFingerprint(&source, "alpha", 20, 0));
}

test "[edge] - [source names]: retain offsets above four GiB" {
    var temporary = std.testing.tmpDir(.{});
    defer temporary.cleanup();
    const io = std.testing.io;
    const file = try temporary.dir.createFile(io, "names.fa", .{ .read = true });
    defer file.close(io);
    const offset: u64 = @as(u64, std.math.maxInt(u32)) + 64;
    const data = "alpha\nalphabet\nbeta\nalpha\n";
    var source = SourceDedup.init(std.testing.allocator, io, file, offset + data.len);
    defer source.deinit();
    source.cache_full = true;
    source.setChunk(data, offset);
    const set = &source.tables.wide;
    try std.testing.expect(!try set.observeWithFingerprint(&source, "alpha", offset, 0));
    try std.testing.expect(!try set.observeWithFingerprint(&source, "alphabet", offset + 6, 0));
    try std.testing.expect(try set.observeWithFingerprint(&source, "alpha", offset + 20, 0));
    try std.testing.expectEqual(offset, set.keyAt(0).offset);
    try std.testing.expectEqual(@as(u64, 0), (try file.stat(io)).size);
}

test "[property] - [source names]: cached names resolve fingerprint collisions" {
    var temporary = std.testing.tmpDir(.{});
    defer temporary.cleanup();
    const io = std.testing.io;
    const file = try temporary.dir.createFile(io, "names.fa", .{ .read = true });
    defer file.close(io);
    const data = "alpha\nalphabet\nbeta\nalpha\n";
    try file.writePositionalAll(io, data, 0);
    var source = SourceDedup.init(std.testing.allocator, io, file, data.len);
    defer source.deinit();
    const set = &source.tables.narrow;
    try std.testing.expect(!try set.observeWithFingerprint(&source, "alpha", 0, 0));
    try std.testing.expect(!try set.observeWithFingerprint(&source, "alphabet", 6, 0));
    try std.testing.expect(!try set.observeWithFingerprint(&source, "beta", 15, 0));
    try std.testing.expect(try set.observeWithFingerprint(&source, "alpha", 20, 0));
}

test "[failure] - [source names]: read errors leave the insertion unfinished" {
    var temporary = std.testing.tmpDir(.{});
    defer temporary.cleanup();
    const io = std.testing.io;
    const file = try temporary.dir.createFile(io, "names.fa", .{});
    defer file.close(io);
    const data = "alpha\nbeta\n";
    try file.writePositionalAll(io, data, 0);
    var source = SourceDedup.init(std.testing.allocator, io, file, data.len);
    defer source.deinit();
    source.cache_full = true;
    const set = &source.tables.narrow;
    try std.testing.expect(!try set.observeWithFingerprint(&source, "alpha", 0, 0));
    try std.testing.expectError(error.SourceReadFailed, set.observeWithFingerprint(&source, "beta", 6, 0));
    try std.testing.expectEqual(@as(u32, 1), set.count());
}

test "[property] - [portable scanner]: whitespace counts match scalar counts" {
    var data: [12289]u8 = undefined;
    const lengths = [_]usize{ 0, 15, 16, 31, 32, 127, 128, 2047, 2048, 2049, 4095, 4096, 4097, 8191, 8192, 8193, data.len };
    for ([_]?u8{ null, ' ', 'X' }) |fill| {
        for (&data, 0..) |*byte, i| byte.* = fill orelse @as(u8, @truncate(i * 17 + i / 31));
        for (lengths) |len| {
            var expected: usize = 0;
            for (data[0..len]) |byte| expected += @intFromBool(byte <= ' ');
            inline for (.{ 1, 8, 16, 32, 64 }) |width| {
                try std.testing.expectEqual(expected, countWhitespace(width, data[0..len]));
            }
        }
    }
}

test "[property] - [portable scanner]: name delimiters match scalar scanning" {
    var data = [_]u8{'A'} ** 65;
    for ([_]u8{ ' ', '\t', '\r', '\n' }) |separator| {
        for ([_]usize{ 0, 1, 15, 16, 17, 31, 32, 63, 64 }) |at| {
            data[at] = separator;
            for ([_]usize{ 0, 1, 15, 16, 32, data.len }) |start| {
                try std.testing.expectEqual(if (start <= at) at else data.len, findNameDelimiter(&data, start));
            }
            data[at] = 'A';
        }
    }
}

test "[failure] - [packed names]: block growth cleans every allocation failure" {
    try std.testing.checkAllAllocationFailures(std.testing.allocator, checkPackedNameAllocations, .{});
}

test "[failure] - [source names]: table and entry growth clean every allocation failure" {
    var data: std.ArrayList(u8) = .empty;
    defer data.deinit(std.testing.allocator);
    for (0..5000) |i| {
        var buffer: [32]u8 = undefined;
        try data.appendSlice(std.testing.allocator, try std.fmt.bufPrint(&buffer, "lookup-{d:0>8}\n", .{i}));
    }
    try std.testing.checkAllAllocationFailures(std.testing.allocator, checkSourceNameAllocations, .{data.items});
}

test "[edge] - [source names]: maximum names cross copied and retained cache blocks" {
    var temporary = std.testing.tmpDir(.{});
    defer temporary.cleanup();
    const io = std.testing.io;
    const file = try temporary.dir.createFile(io, "names.fa", .{ .read = true });
    defer file.close(io);
    var name: [MAX_INDEX_NAME_LEN]u8 = undefined;
    const stride = MAX_INDEX_NAME_LEN + 1;
    for (0..6) |i| {
        @memset(&name, @intCast('A' + i));
        try file.writePositionalAll(io, &name, i * stride);
        try file.writePositionalAll(io, "\n", i * stride + name.len);
    }
    var names = SourceDedup.init(std.testing.allocator, io, file, 6 * stride);
    defer names.deinit();
    for ([_]bool{ false, true, true }) |duplicate| {
        for (0..6) |i| {
            @memset(&name, @intCast('A' + i));
            try std.testing.expectEqual(duplicate, try names.observe(&name, i * stride));
        }
    }
    try std.testing.expect(names.cache_full);
    try std.testing.expect(names.hot_blocks.items.len >= 2);
}

test "[failure] - [source names]: short reads distinguish collisions from truncation" {
    var temporary = std.testing.tmpDir(.{});
    defer temporary.cleanup();
    const io = std.testing.io;
    const file = try temporary.dir.createFile(io, "names.fa", .{ .read = true });
    defer file.close(io);
    try file.writePositionalAll(io, "alpha\n", 0);
    var names = SourceDedup.init(std.testing.allocator, io, file, 6);
    defer names.deinit();
    names.cache_full = true;
    const flag = Tables(u32).CACHE_FLAG;
    try std.testing.expect(!names.sameName("alphabet", 0, flag));
    try std.testing.expect(!names.read_failed);
    try file.setLength(io, 3);
    try std.testing.expect(!names.sameName("alpha", 0, flag));
    try std.testing.expect(names.read_failed);
}

test "[integration] - [catalog staging]: spill replay preserves bytes and removes its file" {
    const allocator = std.testing.allocator;
    const io = std.testing.io;
    const data = try allocator.alloc(u8, CatalogBytes.MEMORY_BYTES + 333);
    defer allocator.free(data);
    for (data, 0..) |*byte, i| byte.* = @truncate(i);
    var spool_dir: ?std.Io.Dir = null;
    defer if (spool_dir) |dir| dir.close(io);
    var spool_name: [FAI_SPOOL_NAME_BUFFER_SIZE]u8 = undefined;
    var spool_name_len: usize = 0;
    {
        var storage = CatalogBytes{ .allocator = allocator, .io = io, .environ = .empty };
        defer storage.deinit();
        try storage.append(data[0..CatalogBytes.MEMORY_BYTES]);
        try std.testing.expectEqual(null, storage.spool);
        try storage.append(data[CatalogBytes.MEMORY_BYTES..]);
        const spool = storage.spool.?;
        spool_dir = try spool.dir.openDir(io, ".", .{});
        spool_name_len = spool.name_len;
        @memcpy(spool_name[0..spool_name_len], spool.name[0..spool_name_len]);
        var output = std.Io.Writer.Allocating.init(allocator);
        defer output.deinit();
        try storage.copyTo(&output.writer);
        try std.testing.expectEqualSlices(u8, data, output.written());
    }
    try std.testing.expectError(error.FileNotFound, spool_dir.?.statFile(io, spool_name[0..spool_name_len], .{}));
}

test "[property] - [source names]: raw and compressed replay retain exact names" {
    const allocator = std.testing.allocator;
    for ([_]bool{ false, true }) |compressed| {
        var data: std.ArrayList(u8) = .empty;
        defer data.deinit(allocator);
        for (0..2) |_| {
            for (0..128) |i| {
                var header: [64]u8 = undefined;
                const name = if (compressed)
                    try std.fmt.bufPrint(&header, ">ENST000000{d:0>6}.1\n", .{i})
                else
                    try std.fmt.bufPrint(&header, ">name-{d}\n", .{i});
                try data.appendSlice(allocator, name);
                try data.appendSlice(allocator, if (i % 2 == 0) "ACGT\nAC\n" else "ACGT\nAC\nTTTT\n");
            }
        }
        try std.testing.checkAllAllocationFailures(allocator, checkStagedNameReplay, .{ data.items, compressed });
    }
}

test "[failure] - [streamed ZFI]: parse errors precede record write failures" {
    const allocator = std.testing.allocator;
    const io = std.testing.io;
    var data: std.ArrayList(u8) = .empty;
    defer data.deinit(allocator);
    for (0..2000) |_| try data.appendSlice(allocator, ">name\nACGT\n");
    try data.append(allocator, '>');
    try data.appendNTimes(allocator, 'A', MAX_INDEX_NAME_LEN + 1);
    try data.appendSlice(allocator, "\nACGT\n");
    var temporary = std.testing.tmpDir(.{});
    defer temporary.cleanup();
    const storage = CatalogBytes{ .allocator = allocator, .io = io, .environ = .empty };
    var catalog = StreamedCatalog{
        .io = io,
        .dir = temporary.dir,
        .path = "missing/out.zfi",
        .names = .{ .raw = storage, .compressed = .{ .encoded = storage } },
        .sides = storage,
    };
    defer catalog.deinit();
    var reader = std.Io.Reader.fixed(data.items);
    var buffer: [8192]u8 = undefined;
    try std.testing.expectError(error.HeaderTooLong, scanFastaReader(
        allocator,
        &reader,
        &buffer,
        .{ .enable_dedup = false },
        null,
        &catalog,
        StreamedCatalog.emit,
    ));
    try std.testing.expect(catalog.write_failed);
    try std.testing.expectError(error.ZfiWriteFailed, catalog.finish(2000, data.items.len, 0));
}

test "[failure] - [source names]: source offsets cannot overlap cache flags" {
    var narrow = SourceDedup.init(std.testing.allocator, std.testing.io, undefined, 0);
    defer narrow.deinit();
    try std.testing.expectError(error.SourceChanged, narrow.observe("name", Tables(u32).CACHE_FLAG));
    try std.testing.expectEqual(@as(u32, 0), narrow.uniqueCount());
    var wide = SourceDedup.init(std.testing.allocator, std.testing.io, undefined, std.math.maxInt(u32));
    defer wide.deinit();
    try std.testing.expectError(error.SourceChanged, wide.observe("name", Tables(u64).CACHE_FLAG));
    try std.testing.expectEqual(@as(u32, 0), wide.uniqueCount());
}

test "[failure] - [packed names]: oversized offsets roll back the insertion" {
    if (comptime @bitSizeOf(usize) <= 32) return error.SkipZigTest;
    var names = PackedNames.init(std.testing.allocator);
    defer names.deinit();
    names.position = @as(usize, std.math.maxInt(u32)) + 1;
    try std.testing.expectError(error.OutOfMemory, names.observeAt("name", 0));
    try std.testing.expectEqual(@as(u32, 0), names.map.count());
}

test "[failure] - [index diagnostics]: retains scanner limit operands and source positions" {
    var failure: diagnostics.Failure = undefined;
    var state = ChunkParseState{ .allocator = std.testing.allocator, .diag = .{ .failure = &failure } };
    defer state.name.deinit(std.testing.allocator);
    defer state.side_table.deinit(std.testing.allocator);
    const max = std.math.maxInt(u32);
    for (0..3) |i| {
        failure.reset(.index);
        state.in_pending_line = true;
        state.pending_line_file_offset = 700;
        state.pending_line_bases = if (i == 0) max else 0;
        state.pending_line_bytes = if (i == 0) 0 else max;
        if (i == 2) {
            try std.testing.expectError(error.InvalidIndex, state.commitPendingLine(false));
        } else {
            try std.testing.expectError(error.InvalidIndex, state.appendPendingFragment("A", 700));
        }
        const primary = failure.primary().?;
        const details = primary.details.@"index.sequence_line_too_long";
        try std.testing.expectEqual(@as(u128, max) + 1, details.observed.exact);
        try std.testing.expectEqual(if (i == 0) diagnostics.IndexField.line_bases else .line_bytes, details.field);
        try std.testing.expectEqual(if (i == 2) diagnostics.LineWidthRole.synthetic_separator else .index_geometry, details.width_role);
        try std.testing.expectEqual(@as(?u64, 700), primary.at.location.offset);
    }
    failure.reset(.index);
    state.in_pending_line = false;
    state.pending_line_bytes = 0;
    state.pending_line_bases = 0;
    try state.startHeader(42);
    try state.name.appendSlice(std.testing.allocator, "name");
    try state.recordSequenceLine(1, 2, true, 1, 48);
    const Ctx = struct {
        fn emit(_: *@This(), _: FastaRecordEmit) !void {
            unreachable;
        }
    };
    var ctx = Ctx{};
    var count: u32 = max;
    try std.testing.expectError(error.InvalidIndex, state.emitRecordIfReady(&ctx, null, Ctx.emit, &count));
    const primary = failure.primary().?;
    try std.testing.expectEqual(@as(u128, max) + 1, primary.details.@"index.record_count_limit_exceeded".record_count.exact);
    try std.testing.expectEqual(@as(?u64, 42), primary.at.location.offset);
    try std.testing.expectEqual(@as(?u64, 1), primary.at.record);
}

test "[failure] - [index diagnostics]: cache representation limits are capacity failures" {
    var failure: diagnostics.Failure = undefined;
    inline for (.{ u32, u64 }) |Offset| {
        failure.reset(.index);
        var names = SourceDedup.init(std.testing.allocator, std.testing.io, undefined, if (Offset == u32) 0 else std.math.maxInt(u32));
        names.diag = .{ .failure = &failure };
        defer names.deinit();
        const limit = Tables(Offset).CACHE_FLAG;
        try std.testing.expectError(error.CapacityExceeded, names.observe("name", limit));
        const details = failure.primary().?.details.@"resource.capacity_limit_exceeded";
        try std.testing.expectEqual(diagnostics.CapacityRole.name_cache_offset, details.capacity_role);
        try std.testing.expectEqual(@as(u128, limit - 1), details.limit);
        try std.testing.expectEqual(@as(u128, limit), details.required.exact);
    }
    if (comptime @bitSizeOf(usize) > 32) {
        failure.reset(.index);
        var names = PackedNames.init(std.testing.allocator);
        names.diag = .{ .failure = &failure };
        defer names.deinit();
        names.position = @as(usize, std.math.maxInt(u32)) + 1;
        try std.testing.expectError(error.CapacityExceeded, names.observeAt("name", 0));
        try std.testing.expectEqual(diagnostics.CapacityRole.name_offset, failure.primary().?.details.@"resource.capacity_limit_exceeded".capacity_role);
        try std.testing.expectEqual(@as(u32, 0), names.map.count());
    }
}

test "[failure] - [index diagnostics]: parser failure promotes over deferred catalog cause" {
    const allocator = std.testing.allocator;
    const io = std.testing.io;
    var data: std.ArrayList(u8) = .empty;
    defer data.deinit(allocator);
    for (0..2000) |_| try data.appendSlice(allocator, ">name\nACGT\n");
    try data.append(allocator, '>');
    try data.appendNTimes(allocator, 'A', MAX_INDEX_NAME_LEN + 1);
    var temporary = std.testing.tmpDir(.{});
    defer temporary.cleanup();
    var failure: diagnostics.Failure = undefined;
    failure.reset(.index);
    const diag = Diagnostics{ .failure = &failure, .deferred = true, .at = .{} };
    const storage = CatalogBytes{ .allocator = allocator, .io = io, .environ = .empty, .diag = diag };
    var catalog = StreamedCatalog{ .io = io, .dir = temporary.dir, .path = "missing/out.zfi", .diag = diag, .names = .{ .raw = storage, .compressed = .{ .encoded = storage } }, .sides = storage };
    defer catalog.deinit();
    var reader = std.Io.Reader.fixed(data.items);
    var buffer: [8192]u8 = undefined;
    try std.testing.expectError(error.InvalidIndex, scanFastaReader(allocator, &reader, &buffer, .{ .enable_dedup = false, .failure = &failure }, null, &catalog, StreamedCatalog.emit));
    try std.testing.expectEqual(diagnostics.Code.@"index.identifier_too_long", failure.primary().?.code());
    try std.testing.expectEqual(@as(?u64, 2001), failure.primary().?.at.record);
    try std.testing.expectEqual(@as(?u64, 22000), failure.primary().?.at.location.offset);
    try std.testing.expectEqual(@as(usize, 1), failure.secondaryItems().len);
    try std.testing.expectEqual(diagnostics.Code.@"io.create_failed", failure.secondaryItems()[0].code());
    try std.testing.expectEqual(error.FileNotFound, failure.secondaryItems()[0].cause.?.zig);
}

test "[failure] - [catalog diagnostics]: spill IO retains causes and exact short reads" {
    const Fault = struct {
        fn write(_: ?*anyopaque, _: std.Io.File, _: []const u8, _: []const []const u8, _: usize, _: u64) std.Io.File.WritePositionalError!usize {
            return error.NoSpaceLeft;
        }
        fn read(_: ?*anyopaque, _: std.Io.File, _: []const []u8, _: u64) std.Io.File.ReadPositionalError!usize {
            return error.InputOutput;
        }
        fn short(_: ?*anyopaque, _: std.Io.File, _: []const []u8, _: u64) std.Io.File.ReadPositionalError!usize {
            return 0;
        }
        fn laterRead(userdata: ?*anyopaque, file: std.Io.File, data: []const []u8, offset: u64) std.Io.File.ReadPositionalError!usize {
            if (offset != 0) return error.InputOutput;
            return std.testing.io.vtable.fileReadPositional(userdata, file, data, offset);
        }
    };
    const allocator = std.testing.allocator;
    const io = std.testing.io;
    const data = try allocator.alloc(u8, CatalogBytes.MEMORY_BYTES + 1);
    defer allocator.free(data);
    @memset(data, 'A');
    for (0..4) |i| {
        var failure: diagnostics.Failure = undefined;
        failure.reset(.index);
        var storage = CatalogBytes{ .allocator = allocator, .io = io, .environ = .empty, .diag = .{ .failure = &failure }, .role = .name_spool };
        defer storage.deinit();
        var vtable = io.vtable.*;
        if (i == 0) {
            vtable.fileWritePositional = Fault.write;
            storage.io = .{ .userdata = io.userdata, .vtable = &vtable };
            try std.testing.expectError(error.Io, storage.append(data));
            try std.testing.expectEqual(diagnostics.Code.@"io.write_failed", failure.primary().?.code());
            try std.testing.expectEqual(error.NoSpaceLeft, failure.primary().?.cause.?.zig);
        } else {
            try storage.append(data);
            try storage.flush();
            vtable.fileReadPositional = switch (i) {
                1 => Fault.read,
                2 => Fault.short,
                else => Fault.laterRead,
            };
            storage.io = .{ .userdata = io.userdata, .vtable = &vtable };
            var output = std.Io.Writer.Allocating.init(allocator);
            defer output.deinit();
            try std.testing.expectError(error.Io, storage.copyTo(&output.writer));
            if (i == 1 or i == 3) {
                try std.testing.expectEqual(diagnostics.Code.@"io.read_failed", failure.primary().?.code());
                try std.testing.expectEqual(error.InputOutput, failure.primary().?.cause.?.zig);
            } else {
                const details = failure.primary().?.details.@"io.unexpected_eof".exact;
                try std.testing.expectEqual(@as(u64, INDEX_OUTPUT_BUFFER_SIZE), details.requested_bytes);
                try std.testing.expectEqual(@as(u64, 0), details.received_bytes);
                try std.testing.expectEqual(@as(?diagnostics.Cause, null), failure.primary().?.cause);
            }
            const expected_offset: u64 = if (i == 3) INDEX_OUTPUT_BUFFER_SIZE else 0;
            try std.testing.expectEqual(expected_offset, output.written().len);
            try std.testing.expectEqual(@as(?u64, expected_offset), failure.primary().?.at.location.offset);
        }
        try std.testing.expectEqual(diagnostics.SourceId.name_spool, failure.primary().?.at.source.?);
    }
}

test "[regression] - [FAI replay diagnostics]: missing delimiters retain EOF without a field length" {
    const io = std.testing.io;
    var temporary = std.testing.tmpDir(.{});
    defer temporary.cleanup();
    const file = try temporary.dir.createFile(io, "replay.fai", .{ .read = true });
    defer file.close(io);
    var names = SourceDedup.init(std.testing.allocator, io, undefined, 10);
    defer names.deinit();
    try std.testing.expect(!try names.observe("name", 0));
    for ([_][]const u8{ "name", "name\t4\t5\t4\t5" }, [_]u8{ '\t', '\n' }) |bytes, delimiter| {
        try file.setLength(io, 0);
        try file.writePositionalAll(io, bytes, 0);
        var failure: diagnostics.Failure = undefined;
        failure.reset(.index);
        var output = std.Io.Writer.Allocating.init(std.testing.allocator);
        defer output.deinit();
        var buffer = FaiEmitBuffer{ .writer = &output.writer, .diag = .{ .failure = &failure } };
        var replay = FaiNameReplay{ .buffer = &buffer, .file = file, .io = io };
        try std.testing.expectError(error.Io, FaiNameReplay.replay(&replay, &names));
        const details = failure.primary().?.details.@"io.unexpected_eof".delimiter;
        try std.testing.expectEqual(delimiter, details.delimiter);
        try std.testing.expectEqual(@as(?diagnostics.Measurement, null), details.observed_bytes);
    }
}

test "[failure] - [index diagnostics]: collected catalog and serialization clean every OOM" {
    try std.testing.checkAllAllocationFailures(std.testing.allocator, checkCollectedIndexAllocations, .{});
}

test "[regression] - [index abort]: removed staging does not remain staged" {
    const io = std.testing.io;
    var temporary = std.testing.tmpDir(.{});
    defer temporary.cleanup();
    var failure: diagnostics.Failure = undefined;
    failure.reset(.index);
    var atomic = try temporary.dir.createFileAtomic(io, "out.zfi", .{ .replace = true });
    const basename = std.fmt.hex(atomic.file_basename_hex);
    failure.outputs.file.state = .staged;
    try std.testing.expectEqual(error.Io, failure.fail(.{ .details = .{ .@"io.write_failed" = {} }, .stage = .write, .cause = .{ .zig = error.NoSpaceLeft } }));

    cleanupAtomic(&atomic, io, &failure);

    try std.testing.expectError(error.FileNotFound, temporary.dir.statFile(io, &basename, .{}));
    try std.testing.expectEqual(.aborted, failure.outputs.file.state);
    try std.testing.expectEqual(error.NoSpaceLeft, failure.primary().?.cause.?.zig);
}

test "[failure] - [index abort]: cleanup uncertainty preserves primary and published state" {
    const Fault = struct {
        fn remove(_: ?*anyopaque, _: std.Io.Dir, _: []const u8) std.Io.Dir.DeleteFileError!void {
            return error.AccessDenied;
        }
    };
    const io = std.testing.io;
    var temporary = std.testing.tmpDir(.{});
    defer temporary.cleanup();
    var failure: diagnostics.Failure = undefined;
    failure.reset(.index);
    var atomic = try temporary.dir.createFileAtomic(io, "out.zfi", .{ .replace = true });
    const basename = std.fmt.hex(atomic.file_basename_hex);
    defer temporary.dir.deleteFile(io, &basename) catch {};
    failure.outputs.file.state = .staged;
    try std.testing.expectEqual(error.Io, failure.fail(.{ .details = .{ .@"io.write_failed" = {} }, .stage = .write, .cause = .{ .zig = error.NoSpaceLeft } }));
    var vtable = io.vtable.*;
    vtable.dirDeleteFile = Fault.remove;

    cleanupAtomic(&atomic, .{ .userdata = io.userdata, .vtable = &vtable }, &failure);

    try std.testing.expectEqual(.unknown, failure.outputs.file.state);
    try std.testing.expectEqual(error.NoSpaceLeft, failure.primary().?.cause.?.zig);
    try std.testing.expectEqual(@as(usize, 1), failure.secondaryItems().len);
    try std.testing.expectEqual(error.AccessDenied, failure.secondaryItems()[0].cause.?.zig);
    _ = try temporary.dir.statFile(io, &basename, .{});

    failure.reset(.index);
    atomic = try temporary.dir.createFileAtomic(io, "published.zfi", .{ .replace = true });
    try atomic.file.writePositionalAll(io, "index", 0);
    try atomic.replace(io);
    failure.outputs.file.state = .published;
    cleanupAtomic(&atomic, io, &failure);
    try std.testing.expectEqual(.published, failure.outputs.file.state);
    try std.testing.expect(failure.primary() == null);
    try std.testing.expectEqual(@as(u64, 5), (try temporary.dir.statFile(io, "published.zfi", .{})).size);
}

test "[edge] - [index preflight]: buffer width and ZFI count reject before narrowing" {
    var failure: diagnostics.Failure = undefined;
    failure.reset(.index);
    var reader = std.Io.Reader.fixed(">name\nACGT\n");
    var empty_buffer: [0]u8 = .{};
    var output_bytes: [16]u8 = undefined;
    var output = std.Io.Writer.fixed(&output_bytes);
    var failing = std.testing.FailingAllocator.init(std.testing.allocator, .{ .fail_index = 0 });
    try std.testing.expectError(error.InvalidArgument, scanFaiReader(failing.allocator(), &reader, &empty_buffer, &output, .{ .enable_dedup = true, .failure = &failure }));
    try std.testing.expectEqual(@as(usize, 0), output.buffered().len);
    try std.testing.expectEqual(@as(usize, 0), reader.seek);
    try std.testing.expectEqual(@as(usize, 0), failing.allocations);
    try std.testing.expectEqual(diagnostics.Code.@"usage.invalid_option_value", failure.primary().?.code());
    failure.reset(.index);
    try checkReadBufferSize(1, &failure);
    try checkReadBufferSize(std.math.maxInt(u32), &failure);
    try std.testing.expectEqual(@as(u32, 0), try checkedZfiRecordCount(0, .{ .failure = &failure }));
    try std.testing.expectEqual(std.math.maxInt(u32), try checkedZfiRecordCount(std.math.maxInt(u32), .{ .failure = &failure }));
    try std.testing.expect(failure.primary() == null);
    if (comptime @sizeOf(usize) > @sizeOf(u32)) {
        const count: usize = @as(usize, std.math.maxInt(u32)) + 1;
        try std.testing.expectError(error.CapacityExceeded, checkReadBufferSize(count, &failure));
        const capacity = failure.primary().?.details.@"resource.capacity_limit_exceeded";
        try std.testing.expectEqual(diagnostics.CapacityRole.reader_buffer, capacity.capacity_role);
        try std.testing.expectEqual(@as(u128, std.math.maxInt(u32)), capacity.limit);
        try std.testing.expectEqual(@as(u128, count), capacity.required.exact);
        failure.reset(.index);
        try std.testing.expectError(error.InvalidIndex, checkedZfiRecordCount(count, .{ .failure = &failure }));
        const primary = failure.primary().?;
        const limit = primary.details.@"index.record_count_limit_exceeded";
        try std.testing.expectEqual(diagnostics.CountLimitRole.zfi_header, limit.limit_role);
        try std.testing.expectEqual(@as(u128, count), limit.record_count.exact);
        try std.testing.expectEqual(@as(u64, std.math.maxInt(u32)), limit.max_records);
        try std.testing.expectEqual(diagnostics.SourceId.zfi, primary.at.source.?);
        try std.testing.expectEqual(diagnostics.Stage.write, primary.stage);
        try std.testing.expectEqual(@as(?diagnostics.Cause, null), primary.cause);
        try std.testing.expectError(error.TooManySequences, checkedZfiRecordCount(count, .{}));
    }
}

test "[failure] - [index diagnostics]: reports the known ZFI path buffer limit" {
    var failure: diagnostics.Failure = undefined;
    failure.reset(.index);
    var buffer: [8]u8 = undefined;
    try std.testing.expectError(error.Io, zfiOutputPath(&buffer, "input", .{ .failure = &failure }));
    const primary = failure.primary().?;
    try std.testing.expectEqual(diagnostics.Code.@"io.invalid_path", primary.code());
    try std.testing.expectEqual(diagnostics.Stage.create, primary.stage);
    try std.testing.expectEqual(diagnostics.SourceId.zfi, primary.at.source);
    try std.testing.expect(primary.cause == null);
    const details = primary.details.@"io.invalid_path";
    try std.testing.expectEqual(diagnostics.PathReason.too_long, details.reason);
    try std.testing.expectEqual(@as(?u64, buffer.len), details.limit_bytes);
}
