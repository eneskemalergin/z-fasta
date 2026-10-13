//! Builds `.zfi` and `.fai` indexes in one pass through a 256 KiB reader.
//! A header starts with `>` at line start; names and bases follow `fasta.zig`. Empty records are
//! skipped. Duplicate names are first-wins by exact string equality, even when hashes collide.
//! A record is uniform when all lines share one width (a shorter last line is allowed); any other
//! layout gets a side table, and `--emit-fai` rejects it. Dense runs of lines are checked in
//! vector blocks; the first odd byte drops back to the line parser at that byte.
//! Inputs over 32 MiB stream records to the output and spill names and side tables to disk after
//! 1 MiB. Output is published only after re-checking the FASTA's size and mtime, and `.fai`
//! reaches stdout only after the whole file is accepted.

const std = @import("std");
const builtin = @import("builtin");
const index_format = @import("index_format.zig");
const fasta = @import("fasta.zig");
const diagnostics = @import("diagnostics");

pub const INDEX_READ_BUFFER_SIZE = 256 * 1024;

pub const MAX_INDEX_NAME_LEN = std.math.maxInt(u16);

pub const ScanOptions = struct {
    enable_dedup: bool,
    max_name_len: ?usize = null,
    collect_side_tables: bool = false,
    require_initial_header: bool = false,
    failure: ?*diagnostics.Failure = null,
};

pub const Options = struct {
    emit_fai: bool = false,
    enable_dedup: bool = true,
};

pub const Result = struct { record_count: u32 };

pub const NameDedup = NameDedupWith(std.hash_map.StringContext);

/// Slices are borrowed only for the duration of the callback.
pub const FastaRecordEmit = struct {
    record: IndexRecord,
    name: []const u8,
    source_reference: ?*const diagnostics.Reference = null,
    nonuniform_reason: diagnostics.NonuniformReason = .nonuniform_layout,
    uses_uniform_formula: bool,
    side_table: []const u8 = &.{},
    name_embedded: bool = false,
};

pub const ZfiIndex = struct {
    records: std.ArrayList(IndexRecord),
    side_tables: std.ArrayList(u8),
    name_blob: std.ArrayList(u8),

    pub fn deinit(self: *ZfiIndex, allocator: std.mem.Allocator) void {
        self.records.deinit(allocator);
        self.side_tables.deinit(allocator);
        self.name_blob.deinit(allocator);
    }
};

pub fn scanZfiReader(
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

pub fn scanZfiData(
    allocator: std.mem.Allocator,
    data: []const u8,
    enable_dedup: bool,
) !ZfiIndex {
    var r = std.Io.Reader.fixed(data);
    var read_buf: [INDEX_READ_BUFFER_SIZE]u8 = undefined;
    return scanZfiReader(allocator, &r, &read_buf, .{ .enable_dedup = enable_dedup });
}

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

pub fn zfiIndexToBytes(allocator: std.mem.Allocator, index: *const ZfiIndex, source_size: u64, source_mtime_ns: u64, failure: *diagnostics.Failure) diagnostics.Error![]u8 {
    return serializeZfiIndex(allocator, index, source_size, source_mtime_ns, .{ .failure = failure }) catch return diagnostics.metadata(failure.primary().?.code()).returned_error;
}

pub fn scanFastaReader(
    allocator: std.mem.Allocator,
    reader: *std.Io.Reader,
    read_buf: []u8,
    options: ScanOptions,
    external_dedup: anytype,
    ctx: anytype,
    comptime emitRecord: fn (@TypeOf(ctx), FastaRecordEmit) anyerror!void,
) !u32 {
    try checkReadBufferSize(read_buf.len, options.failure);
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

    var state = RecordScanner{
        .allocator = allocator,
        .collect_side_tables = options.collect_side_tables,
        .diag = .{ .failure = options.failure, .format = if (options.collect_side_tables) .zfi else .fai },
    };
    defer state.name.deinit(allocator);
    defer if (options.collect_side_tables) state.side_table.deinit(allocator);
    var record_count: u32 = 0;
    var emitter = IndexEmitter(@TypeOf(ctx), emitRecord, ?*Dedup){ .ctx = ctx, .seen_names = dedup_ptr, .record_count = &record_count, .diag = state.diag };
    const max_name_len = @min(options.max_name_len orelse MAX_INDEX_NAME_LEN, MAX_INDEX_NAME_LEN);
    var file_offset: u64 = 0;

    while (true) {
        var read_diag = state.diag;
        read_diag.at = .{ .source = .fasta, .location = .{ .offset = file_offset } };
        const n = reader.readSliceShort(read_buf) catch |err| return read_diag.fail(diagnostics.readFailure(read_buf.len, null), .read, .fasta, err, error.SourceReadFailed);
        if (n == 0) break;
        if (file_offset == 0 and options.require_initial_header and read_buf[0] != '>') {
            return read_diag.fail(.{ .@"index.initial_header_required" = .{ .byte = read_buf[0] } }, .scan, .fasta, null, error.NotFasta);
        }
        const next_offset = std.math.add(u64, file_offset, n) catch return read_diag.fail(.{ .@"resource.count_overflow" = .{ .counter = .source_bytes, .max_value = std.math.maxInt(u64), .current = file_offset, .increment = n } }, .scan, .fasta, null, error.SourceTooLong);
        if (comptime Dedup == SourceAdapter) {
            if (dedup_ptr) |seen| seen.source_names.setChunk(read_buf[0..n], file_offset);
        }
        try state.processChunk(read_buf[0..n], file_offset, max_name_len, &emitter, @TypeOf(emitter).emit);
        file_offset = next_offset;
    }

    if (state.active) {
        try state.emitRecord(&emitter, @TypeOf(emitter).emit);
    }

    return record_count;
}

pub fn scanFaiData(
    allocator: std.mem.Allocator,
    data: []const u8,
    writer: anytype,
    enable_dedup: bool,
) !u32 {
    var r = std.Io.Reader.fixed(data);
    var read_buf: [INDEX_READ_BUFFER_SIZE]u8 = undefined;
    return scanFaiReader(allocator, &r, &read_buf, writer, .{ .enable_dedup = enable_dedup }, null);
}

pub fn runIndex(
    allocator: std.mem.Allocator,
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
        .details = diagnostics.openFailure(err),
        .stage = .open,
        .at = .{ .source = .fasta },
        .cause = .{ .zig = err },
    });
    defer file.close(io);
    const stat = file.stat(io) catch |err| return failure.fail(.{ .details = .{ .@"io.stat_failed" = {} }, .stage = .stat, .at = .{ .source = .fasta }, .cause = .{ .zig = err } });
    if (stat.size == 0) return failure.fail(.{ .details = .{ .@"fasta.no_records" = {} }, .stage = .stat, .at = .{ .source = .fasta } });
    const count = (if (stat.size <= SMALL_INDEX_MAX_BYTES)
        runIndexSmall(allocator, io, environ, path, options, file, stat, diag)
    else
        runIndexLarge(allocator, io, environ, path, options, file, stat, diag)) catch |err| {
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
        return diag.fail(diagnostics.readFailure(requested, received), .read, source, cause, legacy);
    }

    pub fn fail(self: @This(), details: diagnostics.Details, stage: diagnostics.Stage, source: diagnostics.SourceId, cause: ?anyerror, legacy: anyerror) anyerror {
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

    pub fn oom(self: @This(), role: diagnostics.AllocationRole, bytes: ?u64) error{OutOfMemory} {
        @branchHint(.cold);
        if (self.failure) |failure| {
            switch (failure.fail(.{ .details = .{ .@"resource.out_of_memory" = .{ .allocation_role = role, .requested_bytes = bytes } }, .stage = .allocate, .at = self.at.reference(), .cause = .{ .zig = error.OutOfMemory } })) {
                else => {},
            }
        }
        return error.OutOfMemory;
    }

    pub fn capacity(self: @This(), role: diagnostics.CapacityRole, limit: u128, required: u128, legacy: anyerror) anyerror {
        return self.fail(.{ .@"resource.capacity_limit_exceeded" = .{ .capacity_role = role, .limit = limit, .required = .{ .exact = required } } }, .allocate, self.at.source orelse .fasta, null, legacy);
    }

    pub fn consistency(self: @This(), check: diagnostics.Check, legacy: anyerror) anyerror {
        return self.fail(.{ .@"internal.consistency_failed" = .{ .component = .indexer, .check = check } }, .scan, self.at.source orelse .fasta, null, legacy);
    }

    pub fn identifierTooLong(self: @This(), observed: u128, limit: usize) anyerror {
        return self.fail(.{ .@"index.identifier_too_long" = .{ .observed_bytes = .{ .at_least = observed }, .limit_bytes = limit, .format = self.format } }, .scan, .fasta, null, error.HeaderTooLong);
    }

    pub fn lineLimit(self: @This(), field: diagnostics.IndexField, observed: u64, width_role: diagnostics.LineWidthRole) anyerror {
        return self.fail(.{ .@"index.sequence_line_too_long" = .{ .field = field, .observed = .{ .exact = observed }, .limit = std.math.maxInt(u32), .width_role = width_role, .format = self.format } }, .scan, .fasta, null, error.SequenceLineTooLong);
    }

    pub fn recordLimit(self: @This(), count: u64) anyerror {
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

        /// True when the exact name was already seen; otherwise keeps a copy.
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

const RecordScanner = fasta.RecordScannerWith(Diagnostics);

comptime {
    if (@sizeOf(fasta.LineEntry) != @sizeOf(index_format.SideTableLine)) @compileError("side-table line layout differs");
    for (.{ "base_start", "byte_offset", "line_bytes", "line_bases" }) |field| {
        if (@offsetOf(fasta.LineEntry, field) != @offsetOf(index_format.SideTableLine, field)) @compileError("side-table line layout differs");
    }
}

/// Applies the duplicate policy and the record limit to each scanned record, then emits it as an
/// index record.
fn IndexEmitter(comptime Ctx: type, comptime emitRecord: fn (Ctx, FastaRecordEmit) anyerror!void, comptime Seen: type) type {
    return struct {
        ctx: Ctx,
        seen_names: Seen,
        record_count: *u32,
        diag: Diagnostics,

        fn emit(self: *@This(), info: fasta.ScannedRecord) anyerror!void {
            var name_embedded = false;
            var name_offset: u64 = info.header_offset + 1;
            var out_name_len: u16 = @intCast(info.name.len);
            var name_pad: [6]u8 = .{0} ** 6;
            if (self.seen_names) |seen| {
                if (try seen.observeAt(info.name, info.header_offset + 1)) return;
                if (comptime @hasField(@TypeOf(seen.*), "last_offset")) {
                    if (seen.bindsToBlob()) {
                        name_offset = seen.last_offset;
                        out_name_len = seen.last_len;
                        name_pad[0] = index_format.NAME_IN_ZFI_FLAG;
                        name_embedded = true;
                    }
                }
            }
            if (self.record_count.* == std.math.maxInt(u32)) {
                @branchHint(.cold);
                var diag = self.diag;
                diag.at.location.offset = info.header_offset;
                diag.at.record = info.reference.record;
                return diag.recordLimit(@as(u64, self.record_count.*) + 1);
            }
            try emitRecord(self.ctx, .{
                .record = .{
                    .name_offset = name_offset,
                    .name_len = out_name_len,
                    .seq_offset = info.seq_offset,
                    .seq_len = info.seq_len,
                    .line_bases = info.line_bases,
                    .line_bytes = info.line_bytes,
                    .flags_and_offset = name_pad,
                },
                .name = info.name,
                .source_reference = info.reference,
                .nonuniform_reason = info.nonuniform_reason,
                .uses_uniform_formula = info.uses_uniform_formula,
                .side_table = info.side_table,
                .name_embedded = name_embedded,
            });
            self.record_count.* += 1;
        }
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
    const count = try scanFastaReader(allocator, reader, read_buf, options, dedup, &ctx, Ctx.emit);
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
    const count = scanFaiReader(allocator, &reader.interface, read_buf, &tmp_fw.interface, .{ .enable_dedup = options.enable_dedup, .require_initial_header = true, .failure = diag.failure }, dedup) catch |err| {
        retainScanCause(diag.failure, reader, &tmp_fw);
        return err;
    };
    if (diagnostics.flushOutput(&tmp_fw, null)) |cause| return diag.fail(.{ .@"io.flush_failed" = {} }, .flush, .fai_spool, cause, error.FaiSpoolWriteFailed);
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
    if (diagnostics.flushOutput(&stdout_fw, if (diag.failure) |failure| &failure.outputs.stdout else null)) |cause| return diag.fail(.{ .@"io.flush_failed" = {} }, .flush, .stdout, cause, error.StdoutFlushFailed);
    if (diag.failure) |failure| {
        failure.outputs.stdout.records = count;
        failure.outputs.stdout.may_be_partial = false;
        failure.work_complete = true;
    }
    return count;
}

fn runIndexSmall(allocator: std.mem.Allocator, io: std.Io, environ: std.process.Environ, path: []const u8, options: Options, file: std.Io.File, stat: std.Io.File.Stat, diag: Diagnostics) !u32 {
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
    var index = scanZfiReader(allocator, &file_reader.interface, &read_buf, .{ .enable_dedup = options.enable_dedup, .require_initial_header = true, .failure = diag.failure }) catch |err| {
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
    if (diagnostics.flushOutput(&output, if (diag.failure) |failure| &failure.outputs.file else null)) |cause| return diag.fail(.{ .@"io.flush_failed" = {} }, .flush, .zfi, cause, error.ZfiWriteFailed);
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

fn runIndexLarge(allocator: std.mem.Allocator, io: std.Io, environ: std.process.Environ, path: []const u8, options: Options, file: std.Io.File, stat: std.Io.File.Stat, diag: Diagnostics) !u32 {
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
    const count = scanFastaReader(allocator, &reader.interface, &read_buf, .{ .enable_dedup = options.enable_dedup, .require_initial_header = true, .collect_side_tables = true, .failure = diag.failure }, dedup, &catalog, StreamedCatalog.emit) catch |err| {
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

const PackedNames = fasta.PackedNamesWith(Diagnostics);
const SourceDedup = fasta.SourceDedupWith(Diagnostics);

fn Tables(comptime Offset: type) type {
    return SourceDedup.Tables(Offset);
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
            if (err == error.EndOfStream) return self.diag.fail(diagnostics.readFailure(count, reader.bufferedLen()), .read, self.role, null, err);
            return self.diag.fail(diagnostics.readFailure(count, null), .read, self.role, cause, err);
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
            prefix = fasta.commonPrefix(name, self.anchor[0..self.anchor_len]);
            suffix = fasta.commonSuffix(name, self.anchor[0..self.anchor_len], prefix);
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
        if (diagnostics.flushOutput(&output, if (self.diag.failure) |failure| &failure.outputs.file else null)) |cause| return self.diag.fail(.{ .@"io.flush_failed" = {} }, .flush, .zfi, cause, error.ZfiWriteFailed);
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
                    const prefix = fasta.commonPrefix(name, anchor);
                    const suffix = fasta.commonSuffix(name, anchor, prefix);
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
                if (err == error.EndOfStream) return self.raw.diag.fail(diagnostics.readFailure(@sizeOf(IndexRecord), records.interface.bufferedLen()), .read, .catalog_spool, null, err);
                return self.raw.diag.fail(diagnostics.readFailure(@sizeOf(IndexRecord), null), .read, .catalog_spool, if (records.err) |cause| cause else err, err);
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
    const count = try scanFastaReader(
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
    var index = try scanZfiReader(allocator, &reader, &buffer, .{ .enable_dedup = true, .failure = &failure });
    defer index.deinit(allocator);
    const bytes = try zfiIndexToBytes(allocator, &index, data.len, 0, &failure);
    defer allocator.free(bytes);
    try std.testing.expectEqual(@as(usize, 2), index.records.items.len);
    try std.testing.expect(failure.primary() == null);
}

test "[property] - [fai output]: matches standard formatting at integer boundaries" {
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

test "[unit] - [fasta scanner]: gates initial-header validation by option" {
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
            null,
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
            null,
        ),
    );
}

test "[failure] - [fasta scanner]: rejects sequence lines above u32 geometry" {
    var state = RecordScanner{
        .allocator = std.testing.allocator,
        .in_pending_line = true,
        .pending_line_bases = std.math.maxInt(u32),
    };
    const Ctx = struct {
        fn emit(_: *@This(), _: fasta.ScannedRecord) !void {
            unreachable;
        }
    };
    var ctx = Ctx{};
    try std.testing.expectError(error.SequenceLineTooLong, state.processChunk("A", 0, null, &ctx, Ctx.emit));

    state.pending_line_bases = 0;
    state.pending_line_bytes = std.math.maxInt(u32);
    try std.testing.expectError(error.SequenceLineTooLong, state.processChunk("\n", 0, null, &ctx, Ctx.emit));

    state.pending_line_bases = 1;
    try std.testing.expectError(error.SequenceLineTooLong, state.emitRecord(&ctx, Ctx.emit));
}

test "[edge] - [fasta scanner]: rejects excess records before emitting" {
    const Ctx = struct {
        count: usize = 0,
        fn emit(self: *@This(), _: FastaRecordEmit) !void {
            self.count += 1;
        }
    };
    const reference = diagnostics.Reference{ .source = .fasta, .location = .{ .offset = 0 }, .record = 1 };
    const record = fasta.ScannedRecord{
        .name = "limit",
        .header_offset = 0,
        .reference = &reference,
        .seq_offset = 7,
        .seq_len = 1,
        .line_bases = 1,
        .line_bytes = 2,
        .uses_uniform_formula = true,
        .side_table = &.{},
        .nonuniform_reason = .nonuniform_layout,
    };
    var ctx = Ctx{};
    var count: u32 = std.math.maxInt(u32) - 1;
    var plain = IndexEmitter(*Ctx, Ctx.emit, ?*NameDedup){ .ctx = &ctx, .seen_names = null, .record_count = &count, .diag = .{} };
    try plain.emit(record);
    try std.testing.expectEqual(std.math.maxInt(u32), count);
    try std.testing.expectError(error.TooManySequences, plain.emit(record));
    try std.testing.expectEqual(@as(usize, 1), ctx.count);
    try std.testing.expectEqual(std.math.maxInt(u32), count);

    var seen = NameDedup.init(std.testing.allocator);
    defer seen.deinit();
    try std.testing.expect(!try seen.observe("limit"));
    var deduped = IndexEmitter(*Ctx, Ctx.emit, ?*NameDedup){ .ctx = &ctx, .seen_names = &seen, .record_count = &count, .diag = .{} };
    try deduped.emit(record);
    try std.testing.expectEqual(@as(usize, 1), ctx.count);
    try std.testing.expectEqual(std.math.maxInt(u32), count);
}

test "[edge] - [name catalog]: full source tables retain duplicate lookups" {
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

test "[property] - [name catalog]: remains exact under hash collisions" {
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

test "[failure] - [fasta scanner]: maps reader failures to SourceReadFailed" {
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
            null,
        ),
    );
}

test "[failure] - [fai scanner]: maps writer failures to FaiSpoolWriteFailed" {
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
            null,
        ),
    );
}

test "[integration] - [fai output]: uses exclusive files and cleans them" {
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

test "[integration] - [index]: detects path replacement" {
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

test "[property] - [name catalog]: complete names resolve fingerprint collisions" {
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

test "[unit] - [name catalog]: exact comparisons use the current input chunk" {
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

test "[edge] - [name catalog]: retain offsets above four GiB" {
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

test "[property] - [name catalog]: cached names resolve fingerprint collisions" {
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

test "[failure] - [name catalog]: read errors leave the insertion unfinished" {
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

test "[failure] - [name catalog]: block growth cleans every allocation failure" {
    try std.testing.checkAllAllocationFailures(std.testing.allocator, checkPackedNameAllocations, .{});
}

test "[failure] - [name catalog]: table and entry growth clean every allocation failure" {
    var data: std.ArrayList(u8) = .empty;
    defer data.deinit(std.testing.allocator);
    for (0..5000) |i| {
        var buffer: [32]u8 = undefined;
        try data.appendSlice(std.testing.allocator, try std.fmt.bufPrint(&buffer, "lookup-{d:0>8}\n", .{i}));
    }
    try std.testing.checkAllAllocationFailures(std.testing.allocator, checkSourceNameAllocations, .{data.items});
}

test "[edge] - [name catalog]: maximum names cross copied and retained cache blocks" {
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

test "[failure] - [name catalog]: short reads distinguish collisions from truncation" {
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

test "[integration] - [name catalog]: spill replay preserves bytes and removes its file" {
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

test "[property] - [name catalog]: raw and compressed replay retain exact names" {
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

test "[failure] - [zfi writer]: parse errors precede record write failures" {
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

test "[failure] - [name catalog]: source offsets cannot overlap cache flags" {
    var narrow = SourceDedup.init(std.testing.allocator, std.testing.io, undefined, 0);
    defer narrow.deinit();
    try std.testing.expectError(error.SourceChanged, narrow.observe("name", Tables(u32).CACHE_FLAG));
    try std.testing.expectEqual(@as(u32, 0), narrow.uniqueCount());
    var wide = SourceDedup.init(std.testing.allocator, std.testing.io, undefined, std.math.maxInt(u32));
    defer wide.deinit();
    try std.testing.expectError(error.SourceChanged, wide.observe("name", Tables(u64).CACHE_FLAG));
    try std.testing.expectEqual(@as(u32, 0), wide.uniqueCount());
}

test "[failure] - [name catalog]: oversized offsets roll back the insertion" {
    if (comptime @bitSizeOf(usize) <= 32) return error.SkipZigTest;
    var names = PackedNames.init(std.testing.allocator);
    defer names.deinit();
    names.position = @as(usize, std.math.maxInt(u32)) + 1;
    try std.testing.expectError(error.OutOfMemory, names.observeAt("name", 0));
    try std.testing.expectEqual(@as(u32, 0), names.map.count());
}

test "[failure] - [index diagnostics]: retains scanner limit operands and source positions" {
    var failure: diagnostics.Failure = undefined;
    var state = RecordScanner{ .allocator = std.testing.allocator, .diag = .{ .failure = &failure } };
    defer state.name.deinit(std.testing.allocator);
    defer state.side_table.deinit(std.testing.allocator);
    const max = std.math.maxInt(u32);
    const Ctx = struct {
        fn emit(_: *@This(), _: FastaRecordEmit) !void {
            unreachable;
        }
    };
    var ctx = Ctx{};
    var count: u32 = max;
    var emitter = IndexEmitter(*Ctx, Ctx.emit, ?*NameDedup){ .ctx = &ctx, .seen_names = null, .record_count = &count, .diag = state.diag };
    const Emitter = @TypeOf(emitter);
    for (0..3) |i| {
        failure.reset(.index);
        state.in_pending_line = true;
        state.pending_line_file_offset = 700;
        state.pending_line_bases = if (i == 0) max else 0;
        state.pending_line_bytes = if (i == 0) 0 else max;
        if (i == 2) {
            try std.testing.expectError(error.InvalidIndex, state.emitRecord(&emitter, Emitter.emit));
        } else {
            try std.testing.expectError(error.InvalidIndex, state.processChunk("A", 700, null, &emitter, Emitter.emit));
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
    try state.processChunk(">name\nA\n", 42, null, &emitter, Emitter.emit);
    try std.testing.expectError(error.InvalidIndex, state.emitRecord(&emitter, Emitter.emit));
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

test "[failure] - [index diagnostics]: spill IO retains causes and exact short reads" {
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

test "[regression] - [index diagnostics]: missing delimiters retain EOF without a field length" {
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

test "[regression] - [index diagnostics]: removed staging does not remain staged" {
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

test "[failure] - [index diagnostics]: cleanup uncertainty preserves primary and published state" {
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

test "[edge] - [index diagnostics]: buffer width and ZFI count reject before narrowing" {
    var failure: diagnostics.Failure = undefined;
    failure.reset(.index);
    var reader = std.Io.Reader.fixed(">name\nACGT\n");
    var empty_buffer: [0]u8 = .{};
    var output_bytes: [16]u8 = undefined;
    var output = std.Io.Writer.fixed(&output_bytes);
    var failing = std.testing.FailingAllocator.init(std.testing.allocator, .{ .fail_index = 0 });
    try std.testing.expectError(error.InvalidArgument, scanFaiReader(failing.allocator(), &reader, &empty_buffer, &output, .{ .enable_dedup = true, .failure = &failure }, null));
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
