//! Loads and checks `.zfi` and `.fai` sidecars; owns the `.zfi` byte layout and `LoadedIndex`.
//! A present `.zfi` wins; `.fai` is tried only when `.zfi` is absent, unless one is selected.
//! An index older than its FASTA is rejected. `.zfi` also stores the FASTA size and mtime.
//! `.zfi` version 2, little-endian: a 16-byte header, 24-byte records, side tables, the name
//! blob, the name table, then the FASTA's mtime and the blob's length. Everything is read in
//! place. A load for `get` checks the file's frame and each record a request uses; `stats`, the
//! package, and `validate --index` check every record, and the latter two the name table.
//! `LoadedIndex.deinit` releases the FASTA descriptor, mappings, and arena; names and records
//! borrow them until then. Loaders never reset the caller's `Failure`.

const std = @import("std");
const builtin = @import("builtin");
const platform = @import("platform.zig");
pub const diagnostics = @import("diagnostics");

pub const ZFI_MAGIC: [4]u8 = .{ 'Z', 'F', 'I', 0x02 };

pub const MAX_SIDE_TABLE_OFFSET = (1 << 39) - 1;
/// The largest FASTA byte offset or sequence length a `.zfi` record holds.
pub const MAX_ZFI_SOURCE_OFFSET = (1 << 48) - 1;
pub const MAX_ZFI_NAME_OFFSET = (1 << 40) - 1;

pub const ZFI_NAME_FOOTER_MAGIC: [4]u8 = .{ 'Z', 'F', 'N', 'M' };

pub const ZFI_HEADER_BYTES: usize = 16;
pub const ZFI_INDEX_RECORD_BYTES: usize = 24;
pub const ZFI_SIDE_TABLE_LINE_BYTES: usize = 24;

pub const ZFI_NAME_FOOTER_BYTES: usize = 12;

pub const ZFI_SOURCE_ID_MAGIC: [4]u8 = .{ 'Z', 'F', 'I', 'D' };
pub const ZFI_SOURCE_ID_BYTES: usize = 12;

pub fn encodeZfiNameFooter(name_blob_len: u64) [ZFI_NAME_FOOTER_BYTES]u8 {
    var out: [ZFI_NAME_FOOTER_BYTES]u8 = undefined;
    @memcpy(out[0..4], &ZFI_NAME_FOOTER_MAGIC);
    std.mem.writeInt(u64, out[4..12], name_blob_len, .little);
    return out;
}

pub fn encodeZfiSourceId(source_mtime_ns: u64) [ZFI_SOURCE_ID_BYTES]u8 {
    var out: [ZFI_SOURCE_ID_BYTES]u8 = undefined;
    @memcpy(out[0..4], &ZFI_SOURCE_ID_MAGIC);
    std.mem.writeInt(u64, out[4..12], source_mtime_ns, .little);
    return out;
}

pub fn timestampToNs(ts: std.Io.Timestamp) error{TimestampOutOfRange}!u64 {
    return std.math.cast(u64, ts.nanoseconds) orelse error.TimestampOutOfRange;
}

pub const ZfiHeader = extern struct {
    magic: [4]u8,
    record_count: u32,
    source_size: u64,
};

pub const SideTableLine = extern struct {
    base_start: u64,
    byte_offset: u64,
    line_bytes: u32,
    line_bases: u32,
};

/// Name offsets point into the `.zfi` blob or the mapped `.fai`; prefer `LoadedIndex.recordName`.
/// A non-uniform record keeps no line geometry; its side table holds it.
pub const IndexRecord = extern struct {
    name_offset: u64,
    name_len: u16,
    flags_and_offset: [6]u8 = .{0} ** 6,
    seq_offset: u64,
    seq_len: u64,
    line_bases: u32,
    line_bytes: u32,

    /// Caller assumes the name range was validated against `data`.
    pub fn getName(self: IndexRecord, data: []const u8) []const u8 {
        return data[self.name_offset..][0..self.name_len];
    }

    pub fn isUniformWidth(self: IndexRecord) bool {
        return (self.flags_and_offset[0] & NON_UNIFORM_WIDTH_FLAG) == 0;
    }

    pub fn sideTableOffset(self: IndexRecord) u64 {
        var offset: u64 = 0;
        inline for (0..SIDE_TABLE_OFFSET_BYTES) |i| {
            offset |= @as(u64, self.flags_and_offset[i + 1]) << (i * 8);
        }
        return offset;
    }

    pub fn markNonUniform(self: *IndexRecord, offset: u64) !void {
        if (offset > MAX_SIDE_TABLE_OFFSET) return error.SideTableOffsetTooLarge;
        self.flags_and_offset = .{0} ** 6;
        self.flags_and_offset[0] = NON_UNIFORM_WIDTH_FLAG;
        inline for (0..SIDE_TABLE_OFFSET_BYTES) |i| {
            self.flags_and_offset[i + 1] = @intCast((offset >> (i * 8)) & 0xff);
        }
    }
};

pub fn encodeZfiHeader(header: ZfiHeader) [ZFI_HEADER_BYTES]u8 {
    var out: [ZFI_HEADER_BYTES]u8 = undefined;
    @memcpy(out[0..4], &header.magic);
    std.mem.writeInt(u32, out[4..8], header.record_count, .little);
    std.mem.writeInt(u64, out[8..16], header.source_size, .little);
    return out;
}

/// A `.zfi` record. Bits 0-47 hold the sequence offset, 48-63 the name length, 64-111 the
/// sequence length, and 112-151 the name's offset in the blob. Bit 191 marks a non-uniform record,
/// whose side-table offset fills bits 152-190; a uniform record keeps its symbols per line in bits
/// 152-183 and its line ending's bytes in 184-190.
pub const ZfiRecord = extern struct {
    words: [3]u64,

    /// Caller assumes the offsets and lengths fit their fields and a uniform line ends in at most
    /// 127 bytes; only Debug checks it.
    pub fn encode(rec: IndexRecord) ZfiRecord {
        std.debug.assert(rec.seq_offset <= MAX_ZFI_SOURCE_OFFSET and rec.seq_len <= MAX_ZFI_SOURCE_OFFSET and rec.name_offset <= MAX_ZFI_NAME_OFFSET);
        const layout: u64 = if (rec.isUniformWidth()) blk: {
            std.debug.assert(rec.line_bytes >= rec.line_bases and rec.line_bytes - rec.line_bases <= 127);
            break :blk @as(u64, rec.line_bases) | @as(u64, rec.line_bytes - rec.line_bases) << 32;
        } else rec.sideTableOffset() | 1 << 39;
        return .{ .words = .{
            rec.seq_offset | @as(u64, rec.name_len) << 48,
            rec.seq_len | rec.name_offset << 48,
            rec.name_offset >> 16 | layout << 24,
        } };
    }

    pub fn decode(self: ZfiRecord) IndexRecord {
        const w = self.words;
        var rec = IndexRecord{
            .name_offset = w[1] >> 48 | (w[2] & 0xFF_FFFF) << 16,
            .name_len = @truncate(w[0] >> 48),
            .seq_offset = w[0] & MAX_ZFI_SOURCE_OFFSET,
            .seq_len = w[1] & MAX_ZFI_SOURCE_OFFSET,
            .line_bases = 0,
            .line_bytes = 0,
        };
        const layout = w[2] >> 24;
        if (layout >> 39 == 0) {
            rec.line_bases = @truncate(layout);
            rec.line_bytes = rec.line_bases +| @as(u32, @intCast(layout >> 32 & 0x7F));
        } else {
            rec.markNonUniform(layout & MAX_SIDE_TABLE_OFFSET) catch unreachable;
        }
        return rec;
    }
};

/// Index entries: mapped in place from a `.zfi`, or held in memory.
pub const Records = union(enum) {
    mapped: []const ZfiRecord,
    owned: []const IndexRecord,
    /// A `.zfi` read on demand: only `LoadedIndex.checkedRecord` reaches its records.
    on_demand: usize,

    pub fn len(self: Records) usize {
        return switch (self) {
            .on_demand => |count| count,
            inline else => |items| items.len,
        };
    }

    pub fn get(self: Records, index: usize) IndexRecord {
        return switch (self) {
            .mapped => |items| items[index].decode(),
            .owned => |items| items[index],
            .on_demand => unreachable,
        };
    }
};

/// The name table after the name blob: `slotCount(record_count)` little-endian `u32` slots. A
/// slot is 0 when empty; otherwise its low `ordinalBits` bits hold a record's ordinal plus 1 and the
/// bits above hold the low bits of the name's XXH3-64 hash (seed 0). A lookup starts at `home` and
/// moves forward, wrapping. Records go in by ordinal, so the first match is the first record with
/// that name.
pub const NameTable = struct {
    pub fn slotCount(record_count: u64) u64 {
        return record_count + record_count / 4 + 1;
    }

    fn ordinalBits(record_count: u64) u6 {
        return @intCast(64 - @clz(record_count));
    }

    pub fn hash(name: []const u8) u64 {
        return std.hash.XxHash3.hash(0, name);
    }

    fn home(name_hash: u64, slots: u64) u64 {
        return @intCast((@as(u128, name_hash) * slots) >> 64);
    }

    fn fingerprint(name_hash: u64, ordinal_bits: u6) u64 {
        return name_hash & ((@as(u64, 1) << @intCast(32 - @as(u7, ordinal_bits))) - 1);
    }

    pub const Builder = struct {
        slots: []u32,
        ordinal_bits: u6,
        count: u64 = 0,
        /// Hashes not placed yet, oldest at `head`; their home slots were prefetched on arrival, so
        /// the random slot reads overlap instead of stalling one at a time.
        queue: [QUEUE]u64 = undefined,
        head: usize = 0,
        queued: usize = 0,
        const QUEUE = 16;

        pub fn init(allocator: std.mem.Allocator, record_count: u64) error{OutOfMemory}!Builder {
            const slots = try allocator.alloc(u32, std.math.cast(usize, slotCount(record_count)) orelse return error.OutOfMemory);
            // On freshly faulted pages these 64-byte stores clear the table faster than `@memset`.
            const blocks: [*]align(4) volatile @Vector(16, u32) = @ptrCast(slots.ptr);
            for (0..slots.len / 16) |i| blocks[i] = @splat(0);
            @memset(slots[slots.len / 16 * 16 ..], 0);
            return .{ .slots = slots, .ordinal_bits = ordinalBits(record_count) };
        }

        pub fn deinit(self: *Builder, allocator: std.mem.Allocator) void {
            allocator.free(self.slots);
        }

        /// Caller assumes fewer names than the record count it was built for.
        pub fn insert(self: *Builder, name: []const u8) void {
            self.insertHash(hash(name));
        }

        /// `insert` for a name whose `hash` the caller already has.
        pub fn insertHash(self: *Builder, name_hash: u64) void {
            @prefetch(&self.slots[@intCast(home(name_hash, self.slots.len))], .{ .rw = .write });
            if (self.queued == QUEUE) {
                self.place(self.queue[self.head]);
                self.queue[self.head] = name_hash;
                self.head = (self.head + 1) % QUEUE;
            } else {
                self.queue[self.queued] = name_hash;
                self.queued += 1;
            }
        }

        /// The finished table, every inserted name placed.
        pub fn finish(self: *Builder) []u32 {
            for (0..self.queued) |i| self.place(self.queue[(self.head + i) % QUEUE]);
            self.head = 0;
            self.queued = 0;
            return self.slots;
        }

        pub fn bytes(self: *Builder) []const u8 {
            return std.mem.sliceAsBytes(self.finish());
        }

        fn place(self: *Builder, name_hash: u64) void {
            std.debug.assert(self.count + 1 < self.slots.len);
            self.count += 1;
            const entry: u32 = @intCast(fingerprint(name_hash, self.ordinal_bits) << self.ordinal_bits | self.count);
            var slot: usize = @intCast(home(name_hash, self.slots.len));
            // Sixteen slots per step while they fit before the end; the single-slot walk wraps.
            while (slot + 16 <= self.slots.len) : (slot += 16) {
                const window: @Vector(16, u32) = self.slots[slot..][0..16].*;
                const empty: u16 = @bitCast(window == @as(@Vector(16, u32), @splat(0)));
                if (empty != 0) {
                    self.slots[slot + @ctz(empty)] = entry;
                    return;
                }
            }
            while (self.slots[slot] != 0) slot = if (slot + 1 == self.slots.len) 0 else slot + 1;
            self.slots[slot] = entry;
        }
    };
};

pub const IndexSelection = union(enum) {
    automatic,
    zfi: []const u8,
    fai: []const u8,
};

pub const LoadMode = union(enum) {
    /// Finds any name: a `.zfi` through its name table, a `.fai` through a full map.
    lookup_full_map,
    /// Checks every record and the `.zfi` name table, without a `.fai` lookup map.
    catalog,
    /// Checks every record and keeps no names; a `.fai` is streamed.
    stats_scan,
    /// A `.zfi` read on demand; from a `.fai`, only the first line for each name, checked. The
    /// index borrows `names`, which must outlive it.
    requested: []const []const u8,
    /// Makes the file checks and leaves every entry to the stream, which checks each as it reads it.
    stream: StreamTarget,
    /// A mapped `.zfi` with only its file checks; `checkedRecord` checks each record and side
    /// table the first time it reads it, and lookups check what they probe. A `.fai` loads as for
    /// `lookup_full_map`.
    on_access,
};

/// The stream a `.stream` load opens, and whether it rebuilds the `.zfi` name table to compare it.
const StreamTarget = struct { entries: *EntryStream, check_name_table: bool = true };

pub const LoadedIndex = struct {
    io: std.Io,
    fasta_file: std.Io.File,
    zfi_map: ?std.Io.File.MemoryMap = null,
    fai_map: ?std.Io.File.MemoryMap = null,
    records: Records,
    name_map: ?std.StringHashMapUnmanaged(u32) = null,
    name_slices: []const []const u8 = &.{},
    name_blob: ?[]const u8 = null,
    fai_data: ?platform.MappedBytes = null,
    fasta_size: u64,
    zfi_data: ?[]align(@alignOf(SideTableLine)) const u8,
    zfi_side_start: usize = 0,
    zfi_side_end: usize = 0,
    name_table: []const u8 = &.{},
    name_table_offset: usize = 0,
    /// Set for a `.zfi` loaded for `get`; every other load but `on_access` checks all records first.
    on_demand: ?*ZfiOnDemand = null,
    /// One bit per record of an `on_access` `.zfi` or a lazily loaded `.fai`, set once
    /// `checkedRecord` checked it.
    checked_records: []usize = &.{},
    /// Where each line of a lazily loaded `.fai` starts; `checkedRecord` parses it when asked.
    fai_lines: []const FaiLine = &.{},
    source: IndexSource,
    arena: std.heap.ArenaAllocator,

    pub const IndexSource = enum { zfi, fai };

    /// Record `rec_idx`; a `.zfi` loaded for `get` or `on_access` checks it and its side table
    /// first.
    pub fn checkedRecord(self: *const LoadedIndex, rec_idx: usize, failure: *diagnostics.Failure) diagnostics.Error!IndexRecord {
        var context = Diagnostics{ .failure = failure, .at = .{ .source = .zfi } };
        if (self.on_demand) |zfi| return zfi.checked(rec_idx, &context);
        if (self.fai_lines.len != 0) return self.checkedFaiRecord(rec_idx, &context);
        const rec = self.records.get(rec_idx);
        if (self.checked_records.len == 0) return rec;
        const word = &self.checked_records[rec_idx / @bitSizeOf(usize)];
        const bit = @as(usize, 1) << @intCast(rec_idx % @bitSizeOf(usize));
        if (word.* & bit != 0) return rec;
        _ = try checkZfiRecord(self.zfi_data.?, rec, rec_idx, self.fasta_size, self.name_blob.?.len, self.zfi_side_start, self.zfi_side_end, &context);
        word.* |= bit;
        return rec;
    }

    /// The first record named `name`. A `.zfi` looks it up in its name table, and a damaged table
    /// entry or record name fails with its index diagnostic.
    pub fn findName(self: *const LoadedIndex, name: []const u8, failure: *diagnostics.Failure) diagnostics.Error!?usize {
        if (self.on_demand) |zfi| {
            var context = Diagnostics{ .failure = failure, .at = .{ .source = .zfi } };
            return zfi.find(name, &context);
        }
        if (self.name_table.len == 0) return self.lookupName(name);
        var context = Diagnostics{ .failure = failure, .at = .{ .source = .zfi } };
        switch (self.probeName(name)) {
            .found => |ordinal| return ordinal,
            .missing => return null,
            .invalid_slot => |slot| {
                context.at.location.offset = self.name_table_offset + slot * 4;
                return context.fail(.{ .@"index.name_table_invalid" = .{ .slot = slot } });
            },
            .invalid_name => |ordinal| {
                context.at = recordReference(ordinal);
                return context.fail(nameExtentDefect(self.records.get(ordinal), self.name_blob.?.len).?);
            },
        }
    }

    /// Resolves a request batch's names with few reads through `buffer`, for `findNameAt` with the
    /// same slices; the results borrow the bytes of `names` until `clearPrefetch`.
    pub fn prefetchNames(self: *const LoadedIndex, names: []const []const u8, buffer: []u8) void {
        const zfi = self.on_demand orelse return;
        zfi.prefetch(names, buffer) catch zfi.clearPrefetch();
    }

    /// `findName` for request `request` of the batch given to `prefetchNames`; `name` must be that
    /// request's slice. The record it finds is the one `checkedRecord` reads next.
    pub fn findNameAt(self: *const LoadedIndex, request: usize, name: []const u8, failure: *diagnostics.Failure) diagnostics.Error!?usize {
        const zfi = self.on_demand orelse return self.findName(name, failure);
        var context = Diagnostics{ .failure = failure, .at = .{ .source = .zfi } };
        return zfi.findAt(request, name, &context);
    }

    pub fn clearPrefetch(self: *const LoadedIndex) void {
        if (self.on_demand) |zfi| zfi.clearPrefetch();
    }

    const Probe = union(enum) { found: usize, missing, invalid_slot: u64, invalid_name: usize };

    fn probeName(self: *const LoadedIndex, name: []const u8) Probe {
        const count: u64 = self.records.len();
        const slots = self.name_table.len / 4;
        const ordinal_bits = NameTable.ordinalBits(count);
        const name_hash = NameTable.hash(name);
        const wanted = NameTable.fingerprint(name_hash, ordinal_bits);
        var slot = NameTable.home(name_hash, slots);
        for (0..slots) |_| {
            const value: u64 = std.mem.readInt(u32, self.name_table[@intCast(slot * 4)..][0..4], .little);
            if (value == 0) return .missing;
            const ordinal_plus_one = value & ((@as(u64, 1) << ordinal_bits) - 1);
            if (ordinal_plus_one == 0 or ordinal_plus_one > count) return .{ .invalid_slot = slot };
            if (value >> ordinal_bits == wanted) {
                const ordinal: usize = @intCast(ordinal_plus_one - 1);
                const rec = self.records.get(ordinal);
                if (nameExtentDefect(rec, self.name_blob.?.len) != null) return .{ .invalid_name = ordinal };
                if (std.mem.eql(u8, rec.getName(self.name_blob.?), name)) return .{ .found = ordinal };
            }
            slot = if (slot + 1 == slots) 0 else slot + 1;
        }
        return .missing;
    }

    fn recordNameData(self: *const LoadedIndex) []const u8 {
        return if (self.source == .zfi) self.name_blob.? else self.fai_data.?;
    }

    /// Null when the load mode did not keep this record's name.
    pub fn recordName(self: *const LoadedIndex, rec_idx: usize) ?[]const u8 {
        if (rec_idx >= self.records.len() or self.on_demand != null) return null;
        if (self.name_slices.len > 0) {
            return self.name_slices[rec_idx];
        }
        if (self.fai_lines.len != 0) return faiLineName(self.fai_data.?[self.fai_lines[rec_idx].offset..]);
        const rec = self.records.get(rec_idx);
        if (self.source == .zfi) {
            if (nameExtentDefect(rec, self.name_blob.?.len) != null) return null;
            return rec.getName(self.name_blob.?);
        }
        if (self.fai_data) |data| return rec.getName(data);
        return null;
    }

    fn checkedFaiRecord(self: *const LoadedIndex, rec_idx: usize, context: *Diagnostics) diagnostics.Error!IndexRecord {
        const data = self.fai_data.?;
        const line = self.fai_lines[rec_idx];
        const rest = data[line.offset..];
        context.at = .{ .source = .fai, .record = rec_idx + 1, .location = .{ .offset = line.offset, .line = line.row, .column = 1 } };
        const fields = try parseFaiIndexLine(rest[0 .. std.mem.findScalar(u8, rest, '\n') orelse rest.len], context);
        const rec = IndexRecord{ .name_offset = line.offset, .name_len = fields.name_len, .seq_offset = fields.seq_offset, .seq_len = fields.seq_len, .line_bases = fields.line_bases, .line_bytes = fields.line_bytes };
        const word = &self.checked_records[rec_idx / @bitSizeOf(usize)];
        const bit = @as(usize, 1) << @intCast(rec_idx % @bitSizeOf(usize));
        if (word.* & bit == 0) {
            try checkUniformGeometry(rec, self.fasta_size, context);
            word.* |= bit;
        }
        return rec;
    }

    pub fn deinit(self: *LoadedIndex) void {
        if (self.on_demand) |zfi| {
            zfi.file.close(self.io);
            zfi.arena.deinit();
            zfi.prefetch_arena.deinit();
            zfi.backing_allocator.destroy(zfi);
        }
        if (self.zfi_map) |*m| m.destroy(self.io);
        if (self.fai_map) |*m| m.destroy(self.io);
        self.arena.deinit();
        self.fasta_file.close(self.io);
    }

    /// The first record named `name`; null also when a `.zfi` name table entry is damaged.
    pub fn lookupName(self: *const LoadedIndex, name: []const u8) ?usize {
        if (self.on_demand != null) {
            var failure: diagnostics.Failure = undefined;
            failure.reset(.get);
            return self.findName(name, &failure) catch null;
        }
        if (self.name_table.len != 0) return switch (self.probeName(name)) {
            .found => |ordinal| ordinal,
            else => null,
        };
        if (self.name_map) |name_map| {
            const rec_idx = name_map.get(name) orelse return null;
            return rec_idx;
        }

        if (self.name_slices.len > 0) {
            for (self.name_slices, 0..) |rec_name, i| {
                if (std.mem.eql(u8, rec_name, name)) {
                    return i;
                }
            }
            return null;
        }

        for (0..self.records.len()) |i| {
            if (std.mem.eql(u8, self.records.get(i).getName(self.recordNameData()), name)) return i;
        }
        return null;
    }

    pub fn sideTableLines(self: *const LoadedIndex, rec: IndexRecord) []const SideTableLine {
        if (rec.isUniformWidth()) return &.{};
        if (self.on_demand) |zfi| return zfi.sideTable(rec) orelse &.{};
        const zfi = self.zfi_data orelse return &.{};
        const offset = std.math.cast(usize, rec.sideTableOffset()) orelse return &.{};
        if (offset < self.zfi_side_start or offset >= self.zfi_side_end) return &.{};
        if (offset % @alignOf(u64) != 0) return &.{};
        const count_end = std.math.add(usize, offset, @sizeOf(u64)) catch return &.{};
        if (count_end > self.zfi_side_end) return &.{};
        const count_ptr: *const u64 = @ptrCast(@alignCast(zfi[offset..].ptr));
        const line_count = std.math.cast(usize, count_ptr.*) orelse return &.{};
        const table_bytes = std.math.mul(usize, line_count, @sizeOf(SideTableLine)) catch return &.{};
        const table_end = std.math.add(usize, count_end, table_bytes) catch return &.{};
        if (line_count == 0 or table_end > self.zfi_side_end) return &.{};
        if (count_end % @alignOf(SideTableLine) != 0) return &.{};
        return @as(
            [*]const SideTableLine,
            @ptrCast(@alignCast(zfi[count_end..].ptr)),
        )[0..line_count];
    }
};

/// A `.zfi` loaded for `get`, read through small positional reads: memory stays with what the
/// requests reach. Each record and side table is checked on first use, and each name-table entry
/// and name when a lookup reaches it, with the diagnostics a whole-file check would report.
const ZfiOnDemand = struct {
    io: std.Io,
    file: std.Io.File,
    backing_allocator: std.mem.Allocator,
    /// Side tables and the name table copy; lives at the struct's fixed heap address.
    arena: std.heap.ArenaAllocator,
    record_count: u32,
    fasta_size: u64,
    side_start: usize,
    side_end: usize,
    blob_start: usize,
    blob_len: usize,
    table_start: usize,
    table_slots: u64,
    side_tables: std.AutoHashMapUnmanaged(u64, CheckedSides) = .empty,
    /// The record a lookup last read, which a request checks next.
    last: ?struct { index: usize, record: IndexRecord } = null,
    /// What `prefetch` found for each request of the current batch, with the start of its name.
    prefetched: []const Prefetched = &.{},
    prefetched_names: []const [*]const u8 = &.{},
    prefetched_candidates: []const Candidate = &.{},
    prefetch_arena: std.heap.ArenaAllocator,
    name: [4096]u8 = undefined,

    // Reading the bytes between two wanted entries costs about what a second read does up to here.
    const PREFETCH_GAP_BYTES = 4096;
    const PREFETCH_PROBE_SLOTS = 16;

    const CheckedSides = struct { seq_offset: u64, seq_len: u64, lines: []const SideTableLine };

    /// `unknown` leaves the request to `find`, which reports any damage on its path; `found` names
    /// the candidate that matched.
    const Prefetched = struct {
        state: enum(u8) { unknown, missing, found } = .missing,
        candidate: u32 = 0,
    };

    /// A table entry whose stored hash bits match a requested name.
    const Candidate = struct { request: u32, ordinal: u32, record: ZfiRecord = undefined };

    fn read(self: *ZfiOnDemand, offset: u64, buffer: []u8, context: *Diagnostics) diagnostics.Error![]u8 {
        const got = std.Io.File.readPositionalAll(self.file, self.io, buffer, offset) catch |err| return context.read(offset, buffer.len, null, err);
        if (got != buffer.len) return context.read(offset, buffer.len, got, null);
        return buffer;
    }

    fn record(self: *ZfiOnDemand, index: usize, context: *Diagnostics) diagnostics.Error!IndexRecord {
        context.at = recordReference(index);
        if (self.last) |last| if (last.index == index) return last.record;
        var bytes: [ZFI_INDEX_RECORD_BYTES]u8 = undefined;
        const rec = std.mem.bytesToValue(ZfiRecord, try self.read(ZFI_HEADER_BYTES + index * ZFI_INDEX_RECORD_BYTES, &bytes, context)).decode();
        self.last = .{ .index = index, .record = rec };
        return rec;
    }

    fn checked(self: *ZfiOnDemand, index: usize, context: *Diagnostics) diagnostics.Error!IndexRecord {
        const rec = try self.record(index, context);
        try checkNameExtent(rec, self.blob_len, context);
        if (rec.isUniformWidth()) {
            try checkUniformGeometry(rec, self.fasta_size, context);
            return rec;
        }
        try checkSequenceStart(rec, self.fasta_size, context);
        const offset = rec.sideTableOffset();
        if (self.side_tables.get(offset)) |sides| {
            if (sides.seq_offset == rec.seq_offset and sides.seq_len == rec.seq_len) return rec;
        }
        const count_end = try sideTableExtent(rec, self.side_start, self.side_end, null, context);
        context.at.location.offset = offset;
        var count_bytes: [8]u8 = undefined;
        const count = std.mem.readInt(u64, (try self.read(offset, &count_bytes, context))[0..8], .little);
        const table_end = try sideTableExtent(rec, self.side_start, self.side_end, count, context);
        const lines = self.arena.allocator().alloc(SideTableLine, @intCast(count)) catch return context.allocation(.side_tables, table_end - count_end);
        context.at.location.offset = count_end;
        _ = try self.read(count_end, std.mem.sliceAsBytes(lines), context);
        try self.keepSideTable(rec, count_end, lines, context);
        return rec;
    }

    // Checks `lines`, the side table of `rec` that starts at `count_end`, and keeps it for the
    // record's later requests.
    fn keepSideTable(self: *ZfiOnDemand, rec: IndexRecord, count_end: usize, lines: []const SideTableLine, context: *Diagnostics) diagnostics.Error!void {
        var validation = SideTableValidation{};
        for (lines, 0..) |line, i| {
            context.at.location.offset = count_end + i * @sizeOf(SideTableLine);
            try validation.accept(line, rec, self.fasta_size, i + 1, context);
        }
        if (validation.expected_base_start != rec.seq_len) return context.fail(.{ .@"index.side_sequence_length_mismatch" = .{ .indexed_length = rec.seq_len, .side_table_bases = validation.expected_base_start } });
        self.side_tables.put(self.arena.allocator(), rec.sideTableOffset(), .{ .seq_offset = rec.seq_offset, .seq_len = rec.seq_len, .lines = lines }) catch return context.allocation(.side_tables, null);
    }

    // The side tables of the records `prefetch` found, in file order through merged reads; a table
    // that is damaged or does not fit a read is left to `checked`, which reports why.
    fn prefetchSideTables(self: *ZfiOnDemand, keys: []u64, scratch: []u64, candidates: []const Candidate, buffer: []u8) error{OutOfMemory}!void {
        radixSort(keys, scratch[0..keys.len], 24);
        var failure: diagnostics.Failure = undefined;
        failure.reset(.get);
        var context = Diagnostics{ .failure = &failure, .at = .{ .source = .zfi } };
        var i: usize = 0;
        while (i < keys.len) {
            const window_start: usize = @intCast(keys[i] >> 24);
            var last = window_start;
            var k = i + 1;
            while (k < keys.len and (keys[k] >> 24) - last <= 4 * PREFETCH_GAP_BYTES and (keys[k] >> 24) - window_start < buffer.len) : (k += 1) last = @intCast(keys[k] >> 24);
            const window_len = @min(buffer.len, self.side_end - window_start, last - window_start + 4 * PREFETCH_GAP_BYTES);
            const got = std.Io.File.readPositionalAll(self.file, self.io, buffer[0..window_len], window_start) catch 0;
            if (got != window_len) return;
            var j = i;
            while (j < k) : (j += 1) {
                const rec = candidates[@intCast(keys[j] & 0xFF_FFFF)].record.decode();
                const offset: usize = @intCast(keys[j] >> 24);
                const count_end = sideTableExtent(rec, self.side_start, self.side_end, null, &context) catch continue;
                if (count_end - window_start > window_len) break;
                const count = std.mem.readInt(u64, buffer[offset - window_start ..][0..8], .little);
                const table_end = sideTableExtent(rec, self.side_start, self.side_end, count, &context) catch continue;
                if (table_end - window_start > window_len) break;
                const lines = try self.arena.allocator().alloc(SideTableLine, @intCast(count));
                @memcpy(std.mem.sliceAsBytes(lines), buffer[count_end - window_start .. table_end - window_start]);
                self.keepSideTable(rec, count_end, lines, &context) catch |err| if (err == error.OutOfMemory) return error.OutOfMemory;
            }
            i = @max(j, i + 1);
        }
    }

    /// Resolves `names` with reads sorted by position and merged across gaps of at most
    /// `PREFETCH_GAP_BYTES`, through `buffer`, so a large request batch costs few reads. The results
    /// borrow the bytes of `names`. Best effort: a name whose table entries, record, or bytes are
    /// damaged or unreadable stays `unknown`, and its sequential lookup reports why.
    fn prefetch(self: *ZfiOnDemand, names: []const []const u8, buffer: []u8) error{OutOfMemory}!void {
        self.clearPrefetch();
        // Each sort key holds a position (home slot, ordinal, name offset) above an item index.
        const request_bits = 30;
        if (names.len >= 1 << request_bits or buffer.len < PREFETCH_PROBE_SLOTS * 4) return;
        const allocator = self.prefetch_arena.allocator();
        const count: u64 = self.record_count;
        const ordinal_bits = NameTable.ordinalBits(count);
        const results = try allocator.alloc(Prefetched, names.len);
        @memset(results, .{});
        const wanted = try allocator.alloc(u32, names.len);
        const keys = try allocator.alloc(u64, names.len);
        var scratch = try allocator.alloc(u64, names.len);
        for (names, wanted, keys, 0..) |name, *fingerprint, *key, i| {
            const name_hash = NameTable.hash(name);
            fingerprint.* = @intCast(NameTable.fingerprint(name_hash, ordinal_bits));
            key.* = NameTable.home(name_hash, self.table_slots) << request_bits | i;
        }
        radixSort(keys, scratch, request_bits);

        // Name table. A window starts at the first slot a probe needs and reaches the homes that
        // follow it closely; a probe that runs past its end reads a new one.
        var candidates: std.ArrayList(Candidate) = .empty;
        const buffer_slots = buffer.len / 4;
        var window_start: u64 = 0;
        var window_len: u64 = 0;
        // A repeated name shares the lookup of its first request; repeats sort next to each other.
        const leader = try allocator.alloc(u32, names.len);
        var previous: ?u32 = null;
        var previous_home: u64 = 0;
        for (keys, 0..) |key, k| {
            const request: u32 = @intCast(key & ((1 << request_bits) - 1));
            var slot = key >> request_bits;
            if (previous) |first| {
                if (previous_home == slot and std.mem.eql(u8, names[first], names[request])) {
                    leader[request] = first;
                    continue;
                }
            }
            leader[request] = request;
            previous = request;
            previous_home = slot;
            var steps: u64 = 0;
            while (true) : (steps += 1) {
                if (slot == self.table_slots) slot = 0;
                if (steps == self.table_slots) {
                    results[request].state = .unknown;
                    break;
                }
                if (slot < window_start or slot >= window_start + window_len) {
                    var end = slot + PREFETCH_PROBE_SLOTS;
                    var next = k + 1;
                    while (next < keys.len) : (next += 1) {
                        const home = keys[next] >> request_bits;
                        if (home < slot or home > end + PREFETCH_GAP_BYTES / 4 or home + PREFETCH_PROBE_SLOTS - slot > buffer_slots) break;
                        end = @max(end, home + PREFETCH_PROBE_SLOTS);
                    }
                    window_start = slot;
                    window_len = @min(end, self.table_slots) - slot;
                    const bytes = buffer[0..@intCast(window_len * 4)];
                    const got = std.Io.File.readPositionalAll(self.file, self.io, bytes, self.table_start + slot * 4) catch 0;
                    if (got != bytes.len) {
                        window_len = 0;
                        results[request].state = .unknown;
                        break;
                    }
                }
                const value: u64 = std.mem.readInt(u32, buffer[@intCast((slot - window_start) * 4)..][0..4], .little);
                if (value == 0) break;
                const ordinal_plus_one = value & ((@as(u64, 1) << ordinal_bits) - 1);
                if (ordinal_plus_one == 0 or ordinal_plus_one > count) {
                    results[request].state = .unknown;
                    break;
                }
                // The first entry with the name's hash bits is almost always the name: other names
                // share those bits once in 2^(32 - ordinal bits) entries. When it is not, the name
                // goes to `find`, which walks on.
                if (value >> ordinal_bits == wanted[request]) {
                    try candidates.append(allocator, .{ .request = request, .ordinal = @intCast(ordinal_plus_one - 1) });
                    break;
                }
                slot += 1;
            }
        }
        if (candidates.items.len >= 1 << 24) return self.clearPrefetch();

        // Records, in ordinal order; a candidate whose name length matches moves on to the names.
        const order = try allocator.alloc(u64, candidates.items.len);
        if (scratch.len < order.len) scratch = try allocator.alloc(u64, order.len);
        for (order, candidates.items, 0..) |*key, candidate, c| key.* = @as(u64, candidate.ordinal) << 32 | c;
        radixSort(order, scratch[0..order.len], 32);
        var named: usize = 0;
        var i: usize = 0;
        while (i < order.len) {
            const first: u64 = order[i] >> 32;
            var j = i;
            while (j + 1 < order.len) : (j += 1) {
                const next = order[j + 1] >> 32;
                if ((next - (order[j] >> 32)) * ZFI_INDEX_RECORD_BYTES > PREFETCH_GAP_BYTES or (next - first + 1) * ZFI_INDEX_RECORD_BYTES > buffer.len) break;
            }
            const span: usize = @intCast(((order[j] >> 32) - first + 1) * ZFI_INDEX_RECORD_BYTES);
            const got = std.Io.File.readPositionalAll(self.file, self.io, buffer[0..span], ZFI_HEADER_BYTES + first * ZFI_INDEX_RECORD_BYTES) catch 0;
            for (order[i .. j + 1]) |key| {
                const candidate = &candidates.items[@intCast(key & 0xFFFF_FFFF)];
                if (got != span) {
                    results[candidate.request].state = .unknown;
                    continue;
                }
                const at: usize = @intCast((candidate.ordinal - first) * ZFI_INDEX_RECORD_BYTES);
                candidate.record = std.mem.bytesToValue(ZfiRecord, buffer[at..][0..ZFI_INDEX_RECORD_BYTES]);
                const rec = candidate.record.decode();
                if (nameExtentDefect(rec, self.blob_len) != null) {
                    results[candidate.request].state = .unknown;
                    continue;
                }
                if (rec.name_len != names[candidate.request].len) {
                    results[candidate.request].state = .unknown;
                    continue;
                }
                // Name offsets are below 2^40, so 24 bits remain for the candidate.
                order[named] = rec.name_offset << 24 | (key & 0xFFFF_FFFF);
                named += 1;
            }
            i = j + 1;
        }

        // Names, in blob order; a stored name that differs leaves its request to `find`.
        const matches = order[0..named];
        radixSort(matches, scratch[0..named], 24);
        const nameEnd = struct {
            fn get(key: u64, all: []const Candidate, requested: []const []const u8) u64 {
                return (key >> 24) + requested[all[@intCast(key & 0xFF_FFFF)].request].len;
            }
        }.get;
        i = 0;
        while (i < matches.len) {
            const first = matches[i] >> 24;
            var end = nameEnd(matches[i], candidates.items, names);
            var j = i;
            while (j + 1 < matches.len) : (j += 1) {
                const next = matches[j + 1] >> 24;
                const next_end = @max(end, nameEnd(matches[j + 1], candidates.items, names));
                if (next > end + PREFETCH_GAP_BYTES or next_end - first > buffer.len) break;
                end = next_end;
            }
            const span: usize = @intCast(end - first);
            const got = std.Io.File.readPositionalAll(self.file, self.io, buffer[0..span], self.blob_start + first) catch 0;
            for (matches[i .. j + 1]) |key| {
                const candidate = candidates.items[@intCast(key & 0xFF_FFFF)];
                const result = &results[candidate.request];
                if (got != span) {
                    result.state = .unknown;
                    continue;
                }
                const name = names[candidate.request];
                const stored = buffer[@intCast((key >> 24) - first)..][0..name.len];
                if (result.state == .unknown) continue;
                result.* = if (std.mem.eql(u8, stored, name)) .{ .state = .found, .candidate = @intCast(key & 0xFF_FFFF) } else .{ .state = .unknown };
            }
            i = j + 1;
        }
        var sides: std.ArrayList(u64) = .empty;
        for (results, leader, 0..) |result, first, request| {
            if (first != request or result.state != .found) continue;
            const rec = candidates.items[result.candidate].record.decode();
            if (rec.isUniformWidth() or self.side_tables.contains(rec.sideTableOffset())) continue;
            // Side-table offsets are below 2^39, so 24 bits remain for the candidate.
            try sides.append(allocator, rec.sideTableOffset() << 24 | result.candidate);
        }
        if (sides.items.len != 0) {
            if (scratch.len < sides.items.len) scratch = try allocator.alloc(u64, sides.items.len);
            try self.prefetchSideTables(sides.items, scratch, candidates.items, buffer);
        }
        for (results, leader, 0..) |*result, first, request| {
            if (first != request) result.* = results[first];
        }
        const starts = try allocator.alloc([*]const u8, names.len);
        for (starts, names) |*start, name| start.* = name.ptr;
        self.prefetched_names = starts;
        self.prefetched_candidates = candidates.items;
        self.prefetched = results;
    }

    fn clearPrefetch(self: *ZfiOnDemand) void {
        self.prefetched = &.{};
        self.prefetched_names = &.{};
        self.prefetched_candidates = &.{};
        _ = self.prefetch_arena.reset(.retain_capacity);
        self.last = null;
    }

    fn findAt(self: *ZfiOnDemand, request: usize, name: []const u8, context: *Diagnostics) diagnostics.Error!?usize {
        if (request < self.prefetched.len and self.prefetched_names[request] == name.ptr) {
            const result = self.prefetched[request];
            switch (result.state) {
                .found => {
                    const candidate = self.prefetched_candidates[result.candidate];
                    self.last = .{ .index = candidate.ordinal, .record = candidate.record.decode() };
                    return candidate.ordinal;
                },
                .missing => return null,
                .unknown => {},
            }
        }
        return self.find(name, context);
    }

    // Compares the stored name at `offset` with `name`, a buffer at a time.
    fn nameEquals(self: *ZfiOnDemand, offset: u64, name: []const u8, context: *Diagnostics) diagnostics.Error!bool {
        var done: usize = 0;
        while (done < name.len) {
            const len = @min(self.name.len, name.len - done);
            const stored = try self.read(self.blob_start + offset + done, self.name[0..len], context);
            if (!std.mem.eql(u8, stored, name[done..][0..len])) return false;
            done += len;
        }
        return true;
    }

    fn sideTable(self: *ZfiOnDemand, rec: IndexRecord) ?[]const SideTableLine {
        const sides = self.side_tables.get(rec.sideTableOffset()) orelse return null;
        if (sides.seq_offset != rec.seq_offset or sides.seq_len != rec.seq_len) return null;
        return sides.lines;
    }

    fn find(self: *ZfiOnDemand, name: []const u8, context: *Diagnostics) diagnostics.Error!?usize {
        const count: u64 = self.record_count;
        const ordinal_bits = NameTable.ordinalBits(count);
        const name_hash = NameTable.hash(name);
        const wanted = NameTable.fingerprint(name_hash, ordinal_bits);
        var slot = NameTable.home(name_hash, self.table_slots);
        var window_buffer: [16]u32 = undefined;
        var window: []const u32 = &.{};
        var window_start: u64 = 0;
        for (0..self.table_slots) |_| {
            if (slot < window_start or slot >= window_start + window.len) {
                const len = @min(window_buffer.len, self.table_slots - slot);
                context.at = .{ .source = .zfi, .location = .{ .offset = self.table_start + slot * 4 } };
                _ = try self.read(self.table_start + slot * 4, std.mem.sliceAsBytes(window_buffer[0..@intCast(len)]), context);
                window = window_buffer[0..@intCast(len)];
                window_start = slot;
            }
            const value: u64 = std.mem.littleToNative(u32, window[@intCast(slot - window_start)]);
            if (value == 0) return null;
            const ordinal_plus_one = value & ((@as(u64, 1) << ordinal_bits) - 1);
            if (ordinal_plus_one == 0 or ordinal_plus_one > count) {
                context.at = .{ .source = .zfi, .location = .{ .offset = self.table_start + slot * 4 } };
                return context.fail(.{ .@"index.name_table_invalid" = .{ .slot = slot } });
            }
            if (value >> ordinal_bits == wanted) {
                const ordinal: usize = @intCast(ordinal_plus_one - 1);
                const rec = try self.record(ordinal, context);
                try checkNameExtent(rec, self.blob_len, context);
                if (rec.name_len == name.len and try self.nameEquals(rec.name_offset, name, context)) return ordinal;
            }
            slot = if (slot + 1 == self.table_slots) 0 else slot + 1;
        }
        return null;
    }
};

/// Sorts `keys` ascending by their bits from `low_bit` up, 11 bits per pass, or 8 below 4,096 keys,
/// where clearing the larger count table would cost more than the extra passes. Each pass is
/// stable, so keys equal above `low_bit` keep their order: callers put an item index below it.
/// For the key counts a request batch makes, this beats `std.mem.sortUnstable` several times over.
fn radixSort(keys: []u64, scratch: []u64, low_bit: u7) void {
    std.debug.assert(scratch.len == keys.len);
    var used: u64 = 0;
    for (keys) |key| used |= key;
    const digit_bits: u7 = if (keys.len < 4096) 8 else 11;
    const digits: usize = @as(usize, 1) << @intCast(digit_bits);
    var source = keys;
    var target = scratch;
    var shift: u7 = low_bit;
    while (shift < 64 - @clz(used)) : (shift += digit_bits) {
        var starts: [2048]u32 = undefined;
        @memset(starts[0..digits], 0);
        for (source) |key| starts[@intCast(key >> @intCast(shift) & (digits - 1))] += 1;
        var total: u32 = 0;
        for (starts[0..digits]) |*start| {
            const n = start.*;
            start.* = total;
            total += n;
        }
        for (source) |key| {
            const digit: usize = @intCast(key >> @intCast(shift) & (digits - 1));
            target[starts[digit]] = key;
            starts[digit] += 1;
        }
        std.mem.swap([]u64, &source, &target);
    }
    if (source.ptr != keys.ptr) @memcpy(keys, source);
}

/// Every index entry in file order, through fixed buffers. `open` makes `loadIndex`'s file checks;
/// `next` checks each entry as `loadIndex` does, so a damaged one fails where it is read. An
/// `EntryStream` must not move after `open`.
pub const EntryStream = struct {
    io: std.Io,
    file: std.Io.File,
    source: LoadedIndex.IndexSource,
    context: Diagnostics,
    fasta_size: u64,
    reader: std.Io.File.Reader,
    names: ZfiReadWindow = undefined,
    sides: ZfiReadWindow = undefined,
    buffers: [3][FAI_READER_BUFFER_BYTES]u8 = undefined,
    records_end: usize = 0,
    blob_start: usize = 0,
    blob_len: usize = 0,
    table_start: usize = 0,
    /// The name table rebuilt from the entries read so far; compared with the stored one at the end.
    table: ?NameTable.Builder = null,
    allocator: std.mem.Allocator = undefined,
    previous_side_end: usize = 0,
    remaining: u32 = 0,
    ordinal: u64 = 0,
    row: u64 = 0,
    entry: Entry = undefined,

    /// `next` returns the stream's copy, valid until the next call. `name` borrows the stream;
    /// `offset` is where the entry starts in the index file, and `ordinal` counts entries from 1.
    pub const Entry = struct { record: IndexRecord, name: []const u8, offset: u64, ordinal: u64 };

    pub fn open(self: *EntryStream, allocator: std.mem.Allocator, io: std.Io, fasta_path: []const u8, selection: IndexSelection, failure: *diagnostics.Failure) diagnostics.Error!void {
        var shell = try loadIndex(allocator, io, fasta_path, selection, .{ .stream = .{ .entries = self } }, failure);
        shell.deinit();
    }

    /// `open` without the `.zfi` name-table check, keeping the FASTA the index was checked against
    /// open for the caller, who closes it. `.zfi` entries then come without names; `entryName`
    /// reads one.
    pub fn openWithFasta(self: *EntryStream, allocator: std.mem.Allocator, io: std.Io, fasta_path: []const u8, selection: IndexSelection, failure: *diagnostics.Failure) diagnostics.Error!std.Io.File {
        var shell = try loadIndex(allocator, io, fasta_path, selection, .{ .stream = .{ .entries = self, .check_name_table = false } }, failure);
        shell.arena.deinit();
        return shell.fasta_file;
    }

    pub fn close(self: *EntryStream) void {
        if (self.table) |*table| table.deinit(self.allocator);
        self.file.close(self.io);
        self.* = undefined;
    }

    pub fn next(self: *EntryStream) diagnostics.Error!?*const Entry {
        return switch (self.source) {
            .zfi => self.nextZfi(),
            .fai => self.nextFai(),
        };
    }

    /// The entry's name, read from a `.zfi` opened by `openWithFasta`; it borrows the stream until
    /// the next call.
    pub fn entryName(self: *EntryStream, entry: *const Entry) diagnostics.Error![]const u8 {
        if (self.source == .fai or self.table != null) return entry.name;
        self.context.at = .{ .source = .zfi, .record = entry.ordinal, .location = .{ .offset = entry.offset } };
        return self.names.take(self.blob_start + entry.record.name_offset, entry.record.name_len, &self.context);
    }

    /// The side table of a `.zfi` entry that does not use the uniform formula: its line count.
    pub fn sideLineCount(self: *EntryStream, record: IndexRecord) diagnostics.Error!u64 {
        return std.mem.readInt(u64, (try self.sides.take(record.sideTableOffset(), 8, &self.context))[0..8], .little);
    }

    /// Line `index`, from 0, of the entry's side table; lines read in order cost one read per window.
    pub fn sideLine(self: *EntryStream, record: IndexRecord, index: u64) diagnostics.Error!SideTableLine {
        const offset = record.sideTableOffset() + 8 + index * @sizeOf(SideTableLine);
        return std.mem.bytesToValue(SideTableLine, try self.sides.take(offset, @sizeOf(SideTableLine), &self.context));
    }

    fn start(self: *EntryStream, io: std.Io, path: []const u8, source: LoadedIndex.IndexSource, fasta_size: u64, context: *Diagnostics) diagnostics.Error!void {
        const file = platform.openSidecar(io, path) catch |err| return context.ioFailure(.open, err);
        self.* = .{ .io = io, .file = file, .source = source, .context = context.*, .fasta_size = fasta_size, .reader = undefined };
        self.reader = file.reader(io, &self.buffers[0]);
    }

    // The checks `tryLoadZfi` makes before its records, in its order.
    fn startZfi(self: *EntryStream, allocator: std.mem.Allocator, io: std.Io, path: []const u8, zfi_size: u64, fasta_size: u64, fasta_mtime_ns: i96, check_name_table: bool, context: *Diagnostics) diagnostics.Error!void {
        try self.start(io, path, .zfi, fasta_size, context);
        errdefer self.file.close(io);
        const size = std.math.cast(usize, zfi_size) orelse return context.capacity(.map_bytes, std.math.maxInt(usize), zfi_size);
        context.at.location.offset = 0;
        if (size < @sizeOf(ZfiHeader)) return context.truncated(.header, ZFI_HEADER_BYTES, size);
        const header_bytes = self.reader.interface.takeArray(ZFI_HEADER_BYTES) catch |err| switch (err) {
            error.ReadFailed => return context.ioFailure(.read, self.reader.err orelse err),
            error.EndOfStream => return context.read(0, ZFI_HEADER_BYTES, self.reader.interface.buffered().len, null),
        };
        const header = std.mem.bytesToValue(ZfiHeader, header_bytes);
        self.records_end = try checkZfiHeader(header, fasta_size, size, context);
        self.names = .{ .file = self.file, .io = io, .buffer = &self.buffers[1] };
        self.sides = .{ .file = self.file, .io = io, .buffer = &self.buffers[2] };
        var tail: [ZFI_TRAILER_BYTES]u8 = undefined;
        const tail_start = size - @min(tail.len, size);
        const tail_bytes = tail[0 .. size - tail_start];
        @memcpy(tail_bytes, try self.sides.take(tail_start, tail_bytes.len, context));
        context.at = .{ .source = .zfi, .location = .{ .offset = size } };
        const blob = switch (locateZfiNameBlob(tail_bytes, tail_start, size, header.record_count, self.records_end, context)) {
            .invalid => |details| return context.fail(details),
            .blob => |blob| blob,
        };
        context.at.location.offset = blob.trailing.source_mtime_offset;
        try checkSourceTime(blob.trailing.source_mtime_ns, fasta_mtime_ns, context);
        self.blob_start = blob.regions.blob_start;
        self.blob_len = blob.regions.blob_len;
        self.table_start = blob.regions.table_start;
        self.allocator = allocator;
        if (check_name_table) self.table = NameTable.Builder.init(allocator, header.record_count) catch return context.allocation(.name_map, blob.regions.table_len);
        self.previous_side_end = self.records_end;
        self.reader.seekTo(ZFI_HEADER_BYTES) catch |err| return context.ioFailure(.read, err);
        self.remaining = header.record_count;
    }

    // The checks `validateZfiRecords` makes for one record, in its order.
    fn nextZfi(self: *EntryStream) diagnostics.Error!?*const Entry {
        if (self.remaining == 0) {
            if (self.table) |*table| {
                defer {
                    table.deinit(self.allocator);
                    self.table = null;
                }
                try self.compareTable(table.bytes());
            }
            return null;
        }
        const context = &self.context;
        const offset = self.reader.logicalPos();
        const bytes = self.reader.interface.takeArray(ZFI_INDEX_RECORD_BYTES) catch |err| {
            context.at = .{ .source = .zfi, .record = self.ordinal + 1, .location = .{ .offset = offset } };
            switch (err) {
                error.ReadFailed => return context.ioFailure(.read, self.reader.err orelse err),
                error.EndOfStream => return context.read(offset, ZFI_INDEX_RECORD_BYTES, self.reader.interface.buffered().len, null),
            }
        };
        const record = std.mem.bytesToValue(ZfiRecord, bytes).decode();
        // A sound uniform record, nearly every one, is checked without the diagnostic context.
        const plain = record.isUniformWidth() and nameExtentDefect(record, self.blob_len) == null and uniformGeometryDefect(record, self.fasta_size) == null;
        if (!plain) {
            context.at = .{ .source = .zfi, .record = self.ordinal + 1, .location = .{ .offset = offset } };
            try checkNameExtent(record, self.blob_len, context);
            if (record.isUniformWidth()) {
                try checkUniformGeometry(record, self.fasta_size, context);
            } else {
                try checkSequenceStart(record, self.fasta_size, context);
                const end = try checkSideTable(&self.sides, record, self.records_end, self.blob_start, self.fasta_size, context);
                const table = record.sideTableOffset();
                if (table < self.previous_side_end) {
                    context.at.location.offset = table;
                    return context.fail(.{ .@"index.side_table_order_invalid" = .{ .table_start = table, .minimum_start = self.previous_side_end } });
                }
                self.previous_side_end = end;
            }
        }
        const name = if (self.table) |*table| named: {
            context.at = .{ .source = .zfi, .record = self.ordinal + 1, .location = .{ .offset = offset } };
            const name = try self.names.take(self.blob_start + record.name_offset, record.name_len, context);
            table.insert(name);
            break :named name;
        } else "";
        self.remaining -= 1;
        self.ordinal += 1;
        self.entry = .{ .record = record, .name = name, .offset = offset, .ordinal = self.ordinal };
        return &self.entry;
    }

    fn compareTable(self: *EntryStream, expected: []const u8) diagnostics.Error!void {
        var done: usize = 0;
        while (done < expected.len) {
            const len = @min(expected.len - done, self.buffers[2].len);
            const stored = try self.sides.take(self.table_start + done, len, &self.context);
            try compareNameTable(expected[done..][0..len], stored, self.table_start, done, &self.context);
            done += len;
        }
    }

    // The checks `tryLoadFai` makes for one row, in its order.
    fn nextFai(self: *EntryStream) diagnostics.Error!?*const Entry {
        const context = &self.context;
        while (true) {
            const offset = self.reader.logicalPos();
            const available = self.reader.interface.buffered();
            var long: ?FaiLineFields = null;
            const maybe_line = if (std.mem.findScalar(u8, available, '\n')) |end| buffered: {
                self.reader.interface.toss(end + 1);
                break :buffered available[0..end];
            } else self.reader.interface.takeDelimiter('\n') catch |err| switch (err) {
                error.ReadFailed => return context.ioFailure(.read, self.reader.err orelse err),
                error.StreamTooLong => long_line: {
                    self.row = try context.nextOrdinal(self.row, .physical_lines);
                    context.at = .{ .source = .fai, .record = self.ordinal + 1, .location = .{ .offset = offset, .line = self.row, .column = 1 } };
                    const parsed = try takeLongFaiLine(&self.reader, &self.buffers[1], context);
                    long = parsed.fields;
                    break :long_line parsed.name;
                },
            };
            const line = maybe_line orelse {
                if (self.ordinal == 0) return context.fail(.{ .@"index.no_records" = .{ .format = .fai } });
                return null;
            };
            if (long == null) {
                self.row = try context.nextOrdinal(self.row, .physical_lines);
                context.at = .{ .source = .fai, .record = self.ordinal + 1, .location = .{ .offset = offset, .line = self.row, .column = 1 } };
                if (line.len == 0) continue;
            }
            const fields = long orelse try parseFaiIndexLine(line, context);
            const record = IndexRecord{ .name_offset = offset, .name_len = fields.name_len, .seq_offset = fields.seq_offset, .seq_len = fields.seq_len, .line_bases = fields.line_bases, .line_bytes = fields.line_bytes };
            try checkUniformGeometry(record, self.fasta_size, context);
            self.ordinal += 1;
            self.entry = .{ .record = record, .name = line[0..fields.name_end], .offset = offset, .ordinal = self.ordinal };
            return &self.entry;
        }
    }
};

// What `loadIndex` returns in `.stream` mode: the open FASTA, with every entry left to the stream.
fn streamShell(allocator: std.mem.Allocator, io: std.Io, fasta_file: std.Io.File, fasta_size: u64, source: LoadedIndex.IndexSource) LoadedIndex {
    return .{ .io = io, .fasta_file = fasta_file, .records = .{ .owned = &.{} }, .fasta_size = fasta_size, .zfi_data = null, .source = source, .arena = .init(allocator) };
}

pub fn loadIndex(
    allocator: std.mem.Allocator,
    io: std.Io,
    fasta_path: []const u8,
    selection: IndexSelection,
    mode: LoadMode,
    failure: *diagnostics.Failure,
) diagnostics.Error!LoadedIndex {
    // A second load of the same index in one operation keeps the first registration.
    switch (selection) {
        .automatic => {},
        .zfi => |path| if (failure.sources.get(.zfi) == null) failure.sources.set(.zfi, .{ .label = .{ .borrowed = path }, .kind = .path }),
        .fai => |path| if (failure.sources.get(.fai) == null) failure.sources.set(.fai, .{ .label = .{ .borrowed = path }, .kind = .path }),
    }
    inline for (.{ diagnostics.SourceId.fasta, diagnostics.SourceId.zfi, diagnostics.SourceId.fai }) |source| {
        if (failure.sources.get(source) == null) failure.sources.set(source, .{
            .label = .{ .borrowed = fasta_path },
            .kind = .path,
            .suffix = switch (source) {
                .zfi => .zfi,
                .fai => .fai,
                else => .none,
            },
        });
    }
    var context = Diagnostics{ .failure = failure, .at = .{ .source = .fasta } };
    const fasta_file = std.Io.Dir.cwd().openFile(io, fasta_path, .{}) catch |err| return context.ioFailure(.open, err);
    errdefer fasta_file.close(io);
    const fasta_stat = fasta_file.stat(io) catch |err| return context.ioFailure(.stat, err);
    if (fasta_stat.size == 0) return context.fail(.{ .@"fasta.no_records" = {} });

    switch (selection) {
        .automatic => {},
        .zfi => |path| {
            context.at.source = .zfi;
            return (try tryLoadZfi(allocator, io, path, fasta_file, fasta_stat, mode, &context)) orelse context.ioFailure(.open, error.FileNotFound);
        },
        .fai => |path| {
            context.at.source = .fai;
            return (try tryLoadFai(allocator, io, path, fasta_file, fasta_stat, mode, &context)) orelse context.ioFailure(.open, error.FileNotFound);
        },
    }

    // One buffer serves both sidecar paths: neither load keeps its path.
    var path_buf: [4096]u8 = undefined;
    context.at.source = .zfi;
    const zfi_path = std.fmt.bufPrint(&path_buf, "{s}.zfi", .{fasta_path}) catch return context.pathLimit(path_buf.len);
    if (try tryLoadZfi(allocator, io, zfi_path, fasta_file, fasta_stat, mode, &context)) |loaded| return loaded;

    context.at.source = .fai;
    const fai_path = std.fmt.bufPrint(&path_buf, "{s}.fai", .{fasta_path}) catch return context.pathLimit(path_buf.len);
    return (try tryLoadFai(allocator, io, fai_path, fasta_file, fasta_stat, mode, &context)) orelse context.fail(.{ .@"index.not_found" = .{ .zfi_source = .zfi, .fai_source = .fai } });
}

const NON_UNIFORM_WIDTH_FLAG: u8 = 1;

const SIDE_TABLE_OFFSET_BYTES = 5;

const ZFI_TRAILER_BYTES = ZFI_SOURCE_ID_BYTES + ZFI_NAME_FOOTER_BYTES;

const FAI_READER_BUFFER_BYTES: usize = 65_600;

fn decodeZfiNameFooterBytes(footer_bytes: []const u8) ?u64 {
    if (footer_bytes.len < ZFI_NAME_FOOTER_BYTES) return null;
    if (!std.mem.eql(u8, footer_bytes[0..4], &ZFI_NAME_FOOTER_MAGIC)) return null;
    return std.mem.readInt(u64, footer_bytes[4..12], .little);
}

fn decodeZfiSourceIdBytes(id_bytes: []const u8) ?u64 {
    if (id_bytes.len < ZFI_SOURCE_ID_BYTES) return null;
    if (!std.mem.eql(u8, id_bytes[0..4], &ZFI_SOURCE_ID_MAGIC)) return null;
    return std.mem.readInt(u64, id_bytes[4..12], .little);
}

const ZfiTrailingMeta = struct {
    payload_end: usize,
    name_blob_len: u64,
    source_mtime_ns: u64,
    source_mtime_offset: usize,
};

// The source identity and name footer that end `zfi_data`.
fn parseZfiTrailingMeta(zfi_data: []const u8) ?ZfiTrailingMeta {
    if (zfi_data.len < ZFI_TRAILER_BYTES) return null;
    const trailer = zfi_data[zfi_data.len - ZFI_TRAILER_BYTES ..];
    const mtime_ns = decodeZfiSourceIdBytes(trailer[0..ZFI_SOURCE_ID_BYTES]) orelse return null;
    const name_blob_len = decodeZfiNameFooterBytes(trailer[ZFI_SOURCE_ID_BYTES..]) orelse return null;
    return .{
        .payload_end = zfi_data.len - ZFI_TRAILER_BYTES,
        .name_blob_len = name_blob_len,
        .source_mtime_ns = mtime_ns,
        .source_mtime_offset = zfi_data.len - ZFI_TRAILER_BYTES + 4,
    };
}

comptime {
    if (builtin.cpu.arch.endian() != .little) {
        @compileError(".zfi records are read in place, which needs a little-endian target");
    }
    if (@sizeOf(ZfiHeader) != ZFI_HEADER_BYTES) @compileError("ZfiHeader size drifted from wire contract");
    if (@sizeOf(ZfiRecord) != ZFI_INDEX_RECORD_BYTES or @alignOf(ZfiRecord) != 8) @compileError("ZfiRecord size drifted from wire contract");
    if (@sizeOf(SideTableLine) != ZFI_SIDE_TABLE_LINE_BYTES) @compileError("SideTableLine size drifted from wire contract");
    if (@offsetOf(ZfiHeader, "magic") != 0) @compileError("ZfiHeader.magic offset");
    if (@offsetOf(ZfiHeader, "record_count") != 4) @compileError("ZfiHeader.record_count offset");
    if (@offsetOf(ZfiHeader, "source_size") != 8) @compileError("ZfiHeader.source_size offset");
    if (@offsetOf(SideTableLine, "base_start") != 0) @compileError("SideTableLine.base_start offset");
    if (@offsetOf(SideTableLine, "byte_offset") != 8) @compileError("SideTableLine.byte_offset offset");
    if (@offsetOf(SideTableLine, "line_bytes") != 16) @compileError("SideTableLine.line_bytes offset");
    if (@offsetOf(SideTableLine, "line_bases") != 20) @compileError("SideTableLine.line_bases offset");
}

const Diagnostics = struct {
    failure: *diagnostics.Failure,
    at: diagnostics.Reference,
    deferred: bool = false,

    fn failureStage(details: diagnostics.Details) diagnostics.Stage {
        return switch (std.meta.activeTag(details)) {
            .@"fasta.no_records" => .stat,
            .@"index.not_found" => .open,
            .@"index.older_than_source", .@"index.source_size_mismatch", .@"index.source_mtime_mismatch" => .check_identity,
            else => .parse,
        };
    }

    fn fail(self: *const Diagnostics, details: diagnostics.Details) diagnostics.Error {
        @branchHint(.cold);
        const diagnostic = diagnostics.Diagnostic{ .details = details, .stage = failureStage(details), .at = self.at };
        if (self.deferred) {
            self.failure.deferFailure(diagnostic);
            return diagnostics.metadata(diagnostic.code()).returned_error;
        }
        return self.failure.fail(diagnostic);
    }

    const IoStage = enum { open, stat, read };

    fn ioFailure(self: *const Diagnostics, stage: IoStage, err: anyerror) diagnostics.Error {
        @branchHint(.cold);
        const details: diagnostics.Details = switch (stage) {
            .open => diagnostics.openFailure(err),
            .stat => .{ .@"io.stat_failed" = {} },
            .read => .{ .@"io.read_failed" = .{} },
        };
        return self.failure.fail(.{ .details = details, .stage = switch (stage) {
            .open => .open,
            .stat => .stat,
            .read => .read,
        }, .at = self.at, .cause = .{ .zig = err } });
    }

    fn read(self: *const Diagnostics, offset: u64, requested: u64, received: ?u64, cause: ?anyerror) diagnostics.Error {
        @branchHint(.cold);
        var at = self.at;
        at.location.offset = offset;
        return self.failure.fail(.{ .details = diagnostics.readFailure(requested, received), .stage = .read, .at = at, .cause = if (cause) |err| .{ .zig = err } else null });
    }

    fn allocation(self: *const Diagnostics, role: diagnostics.AllocationRole, bytes: ?u64) diagnostics.Error {
        @branchHint(.cold);
        return self.failure.fail(.{ .details = .{ .@"resource.out_of_memory" = .{ .allocation_role = role, .requested_bytes = bytes } }, .stage = .allocate, .at = self.at, .cause = .{ .zig = error.OutOfMemory } });
    }

    fn mapping(self: *const Diagnostics, bytes: u64, err: ?anyerror) diagnostics.Error {
        @branchHint(.cold);
        return self.failure.fail(.{ .details = .{ .@"io.map_failed" = .{ .map_bytes = bytes } }, .stage = .map, .at = self.at, .cause = if (err) |cause| .{ .zig = cause } else null });
    }

    fn pathLimit(self: *const Diagnostics, limit: u64) diagnostics.Error {
        @branchHint(.cold);
        return self.failure.fail(.{ .details = .{ .@"io.invalid_path" = .{ .reason = .too_long, .limit_bytes = limit } }, .stage = .open, .at = self.at });
    }

    fn capacity(self: *const Diagnostics, role: diagnostics.CapacityRole, limit: u128, required: u128) diagnostics.Error {
        @branchHint(.cold);
        return self.failure.fail(.{ .details = .{ .@"resource.capacity_limit_exceeded" = .{ .capacity_role = role, .limit = limit, .required = .{ .exact = required } } }, .stage = switch (role) {
            .map_bytes => .map,
            .name_map, .records => .allocate,
            .input_bytes => .read,
            else => .parse,
        }, .at = self.at });
    }

    fn nextOrdinal(self: *const Diagnostics, value: u64, counter: diagnostics.Counter) diagnostics.Error!u64 {
        return std.math.add(u64, value, 1) catch return self.failure.fail(.{ .details = .{ .@"resource.count_overflow" = .{ .counter = counter, .max_value = std.math.maxInt(u64), .current = value, .increment = 1 } }, .stage = .parse, .at = self.at });
    }

    fn overflow(self: *const Diagnostics, field: diagnostics.IndexField, arithmetic: diagnostics.Arithmetic, lhs: u128, rhs: u128, max: u128) diagnostics.Error {
        @branchHint(.cold);
        return self.fail(.{ .@"index.arithmetic_overflow" = .{ .field = field, .arithmetic = arithmetic, .lhs = lhs, .rhs = rhs, .max_value = max } });
    }

    fn truncated(self: *const Diagnostics, section: diagnostics.Section, required: u128, available: u64) diagnostics.Error {
        return self.fail(.{ .@"index.truncated" = .{ .section = section, .required_end = required, .available_end = available } });
    }
};

const ZfiNameRegions = struct { blob_start: usize, blob_len: usize, table_start: usize, table_len: usize };

// The name blob and name table that end at `payload_end`: the table's size follows from the record
// count, and the blob ends where the table starts.
fn zfiNameRegions(payload_end: usize, name_blob_len: u64, record_count: u32, records_end: usize) union(enum) { regions: ZfiNameRegions, invalid: diagnostics.Details } {
    const table_len: usize = @intCast(NameTable.slotCount(record_count) * 4);
    if (payload_end < records_end or payload_end - records_end < table_len) return .{ .invalid = .{ .@"index.truncated" = .{ .section = .name_table, .required_end = @as(u128, records_end) + table_len, .available_end = payload_end } } };
    const table_start = payload_end - table_len;
    const out_of_bounds: diagnostics.Details = .{ .@"index.name_blob_out_of_bounds" = .{ .name_blob_bytes = name_blob_len, .payload_end = table_start, .records_end = records_end } };
    const blob_len = std.math.cast(usize, name_blob_len) orelse return .{ .invalid = out_of_bounds };
    const blob_start = std.math.sub(usize, table_start, blob_len) catch return .{ .invalid = out_of_bounds };
    if (blob_start < records_end) return .{ .invalid = out_of_bounds };
    return .{ .regions = .{ .blob_start = blob_start, .blob_len = blob_len, .table_start = table_start, .table_len = table_len } };
}

fn zfiRecordsEnd(record_count: u32, context: *Diagnostics) diagnostics.Error!usize {
    const records_bytes = std.math.mul(usize, @as(usize, record_count), ZFI_INDEX_RECORD_BYTES) catch return context.overflow(.record_count, .multiply, record_count, ZFI_INDEX_RECORD_BYTES, std.math.maxInt(usize));
    return std.math.add(usize, @sizeOf(ZfiHeader), records_bytes) catch return context.overflow(.record_count, .add, @sizeOf(ZfiHeader), records_bytes, std.math.maxInt(usize));
}

fn tryLoadZfi(
    backing_allocator: std.mem.Allocator,
    io: std.Io,
    zfi_path: []const u8,
    fasta_file: std.Io.File,
    fasta_stat: std.Io.File.Stat,
    mode: LoadMode,
    context: *Diagnostics,
) diagnostics.Error!?LoadedIndex {
    const zfi_file = platform.openSidecar(io, zfi_path) catch |err| switch (err) {
        error.FileNotFound => return null,
        else => return context.ioFailure(.open, err),
    };
    // An index read on demand keeps the file.
    var keep_file = false;
    defer if (!keep_file) zfi_file.close(io);

    const zfi_stat = zfi_file.stat(io) catch |err| return context.ioFailure(.stat, err);

    if (zfi_stat.kind == .directory) return context.mapping(zfi_stat.size, null);
    if (zfi_stat.size == 0) {
        return context.truncated(.header, ZFI_HEADER_BYTES, zfi_stat.size);
    }

    if (zfi_stat.mtime.nanoseconds < fasta_stat.mtime.nanoseconds) {
        return context.fail(.{ .@"index.older_than_source" = .{ .index_mtime_ns = zfi_stat.mtime.nanoseconds, .source_mtime_ns = fasta_stat.mtime.nanoseconds } });
    }

    if (mode == .stream) {
        try mode.stream.entries.startZfi(backing_allocator, io, zfi_path, zfi_stat.size, fasta_stat.size, fasta_stat.mtime.nanoseconds, mode.stream.check_name_table, context);
        return streamShell(backing_allocator, io, fasta_file, fasta_stat.size, .zfi);
    }

    switch (mode) {
        .lookup_full_map, .requested => {
            const loaded = try openZfiOnDemand(backing_allocator, io, zfi_file, zfi_stat.size, fasta_file, fasta_stat, context);
            keep_file = true;
            return loaded;
        },
        .catalog, .stats_scan, .stream, .on_access => {},
    }
    const zfi_size = std.math.cast(usize, zfi_stat.size) orelse return context.capacity(.map_bytes, std.math.maxInt(usize), zfi_stat.size);
    var zfi_map = platform.mapFileReadOnly(io, zfi_file, zfi_size) catch |err| return context.mapping(zfi_stat.size, err);
    errdefer zfi_map.destroy(io);

    const zfi_data: platform.MappedBytes = zfi_map.memory;

    context.at.location.offset = 0;
    if (zfi_data.len < @sizeOf(ZfiHeader)) {
        return context.truncated(.header, ZFI_HEADER_BYTES, zfi_data.len);
    }

    const header: *const ZfiHeader = @ptrCast(@alignCast(zfi_data.ptr));
    try checkZfiMagic(header.magic, context);

    context.at.location.offset = 8;
    if (header.source_size != fasta_stat.size) {
        return context.fail(.{ .@"index.source_size_mismatch" = .{ .stored_bytes = header.source_size, .actual_bytes = fasta_stat.size } });
    }

    context.at.location.offset = 4;
    if (header.record_count == 0) return context.fail(.{ .@"index.no_records" = .{ .format = .zfi } });

    const records_end = try zfiRecordsEnd(header.record_count, context);
    if (zfi_data.len < records_end) return context.truncated(.records, records_end, zfi_data.len);

    const record_bytes = zfi_data[@sizeOf(ZfiHeader)..records_end];
    const records: []const ZfiRecord = @as(
        [*]const ZfiRecord,
        @ptrCast(@alignCast(record_bytes.ptr)),
    )[0..header.record_count];

    context.at = .{ .source = .zfi, .location = .{ .offset = zfi_data.len } };
    const trailing = parseZfiTrailingMeta(zfi_data) orelse return context.fail(.{ .@"index.name_footer_missing" = {} });
    const regions = switch (zfiNameRegions(trailing.payload_end, trailing.name_blob_len, header.record_count, records_end)) {
        .invalid => |details| return context.fail(details),
        .regions => |regions| regions,
    };
    context.at.location.offset = trailing.source_mtime_offset;
    try checkSourceTime(trailing.source_mtime_ns, fasta_stat.mtime.nanoseconds, context);
    const name_blob = zfi_data[regions.blob_start..][0..regions.blob_len];
    const name_table = zfi_data[regions.table_start..][0..regions.table_len];

    const side_region_start = records_end;
    const side_region_end = regions.blob_start;

    if (mode != .on_access) try validateZfiRecords(records, fasta_stat.size, zfi_data, name_blob, side_region_start, side_region_end, context);
    if (mode == .catalog) try verifyNameTable(backing_allocator, .{ .mapped = records }, name_blob, name_table, regions.table_start, context);

    var arena = std.heap.ArenaAllocator.init(backing_allocator);
    errdefer arena.deinit();
    const checked_records: []usize = if (mode == .on_access) checked: {
        const words = arena.allocator().alloc(usize, (header.record_count + @bitSizeOf(usize) - 1) / @bitSizeOf(usize)) catch return context.allocation(.records, null);
        @memset(words, 0);
        break :checked words;
    } else &.{};

    return LoadedIndex{
        .io = io,
        .fasta_file = fasta_file,
        .zfi_map = zfi_map,
        .records = .{ .mapped = records },
        .name_blob = name_blob,
        .fasta_size = fasta_stat.size,
        .zfi_data = zfi_data,
        .zfi_side_start = side_region_start,
        .zfi_side_end = side_region_end,
        .name_table = name_table,
        .name_table_offset = regions.table_start,
        .checked_records = checked_records,
        .source = .zfi,
        .arena = arena,
    };
}

// The file checks of `tryLoadZfi`, in its order, through positional reads; records are left for
// `ZfiOnDemand` to check as requests reach them.
fn openZfiOnDemand(backing_allocator: std.mem.Allocator, io: std.Io, file: std.Io.File, zfi_size: u64, fasta_file: std.Io.File, fasta_stat: std.Io.File.Stat, context: *Diagnostics) diagnostics.Error!LoadedIndex {
    const size = std.math.cast(usize, zfi_size) orelse return context.capacity(.map_bytes, std.math.maxInt(usize), zfi_size);
    context.at.location.offset = 0;
    if (size < ZFI_HEADER_BYTES) return context.truncated(.header, ZFI_HEADER_BYTES, size);
    var header_bytes: [ZFI_HEADER_BYTES]u8 = undefined;
    const got = std.Io.File.readPositionalAll(file, io, &header_bytes, 0) catch |err| return context.read(0, ZFI_HEADER_BYTES, null, err);
    if (got != ZFI_HEADER_BYTES) return context.read(0, ZFI_HEADER_BYTES, got, null);
    const header = std.mem.bytesToValue(ZfiHeader, &header_bytes);
    const records_end = try checkZfiHeader(header, fasta_stat.size, size, context);
    var tail: [ZFI_TRAILER_BYTES]u8 = undefined;
    const tail_start = size - @min(tail.len, size);
    const tail_bytes = tail[0 .. size - tail_start];
    const tail_got = std.Io.File.readPositionalAll(file, io, tail_bytes, tail_start) catch |err| return context.read(tail_start, tail_bytes.len, null, err);
    if (tail_got != tail_bytes.len) return context.read(tail_start, tail_bytes.len, tail_got, null);
    context.at = .{ .source = .zfi, .location = .{ .offset = size } };
    const blob = switch (locateZfiNameBlob(tail_bytes, tail_start, size, header.record_count, records_end, context)) {
        .invalid => |details| return context.fail(details),
        .blob => |blob| blob,
    };
    context.at.location.offset = blob.trailing.source_mtime_offset;
    try checkSourceTime(blob.trailing.source_mtime_ns, fasta_stat.mtime.nanoseconds, context);
    const zfi = backing_allocator.create(ZfiOnDemand) catch return context.allocation(.records, @sizeOf(ZfiOnDemand));
    zfi.* = .{
        .io = io,
        .file = file,
        .backing_allocator = backing_allocator,
        .arena = .init(backing_allocator),
        .prefetch_arena = .init(backing_allocator),
        .record_count = header.record_count,
        .fasta_size = fasta_stat.size,
        .side_start = records_end,
        .side_end = blob.regions.blob_start,
        .blob_start = blob.regions.blob_start,
        .blob_len = blob.regions.blob_len,
        .table_start = blob.regions.table_start,
        .table_slots = NameTable.slotCount(header.record_count),
    };
    return .{
        .io = io,
        .fasta_file = fasta_file,
        .records = .{ .on_demand = header.record_count },
        .fasta_size = fasta_stat.size,
        .zfi_data = null,
        .on_demand = zfi,
        .source = .zfi,
        .arena = .init(backing_allocator),
    };
}

const ZfiReadWindow = struct {
    file: std.Io.File,
    io: std.Io,
    buffer: []u8 = &.{},
    start: u64 = 0,
    len: usize = 0,

    fn take(self: *ZfiReadWindow, offset: u64, len: usize, context: *Diagnostics) diagnostics.Error![]const u8 {
        if (len == 0) return self.buffer[0..0];
        if (offset >= self.start) {
            const relative = offset - self.start;
            if (relative <= self.len and len <= self.len - relative) {
                return self.buffer[@intCast(relative)..][0..len];
            }
        }
        self.len = std.Io.File.readPositionalAll(self.file, self.io, self.buffer, offset) catch |err| return context.read(offset, self.buffer.len, null, err);
        self.start = offset;
        if (self.len < len) return context.read(offset, len, self.len, null);
        return self.buffer[0..len];
    }
};

const UNMATCHED_RECORD_INDEX = std.math.maxInt(u32);

fn requestedNameMap(allocator: std.mem.Allocator, names: []const []const u8, context: *Diagnostics) diagnostics.Error!std.StringHashMapUnmanaged(u32) {
    var map: std.StringHashMapUnmanaged(u32) = .empty;
    const count = std.math.cast(u32, names.len) orelse return context.capacity(.name_map, std.math.maxInt(u32), names.len);
    map.ensureTotalCapacity(allocator, count) catch return context.allocation(.name_map, null);
    for (names) |name| {
        const entry = map.getOrPutAssumeCapacity(name);
        if (!entry.found_existing) entry.value_ptr.* = UNMATCHED_RECORD_INDEX;
    }
    return map;
}

fn discardUnmatchedNames(map: *std.StringHashMapUnmanaged(u32), storage: [][]const u8, matched_count: usize) void {
    var end = matched_count;
    var iterator = map.iterator();
    while (iterator.next()) |entry| {
        if (entry.value_ptr.* != UNMATCHED_RECORD_INDEX) continue;
        storage[end] = entry.key_ptr.*;
        end += 1;
    }
    for (storage[matched_count..end]) |name| _ = map.remove(name);
}

const SideTableValidation = struct {
    expected_base_start: u64 = 0,
    previous_byte_end: u64 = 0,

    fn accept(self: *SideTableValidation, line: SideTableLine, rec: IndexRecord, fasta_len: u64, entry: u64, context: *Diagnostics) diagnostics.Error!void {
        if (line.base_start != self.expected_base_start) return context.fail(.{ .@"index.side_base_start_mismatch" = .{ .entry = entry, .base_start = line.base_start, .expected_base_start = self.expected_base_start } });
        const geometry = diagnostics.Geometry{ .line_bases = line.line_bases, .line_bytes = line.line_bytes };
        if (geometryConstraint(geometry)) |constraint| return context.fail(.{ .@"index.side_line_geometry_invalid" = .{ .entry = entry, .geometry = geometry, .constraint = constraint } });
        const byte_end = std.math.add(u64, line.byte_offset, line.line_bytes) catch return context.overflow(.side_byte_end, .add, line.byte_offset, line.line_bytes, std.math.maxInt(u64));
        if (byte_end > fasta_len) return context.fail(.{ .@"index.sequence_span_out_of_bounds" = .{ .required_end = byte_end, .source_bytes = fasta_len, .entry = entry } });
        if (self.expected_base_start == 0) {
            if (line.byte_offset != rec.seq_offset) return context.fail(.{ .@"index.side_first_offset_mismatch" = .{ .first_byte_offset = line.byte_offset, .sequence_offset = rec.seq_offset } });
        } else if (line.byte_offset < self.previous_byte_end) return context.fail(.{ .@"index.side_line_order_invalid" = .{ .entry = entry, .byte_offset = line.byte_offset, .minimum_offset = self.previous_byte_end } });
        self.expected_base_start = std.math.add(u64, self.expected_base_start, line.line_bases) catch return context.overflow(.base_start, .add, self.expected_base_start, line.line_bases, std.math.maxInt(u64));
        self.previous_byte_end = byte_end;
    }
};

fn geometryConstraint(geometry: diagnostics.Geometry) ?diagnostics.GeometryConstraint {
    if (geometry.line_bases == 0) return .bases_zero;
    if (geometry.line_bytes == 0) return .bytes_zero;
    if (geometry.line_bytes < geometry.line_bases) return .bytes_below_bases;
    return null;
}

fn sequenceStartDefect(rec: IndexRecord, fasta_len: u64) ?diagnostics.Details {
    if (rec.seq_len == 0) return .{ .@"index.zero_sequence_length" = .{} };
    if (rec.seq_offset >= fasta_len) return .{ .@"index.sequence_offset_out_of_bounds" = .{ .sequence_offset = rec.seq_offset, .source_bytes = fasta_len } };
    return null;
}

fn uniformGeometryDefect(rec: IndexRecord, fasta_len: u64) ?diagnostics.Details {
    if (sequenceStartDefect(rec, fasta_len)) |defect| return defect;
    const geometry = diagnostics.Geometry{ .line_bases = rec.line_bases, .line_bytes = rec.line_bytes };
    if (geometryConstraint(geometry)) |constraint| return .{ .@"index.invalid_line_geometry" = .{ .geometry = geometry, .constraint = constraint } };
    const last_index: u128 = rec.seq_len - 1;
    const required_end = rec.seq_offset + (last_index / rec.line_bases) * rec.line_bytes + (last_index % rec.line_bases) + 1;
    if (required_end > fasta_len) return .{ .@"index.sequence_span_out_of_bounds" = .{ .required_end = required_end, .source_bytes = fasta_len } };
    return null;
}

fn nameExtentDefect(rec: IndexRecord, blob_len: usize) ?diagnostics.Details {
    if (!rangeFitsUsize(rec.name_offset, rec.name_len, blob_len)) return .{ .@"index.name_out_of_bounds" = .{ .name_offset = rec.name_offset, .name_bytes = rec.name_len, .blob_bytes = blob_len } };
    return null;
}

fn checkSequenceStart(rec: IndexRecord, fasta_len: u64, context: *Diagnostics) diagnostics.Error!void {
    if (sequenceStartDefect(rec, fasta_len)) |defect| return context.fail(defect);
}

fn checkUniformGeometry(rec: IndexRecord, fasta_len: u64, context: *Diagnostics) diagnostics.Error!void {
    if (uniformGeometryDefect(rec, fasta_len)) |defect| return context.fail(defect);
}

fn checkNameExtent(rec: IndexRecord, blob_len: usize, context: *Diagnostics) diagnostics.Error!void {
    if (nameExtentDefect(rec, blob_len)) |defect| return context.fail(defect);
}

fn sideTableExtent(rec: IndexRecord, side_start: usize, side_end: usize, count: ?u64, context: *Diagnostics) diagnostics.Error!usize {
    const wire_offset = rec.sideTableOffset();
    const offset = std.math.cast(usize, wire_offset) orelse return context.fail(.{ .@"index.side_table_offset_out_of_bounds" = .{ .offset = wire_offset, .region_start = side_start, .region_end = side_end } });
    if (offset < side_start or offset >= side_end) return context.fail(.{ .@"index.side_table_offset_out_of_bounds" = .{ .offset = offset, .region_start = side_start, .region_end = side_end } });
    if (offset % @alignOf(u64) != 0) return context.fail(.{ .@"index.side_table_misaligned" = .{ .offset = offset, .alignment_bytes = @alignOf(u64), .section = .side_table_count } });
    const count_end = std.math.add(usize, offset, 8) catch return context.overflow(.side_offset, .add, offset, 8, std.math.maxInt(usize));
    if (count_end > side_end) return context.truncated(.side_table_count, count_end, side_end);
    const line_count = count orelse return count_end;
    if (line_count == 0) return context.fail(.{ .@"index.side_table_empty" = .{} });
    const n = std.math.cast(usize, line_count) orelse return context.overflow(.side_count, .multiply, line_count, @sizeOf(SideTableLine), std.math.maxInt(usize));
    const table_bytes = std.math.mul(usize, n, @sizeOf(SideTableLine)) catch return context.overflow(.side_count, .multiply, line_count, @sizeOf(SideTableLine), std.math.maxInt(usize));
    const end = std.math.add(usize, count_end, table_bytes) catch return context.overflow(.side_count, .add, count_end, table_bytes, std.math.maxInt(usize));
    if (end > side_end) return context.truncated(.side_table_entries, end, side_end);
    return end;
}

// Checks a side table read through `source`; returns its end, with the location at its start.
fn checkSideTable(source: anytype, rec: IndexRecord, side_start: usize, side_end: usize, fasta_len: u64, context: *Diagnostics) diagnostics.Error!usize {
    const count_end = try sideTableExtent(rec, side_start, side_end, null, context);
    const offset: usize = @intCast(rec.sideTableOffset());
    context.at.location.offset = offset;
    const count = std.mem.readInt(u64, (try source.take(offset, 8, context))[0..8], .little);
    const table_end = try sideTableExtent(rec, side_start, side_end, count, context);
    var validation = SideTableValidation{};
    for (0..@intCast(count)) |i| {
        context.at.location.offset = count_end + i * @sizeOf(SideTableLine);
        const line = std.mem.bytesToValue(SideTableLine, try source.take(context.at.location.offset.?, @sizeOf(SideTableLine), context));
        try validation.accept(line, rec, fasta_len, i + 1, context);
    }
    if (validation.expected_base_start != rec.seq_len) return context.fail(.{ .@"index.side_sequence_length_mismatch" = .{ .indexed_length = rec.seq_len, .side_table_bases = validation.expected_base_start } });
    context.at.location.offset = offset;
    return table_end;
}

// Mapped bytes for `checkSideTable`, whose `sideTableExtent` checks bound every range first.
const MappedSideTables = struct {
    data: []const u8,

    fn take(self: MappedSideTables, offset: u64, len: usize, _: *Diagnostics) diagnostics.Error![]const u8 {
        return self.data[@intCast(offset)..][0..len];
    }
};

// A streamed `.zfi` header: magic, source size, record count, and room for the records.
fn checkZfiHeader(header: ZfiHeader, fasta_size: u64, zfi_size: usize, context: *Diagnostics) diagnostics.Error!usize {
    try checkZfiMagic(header.magic, context);
    context.at.location.offset = 8;
    if (header.source_size != fasta_size) return context.fail(.{ .@"index.source_size_mismatch" = .{ .stored_bytes = header.source_size, .actual_bytes = fasta_size } });
    context.at.location.offset = 4;
    if (header.record_count == 0) return context.fail(.{ .@"index.no_records" = .{ .format = .zfi } });
    const records_end = try zfiRecordsEnd(header.record_count, context);
    if (records_end > zfi_size) return context.truncated(.records, records_end, zfi_size);
    return records_end;
}

fn checkZfiMagic(magic: [4]u8, context: *Diagnostics) diagnostics.Error!void {
    context.at.location.offset = 0;
    if (!std.mem.eql(u8, magic[0..3], ZFI_MAGIC[0..3])) return context.fail(.{ .@"index.invalid_magic" = .{ .expected_magic = ZFI_MAGIC, .actual = context.failure.retain(&magic, 4) } });
    if (magic[3] != ZFI_MAGIC[3]) {
        context.at.location.offset = 3;
        return context.fail(.{ .@"index.unsupported_version" = .{ .found_version = magic[3], .supported_version = ZFI_MAGIC[3] } });
    }
}

fn checkSourceTime(stored: u64, actual: i128, context: *Diagnostics) diagnostics.Error!void {
    if (actual < 0 or actual > std.math.maxInt(u64) or stored != @as(u64, @intCast(actual))) return context.fail(.{ .@"index.source_mtime_mismatch" = .{ .stored_mtime_ns = stored, .actual_mtime_ns = actual } });
}

const ZfiNameBlob = struct { regions: ZfiNameRegions, trailing: ZfiTrailingMeta };

// The name blob and table of a `.zfi` from the file's last bytes (`tail_bytes`, which start at
// `tail_start`). A trailer that cannot be followed comes back as `invalid`.
fn locateZfiNameBlob(tail_bytes: []const u8, tail_start: usize, zfi_size: usize, record_count: u32, records_end: usize, context: *Diagnostics) union(enum) { blob: ZfiNameBlob, invalid: diagnostics.Details } {
    context.at.location.offset = zfi_size;
    var trailing = parseZfiTrailingMeta(tail_bytes) orelse return .{ .invalid = .{ .@"index.name_footer_missing" = {} } };
    trailing.payload_end += tail_start;
    trailing.source_mtime_offset += tail_start;
    return switch (zfiNameRegions(trailing.payload_end, trailing.name_blob_len, record_count, records_end)) {
        .invalid => |details| .{ .invalid = details },
        .regions => |regions| .{ .blob = .{ .regions = regions, .trailing = trailing } },
    };
}

fn tryLoadFai(
    backing_allocator: std.mem.Allocator,
    io: std.Io,
    fai_path: []const u8,
    fasta_file: std.Io.File,
    fasta_stat: std.Io.File.Stat,
    mode: LoadMode,
    context: *Diagnostics,
) diagnostics.Error!?LoadedIndex {
    const fai_file = platform.openSidecar(io, fai_path) catch |err| switch (err) {
        error.FileNotFound => return null,
        else => return context.ioFailure(.open, err),
    };
    defer fai_file.close(io);

    const fai_stat = fai_file.stat(io) catch |err| return context.ioFailure(.stat, err);
    if (fai_stat.kind == .directory) return context.mapping(fai_stat.size, null);

    if (fai_stat.mtime.nanoseconds < fasta_stat.mtime.nanoseconds) {
        return context.fail(.{ .@"index.older_than_source" = .{ .index_mtime_ns = fai_stat.mtime.nanoseconds, .source_mtime_ns = fasta_stat.mtime.nanoseconds } });
    }

    if (fai_stat.size == 0) {
        return context.fail(.{ .@"index.no_records" = .{ .format = .fai } });
    }
    if (mode == .stream) {
        try mode.stream.entries.start(io, fai_path, .fai, fasta_stat.size, context);
        return streamShell(backing_allocator, io, fasta_file, fasta_stat.size, .fai);
    }

    const streamed_mode: ?StreamedFaiMode = switch (mode) {
        .stats_scan => .stats_scan,
        .requested => |names| .{ .requested = names },
        .lookup_full_map, .catalog, .on_access => null,
        .stream => unreachable,
    };
    if (streamed_mode) |selected| return @as(?LoadedIndex, try loadFaiStreamed(backing_allocator, io, fai_file, fasta_file, fasta_stat, selected, context));

    const fai_size = std.math.cast(usize, fai_stat.size) orelse return context.capacity(.map_bytes, std.math.maxInt(usize), fai_stat.size);
    var fai_map = platform.mapFileReadOnly(io, fai_file, fai_size) catch |err| return context.mapping(fai_stat.size, err);

    errdefer fai_map.destroy(io);

    var arena = std.heap.ArenaAllocator.init(backing_allocator);
    errdefer arena.deinit();
    const allocator = arena.allocator();

    const fai_data: platform.MappedBytes = fai_map.memory;

    if (mode != .catalog) {
        const lazy = try loadFaiLazy(allocator, fai_data, context);
        return LoadedIndex{
            .io = io,
            .fasta_file = fasta_file,
            .fai_map = fai_map,
            .records = .{ .on_demand = lazy.lines.len },
            .name_map = lazy.map,
            .fai_data = fai_data,
            .fai_lines = lazy.lines,
            .checked_records = lazy.checked,
            .fasta_size = fasta_stat.size,
            .zfi_data = null,
            .source = .fai,
            .arena = arena,
        };
    }

    var approx_records = std.mem.count(u8, fai_data, "\n");
    if (fai_data.len > 0 and fai_data[fai_data.len - 1] != '\n') {
        approx_records = std.math.add(usize, approx_records, 1) catch return context.overflow(.record_count, .add, approx_records, 1, std.math.maxInt(usize));
    }
    _ = std.math.mul(usize, approx_records, @sizeOf(IndexRecord)) catch return context.overflow(.record_count, .multiply, approx_records, @sizeOf(IndexRecord), std.math.maxInt(usize));
    if (approx_records == 0) return context.fail(.{ .@"index.no_records" = .{ .format = .fai } });
    const records = allocator.alloc(IndexRecord, approx_records) catch return context.allocation(.records, @as(u64, approx_records) * @sizeOf(IndexRecord));

    var record_count: usize = 0;
    var pos: usize = 0;
    var row_number: u64 = 0;
    while (pos < fai_data.len) {
        const line_start = pos;
        const rel_eol = std.mem.findScalar(u8, fai_data[pos..], '\n') orelse fai_data.len - pos;
        const line_len = rel_eol;
        const line_end = std.math.add(usize, line_start, line_len) catch return context.overflow(.offset, .add, line_start, line_len, std.math.maxInt(usize));
        if (line_end < fai_data.len) {
            pos = std.math.add(usize, line_end, 1) catch return context.overflow(.offset, .add, line_end, 1, std.math.maxInt(usize));
        } else {
            pos = line_end;
        }
        row_number = try context.nextOrdinal(row_number, .physical_lines);
        context.at = .{ .source = .fai, .record = record_count + 1, .location = .{ .offset = line_start, .line = row_number, .column = 1 } };
        if (line_len == 0) continue;

        const line = fai_data[line_start..][0..line_len];
        const fields = try parseFaiIndexLine(line, context);

        const rec = IndexRecord{
            .name_offset = line_start,
            .name_len = fields.name_len,
            .seq_offset = fields.seq_offset,
            .seq_len = fields.seq_len,
            .line_bases = fields.line_bases,
            .line_bytes = fields.line_bytes,
        };
        try checkUniformGeometry(rec, fasta_stat.size, context);

        if (record_count >= records.len) return context.fail(.{ .@"internal.consistency_failed" = .{ .component = .index_loader, .check = .record_reservation } });
        records[record_count] = rec;
        record_count += 1;
    }

    if (record_count == 0) {
        return context.fail(.{ .@"index.no_records" = .{ .format = .fai } });
    }

    const loaded_records = records[0..record_count];

    return LoadedIndex{
        .io = io,
        .fasta_file = fasta_file,
        .fai_map = fai_map,
        .records = .{ .owned = loaded_records },
        .name_slices = &.{},
        .fai_data = fai_data,
        .fasta_size = fasta_stat.size,
        .zfi_data = null,
        .source = .fai,
        .arena = arena,
    };
}

const FaiLine = struct { offset: u64, row: u64 };

fn faiLineName(line: []const u8) []const u8 {
    for (line, 0..) |byte, i| if (byte == '\t' or byte == '\n') return line[0..i];
    return line;
}

// The lines of a mapped `.fai`, for any-name lookups: each name goes into the name map (first line
// wins) and each line's start is kept; `checkedRecord` parses and checks a line when it is used.
fn loadFaiLazy(allocator: std.mem.Allocator, data: []const u8, context: *Diagnostics) diagnostics.Error!struct { lines: []const FaiLine, map: std.StringHashMapUnmanaged(u32), checked: []usize } {
    var lines: std.ArrayList(FaiLine) = .empty;
    var map: std.StringHashMapUnmanaged(u32) = .empty;
    const approx = std.mem.count(u8, data, "\n") + 1;
    lines.ensureTotalCapacityPrecise(allocator, approx) catch return context.allocation(.records, @as(u64, approx) * @sizeOf(FaiLine));
    map.ensureTotalCapacity(allocator, std.math.cast(u32, approx) orelse return context.capacity(.name_map, std.math.maxInt(u32), approx)) catch return context.allocation(.name_map, null);
    var pos: usize = 0;
    var row: u64 = 0;
    while (pos < data.len) {
        row += 1;
        const breaks = @call(.always_inline, faiLineBreaks, .{data[pos..]});
        const end = if (breaks) |found| pos + found.end else data.len;
        defer pos = end + 1;
        if (end == pos) continue;
        const name_end = if (breaks) |found| if (found.tab) |tab| pos + tab else end else if (std.mem.findScalar(u8, data[pos..end], '\t')) |tab| pos + tab else end;
        const ordinal = std.math.cast(u32, lines.items.len) orelse return context.capacity(.records, std.math.maxInt(u32), lines.items.len);
        const entry = map.getOrPutAssumeCapacity(data[pos..name_end]);
        if (!entry.found_existing) entry.value_ptr.* = ordinal;
        lines.appendAssumeCapacity(.{ .offset = pos, .row = row });
    }
    if (lines.items.len == 0) return context.fail(.{ .@"index.no_records" = .{ .format = .fai } });
    const checked = allocator.alloc(usize, (lines.items.len + @bitSizeOf(usize) - 1) / @bitSizeOf(usize)) catch return context.allocation(.records, null);
    @memset(checked, 0);
    return .{ .lines = lines.items, .map = map, .checked = checked };
}

const StreamedFaiMode = union(enum) {
    stats_scan,
    requested: []const []const u8,
};

fn loadFaiStreamed(
    backing_allocator: std.mem.Allocator,
    io: std.Io,
    fai_file: std.Io.File,
    fasta_file: std.Io.File,
    fasta_stat: std.Io.File.Stat,
    mode: StreamedFaiMode,
    context: *Diagnostics,
) diagnostics.Error!LoadedIndex {
    var arena = std.heap.ArenaAllocator.init(backing_allocator);
    errdefer arena.deinit();
    const allocator = arena.allocator();

    var records_list: std.ArrayList(IndexRecord) = .empty;
    const records_allocator = if (mode == .stats_scan) backing_allocator else allocator;
    defer if (mode == .stats_scan) records_list.deinit(backing_allocator);
    var slices_list: std.ArrayList([]const u8) = .empty;
    const single_name = switch (mode) {
        .requested => |names| if (names.len == 1) names[0] else null,
        .stats_scan => null,
    };
    var requested: std.StringHashMapUnmanaged(u32) = .empty;
    var maybe_requested: ?*NamePrefilter = null;
    switch (mode) {
        .requested => |names| {
            if (single_name == null) {
                requested = try requestedNameMap(allocator, names, context);
                maybe_requested = allocator.create(NamePrefilter) catch return context.allocation(.name_map, @sizeOf(NamePrefilter));
                maybe_requested.?.* = .{};
                for (names) |name| maybe_requested.?.add(name);
            }
            const count = if (single_name != null) 1 else requested.count();
            records_list.ensureTotalCapacityPrecise(allocator, count) catch return context.allocation(.records, @as(u64, count) * @sizeOf(IndexRecord));
            slices_list.ensureTotalCapacityPrecise(allocator, count) catch return context.allocation(.name_slices, @as(u64, count) * @sizeOf([]const u8));
        },
        .stats_scan => {},
    }

    const buffers = backing_allocator.alloc(u8, 2 * FAI_READER_BUFFER_BYTES) catch return context.allocation(.reader_buffer, 2 * FAI_READER_BUFFER_BYTES);
    defer backing_allocator.free(buffers);
    var file_reader = fai_file.reader(io, buffers[0..FAI_READER_BUFFER_BYTES]);
    const long_line_buffer = buffers[FAI_READER_BUFFER_BYTES..][0..FAI_READER_BUFFER_BYTES];
    var saw_record = false;

    var row_number: u64 = 0;
    var record_number: u64 = 0;
    while (true) {
        const current_offset = file_reader.logicalPos();
        context.at.location.offset = current_offset;
        context.at.location.line = try context.nextOrdinal(row_number, .physical_lines);
        context.at.location.column = 1;
        context.at.record = null;
        var long_fields: ?FaiLineFields = null;
        const available = file_reader.interface.buffered();
        var name_end: ?usize = null;
        // Inlined: as a call per line, the scan of a transcriptome `.fai` took 7 M more instructions.
        const maybe_line = if (@call(.always_inline, faiLineBreaks, .{available})) |breaks| buffered: {
            file_reader.interface.toss(breaks.end + 1);
            name_end = breaks.tab orelse breaks.end;
            break :buffered available[0..breaks.end];
        } else file_reader.interface.takeDelimiter('\n') catch |err| switch (err) {
            error.ReadFailed => return context.ioFailure(.read, file_reader.err orelse err),
            error.StreamTooLong => long: {
                context.at.record = try context.nextOrdinal(record_number, .records);
                const parsed = try takeLongFaiLine(&file_reader, long_line_buffer, context);
                long_fields = parsed.fields;
                break :long parsed.name;
            },
        };
        const line = maybe_line orelse break;
        row_number += 1;
        if (line.len == 0 and long_fields == null) continue;
        record_number = try context.nextOrdinal(record_number, .records);
        context.at.record = record_number;
        saw_record = true;

        // GET parses and checks only the lines it retains, and stops once every name is found.
        var matched: ?std.StringHashMapUnmanaged(u32).Entry = null;
        if (mode != .stats_scan) {
            const name = if (long_fields != null) line else line[0 .. name_end orelse std.mem.findScalar(u8, line, '\t') orelse line.len];
            if (single_name == null and maybe_requested.?.mayHold(name)) matched = requested.getEntry(name);
            const wanted = if (single_name) |only| std.mem.eql(u8, name, only) else if (matched) |entry| entry.value_ptr.* == UNMATCHED_RECORD_INDEX else false;
            if (!wanted) continue;
        }

        const fields = long_fields orelse try parseFaiIndexLine(line, context);

        const rec = IndexRecord{
            .name_offset = 0,
            .name_len = fields.name_len,
            .seq_offset = fields.seq_offset,
            .seq_len = fields.seq_len,
            .line_bases = fields.line_bases,
            .line_bytes = fields.line_bytes,
        };
        try checkUniformGeometry(rec, fasta_stat.size, context);

        const rec_index = records_list.items.len;
        switch (mode) {
            .requested => {
                // The requested slice equals the line's name and outlives the index.
                const name = if (matched) |entry| entry.key_ptr.* else single_name.?;
                slices_list.append(allocator, name) catch return context.allocation(.name_slices, null);
                if (matched) |entry| entry.value_ptr.* = @intCast(rec_index);
            },
            .stats_scan => {},
        }
        records_list.append(records_allocator, rec) catch return context.allocation(.records, null);
        if (mode != .stats_scan and records_list.items.len == if (single_name != null) 1 else requested.count()) break;
    }

    if (!saw_record) return context.fail(.{ .@"index.no_records" = .{ .format = .fai } });

    if (mode != .stats_scan and single_name == null) discardUnmatchedNames(&requested, slices_list.allocatedSlice(), records_list.items.len);

    const owned_records = if (mode == .stats_scan)
        allocator.dupe(IndexRecord, records_list.items) catch return context.allocation(.records, null)
    else
        records_list.toOwnedSlice(allocator) catch return context.allocation(.records, null);
    const owned_slices: []const []const u8 = switch (mode) {
        .requested => slices_list.toOwnedSlice(allocator) catch return context.allocation(.name_slices, null),
        .stats_scan => &.{},
    };

    const name_map: ?std.StringHashMapUnmanaged(u32) = switch (mode) {
        .requested => if (single_name != null) null else requested,
        .stats_scan => null,
    };

    return LoadedIndex{
        .io = io,
        .fasta_file = fasta_file,
        .records = .{ .owned = owned_records },
        .name_map = name_map,
        .name_slices = owned_slices,
        .fai_data = null,
        .fasta_size = fasta_stat.size,
        .zfi_data = null,
        .source = .fai,
        .arena = arena,
    };
}

const FaiLineFields = struct {
    name_end: usize,
    name_len: u16,
    seq_len: u64,
    seq_offset: u64,
    line_bases: u32,
    line_bytes: u32,
};

const FAI_NUMERIC_FIELDS = [_]diagnostics.FaiField{ .length, .offset, .line_bases, .line_bytes };

fn faiFieldContext(context: *Diagnostics, start: u64) Diagnostics {
    var field_context = context.*;
    field_context.at.location.offset = context.at.location.offset.? + start;
    field_context.at.location.column = start + 1;
    return field_context;
}

fn missingFaiField(context: *Diagnostics, field: diagnostics.FaiField) diagnostics.Error {
    @branchHint(.cold);
    return context.fail(.{ .@"index.fai.field_missing" = .{ .field = field, .column_number = @intCast(@intFromEnum(field) + 1) } });
}

fn checkFaiInteger(value: u64, lexical: bool, overflow: bool, bytes: u64, prefix: []const u8, field: diagnostics.FaiField, max: u64, context: *Diagnostics) diagnostics.Error!u64 {
    if (bytes == 0) return missingFaiField(context, field);
    if (lexical) return context.fail(.{ .@"index.fai.invalid_integer" = .{ .field = field, .value = context.failure.retainPrefix(prefix, .{ .exact = bytes }, 32) } });
    if (overflow or value > max) return context.fail(.{ .@"index.fai.integer_out_of_range" = .{ .field = field, .max_value = max, .value = context.failure.retainPrefix(prefix, .{ .exact = bytes }, 32) } });
    return value;
}

fn parseFaiAsciiU64(line: []const u8, start: usize) ?struct { value: u64, end: usize } {
    var end = start;
    var value: u64 = 0;
    while (end < line.len and line[end] != '\t') : (end += 1) {
        const byte = line[end];
        if (byte < '0' or byte > '9') return null;
        // Nineteen digits always fit a u64.
        value = if (end - start < 19) value * 10 + (byte - '0') else appendFaiDigit(value, byte) orelse return null;
    }
    if (end == start) return null;
    return .{ .value = value, .end = end };
}

fn faiIntegerFailure(text: []const u8, field: diagnostics.FaiField, max: u64, context: *Diagnostics) diagnostics.Error {
    @branchHint(.cold);
    if (text.len == 0) return missingFaiField(context, field);
    const preview = context.failure.retain(text, 32);
    for (text) |byte| {
        if (byte < '0' or byte > '9') return context.fail(.{ .@"index.fai.invalid_integer" = .{ .field = field, .value = preview } });
    }
    return context.fail(.{ .@"index.fai.integer_out_of_range" = .{ .field = field, .max_value = max, .value = preview } });
}

// The end of the first line of `bytes` and the tab before it, if any, in one pass.
fn faiLineBreaks(bytes: []const u8) ?struct { tab: ?usize, end: usize } {
    const Block = @Vector(32, u8);
    var tab: ?usize = null;
    var i: usize = 0;
    while (i + 32 <= bytes.len) : (i += 32) {
        const block: Block = bytes[i..][0..32].*;
        const newlines: u32 = @bitCast(block == @as(Block, @splat('\n')));
        if (tab == null) {
            const tabs: u32 = @bitCast(block == @as(Block, @splat('\t')));
            if (tabs != 0 and (newlines == 0 or @ctz(tabs) < @ctz(newlines))) tab = i + @ctz(tabs);
        }
        if (newlines != 0) return .{ .tab = tab, .end = i + @ctz(newlines) };
    }
    while (i < bytes.len) : (i += 1) {
        if (bytes[i] == '\n') return .{ .tab = tab, .end = i };
        if (tab == null and bytes[i] == '\t') tab = i;
    }
    return null;
}

/// A 65,536-bit filter over requested names, keyed by length and the first 16 and last 8 bytes;
/// a `.fai` line whose name it rules out skips the hash lookup.
const NamePrefilter = struct {
    bits: [1024]u64 = .{0} ** 1024,

    fn key(name: []const u8) u16 {
        var mixed: u64 = name.len *% 0x9E3779B97F4A7C15;
        if (name.len >= 16) {
            mixed ^= std.mem.readInt(u64, name[0..8], .little) *% 0xC2B2AE3D27D4EB4F;
            mixed ^= std.mem.readInt(u64, name[8..16], .little) *% 0x165667B19E3779F9;
            mixed ^= std.mem.readInt(u64, name[name.len - 8 ..][0..8], .little) *% 0x27D4EB2F165667C5;
        } else {
            var padded = [_]u8{0} ** 16;
            @memcpy(padded[0..name.len], name);
            mixed ^= std.mem.readInt(u64, padded[0..8], .little) *% 0xC2B2AE3D27D4EB4F;
            mixed ^= std.mem.readInt(u64, padded[8..16], .little) *% 0x165667B19E3779F9;
        }
        mixed ^= mixed >> 29;
        return @truncate((mixed *% 0xBF58476D1CE4E5B9) >> 48);
    }

    fn add(self: *NamePrefilter, name: []const u8) void {
        const k = key(name);
        self.bits[k >> 6] |= @as(u64, 1) << @intCast(k & 63);
    }

    fn mayHold(self: *const NamePrefilter, name: []const u8) bool {
        const k = key(name);
        return self.bits[k >> 6] & @as(u64, 1) << @intCast(k & 63) != 0;
    }
};

fn parseFaiIndexLine(line: []const u8, context: *Diagnostics) diagnostics.Error!FaiLineFields {
    const name_end = std.mem.findScalar(u8, line, '\t') orelse {
        var field_context = faiFieldContext(context, line.len);
        return missingFaiField(&field_context, .length);
    };
    const name_len = std.math.cast(u16, name_end) orelse return context.fail(.{ .@"index.identifier_too_long" = .{ .observed_bytes = .{ .exact = name_end }, .limit_bytes = std.math.maxInt(u16), .format = .fai } });
    var start = name_end;
    var values: [4]u64 = undefined;
    inline for (FAI_NUMERIC_FIELDS, 0..) |field, i| {
        if (start < line.len and line[start] == '\t') start += 1;
        const max = if (i < 2) std.math.maxInt(u64) else std.math.maxInt(u32);
        const parsed = parseFaiAsciiU64(line, start) orelse {
            const end = std.mem.findScalarPos(u8, line, start, '\t') orelse line.len;
            var field_context = faiFieldContext(context, start);
            return faiIntegerFailure(line[start..end], field, max, &field_context);
        };
        values[i] = parsed.value;
        if (values[i] > max) {
            var field_context = faiFieldContext(context, start);
            return faiIntegerFailure(line[start..parsed.end], field, max, &field_context);
        }
        start = parsed.end;
    }
    return .{ .name_end = name_end, .name_len = name_len, .seq_len = values[0], .seq_offset = values[1], .line_bases = @intCast(values[2]), .line_bytes = @intCast(values[3]) };
}

fn takeLongFaiLine(file_reader: *std.Io.File.Reader, buffer: []u8, context: *Diagnostics) diagnostics.Error!struct { name: []const u8, fields: FaiLineFields } {
    const reader = &file_reader.interface;
    var name_bytes: u64 = 0;
    while (true) {
        const byte = reader.takeByte() catch |err| switch (err) {
            error.ReadFailed => return context.ioFailure(.read, file_reader.err orelse err),
            error.EndOfStream => {
                var field_context = faiFieldContext(context, name_bytes);
                return missingFaiField(&field_context, .length);
            },
        };
        if (byte == '\t') break;
        if (byte == '\n') {
            var field_context = faiFieldContext(context, name_bytes);
            return missingFaiField(&field_context, .length);
        }
        if (name_bytes < std.math.maxInt(u16)) buffer[@intCast(name_bytes)] = byte;
        name_bytes = std.math.add(u64, name_bytes, 1) catch return context.capacity(.input_bytes, std.math.maxInt(u64), @as(u128, name_bytes) + 1);
    }
    const name_len = std.math.cast(u16, name_bytes) orelse return context.fail(.{ .@"index.identifier_too_long" = .{ .observed_bytes = .{ .exact = name_bytes }, .limit_bytes = std.math.maxInt(u16), .format = .fai } });
    var values: [4]u64 = .{0} ** 4;
    var field: usize = 0;
    var position = name_bytes + 1;
    var field_start = position;
    var bytes: u64 = 0;
    var prefix: [32]u8 = undefined;
    var lexical = false;
    var overflow = false;
    while (true) {
        const byte = reader.takeByte() catch |err| switch (err) {
            error.ReadFailed => return context.ioFailure(.read, file_reader.err orelse err),
            error.EndOfStream => break,
        };
        if (byte == '\n') break;
        if (field == values.len) continue;
        if (byte == '\t') {
            var field_context = faiFieldContext(context, field_start);
            values[field] = try checkFaiInteger(values[field], lexical, overflow, bytes, prefix[0..@intCast(@min(bytes, prefix.len))], FAI_NUMERIC_FIELDS[field], if (field < 2) std.math.maxInt(u64) else std.math.maxInt(u32), &field_context);
            field += 1;
            bytes = 0;
            lexical = false;
            overflow = false;
            field_start = position + 1;
        } else {
            if (bytes < prefix.len) prefix[@intCast(bytes)] = byte;
            bytes += 1;
            if (byte < '0' or byte > '9') lexical = true else if (!overflow) values[field] = appendFaiDigit(values[field], byte) orelse overflowed: {
                overflow = true;
                break :overflowed values[field];
            };
        }
        position += 1;
    }
    if (field < values.len) {
        var field_context = faiFieldContext(context, field_start);
        values[field] = try checkFaiInteger(values[field], lexical, overflow, bytes, prefix[0..@intCast(@min(bytes, prefix.len))], FAI_NUMERIC_FIELDS[field], if (field < 2) std.math.maxInt(u64) else std.math.maxInt(u32), &field_context);
        if (field < 3) {
            field_context = faiFieldContext(context, position);
            return missingFaiField(&field_context, FAI_NUMERIC_FIELDS[field + 1]);
        }
    }
    return .{ .name = buffer[0..name_len], .fields = .{ .name_end = name_len, .name_len = name_len, .seq_len = values[0], .seq_offset = values[1], .line_bases = @intCast(values[2]), .line_bytes = @intCast(values[3]) } };
}

fn appendFaiDigit(value: u64, byte: u8) ?u64 {
    const shifted = std.math.mul(u64, value, 10) catch return null;
    return std.math.add(u64, shifted, byte - '0') catch null;
}

fn rangeFitsUsize(offset: u64, len: u64, limit: usize) bool {
    const off = std.math.cast(usize, offset) orelse return false;
    const len_usz = std.math.cast(usize, len) orelse return false;
    const end = std.math.add(usize, off, len_usz) catch return false;
    return end <= limit;
}

fn validateZfiRecords(records: []const ZfiRecord, fasta_len: u64, data: []const u8, blob: []const u8, side_start: usize, side_end: usize, context: *Diagnostics) diagnostics.Error!void {
    var previous_end = side_start;
    for (records, 0..) |wire, record_index| {
        const rec = wire.decode();
        if (rec.isUniformWidth() and nameExtentDefect(rec, blob.len) == null and uniformGeometryDefect(rec, fasta_len) == null) continue;
        const table_end = try checkZfiRecord(data, rec, record_index, fasta_len, blob.len, side_start, side_end, context) orelse continue;
        const offset: usize = @intCast(rec.sideTableOffset());
        if (offset < previous_end) return context.fail(.{ .@"index.side_table_order_invalid" = .{ .table_start = offset, .minimum_start = previous_end } });
        previous_end = table_end;
    }
}

fn recordReference(index: usize) diagnostics.Reference {
    return .{ .source = .zfi, .record = index + 1, .location = .{ .offset = ZFI_HEADER_BYTES + index * ZFI_INDEX_RECORD_BYTES } };
}

// Checks record `index` of a mapped `.zfi`: its name, then its geometry or its whole side table.
// Returns the side table's end, or null for a uniform record.
fn checkZfiRecord(data: []const u8, rec: IndexRecord, index: usize, fasta_len: u64, blob_len: usize, side_start: usize, side_end: usize, context: *Diagnostics) diagnostics.Error!?usize {
    context.at = recordReference(index);
    try checkNameExtent(rec, blob_len, context);
    if (rec.isUniformWidth()) {
        try checkUniformGeometry(rec, fasta_len, context);
        return null;
    }
    try checkSequenceStart(rec, fasta_len, context);
    return try checkSideTable(MappedSideTables{ .data = data }, rec, side_start, side_end, fasta_len, context);
}

// Rebuilds the name table from the records and compares it with the stored one, slot by slot.
fn verifyNameTable(allocator: std.mem.Allocator, records: Records, blob: []const u8, table: []const u8, table_offset: usize, context: *Diagnostics) diagnostics.Error!void {
    var expected = NameTable.Builder.init(allocator, records.len()) catch return context.allocation(.name_map, table.len);
    defer expected.deinit(allocator);
    for (0..records.len()) |i| expected.insert(records.get(i).getName(blob));
    try compareNameTable(expected.bytes(), table, table_offset, 0, context);
}

// Compares the table bytes from `first_byte` on; both slices start at the same slot boundary.
fn compareNameTable(expected: []const u8, stored: []const u8, table_start: usize, first_byte: usize, context: *Diagnostics) diagnostics.Error!void {
    std.debug.assert(expected.len == stored.len and first_byte % 4 == 0);
    if (std.mem.indexOfDiff(u8, expected, stored)) |byte| {
        const slot = (first_byte + byte) / 4;
        context.at = .{ .source = .zfi, .location = .{ .offset = table_start + slot * 4 } };
        return context.fail(.{ .@"index.name_table_invalid" = .{ .slot = slot } });
    }
}

test "[property] - [fai parser]: integer fields match a reference parse" {
    var prng = std.Random.DefaultPrng.init(0xfa1);
    const random = prng.random();
    var line: [40]u8 = undefined;
    for (0..20_000) |_| {
        const digits = random.intRangeAtMost(usize, 0, 22);
        const tail = random.intRangeAtMost(usize, 0, 12);
        for (line[0..digits]) |*byte| byte.* = '0' + random.uintLessThan(u8, 10);
        for (line[digits..][0..tail]) |*byte| byte.* = switch (random.uintLessThan(u8, 4)) {
            0 => '\t',
            1 => '0' + random.uintLessThan(u8, 10),
            2 => '/',
            else => ':',
        };
        const text = line[0 .. digits + tail];
        const field_end = std.mem.findScalar(u8, text, '\t') orelse text.len;
        const expected = std.fmt.parseInt(u64, text[0..field_end], 10) catch null;
        const got = parseFaiAsciiU64(text, 0);
        if (expected) |value| {
            try std.testing.expectEqual(value, got.?.value);
            try std.testing.expectEqual(field_end, got.?.end);
        } else {
            try std.testing.expect(got == null);
        }
    }
}
