//! `.zfi` binary and `.fai` text index load/write helpers (`LoadedIndex`, staleness, geometry).
//!
//! On-disk `.zfi` bytes follow the frozen little-endian contract in `plan/docs/zfi-format.md`.
//! Host `extern struct` views are allowed only when comptime checks prove native layout
//! matches that contract (little-endian + fixed sizes/offsets).

const std = @import("std");
const builtin = @import("builtin");
const platform = @import("platform.zig");

// --- Types shared by indexer, getter, and stats ---

/// Four-byte magic and version for the supported `.zfi` wire format.
pub const ZFI_MAGIC: [4]u8 = .{ 'Z', 'F', 'I', 0x01 };

const NON_UNIFORM_WIDTH_FLAG: u8 = 1;
/// Record flag indicating that its name is embedded in the `.zfi` name blob.
pub const NAME_IN_ZFI_FLAG: u8 = 2;
const SIDE_TABLE_OFFSET_BYTES = 5;
const MAX_SIDE_TABLE_OFFSET = (1 << (SIDE_TABLE_OFFSET_BYTES * 8)) - 1;
// Holds the maximum u16 name, four tabs, four maximum-width numeric fields,
// and the line delimiter.
const FAI_READER_BUFFER_BYTES: usize = 65_600;

pub const ZFI_NAME_FOOTER_MAGIC: [4]u8 = .{ 'Z', 'F', 'N', 'M' };

/// On-disk sizes from the wire contract (must match `@sizeOf` under the frozen layout).
pub const ZFI_HEADER_BYTES: usize = 16;
pub const ZFI_INDEX_RECORD_BYTES: usize = 40;
pub const ZFI_SIDE_TABLE_LINE_BYTES: usize = 32;

/// On-disk name-table footer size (4-byte magic + little-endian u64 length).
pub const ZFI_NAME_FOOTER_BYTES: usize = 12;

/// Trailing source-identity block (`ZFID` + FASTA mtime ns) written immediately
/// before the name footer. Production writers always include it. Legacy files
/// without it use the weaker mtime policy.
pub const ZFI_SOURCE_ID_MAGIC: [4]u8 = .{ 'Z', 'F', 'I', 'D' };
pub const ZFI_SOURCE_ID_BYTES: usize = 12;

/// Encodes the trailing name-blob length in the frozen little-endian layout.
pub fn encodeZfiNameFooter(name_blob_len: u64) [ZFI_NAME_FOOTER_BYTES]u8 {
    var out: [ZFI_NAME_FOOTER_BYTES]u8 = undefined;
    @memcpy(out[0..4], &ZFI_NAME_FOOTER_MAGIC);
    std.mem.writeInt(u64, out[4..12], name_blob_len, .little);
    return out;
}

/// Encodes the source modification time in the frozen little-endian layout.
pub fn encodeZfiSourceId(source_mtime_ns: u64) [ZFI_SOURCE_ID_BYTES]u8 {
    var out: [ZFI_SOURCE_ID_BYTES]u8 = undefined;
    @memcpy(out[0..4], &ZFI_SOURCE_ID_MAGIC);
    std.mem.writeInt(u64, out[4..12], source_mtime_ns, .little);
    return out;
}

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

// Legacy on-disk footer written before tight layout (16 bytes: magic + 4 pad + u64).
const ZFI_NAME_FOOTER_LEGACY_BYTES: usize = 16;

/// Converts a non-negative timestamp to the unsigned `.zfi` representation.
pub fn timestampToNs(ts: std.Io.Timestamp) error{TimestampOutOfRange}!u64 {
    return std.math.cast(u64, ts.nanoseconds) orelse error.TimestampOutOfRange;
}

// Trailing layout after records/side tables:
// `[name blob][optional ZFID][ZFNM footer]`.
// `ZFID` sits before the footer so legacy loaders that seek `ZFNM` at EOF still work.
const ZfiTrailingMeta = struct {
    payload_end: usize,
    name_blob_len: u64,
    source_mtime_ns: ?u64,
};

fn parseZfiTrailingMeta(zfi_data: []const u8) ?ZfiTrailingMeta {
    const footer = parseZfiNameFooterAtEnd(zfi_data) orelse return null;
    const footer_start = zfi_data.len - footer.footer_bytes;

    var source_mtime_ns: ?u64 = null;
    var payload_end = footer_start;
    if (footer_start >= ZFI_SOURCE_ID_BYTES) {
        if (decodeZfiSourceIdBytes(zfi_data[footer_start - ZFI_SOURCE_ID_BYTES .. footer_start])) |mtime_ns| {
            source_mtime_ns = mtime_ns;
            payload_end = footer_start - ZFI_SOURCE_ID_BYTES;
        }
    }

    return .{
        .payload_end = payload_end,
        .name_blob_len = footer.name_blob_len,
        .source_mtime_ns = source_mtime_ns,
    };
}

fn parseZfiNameFooterAtEnd(zfi_data: []const u8) ?struct {
    name_blob_len: u64,
    footer_bytes: usize,
} {
    if (zfi_data.len >= ZFI_NAME_FOOTER_BYTES) {
        if (decodeZfiNameFooterBytes(zfi_data[zfi_data.len - ZFI_NAME_FOOTER_BYTES ..])) |name_blob_len| {
            return .{ .name_blob_len = name_blob_len, .footer_bytes = ZFI_NAME_FOOTER_BYTES };
        }
    }

    // Legacy 16-byte footer: magic + 4 pad bytes + little-endian u64 length.
    if (zfi_data.len >= ZFI_NAME_FOOTER_LEGACY_BYTES) {
        const tail = zfi_data[zfi_data.len - ZFI_NAME_FOOTER_LEGACY_BYTES ..];
        if (std.mem.eql(u8, tail[0..4], &ZFI_NAME_FOOTER_MAGIC)) {
            const name_blob_len = std.mem.readInt(u64, tail[8..16], .little);
            return .{ .name_blob_len = name_blob_len, .footer_bytes = ZFI_NAME_FOOTER_LEGACY_BYTES };
        }
    }

    return null;
}

/// Frozen 16-byte `.zfi` header.
pub const ZfiHeader = extern struct {
    magic: [4]u8,
    record_count: u32,
    source_size: u64,
};

/// One line entry for non-uniform FASTA records.
///
/// `base_start` is the 0-based sequence coordinate of the first base on this line.
/// It makes side-table lookup a binary search instead of a linear prefix sum.
pub const SideTableLine = extern struct {
    base_start: u64,
    byte_offset: u64,
    line_bytes: u64,
    line_bases: u64,
};

/// Frozen 40-byte record shared by `.zfi` and normalized in-memory `.fai` loads.
///
/// For `.zfi` indexes, `name_offset` / `name_len` point into the embedded name blob.
/// For full-map `.fai` loads, they point into the mapped `.fai` line. Streamed
/// `.fai` loads keep copied names separately. Use `LoadedIndex.recordName` when
/// the index source is not already known.
///
/// `_pad[0] & 1 == 0` selects the uniform O(1) line formula. Non-uniform records
/// store a 40-bit little-endian absolute side-table offset in `_pad[1..6]`.
pub const IndexRecord = extern struct {
    name_offset: u64,
    name_len: u16,
    _pad: [6]u8 = .{0} ** 6,
    seq_offset: u64,
    seq_len: u64,
    line_bases: u32,
    line_bytes: u32,

    /// Caller assumes the record name range was validated against `data`; violating
    /// this is undefined behavior.
    pub fn getName(self: IndexRecord, data: []const u8) []const u8 {
        return data[self.name_offset..][0..self.name_len];
    }

    pub fn isUniformWidth(self: IndexRecord) bool {
        return (self._pad[0] & NON_UNIFORM_WIDTH_FLAG) == 0;
    }

    pub fn nameInZfi(self: IndexRecord) bool {
        return (self._pad[0] & NAME_IN_ZFI_FLAG) != 0;
    }

    pub fn sideTableOffset(self: IndexRecord) u64 {
        var offset: u64 = 0;
        inline for (0..SIDE_TABLE_OFFSET_BYTES) |i| {
            offset |= @as(u64, self._pad[i + 1]) << (i * 8);
        }
        return offset;
    }

    /// Marks the record non-uniform or returns `SideTableOffsetTooLarge`.
    pub fn markNonUniform(self: *IndexRecord, offset: u64) !void {
        if (offset > MAX_SIDE_TABLE_OFFSET) return error.SideTableOffsetTooLarge;
        const preserve = self._pad[0] & NAME_IN_ZFI_FLAG;
        self._pad = .{0} ** 6;
        self._pad[0] = NON_UNIFORM_WIDTH_FLAG | preserve;
        inline for (0..SIDE_TABLE_OFFSET_BYTES) |i| {
            self._pad[i + 1] = @intCast((offset >> (i * 8)) & 0xff);
        }
    }
};

// Frozen wire layout (`plan/docs/zfi-format.md`). `asBytes` / mmap views are valid only when
// native endianness and struct layout match the little-endian contract.
comptime {
    if (builtin.cpu.arch.endian() != .little) {
        @compileError(".zfi wire format requires little-endian (see plan/docs/zfi-format.md); big-endian codec not implemented");
    }
    if (@sizeOf(ZfiHeader) != ZFI_HEADER_BYTES) @compileError("ZfiHeader size drifted from wire contract");
    if (@sizeOf(IndexRecord) != ZFI_INDEX_RECORD_BYTES) @compileError("IndexRecord size drifted from wire contract");
    if (@sizeOf(SideTableLine) != ZFI_SIDE_TABLE_LINE_BYTES) @compileError("SideTableLine size drifted from wire contract");
    if (@offsetOf(ZfiHeader, "magic") != 0) @compileError("ZfiHeader.magic offset");
    if (@offsetOf(ZfiHeader, "record_count") != 4) @compileError("ZfiHeader.record_count offset");
    if (@offsetOf(ZfiHeader, "source_size") != 8) @compileError("ZfiHeader.source_size offset");
    if (@offsetOf(IndexRecord, "name_offset") != 0) @compileError("IndexRecord.name_offset offset");
    if (@offsetOf(IndexRecord, "name_len") != 8) @compileError("IndexRecord.name_len offset");
    if (@offsetOf(IndexRecord, "_pad") != 10) @compileError("IndexRecord._pad offset");
    if (@offsetOf(IndexRecord, "seq_offset") != 16) @compileError("IndexRecord.seq_offset offset");
    if (@offsetOf(IndexRecord, "seq_len") != 24) @compileError("IndexRecord.seq_len offset");
    if (@offsetOf(IndexRecord, "line_bases") != 32) @compileError("IndexRecord.line_bases offset");
    if (@offsetOf(IndexRecord, "line_bytes") != 36) @compileError("IndexRecord.line_bytes offset");
    if (@offsetOf(SideTableLine, "base_start") != 0) @compileError("SideTableLine.base_start offset");
    if (@offsetOf(SideTableLine, "byte_offset") != 8) @compileError("SideTableLine.byte_offset offset");
    if (@offsetOf(SideTableLine, "line_bytes") != 16) @compileError("SideTableLine.line_bytes offset");
    if (@offsetOf(SideTableLine, "line_bases") != 24) @compileError("SideTableLine.line_bases offset");
}

/// Explicit little-endian header bytes (must match `asBytes` on little-endian hosts).
pub fn encodeZfiHeader(header: ZfiHeader) [ZFI_HEADER_BYTES]u8 {
    var out: [ZFI_HEADER_BYTES]u8 = undefined;
    @memcpy(out[0..4], &header.magic);
    std.mem.writeInt(u32, out[4..8], header.record_count, .little);
    std.mem.writeInt(u64, out[8..16], header.source_size, .little);
    return out;
}

/// Explicit little-endian record bytes (must match `asBytes` on little-endian hosts).
pub fn encodeIndexRecord(rec: IndexRecord) [ZFI_INDEX_RECORD_BYTES]u8 {
    var out: [ZFI_INDEX_RECORD_BYTES]u8 = undefined;
    std.mem.writeInt(u64, out[0..8], rec.name_offset, .little);
    std.mem.writeInt(u16, out[8..10], rec.name_len, .little);
    @memcpy(out[10..16], &rec._pad);
    std.mem.writeInt(u64, out[16..24], rec.seq_offset, .little);
    std.mem.writeInt(u64, out[24..32], rec.seq_len, .little);
    std.mem.writeInt(u32, out[32..36], rec.line_bases, .little);
    std.mem.writeInt(u32, out[36..40], rec.line_bytes, .little);
    return out;
}

/// Explicit little-endian side-table line bytes (must match `asBytes` on little-endian hosts).
pub fn encodeSideTableLine(line: SideTableLine) [ZFI_SIDE_TABLE_LINE_BYTES]u8 {
    var out: [ZFI_SIDE_TABLE_LINE_BYTES]u8 = undefined;
    std.mem.writeInt(u64, out[0..8], line.base_start, .little);
    std.mem.writeInt(u64, out[8..16], line.byte_offset, .little);
    std.mem.writeInt(u64, out[16..24], line.line_bytes, .little);
    std.mem.writeInt(u64, out[24..32], line.line_bases, .little);
    return out;
}

/// Selects the records, names, and lookup state retained by a loader.
pub const LoadMode = union(enum) {
    lookup_full_map,
    /// Stats composition: retain every record and only the two report names for streamed FAI.
    stats_scan,
    /// Positional GET: retain only the first sidecar record for each requested name.
    positional: []const []const u8,
};

/// Loaded FASTA + index state (from `.zfi` or samtools-compatible `.fai`).
///
/// Ownership (one owner each; `deinit()` is the only cleanup entry point):
/// - **FASTA file**: `fasta_file` remains open until `deinit`.
/// - **Sidecar maps**: `deinit` destroys the optional `.zfi` / `.fai` maps. Byte
///   slices and typed views borrow from those maps.
/// - **Arena**: owns heap for streamed `.fai` record arrays and retained names,
///   plus `name_map` table storage. Reclaimed
///   only via `arena.deinit()` (do not `name_map.deinit()` on an arena-backed map).
/// - **`name_map` keys**: borrowed from a sidecar map for full loads, or from
///   arena-owned `name_slices` for positional streamed `.fai` loads.
/// - **`io`**: borrowed for reads and cleanup; caller must keep it alive until `deinit`.
///
/// GET and stats load through this type and only call `deinit()`. Validator maps
/// the FASTA independently.
pub const LoadedIndex = struct {
    io: std.Io,
    fasta_file: std.Io.File,
    zfi_map: ?std.Io.File.MemoryMap = null,
    fai_map: ?std.Io.File.MemoryMap = null,
    records: []const IndexRecord,
    name_map: ?std.StringHashMapUnmanaged(u32) = null,
    /// Arena-owned names: one per retained positional record, or stats extrema only.
    name_slices: []const []const u8 = &.{},
    /// Record indices for the two extrema names retained by streamed FAI stats loads.
    stats_name_indices: ?[2]usize = null,
    /// Borrow into `zfi_map` when `.zfi` embeds a name blob.
    name_blob: ?[]const u8 = null,
    fai_data: ?platform.MappedBytes = null,
    fasta_size: u64,
    zfi_data: ?platform.MappedBytes,
    zfi_side_start: usize = 0,
    zfi_side_end: usize = 0,
    source: IndexSource,
    arena: std.heap.ArenaAllocator,

    pub const IndexSource = enum { zfi, fai };

    fn recordNameData(self: *const LoadedIndex) []const u8 {
        return if (self.source == .zfi) self.name_blob.? else self.fai_data.?;
    }

    /// Returns a retained record name, or `null` when the selected load mode did
    /// not retain that record's name.
    pub fn recordName(self: *const LoadedIndex, rec_idx: usize) ?[]const u8 {
        if (rec_idx >= self.records.len) return null;
        if (self.stats_name_indices) |indices| {
            if (rec_idx == indices[0]) return self.name_slices[0];
            if (rec_idx == indices[1]) return self.name_slices[1];
            return null;
        }
        if (self.name_slices.len > 0) {
            return self.name_slices[rec_idx];
        }
        const rec = self.records[rec_idx];
        if (self.source == .zfi) return rec.getName(self.name_blob.?);
        if (self.fai_data) |data| return rec.getName(data);
        return null;
    }

    /// Releases all mappings, arena allocations, and the FASTA descriptor.
    pub fn deinit(self: *LoadedIndex) void {
        if (self.zfi_map) |*m| m.destroy(self.io);
        if (self.fai_map) |*m| m.destroy(self.io);
        // `name_map` table bytes and copied names are arena-owned. The arena
        // reclaims them together; individual container deinit calls are redundant.
        self.arena.deinit();
        self.fasta_file.close(self.io);
    }

    /// Returns the first retained record whose name exactly matches `name`.
    pub fn lookupName(self: *const LoadedIndex, name: []const u8) ?usize {
        if (self.name_map) |name_map| {
            const rec_idx = name_map.get(name) orelse return null;
            return rec_idx;
        }

        if (self.stats_name_indices) |indices| {
            for (self.name_slices, indices) |retained, rec_idx| {
                if (std.mem.eql(u8, retained, name)) return rec_idx;
            }
            return null;
        }

        if (self.name_slices.len > 0) {
            for (self.name_slices, 0..) |rec_name, i| {
                if (std.mem.eql(u8, rec_name, name)) {
                    return i;
                }
            }
            return null;
        }

        for (self.records, 0..) |rec, i| {
            if (std.mem.eql(u8, rec.getName(self.recordNameData()), name)) {
                return i;
            }
        }
        return null;
    }

    /// Returns borrowed lines for a loaded non-uniform record, or empty if unavailable.
    pub fn sideTableLines(self: *const LoadedIndex, rec: IndexRecord) []const SideTableLine {
        if (rec.isUniformWidth()) return &.{};
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

pub const LoadIndexError = error{
    FileNotFound,
    AccessDenied,
    Io,
    EmptyFile,
    PathTooLong,
    MmapFailed,
    CorruptIndex,
    StaleIndex,
    NoIndexFound,
    OutOfMemory,
};

// --- Shared error helper ---

pub fn printErrorAndExit(comptime fmt: []const u8, args: anytype) noreturn {
    std.debug.print(fmt, args);
    std.process.exit(1);
}

// --- Shared index loader ---

/// Load the index for a FASTA file. A present `.zfi` is authoritative;
/// `.fai` is considered only when `.zfi` is absent.
/// The caller must call `deinit()` on the returned LoadedIndex.
pub fn loadIndex(allocator: std.mem.Allocator, io: std.Io, fasta_path: []const u8) LoadedIndex {
    return loadIndexOrExit(allocator, io, fasta_path, .lookup_full_map);
}

/// Load all stats metadata while keeping the FASTA descriptor-backed.
pub fn loadIndexForStats(allocator: std.mem.Allocator, io: std.Io, fasta_path: []const u8) LoadedIndex {
    return loadIndexOrExit(allocator, io, fasta_path, .stats_scan);
}

/// Keeps the FASTA open without retaining its contents in memory.
pub fn loadIndexForGet(allocator: std.mem.Allocator, io: std.Io, fasta_path: []const u8, mode: LoadMode) LoadedIndex {
    return loadIndexOrExit(allocator, io, fasta_path, mode);
}

fn loadIndexOrExit(allocator: std.mem.Allocator, io: std.Io, fasta_path: []const u8, mode: LoadMode) LoadedIndex {
    return loadIndexCheckedWithMode(allocator, io, fasta_path, mode) catch |err| switch (err) {
        error.FileNotFound => printErrorAndExit("error: file not found: {s}\n", .{fasta_path}),
        error.AccessDenied => printErrorAndExit("error: access denied: {s}\n", .{fasta_path}),
        error.EmptyFile => printErrorAndExit("error: file is empty: {s}\n", .{fasta_path}),
        error.PathTooLong => printErrorAndExit("error: path too long\n", .{}),
        error.MmapFailed => printErrorAndExit("error: failed to mmap file: {s}\n", .{fasta_path}),
        error.CorruptIndex => printErrorAndExit("error: corrupt index file for: {s}\n", .{fasta_path}),
        error.StaleIndex => printErrorAndExit(
            "error: index is stale (FASTA size or mtime does not match the index). Re-run 'z-fasta index <path>'.\n",
            .{},
        ),
        error.NoIndexFound => printErrorAndExit("error: no index found for {s}. Run 'z-fasta index {s}' first.\n", .{ fasta_path, fasta_path }),
        error.OutOfMemory => printErrorAndExit("error: out of memory loading index\n", .{}),
        error.Io => printErrorAndExit("error: failed to load index for: {s}\n", .{fasta_path}),
    };
}

/// Load the index for a FASTA file with typed errors for testing.
pub fn loadIndexChecked(allocator: std.mem.Allocator, io: std.Io, fasta_path: []const u8) LoadIndexError!LoadedIndex {
    return loadIndexCheckedWithMode(allocator, io, fasta_path, .lookup_full_map);
}

/// Loads one index mode while returning typed errors instead of exiting.
pub fn loadIndexCheckedWithMode(
    allocator: std.mem.Allocator,
    io: std.Io,
    fasta_path: []const u8,
    mode: LoadMode,
) LoadIndexError!LoadedIndex {
    const fasta_file = std.Io.Dir.cwd().openFile(io, fasta_path, .{}) catch |err| switch (err) {
        error.FileNotFound => return error.FileNotFound,
        error.AccessDenied => return error.AccessDenied,
        else => return error.Io,
    };
    errdefer fasta_file.close(io);

    const fasta_stat = fasta_file.stat(io) catch return error.Io;

    if (fasta_stat.size == 0) {
        return error.EmptyFile;
    }

    // A present `.zfi` is authoritative, including when it is invalid.
    var zfi_path_buf: [4096]u8 = undefined;
    const zfi_path = std.fmt.bufPrint(&zfi_path_buf, "{s}.zfi", .{fasta_path}) catch return error.PathTooLong;

    if (try tryLoadZfi(allocator, io, zfi_path, fasta_file, fasta_stat, mode)) |loaded| return loaded;

    // Try `.fai` only when `.zfi` is absent.
    var fai_path_buf: [4096]u8 = undefined;
    const fai_path = std.fmt.bufPrint(&fai_path_buf, "{s}.fai", .{fasta_path}) catch return error.PathTooLong;

    return (try tryLoadFai(allocator, io, fai_path, fasta_file, fasta_stat, mode)) orelse error.NoIndexFound;
}

// Current indexes embed every name and end with `[name blob][ZFID][ZFNM footer]`.
// The source-identity block remains optional for supported early embedded-name files.
const ZfiNameLayout = struct {
    bytes: []const u8,
    start: usize,
    source_mtime_ns: ?u64,
};

fn resolveZfiNameLayout(zfi_data: []const u8, records_end: usize, records: []const IndexRecord) ?ZfiNameLayout {
    for (records) |rec| {
        if (!rec.nameInZfi()) return null;
    }

    const t = parseZfiTrailingMeta(zfi_data) orelse return null;
    const blob_len = std.math.cast(usize, t.name_blob_len) orelse return null;
    const blob_start = std.math.sub(usize, t.payload_end, blob_len) catch return null;
    if (blob_start < records_end) return null;
    return .{
        .bytes = zfi_data[blob_start..t.payload_end],
        .start = blob_start,
        .source_mtime_ns = t.source_mtime_ns,
    };
}

fn zfiRecordsEnd(record_count: u32) ?usize {
    const records_bytes = std.math.mul(usize, @as(usize, record_count), @sizeOf(IndexRecord)) catch return null;
    return std.math.add(usize, @sizeOf(ZfiHeader), records_bytes) catch null;
}

fn buildNameMapFast(
    allocator: std.mem.Allocator,
    records: []const IndexRecord,
    name_data: []const u8,
) LoadIndexError!std.StringHashMapUnmanaged(u32) {
    var name_map: std.StringHashMapUnmanaged(u32) = .empty;
    name_map.ensureTotalCapacity(allocator, @intCast(records.len)) catch return error.OutOfMemory;
    for (records, 0..) |rec, i| {
        const name = rec.getName(name_data);
        const entry = name_map.getOrPutAssumeCapacity(name);
        if (!entry.found_existing) entry.value_ptr.* = @intCast(i);
    }
    return name_map;
}

fn tryLoadZfi(
    backing_allocator: std.mem.Allocator,
    io: std.Io,
    zfi_path: []const u8,
    fasta_file: std.Io.File,
    fasta_stat: std.Io.File.Stat,
    mode: LoadMode,
) LoadIndexError!?LoadedIndex {
    const zfi_file = std.Io.Dir.cwd().openFile(io, zfi_path, .{}) catch |err| switch (err) {
        error.FileNotFound => return null,
        error.AccessDenied => return error.AccessDenied,
        else => return error.Io,
    };
    defer zfi_file.close(io);

    const zfi_stat = zfi_file.stat(io) catch return error.Io;

    // Check for 0-byte file before mmap (POSIX returns EINVAL)
    if (zfi_stat.size == 0) {
        return error.CorruptIndex;
    }

    // Staleness: index file older than FASTA (both formats). `.fai` has no embedded
    // source size/mtime, so this plus "index not older" is its entire identity check.
    if (zfi_stat.mtime.nanoseconds < fasta_stat.mtime.nanoseconds) {
        return error.StaleIndex;
    }

    const zfi_size = std.math.cast(usize, zfi_stat.size) orelse return error.MmapFailed;
    var zfi_map = platform.mapFileReadOnly(io, zfi_file, zfi_size) catch return error.MmapFailed;
    errdefer zfi_map.destroy(io);

    const zfi_data: platform.MappedBytes = zfi_map.memory;

    if (zfi_data.len < @sizeOf(ZfiHeader)) {
        return error.CorruptIndex;
    }

    const header: *const ZfiHeader = @ptrCast(@alignCast(zfi_data.ptr));
    if (!std.mem.eql(u8, &header.magic, &ZFI_MAGIC)) {
        return error.CorruptIndex;
    }

    if (header.source_size != fasta_stat.size) {
        return error.StaleIndex;
    }

    // Empty catalogs are not useful and match the `.fai` reject-empty policy.
    if (header.record_count == 0) return error.CorruptIndex;

    // Checked records extent before forming the typed record slice.
    const records_end = zfiRecordsEnd(header.record_count) orelse return error.CorruptIndex;
    if (zfi_data.len < records_end) return error.CorruptIndex;

    const record_bytes = zfi_data[@sizeOf(ZfiHeader)..records_end];
    const records: []const IndexRecord = @as(
        [*]const IndexRecord,
        @ptrCast(@alignCast(record_bytes.ptr)),
    )[0..header.record_count];

    const name_layout = resolveZfiNameLayout(zfi_data, records_end, records) orelse return error.CorruptIndex;
    const name_blob = name_layout.bytes;

    // Strong identity: production trailers store the FASTA mtime at index time.
    // Exact match catches same-size content replacement even when the `.zfi` file
    // mtime was touched forward. Legacy files without a trailer keep the weaker
    // index-mtime check above only.
    if (name_layout.source_mtime_ns) |stored_mtime| {
        const source_mtime = timestampToNs(fasta_stat.mtime) catch return error.StaleIndex;
        if (stored_mtime != source_mtime) {
            return error.StaleIndex;
        }
    }

    // Side tables live after records and before the embedded name blob.
    const side_region_start = records_end;
    const side_region_end = name_layout.start;

    if (!validateZfiRecords(records, fasta_stat.size, zfi_data, name_blob, side_region_start, side_region_end)) {
        return error.CorruptIndex;
    }

    var arena = std.heap.ArenaAllocator.init(backing_allocator);
    errdefer arena.deinit();
    const allocator = arena.allocator();

    const build_name_map = switch (mode) {
        .lookup_full_map => true,
        .positional => |names| names.len > 1,
        .stats_scan => false,
    };

    const name_map = if (build_name_map)
        buildNameMapFast(allocator, records, name_blob) catch return error.OutOfMemory
    else
        null;

    return LoadedIndex{
        .io = io,
        .fasta_file = fasta_file,
        .zfi_map = zfi_map,
        .records = records,
        .name_map = name_map,
        .name_slices = &.{},
        .name_blob = name_blob,
        .fasta_size = fasta_stat.size,
        .zfi_data = zfi_data,
        .zfi_side_start = side_region_start,
        .zfi_side_end = side_region_end,
        .source = .zfi,
        .arena = arena,
    };
}

fn tryLoadFai(
    backing_allocator: std.mem.Allocator,
    io: std.Io,
    fai_path: []const u8,
    fasta_file: std.Io.File,
    fasta_stat: std.Io.File.Stat,
    mode: LoadMode,
) LoadIndexError!?LoadedIndex {
    const fai_file = std.Io.Dir.cwd().openFile(io, fai_path, .{}) catch |err| switch (err) {
        error.FileNotFound => return null,
        error.AccessDenied => return error.AccessDenied,
        else => return error.Io,
    };
    defer fai_file.close(io);

    const fai_stat = fai_file.stat(io) catch return error.Io;

    // `.fai` identity is mtime-only: the text format has no source size or mtime field.
    // Same-size FASTA replacement is detected only when the FASTA mtime moves forward.
    if (fai_stat.mtime.nanoseconds < fasta_stat.mtime.nanoseconds) {
        return error.StaleIndex;
    }

    if (fai_stat.size == 0) {
        return error.CorruptIndex;
    }

    switch (mode) {
        .stats_scan => return @as(?LoadedIndex, try loadFaiStreamed(backing_allocator, io, fai_file, fasta_file, fasta_stat, .stats_scan)),
        .positional => |names| return @as(?LoadedIndex, try loadFaiStreamed(backing_allocator, io, fai_file, fasta_file, fasta_stat, .{ .matched = names })),
        .lookup_full_map => {},
    }

    const fai_size = std.math.cast(usize, fai_stat.size) orelse return error.MmapFailed;
    var fai_map = platform.mapFileReadOnly(io, fai_file, fai_size) catch return error.MmapFailed;

    errdefer fai_map.destroy(io);

    var arena = std.heap.ArenaAllocator.init(backing_allocator);
    errdefer arena.deinit();
    const allocator = arena.allocator();

    const fai_data: platform.MappedBytes = fai_map.memory;

    // Newline count undercounts by one when the final `.fai` line has no trailing `\n`.
    var approx_records = std.mem.count(u8, fai_data, "\n");
    if (fai_data.len > 0 and fai_data[fai_data.len - 1] != '\n') {
        approx_records = std.math.add(usize, approx_records, 1) catch return error.CorruptIndex;
    }
    // Reject sizes that would overflow the records allocation.
    _ = std.math.mul(usize, approx_records, @sizeOf(IndexRecord)) catch return error.CorruptIndex;
    if (approx_records == 0) return error.CorruptIndex;
    const records = allocator.alloc(IndexRecord, approx_records) catch return error.OutOfMemory;

    var record_count: usize = 0;
    var pos: usize = 0;
    while (pos < fai_data.len) {
        const line_start = pos;
        const rel_eol = std.mem.findScalar(u8, fai_data[pos..], '\n') orelse fai_data.len - pos;
        const line_len = rel_eol;
        const line_end = std.math.add(usize, line_start, line_len) catch return error.CorruptIndex;
        // Advance past `\n` only when present; a final line may omit it.
        if (line_end < fai_data.len) {
            pos = std.math.add(usize, line_end, 1) catch return error.CorruptIndex;
        } else {
            pos = line_end;
        }
        if (line_len == 0) continue;

        const line = fai_data[line_start..][0..line_len];
        const fields = parseFaiIndexLine(line) catch return error.CorruptIndex;

        const rec = IndexRecord{
            .name_offset = line_start,
            .name_len = fields.name_len,
            .seq_offset = fields.seq_offset,
            .seq_len = fields.seq_len,
            .line_bases = fields.line_bases,
            .line_bytes = fields.line_bytes,
        };
        if (!isValidUniformSequenceGeometry(rec, fasta_stat.size)) return error.CorruptIndex;

        if (record_count >= records.len) return error.CorruptIndex;
        records[record_count] = rec;
        record_count += 1;
    }

    if (record_count == 0) {
        return error.CorruptIndex;
    }

    const loaded_records = records[0..record_count];
    const name_map = buildNameMapFast(allocator, loaded_records, fai_data) catch return error.OutOfMemory;

    return LoadedIndex{
        .io = io,
        .fasta_file = fasta_file,
        .fai_map = fai_map,
        .records = loaded_records,
        .name_map = name_map,
        .name_slices = &.{},
        .fai_data = fai_data,
        .fasta_size = fasta_stat.size,
        .zfi_data = null,
        .source = .fai,
        .arena = arena,
    };
}

const FaiStreamMode = union(enum) {
    matched: []const []const u8,
    stats_scan,
};

const FaiNameRef = struct {
    offset: u64,
    len: usize,
};

fn loadFaiStreamed(
    backing_allocator: std.mem.Allocator,
    io: std.Io,
    fai_file: std.Io.File,
    fasta_file: std.Io.File,
    fasta_stat: std.Io.File.Stat,
    mode: FaiStreamMode,
) LoadIndexError!LoadedIndex {
    var arena = std.heap.ArenaAllocator.init(backing_allocator);
    errdefer arena.deinit();
    const allocator = arena.allocator();

    var records_list: std.ArrayList(IndexRecord) = .empty;
    var slices_list: std.ArrayList([]const u8) = .empty;
    var requested: std.StringHashMapUnmanaged(bool) = .empty;
    defer requested.deinit(backing_allocator);
    switch (mode) {
        .matched => |names| {
            requested.ensureTotalCapacity(backing_allocator, @intCast(names.len)) catch return error.OutOfMemory;
            for (names) |name| {
                const entry = requested.getOrPutAssumeCapacity(name);
                if (!entry.found_existing) entry.value_ptr.* = false;
            }
        },
        .stats_scan => {},
    }

    var io_buf: [FAI_READER_BUFFER_BYTES]u8 = undefined;
    var file_reader = fai_file.reader(io, &io_buf);
    var saw_record = false;
    var line_offset: u64 = 0;
    var shortest_ref: FaiNameRef = undefined;
    var longest_ref: FaiNameRef = undefined;
    var shortest_index: usize = 0;
    var longest_index: usize = 0;
    var shortest_len: u64 = 0;
    var longest_len: u64 = 0;

    while (true) {
        const maybe_line = file_reader.interface.takeDelimiter('\n') catch |err| switch (err) {
            error.ReadFailed => return error.Io,
            error.StreamTooLong => return error.CorruptIndex,
        };
        const line = maybe_line orelse break;
        const current_offset = line_offset;
        line_offset = std.math.add(u64, line_offset, line.len + 1) catch return error.CorruptIndex;
        if (line.len == 0) continue;

        const fields = parseFaiIndexLine(line) catch return error.CorruptIndex;

        const rec = IndexRecord{
            .name_offset = 0,
            .name_len = fields.name_len,
            .seq_offset = fields.seq_offset,
            .seq_len = fields.seq_len,
            .line_bases = fields.line_bases,
            .line_bytes = fields.line_bytes,
        };
        if (!isValidUniformSequenceGeometry(rec, fasta_stat.size)) return error.CorruptIndex;

        saw_record = true;
        const retain = switch (mode) {
            .stats_scan => true,
            .matched => blk: {
                const matched = requested.getPtr(line[0..fields.name_end]) orelse break :blk false;
                if (matched.*) break :blk false;
                matched.* = true;
                break :blk true;
            },
        };
        if (!retain) continue;

        const rec_index = records_list.items.len;
        switch (mode) {
            .matched => {
                const name = allocator.dupe(u8, line[0..fields.name_end]) catch return error.OutOfMemory;
                try slices_list.append(allocator, name);
            },
            .stats_scan => {
                const name_ref = FaiNameRef{ .offset = current_offset, .len = fields.name_end };
                if (rec_index == 0 or rec.seq_len < shortest_len) {
                    shortest_ref = name_ref;
                    shortest_index = rec_index;
                    shortest_len = rec.seq_len;
                }
                if (rec_index == 0 or rec.seq_len > longest_len) {
                    longest_ref = name_ref;
                    longest_index = rec_index;
                    longest_len = rec.seq_len;
                }
            },
        }
        try records_list.append(allocator, rec);
    }

    if (!saw_record) return error.CorruptIndex;

    const owned_records = records_list.toOwnedSlice(allocator) catch return error.OutOfMemory;
    const owned_slices = switch (mode) {
        .matched => slices_list.toOwnedSlice(allocator) catch return error.OutOfMemory,
        .stats_scan => blk: {
            const names = allocator.alloc([]const u8, 2) catch return error.OutOfMemory;
            names[0] = try readFaiName(allocator, io, fai_file, shortest_ref);
            names[1] = if (shortest_index == longest_index)
                names[0]
            else
                try readFaiName(allocator, io, fai_file, longest_ref);
            break :blk names;
        },
    };

    const name_map: ?std.StringHashMapUnmanaged(u32) = switch (mode) {
        .matched => blk: {
            var map: std.StringHashMapUnmanaged(u32) = .empty;
            map.ensureTotalCapacity(allocator, @intCast(owned_slices.len)) catch return error.OutOfMemory;
            for (owned_slices, 0..) |name, i| {
                map.putAssumeCapacity(name, @intCast(i));
            }
            break :blk map;
        },
        .stats_scan => null,
    };

    return LoadedIndex{
        .io = io,
        .fasta_file = fasta_file,
        .records = owned_records,
        .name_map = name_map,
        .name_slices = owned_slices,
        .stats_name_indices = switch (mode) {
            .matched => null,
            .stats_scan => .{ shortest_index, longest_index },
        },
        .fai_data = null,
        .fasta_size = fasta_stat.size,
        .zfi_data = null,
        .source = .fai,
        .arena = arena,
    };
}

fn readFaiName(
    allocator: std.mem.Allocator,
    io: std.Io,
    fai_file: std.Io.File,
    name_ref: FaiNameRef,
) LoadIndexError![]const u8 {
    const name = allocator.alloc(u8, name_ref.len) catch return error.OutOfMemory;
    if (name.len == 0) return name;
    const got = std.Io.File.readPositionalAll(fai_file, io, name, name_ref.offset) catch return error.Io;
    if (got != name.len) return error.CorruptIndex;
    return name;
}

fn parseFaiFieldU64(line: []const u8, field_start: *usize) LoadIndexError!u64 {
    if (field_start.* < line.len and line[field_start.*] == '\t') {
        field_start.* += 1;
    }
    if (field_start.* >= line.len) return error.CorruptIndex;

    const field_end = std.mem.findScalarPos(u8, line, field_start.*, '\t') orelse line.len;
    const value = parseFaiAsciiU64(line[field_start.*..field_end]) orelse return error.CorruptIndex;
    field_start.* = field_end;
    return value;
}

fn parseFaiFieldU32(line: []const u8, field_start: *usize) LoadIndexError!u32 {
    if (field_start.* < line.len and line[field_start.*] == '\t') {
        field_start.* += 1;
    }
    if (field_start.* >= line.len) return error.CorruptIndex;

    const field_end = std.mem.findScalarPos(u8, line, field_start.*, '\t') orelse line.len;
    const value = parseFaiAsciiU32(line[field_start.*..field_end]) orelse return error.CorruptIndex;
    field_start.* = field_end;
    return value;
}

// --- FAI field parsing and uniform geometry checks ---

// One parsed `.fai` text line. Shared by mmap and streaming loaders so a missing
// terminal newline cannot make the two paths disagree on field values.
const FaiLineFields = struct {
    name_end: usize,
    name_len: u16,
    seq_len: u64,
    seq_offset: u64,
    line_bases: u32,
    line_bytes: u32,
};

fn parseFaiIndexLine(line: []const u8) LoadIndexError!FaiLineFields {
    const name_end = std.mem.findScalar(u8, line, '\t') orelse return error.CorruptIndex;
    const name_len = faiNameLen(name_end) orelse return error.CorruptIndex;
    var field_start: usize = name_end + 1;
    const seq_len = try parseFaiFieldU64(line, &field_start);
    const seq_offset = try parseFaiFieldU64(line, &field_start);
    const line_bases = try parseFaiFieldU32(line, &field_start);
    const line_bytes = try parseFaiFieldU32(line, &field_start);
    return .{
        .name_end = name_end,
        .name_len = name_len,
        .seq_len = seq_len,
        .seq_offset = seq_offset,
        .line_bases = line_bases,
        .line_bytes = line_bytes,
    };
}

// `IndexRecord.name_len` is u16. Names longer than 65535 are not representable;
// empty names (`name_end == 0`, samtools `>\n…` → `.fai` line starting with tab) are.
fn faiNameLen(name_end: usize) ?u16 {
    return std.math.cast(u16, name_end);
}

// Checked decimal parse for `.fai` numeric fields. Wrapping `value * 10 + digit`
// would accept oversized ASCII as a smaller integer and pass geometry checks.
fn parseFaiAsciiU64(text: []const u8) ?u64 {
    if (text.len == 0) return null;
    var value: u64 = 0;
    for (text) |byte| {
        if (byte < '0' or byte > '9') return null;
        const digit: u64 = byte - '0';
        value = std.math.mul(u64, value, 10) catch return null;
        value = std.math.add(u64, value, digit) catch return null;
    }
    return value;
}

fn parseFaiAsciiU32(text: []const u8) ?u32 {
    const wide = parseFaiAsciiU64(text) orelse return null;
    return std.math.cast(u32, wide);
}

// Shared by `.fai` loaders and uniform `.zfi` metadata checks.
//
// Reject zero or impossible `line_bases` / `line_bytes`, `seq_offset` past EOF, and
// records whose last base falls outside the FASTA. Samtools-style `line_bytes`
// includes the line separator even when the FASTA omits a terminal newline on the
// final line, so a full `full_lines * line_bytes` span would demand a phantom
// trailing separator. Require the last base byte instead:
// `seq_offset + ((seq_len - 1) / line_bases) * line_bytes + ((seq_len - 1) % line_bases)`.
fn isValidUniformSequenceGeometry(rec: IndexRecord, fasta_len: u64) bool {
    if (rec.seq_len == 0) return false;
    if (rec.seq_offset >= fasta_len) return false;
    if (rec.line_bases == 0 or rec.line_bytes == 0) return false;
    if (rec.line_bytes < rec.line_bases) return false;

    const last_index: u128 = rec.seq_len - 1;
    const full_lines = last_index / rec.line_bases;
    const col = last_index % rec.line_bases;
    const last_byte = rec.seq_offset + full_lines * rec.line_bytes + col;
    return last_byte < fasta_len;
}

fn rangeFitsUsize(offset: u64, len: u64, limit: usize) bool {
    const off = std.math.cast(usize, offset) orelse return false;
    const len_usz = std.math.cast(usize, len) orelse return false;
    const end = std.math.add(usize, off, len_usz) catch return false;
    return end <= limit;
}

fn isValidZfiRecordMetadata(
    rec: IndexRecord,
    fasta_len: u64,
    zfi_data: platform.MappedBytes,
    name_blob: []const u8,
    side_region_start: usize,
    side_region_end: usize,
    prev_side_end: *usize,
) bool {
    if (!rangeFitsUsize(rec.name_offset, rec.name_len, name_blob.len)) return false;

    if (!rec.isUniformWidth()) {
        if (rec.seq_len == 0) return false;
        if (rec.seq_offset >= fasta_len) return false;
        const parsed = parseSideTable(zfi_data, rec, side_region_start, side_region_end, fasta_len) orelse return false;
        if (parsed.start < prev_side_end.*) return false;
        prev_side_end.* = parsed.end;
        return true;
    }

    return isValidUniformSequenceGeometry(rec, fasta_len);
}

fn validateZfiRecords(
    records: []const IndexRecord,
    fasta_len: u64,
    zfi_data: platform.MappedBytes,
    name_blob: []const u8,
    side_region_start: usize,
    side_region_end: usize,
) bool {
    if (side_region_end < side_region_start) return false;
    var prev_side_end: usize = side_region_start;
    for (records) |rec| {
        if (!isValidZfiRecordMetadata(rec, fasta_len, zfi_data, name_blob, side_region_start, side_region_end, &prev_side_end)) {
            return false;
        }
    }
    return true;
}

const ParsedSideTable = struct {
    start: usize,
    end: usize,
    lines: []const SideTableLine,
};

// Validate one non-uniform record's side table before creating typed line pointers.
//
// The table must sit in `[side_region_start, side_region_end)`, be 8-byte aligned,
// and describe `rec.seq_len` bases starting at `rec.seq_offset`. Line byte ranges
// must be ordered and non-overlapping; checked arithmetic rejects wraps before any slice.
fn parseSideTable(
    zfi_data: platform.MappedBytes,
    rec: IndexRecord,
    side_region_start: usize,
    side_region_end: usize,
    fasta_len: u64,
) ?ParsedSideTable {
    const offset = std.math.cast(usize, rec.sideTableOffset()) orelse return null;
    if (offset < side_region_start or offset >= side_region_end) return null;
    if (offset % @alignOf(u64) != 0) return null;

    const count_end = std.math.add(usize, offset, @sizeOf(u64)) catch return null;
    if (count_end > side_region_end) return null;

    const count_ptr: *const u64 = @ptrCast(@alignCast(zfi_data[offset..].ptr));
    const line_count = count_ptr.*;
    if (line_count == 0) return null;
    if (line_count > std.math.maxInt(usize) / @sizeOf(SideTableLine)) return null;

    const table_bytes = std.math.mul(usize, @as(usize, @intCast(line_count)), @sizeOf(SideTableLine)) catch return null;
    const lines_offset = count_end;
    const table_end = std.math.add(usize, lines_offset, table_bytes) catch return null;
    if (table_end > side_region_end) return null;
    if (lines_offset % @alignOf(SideTableLine) != 0) return null;

    const lines = @as(
        [*]const SideTableLine,
        @ptrCast(@alignCast(zfi_data[lines_offset..].ptr)),
    )[0..@intCast(line_count)];

    var expected_base_start: u64 = 0;
    var previous_byte_end: u64 = 0;
    for (lines, 0..) |line, i| {
        if (line.base_start != expected_base_start) return null;
        if (line.line_bases == 0 or line.line_bytes == 0) return null;
        if (line.line_bytes < line.line_bases) return null;
        const byte_end = std.math.add(u64, line.byte_offset, line.line_bytes) catch return null;
        if (byte_end > fasta_len) return null;
        if (i == 0) {
            // Writer invariant: first base-bearing line begins at seq_offset.
            if (line.byte_offset != rec.seq_offset) return null;
        } else if (line.byte_offset < previous_byte_end) {
            return null;
        }

        expected_base_start = std.math.add(u64, expected_base_start, line.line_bases) catch return null;
        previous_byte_end = byte_end;
    }

    if (expected_base_start != rec.seq_len) return null;
    return .{ .start = offset, .end = table_end, .lines = lines };
}
