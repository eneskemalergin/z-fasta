//! `.zfi` binary and `.fai` text index load/write helpers (`LoadedIndex`, staleness, geometry).
//!
//! On-disk `.zfi` bytes follow the frozen little-endian contract in `plan/docs/zfi-format.md`.
//! Host `extern struct` views are allowed only when comptime checks prove native layout
//! matches that contract (little-endian + fixed sizes/offsets).

const std = @import("std");
const builtin = @import("builtin");
const platform = @import("platform.zig");
pub const diagnostics = @import("diagnostics.zig");

// --- Types shared by indexer, getter, and stats ---

/// Four-byte magic and version for the supported `.zfi` wire format.
pub const ZFI_MAGIC: [4]u8 = .{ 'Z', 'F', 'I', 0x01 };

const NON_UNIFORM_WIDTH_FLAG: u8 = 1;
/// Record flag indicating that its name is embedded in the `.zfi` name blob.
pub const NAME_IN_ZFI_FLAG: u8 = 2;
const SIDE_TABLE_OFFSET_BYTES = 5;
pub const MAX_SIDE_TABLE_OFFSET = (1 << (SIDE_TABLE_OFFSET_BYTES * 8)) - 1;
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
    source_mtime_offset: ?usize,
};

fn parseZfiTrailingMeta(zfi_data: []const u8) ?ZfiTrailingMeta {
    const footer = parseZfiNameFooterAtEnd(zfi_data) orelse return null;
    const footer_start = zfi_data.len - footer.footer_bytes;

    var source_mtime_ns: ?u64 = null;
    var source_mtime_offset: ?usize = null;
    var payload_end = footer_start;
    if (footer_start >= ZFI_SOURCE_ID_BYTES) {
        if (decodeZfiSourceIdBytes(zfi_data[footer_start - ZFI_SOURCE_ID_BYTES .. footer_start])) |mtime_ns| {
            source_mtime_ns = mtime_ns;
            source_mtime_offset = footer_start - ZFI_SOURCE_ID_BYTES + 4;
            payload_end = footer_start - ZFI_SOURCE_ID_BYTES;
        }
    }

    return .{
        .payload_end = payload_end,
        .name_blob_len = footer.name_blob_len,
        .source_mtime_ns = source_mtime_ns,
        .source_mtime_offset = source_mtime_offset,
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
/// For full-map `.fai` loads, they point into the mapped `.fai` line. Matched
/// `.zfi` and streamed `.fai` loads keep copied names separately. Use `LoadedIndex.recordName` when
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
    /// Positional GET: retain first FAI matches; keep the complete mapped ZFI catalog.
    positional: []const []const u8,
    /// Seekable batch GET: validate the complete sidecar and retain only first matches.
    matched: []const []const u8,
};

/// Loaded FASTA + index state (from `.zfi` or samtools-compatible `.fai`).
///
/// Ownership (one owner each; `deinit()` is the only cleanup entry point):
/// - **FASTA file**: `fasta_file` remains open until `deinit`.
/// - **Sidecar maps**: `deinit` destroys the optional `.zfi` / `.fai` maps. Byte
///   slices and typed views borrow from those maps.
/// - **Arena**: owns heap for matched `.zfi` and streamed `.fai` records, names,
///   and retained side-table bytes,
///   plus `name_map` table storage. Reclaimed
///   only via `arena.deinit()` (do not `name_map.deinit()` on an arena-backed map).
/// - **`name_map` keys**: borrowed from a sidecar map for full loads, or from
///   arena-owned `name_slices` for matched loads.
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
    /// Arena-owned names: one per retained matched record, or stats extrema only.
    name_slices: []const []const u8 = &.{},
    /// Record indices for the two extrema names retained by streamed FAI stats loads.
    stats_name_indices: ?[2]usize = null,
    /// Borrow into `zfi_map` for a mapped `.zfi` load.
    name_blob: ?[]const u8 = null,
    fai_data: ?platform.MappedBytes = null,
    fasta_size: u64,
    /// Complete mapped ZFI bytes, or arena-owned side-table bytes for matched loads.
    zfi_data: ?[]align(@alignOf(SideTableLine)) const u8,
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
    return loadIndexCheckedWithModeLegacy(allocator, io, fasta_path, mode) catch |err| switch (err) {
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

/// Loads a complete catalog. Reset Failure once in the owning operation before
/// calling. Path labels borrow fasta_path through handling or Failure.reset.
pub fn loadIndexChecked(allocator: std.mem.Allocator, io: std.Io, fasta_path: []const u8, failure: *diagnostics.Failure) diagnostics.Error!LoadedIndex {
    return loadIndexCheckedWithMode(allocator, io, fasta_path, .lookup_full_map, failure);
}

/// Loads one mode without resetting the owning operation's diagnostics.
/// Missing optional sidecars are silent; a present ZFI remains authoritative.
pub fn loadIndexCheckedWithMode(
    allocator: std.mem.Allocator,
    io: std.Io,
    fasta_path: []const u8,
    mode: LoadMode,
    failure: *diagnostics.Failure,
) diagnostics.Error!LoadedIndex {
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
    var context = LoadDiagnostics{ .failure = failure, .at = .{ .source = .fasta } };
    const fasta_file = std.Io.Dir.cwd().openFile(io, fasta_path, .{}) catch |err| return context.ioFailure(.open, err);
    errdefer fasta_file.close(io);
    const fasta_stat = fasta_file.stat(io) catch |err| return context.ioFailure(.stat, err);
    if (fasta_stat.size == 0) return context.fail(.{ .@"fasta.no_records" = {} });

    var zfi_path_buf: [4096]u8 = undefined;
    context.at.source = .zfi;
    const zfi_path = std.fmt.bufPrint(&zfi_path_buf, "{s}.zfi", .{fasta_path}) catch return context.pathLimit(zfi_path_buf.len);
    if (try tryLoadZfi(allocator, io, zfi_path, fasta_file, fasta_stat, mode, &context)) |loaded| return loaded;

    var fai_path_buf: [4096]u8 = undefined;
    context.at.source = .fai;
    const fai_path = std.fmt.bufPrint(&fai_path_buf, "{s}.fai", .{fasta_path}) catch return context.pathLimit(fai_path_buf.len);
    return (try tryLoadFai(allocator, io, fai_path, fasta_file, fasta_stat, mode, &context)) orelse context.fail(.{ .@"index.not_found" = .{ .zfi_source = .zfi, .fai_source = .fai } });
}

/// Temporary old-error adapter for commands awaiting their diagnostic conversion.
pub fn loadIndexCheckedLegacy(allocator: std.mem.Allocator, io: std.Io, fasta_path: []const u8) LoadIndexError!LoadedIndex {
    return loadIndexCheckedWithModeLegacy(allocator, io, fasta_path, .lookup_full_map);
}

pub fn loadIndexCheckedWithModeLegacy(allocator: std.mem.Allocator, io: std.Io, fasta_path: []const u8, mode: LoadMode) LoadIndexError!LoadedIndex {
    var failure: diagnostics.Failure = undefined;
    failure.reset(if (mode == .stats_scan) .stats else .get);
    return loadIndexCheckedWithMode(allocator, io, fasta_path, mode, &failure) catch |err| {
        return switch (failure.primary().?.code()) {
            .@"io.file_not_found" => error.FileNotFound,
            .@"io.access_denied" => error.AccessDenied,
            .@"io.invalid_path" => error.PathTooLong,
            .@"io.map_failed" => error.MmapFailed,
            .@"fasta.no_records" => error.EmptyFile,
            .@"index.older_than_source", .@"index.source_size_mismatch", .@"index.source_mtime_mismatch" => error.StaleIndex,
            .@"index.not_found" => error.NoIndexFound,
            .@"resource.line_buffer_exceeded" => error.CorruptIndex,
            .@"resource.capacity_limit_exceeded" => if (failure.primary().?.details.@"resource.capacity_limit_exceeded".capacity_role == .map_bytes) error.MmapFailed else error.OutOfMemory,
            else => switch (err) {
                error.OutOfMemory => error.OutOfMemory,
                error.InvalidIndex => error.CorruptIndex,
                else => error.Io,
            },
        };
    };
}

const LoadDiagnostics = struct {
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

    fn fail(self: *const LoadDiagnostics, details: diagnostics.Details) diagnostics.Error {
        @branchHint(.cold);
        const diagnostic = diagnostics.Diagnostic{ .details = details, .stage = failureStage(details), .at = self.at };
        if (self.deferred) {
            self.failure.deferFailure(diagnostic);
            return diagnostics.metadata(diagnostic.code()).returned_error;
        }
        return self.failure.fail(diagnostic);
    }

    fn deferFailure(self: *const LoadDiagnostics, details: diagnostics.Details) void {
        self.failure.deferFailure(.{ .details = details, .stage = failureStage(details), .at = self.at });
    }

    fn ioFailure(self: *const LoadDiagnostics, stage: diagnostics.Stage, err: anyerror) diagnostics.Error {
        @branchHint(.cold);
        const details: diagnostics.Details = switch (stage) {
            .open => switch (err) {
                error.FileNotFound => .{ .@"io.file_not_found" = {} },
                error.AccessDenied => .{ .@"io.access_denied" = {} },
                else => .{ .@"io.open_failed" = {} },
            },
            .stat => .{ .@"io.stat_failed" = {} },
            .read => .{ .@"io.read_failed" = .{} },
            else => unreachable,
        };
        return self.failure.fail(.{ .details = details, .stage = stage, .at = self.at, .cause = .{ .zig = err } });
    }

    fn readFailure(self: *const LoadDiagnostics, offset: u64, requested: u64, err: anyerror) diagnostics.Error {
        var at = self.at;
        at.location.offset = offset;
        return self.failure.fail(.{ .details = .{ .@"io.read_failed" = .{ .requested_bytes = requested } }, .stage = .read, .at = at, .cause = .{ .zig = err } });
    }

    fn lineBuffer(self: *const LoadDiagnostics, observed: usize) diagnostics.Error {
        return self.failure.fail(.{ .details = .{ .@"resource.line_buffer_exceeded" = .{ .buffer_bytes = FAI_READER_BUFFER_BYTES, .observed_bytes = .{ .at_least = observed } } }, .stage = .read, .at = self.at, .cause = .{ .zig = error.StreamTooLong } });
    }

    fn allocation(self: *const LoadDiagnostics, role: diagnostics.AllocationRole, bytes: ?u64) diagnostics.Error {
        @branchHint(.cold);
        return self.failure.fail(.{ .details = .{ .@"resource.out_of_memory" = .{ .allocation_role = role, .requested_bytes = bytes } }, .stage = .allocate, .at = self.at, .cause = .{ .zig = error.OutOfMemory } });
    }

    fn mapping(self: *const LoadDiagnostics, bytes: u64, err: ?anyerror) diagnostics.Error {
        @branchHint(.cold);
        return self.failure.fail(.{ .details = .{ .@"io.map_failed" = .{ .map_bytes = bytes } }, .stage = .map, .at = self.at, .cause = if (err) |cause| .{ .zig = cause } else null });
    }

    fn pathLimit(self: *const LoadDiagnostics, limit: u64) diagnostics.Error {
        @branchHint(.cold);
        return self.failure.fail(.{ .details = .{ .@"io.invalid_path" = .{ .reason = .too_long, .limit_bytes = limit } }, .stage = .open, .at = self.at });
    }

    fn capacity(self: *const LoadDiagnostics, role: diagnostics.CapacityRole, limit: u128, required: u128) diagnostics.Error {
        return self.failure.fail(.{ .details = .{ .@"resource.capacity_limit_exceeded" = .{ .capacity_role = role, .limit = limit, .required = .{ .exact = required } } }, .stage = switch (role) {
            .map_bytes => .map,
            .name_map, .records => .allocate,
            .input_bytes => .read,
            else => .parse,
        }, .at = self.at });
    }

    inline fn nextOrdinal(self: *const LoadDiagnostics, value: u64, counter: diagnostics.Counter) diagnostics.Error!u64 {
        return std.math.add(u64, value, 1) catch return self.failure.fail(.{ .details = .{ .@"resource.count_overflow" = .{ .counter = counter, .max_value = std.math.maxInt(u64), .current = value, .increment = 1 } }, .stage = .parse, .at = self.at });
    }

    fn overflow(self: *const LoadDiagnostics, field: diagnostics.IndexField, arithmetic: diagnostics.Arithmetic, lhs: u128, rhs: u128, max: u128) diagnostics.Error {
        return self.fail(.{ .@"index.arithmetic_overflow" = .{ .field = field, .arithmetic = arithmetic, .lhs = lhs, .rhs = rhs, .max_value = max } });
    }

    fn truncated(self: *const LoadDiagnostics, section: diagnostics.Section, required: u128, available: u64) diagnostics.Error {
        return self.fail(.{ .@"index.truncated" = .{ .section = section, .required_end = required, .available_end = available } });
    }

    fn shortRead(self: *const LoadDiagnostics, offset: u64, requested: u64, received: u64) diagnostics.Error {
        var at = self.at;
        at.location.offset = offset;
        return self.failure.fail(.{ .details = .{ .@"io.unexpected_eof" = .{ .exact = .{ .requested_bytes = requested, .received_bytes = received } } }, .stage = .read, .at = at });
    }
};

// Current indexes embed every name and end with `[name blob][ZFID][ZFNM footer]`.
// The source-identity block remains optional for supported early embedded-name files.
const ZfiNameLayout = struct {
    bytes: []const u8,
    start: usize,
    source_mtime_ns: ?u64,
    source_mtime_offset: ?usize,
};

fn resolveZfiNameLayout(zfi_data: []const u8, records_end: usize, records: []const IndexRecord, context: *LoadDiagnostics) diagnostics.Error!ZfiNameLayout {
    for (records, 0..) |rec, i| {
        if (!rec.nameInZfi()) {
            context.at = .{ .source = .zfi, .record = i + 1, .location = .{ .offset = ZFI_HEADER_BYTES + i * ZFI_INDEX_RECORD_BYTES + 10 } };
            return context.fail(.{ .@"index.name_not_embedded" = .{ .flags = rec._pad[0], .required_flag = NAME_IN_ZFI_FLAG } });
        }
    }
    context.at = .{ .source = .zfi, .location = .{ .offset = zfi_data.len } };
    const t = parseZfiTrailingMeta(zfi_data) orelse return context.fail(.{ .@"index.name_footer_missing" = .{ .tight = true, .legacy = true } });
    return resolveNameBlob(t, records_end, zfi_data, context);
}

fn resolveNameBlob(t: ZfiTrailingMeta, records_end: usize, data: []const u8, context: *LoadDiagnostics) diagnostics.Error!ZfiNameLayout {
    const blob_start = std.math.sub(u64, t.payload_end, t.name_blob_len) catch return context.fail(.{ .@"index.name_blob_out_of_bounds" = .{ .name_blob_bytes = t.name_blob_len, .payload_end = t.payload_end, .records_end = records_end } });
    if (blob_start < records_end) return context.fail(.{ .@"index.name_blob_out_of_bounds" = .{ .name_blob_bytes = t.name_blob_len, .payload_end = t.payload_end, .records_end = records_end } });
    return .{ .bytes = data[@intCast(blob_start)..t.payload_end], .start = @intCast(blob_start), .source_mtime_ns = t.source_mtime_ns, .source_mtime_offset = t.source_mtime_offset };
}

fn zfiRecordsEnd(record_count: u32, context: *LoadDiagnostics) diagnostics.Error!usize {
    const records_bytes = std.math.mul(usize, @as(usize, record_count), @sizeOf(IndexRecord)) catch return context.overflow(.record_count, .multiply, record_count, @sizeOf(IndexRecord), std.math.maxInt(usize));
    return std.math.add(usize, @sizeOf(ZfiHeader), records_bytes) catch return context.overflow(.record_count, .add, @sizeOf(ZfiHeader), records_bytes, std.math.maxInt(usize));
}

fn buildNameMapFast(
    allocator: std.mem.Allocator,
    records: []const IndexRecord,
    name_data: []const u8,
    context: *LoadDiagnostics,
) diagnostics.Error!std.StringHashMapUnmanaged(u32) {
    if (records.len > std.math.maxInt(u32)) return context.capacity(.records, std.math.maxInt(u32), records.len);
    var name_map: std.StringHashMapUnmanaged(u32) = .empty;
    name_map.ensureTotalCapacity(allocator, @intCast(records.len)) catch return context.allocation(.name_map, null);
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
    context: *LoadDiagnostics,
) diagnostics.Error!?LoadedIndex {
    const zfi_file = std.Io.Dir.cwd().openFile(io, zfi_path, .{}) catch |err| switch (err) {
        error.FileNotFound => return null,
        else => return context.ioFailure(.open, err),
    };
    defer zfi_file.close(io);

    const zfi_stat = zfi_file.stat(io) catch |err| return context.ioFailure(.stat, err);

    // Check for 0-byte file before mmap (POSIX returns EINVAL)
    if (zfi_stat.size == 0) {
        return context.truncated(.header, ZFI_HEADER_BYTES, zfi_stat.size);
    }

    // Staleness: index file older than FASTA (both formats). `.fai` has no embedded
    // source size/mtime, so this plus "index not older" is its entire identity check.
    if (zfi_stat.mtime.nanoseconds < fasta_stat.mtime.nanoseconds) {
        return context.fail(.{ .@"index.older_than_source" = .{ .index_mtime_ns = zfi_stat.mtime.nanoseconds, .source_mtime_ns = fasta_stat.mtime.nanoseconds } });
    }

    if (mode == .matched) {
        if (zfi_stat.kind == .directory) return context.mapping(zfi_stat.size, null);
        return try loadZfiMatched(backing_allocator, io, zfi_file, zfi_stat.size, fasta_file, fasta_stat, mode.matched, context);
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

    // Empty catalogs are not useful and match the `.fai` reject-empty policy.
    context.at.location.offset = 4;
    if (header.record_count == 0) return context.fail(.{ .@"index.no_records" = .{ .format = .zfi } });

    // Checked records extent before forming the typed record slice.
    const records_end = try zfiRecordsEnd(header.record_count, context);
    if (zfi_data.len < records_end) return context.truncated(.records, records_end, zfi_data.len);

    const record_bytes = zfi_data[@sizeOf(ZfiHeader)..records_end];
    const records: []const IndexRecord = @as(
        [*]const IndexRecord,
        @ptrCast(@alignCast(record_bytes.ptr)),
    )[0..header.record_count];

    const name_layout = try resolveZfiNameLayout(zfi_data, records_end, records, context);
    const name_blob = name_layout.bytes;

    // Strong identity: production trailers store the FASTA mtime at index time.
    // Exact match catches same-size content replacement even when the `.zfi` file
    // mtime was touched forward. Legacy files without a trailer keep the weaker
    // index-mtime check above only.
    if (name_layout.source_mtime_ns) |stored_mtime| {
        context.at.location.offset = name_layout.source_mtime_offset;
        try checkSourceTime(stored_mtime, fasta_stat.mtime.nanoseconds, context);
    }

    // Side tables live after records and before the embedded name blob.
    const side_region_start = records_end;
    const side_region_end = name_layout.start;

    try validateZfiRecords(records, fasta_stat.size, zfi_data, name_blob, side_region_start, side_region_end, context);

    var arena = std.heap.ArenaAllocator.init(backing_allocator);
    errdefer arena.deinit();
    const allocator = arena.allocator();

    const build_name_map = switch (mode) {
        .lookup_full_map => true,
        .positional => |names| names.len > 1,
        .matched => true,
        .stats_scan => false,
    };

    const name_map = if (build_name_map)
        try buildNameMapFast(allocator, records, name_blob, context)
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

const ZfiReadWindow = struct {
    file: std.Io.File,
    io: std.Io,
    buffer: []u8,
    start: u64 = 0,
    len: usize = 0,

    fn take(self: *ZfiReadWindow, offset: u64, len: usize, context: *LoadDiagnostics) diagnostics.Error![]const u8 {
        if (len == 0) return self.buffer[0..0];
        if (offset >= self.start) {
            const relative = offset - self.start;
            if (relative <= self.len and len <= self.len - relative) {
                return self.buffer[@intCast(relative)..][0..len];
            }
        }
        self.len = std.Io.File.readPositionalAll(self.file, self.io, self.buffer, offset) catch |err| return context.readFailure(offset, self.buffer.len, err);
        self.start = offset;
        if (self.len < len) return context.shortRead(offset, len, self.len);
        return self.buffer[0..len];
    }
};

const UNMATCHED_RECORD_INDEX = std.math.maxInt(u32);

fn requestedNameMap(allocator: std.mem.Allocator, names: []const []const u8, reserve_multiplier: u32, context: *LoadDiagnostics) diagnostics.Error!std.StringHashMapUnmanaged(u32) {
    var map: std.StringHashMapUnmanaged(u32) = .empty;
    const count = std.math.cast(u32, names.len) orelse return context.capacity(.name_map, std.math.maxInt(u32), names.len);
    // Most catalog names miss this table. Leave spare slots to shorten those probes.
    const capacity = std.math.mul(u32, count, reserve_multiplier) catch return context.capacity(.name_map, std.math.maxInt(u32), @as(u128, count) * reserve_multiplier);
    map.ensureTotalCapacity(allocator, capacity) catch return context.allocation(.name_map, null);
    for (names) |name| {
        const entry = map.getOrPutAssumeCapacity(name);
        if (!entry.found_existing) entry.value_ptr.* = UNMATCHED_RECORD_INDEX;
    }
    return map;
}

fn discardUnmatchedNames(map: *std.StringHashMapUnmanaged(u32), storage: [][]const u8, matched_count: usize) void {
    // Unused name slots hold borrowed keys until iteration finishes.
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

    fn accept(self: *SideTableValidation, line: SideTableLine, rec: IndexRecord, fasta_len: u64, entry: u64, context: *LoadDiagnostics) diagnostics.Error!void {
        if (line.base_start != self.expected_base_start) return context.fail(.{ .@"index.side_base_start_mismatch" = .{ .entry = entry, .base_start = line.base_start, .expected_base_start = self.expected_base_start } });
        const geometry = diagnostics.Geometry{ .line_bases = line.line_bases, .line_bytes = line.line_bytes };
        if (geometryConstraint(geometry)) |constraint| return context.fail(.{ .@"index.side_line_geometry_invalid" = .{ .entry = entry, .geometry = geometry, .constraint = constraint } });
        const byte_end = std.math.add(u64, line.byte_offset, line.line_bytes) catch return context.overflow(.side_byte_end, .add, line.byte_offset, line.line_bytes, std.math.maxInt(u64));
        if (byte_end > fasta_len) return context.fail(.{ .@"index.sequence_span_out_of_bounds" = .{ .required_end = byte_end, .source_bytes = fasta_len, .bound_source = .actual, .entry = entry } });
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

inline fn checkSequenceStart(rec: IndexRecord, fasta_len: u64, context: *LoadDiagnostics) diagnostics.Error!void {
    if (rec.seq_len == 0) return context.fail(.{ .@"index.zero_sequence_length" = .{} });
    if (rec.seq_offset >= fasta_len) return context.fail(.{ .@"index.sequence_offset_out_of_bounds" = .{ .sequence_offset = rec.seq_offset, .source_bytes = fasta_len, .bound_source = .actual } });
}

inline fn checkUniformGeometry(rec: IndexRecord, fasta_len: u64, context: *LoadDiagnostics) diagnostics.Error!void {
    try checkSequenceStart(rec, fasta_len, context);
    const geometry = diagnostics.Geometry{ .line_bases = rec.line_bases, .line_bytes = rec.line_bytes };
    if (geometryConstraint(geometry)) |constraint| return context.fail(.{ .@"index.invalid_line_geometry" = .{ .geometry = geometry, .constraint = constraint } });
    const last_index: u128 = rec.seq_len - 1;
    const required_end = rec.seq_offset + (last_index / rec.line_bases) * rec.line_bytes + (last_index % rec.line_bases) + 1;
    if (required_end > fasta_len) return context.fail(.{ .@"index.sequence_span_out_of_bounds" = .{ .required_end = required_end, .source_bytes = fasta_len, .bound_source = .actual } });
}

fn checkNameExtent(rec: IndexRecord, blob_len: usize, context: *LoadDiagnostics) diagnostics.Error!void {
    if (!rangeFitsUsize(rec.name_offset, rec.name_len, blob_len)) return context.fail(.{ .@"index.name_out_of_bounds" = .{ .name_offset = rec.name_offset, .name_bytes = rec.name_len, .blob_bytes = blob_len } });
}

fn sideTableExtent(rec: IndexRecord, side_start: usize, side_end: usize, count: ?u64, context: *LoadDiagnostics) diagnostics.Error!usize {
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

fn loadZfiSideTable(
    allocator: std.mem.Allocator,
    window: *ZfiReadWindow,
    rec: IndexRecord,
    side_start: usize,
    side_end: usize,
    fasta_len: u64,
    retained: ?*std.ArrayList(u64),
    context: *LoadDiagnostics,
) diagnostics.Error!usize {
    const record_at = context.at;
    defer context.at = record_at;
    const count_end = try sideTableExtent(rec, side_start, side_end, null, context);
    const offset: usize = @intCast(rec.sideTableOffset());
    context.at.location.offset = offset;
    const count = std.mem.readInt(u64, (try window.take(offset, 8, context))[0..8], .little);
    const table_end = try sideTableExtent(rec, side_start, side_end, count, context);
    const table_bytes = table_end - count_end;
    const destination = if (retained) |words| blk: {
        const new_words = words.addManyAsSlice(allocator, 1 + table_bytes / @sizeOf(u64)) catch return context.allocation(.side_tables, table_bytes + 8);
        new_words[0] = count;
        break :blk std.mem.sliceAsBytes(new_words[1..]);
    } else null;
    var validation = SideTableValidation{};
    for (0..@intCast(count)) |i| {
        context.at.location.offset = count_end + i * @sizeOf(SideTableLine);
        const line_bytes = try window.take(context.at.location.offset.?, @sizeOf(SideTableLine), context);
        const line = std.mem.bytesToValue(SideTableLine, line_bytes);
        try validation.accept(line, rec, fasta_len, i + 1, context);
        if (destination) |bytes| @memcpy(bytes[i * @sizeOf(SideTableLine) ..][0..@sizeOf(SideTableLine)], line_bytes);
    }
    if (validation.expected_base_start != rec.seq_len) return context.fail(.{ .@"index.side_sequence_length_mismatch" = .{ .indexed_length = rec.seq_len, .side_table_bases = validation.expected_base_start } });
    return table_end;
}

fn checkZfiMagic(magic: [4]u8, context: *LoadDiagnostics) diagnostics.Error!void {
    context.at.location.offset = 0;
    if (!std.mem.eql(u8, magic[0..3], ZFI_MAGIC[0..3])) return context.fail(.{ .@"index.invalid_magic" = .{ .expected_magic = ZFI_MAGIC, .actual = context.failure.retain(&magic, 4) } });
    if (magic[3] != ZFI_MAGIC[3]) {
        context.at.location.offset = 3;
        return context.fail(.{ .@"index.unsupported_version" = .{ .found_version = magic[3], .supported_version = ZFI_MAGIC[3] } });
    }
}

fn checkSourceTime(stored: u64, actual: i128, context: *LoadDiagnostics) diagnostics.Error!void {
    if (actual < 0 or actual > std.math.maxInt(u64) or stored != @as(u64, @intCast(actual))) return context.fail(.{ .@"index.source_mtime_mismatch" = .{ .stored_mtime_ns = stored, .actual_mtime_ns = actual } });
}

// The mapped loader checks flags before trailers. Scan flags on this error path
// only, so matched success still reads each record once with bounded storage.
fn rejectMatchedTrailer(reader: *std.Io.File.Reader, count: u32, context: *LoadDiagnostics) diagnostics.Error {
    for (0..count) |i| {
        context.at = .{ .source = .zfi, .record = i + 1, .location = .{ .offset = ZFI_HEADER_BYTES + i * ZFI_INDEX_RECORD_BYTES } };
        const bytes = reader.interface.takeArray(ZFI_INDEX_RECORD_BYTES) catch |err| switch (err) {
            error.ReadFailed => return context.ioFailure(.read, reader.err orelse err),
            error.EndOfStream => return context.shortRead(reader.logicalPos(), ZFI_INDEX_RECORD_BYTES, reader.interface.buffered().len),
        };
        const rec = std.mem.bytesToValue(IndexRecord, bytes);
        if (!rec.nameInZfi()) {
            context.at.location.offset.? += 10;
            return context.fail(.{ .@"index.name_not_embedded" = .{ .flags = rec._pad[0], .required_flag = NAME_IN_ZFI_FLAG } });
        }
    }
    context.failure.finishDeferred() catch |err| return err;
    unreachable;
}

fn loadZfiMatched(backing_allocator: std.mem.Allocator, io: std.Io, zfi_file: std.Io.File, file_size: u64, fasta_file: std.Io.File, fasta_stat: std.Io.File.Stat, names: []const []const u8, context: *LoadDiagnostics) diagnostics.Error!LoadedIndex {
    const zfi_size = std.math.cast(usize, file_size) orelse return context.capacity(.map_bytes, std.math.maxInt(usize), file_size);
    if (zfi_size < ZFI_HEADER_BYTES) return context.truncated(.header, ZFI_HEADER_BYTES, zfi_size);
    var record_buffer: [64 * 1024]u8 = undefined;
    var reader = zfi_file.reader(io, &record_buffer);
    const header_bytes = reader.interface.takeArray(ZFI_HEADER_BYTES) catch |err| switch (err) {
        error.ReadFailed => return context.ioFailure(.read, reader.err orelse err),
        error.EndOfStream => return context.shortRead(0, ZFI_HEADER_BYTES, reader.interface.buffered().len),
    };
    const header = std.mem.bytesToValue(ZfiHeader, header_bytes);
    try checkZfiMagic(header.magic, context);
    context.at.location.offset = 8;
    if (header.source_size != fasta_stat.size) return context.fail(.{ .@"index.source_size_mismatch" = .{ .stored_bytes = header.source_size, .actual_bytes = fasta_stat.size } });
    context.at.location.offset = 4;
    if (header.record_count == 0) return context.fail(.{ .@"index.no_records" = .{ .format = .zfi } });
    const records_end = try zfiRecordsEnd(header.record_count, context);
    if (records_end > zfi_size) return context.truncated(.records, records_end, zfi_size);
    var tail: [ZFI_NAME_FOOTER_LEGACY_BYTES + ZFI_SOURCE_ID_BYTES]u8 = undefined;
    const tail_len = @min(tail.len, zfi_size);
    const tail_start = zfi_size - tail_len;
    const buffered_body = if (zfi_size <= record_buffer.len)
        reader.interface.peek(zfi_size - ZFI_HEADER_BYTES) catch |err| switch (err) {
            error.ReadFailed => return context.ioFailure(.read, reader.err orelse err),
            error.EndOfStream => return context.shortRead(ZFI_HEADER_BYTES, zfi_size - ZFI_HEADER_BYTES, reader.interface.buffered().len),
        }
    else
        null;
    const tail_bytes = if (buffered_body) |body| body[body.len - tail_len ..] else tail_read: {
        const got = std.Io.File.readPositionalAll(zfi_file, io, tail[0..tail_len], tail_start) catch |err| return context.readFailure(tail_start, tail_len, err);
        if (got != tail_len) return context.shortRead(tail_start, tail_len, got);
        break :tail_read tail[0..tail_len];
    };
    context.at.location.offset = zfi_size;
    var trailing = parseZfiTrailingMeta(tail_bytes) orelse {
        context.deferFailure(.{ .@"index.name_footer_missing" = .{ .tight = true, .legacy = true } });
        return rejectMatchedTrailer(&reader, header.record_count, context);
    };
    trailing.payload_end += tail_start;
    const blob_len = std.math.cast(usize, trailing.name_blob_len) orelse {
        context.deferFailure(.{ .@"index.name_blob_out_of_bounds" = .{ .name_blob_bytes = trailing.name_blob_len, .payload_end = trailing.payload_end, .records_end = records_end } });
        return rejectMatchedTrailer(&reader, header.record_count, context);
    };
    const blob_start = std.math.sub(usize, trailing.payload_end, blob_len) catch {
        context.deferFailure(.{ .@"index.name_blob_out_of_bounds" = .{ .name_blob_bytes = trailing.name_blob_len, .payload_end = trailing.payload_end, .records_end = records_end } });
        return rejectMatchedTrailer(&reader, header.record_count, context);
    };
    if (blob_start < records_end) {
        context.deferFailure(.{ .@"index.name_blob_out_of_bounds" = .{ .name_blob_bytes = trailing.name_blob_len, .payload_end = trailing.payload_end, .records_end = records_end } });
        return rejectMatchedTrailer(&reader, header.record_count, context);
    }
    const stale_identity = if (trailing.source_mtime_ns) |stored| stale: {
        context.deferred = true;
        defer context.deferred = false;
        context.at.location.offset = tail_start + trailing.source_mtime_offset.?;
        checkSourceTime(stored, fasta_stat.mtime.nanoseconds, context) catch break :stale true;
        break :stale false;
    } else false;
    var arena = std.heap.ArenaAllocator.init(backing_allocator);
    errdefer arena.deinit();
    const allocator = arena.allocator();
    const single_name = if (names.len == 1) names[0] else null;
    var requested: std.StringHashMapUnmanaged(u32) = if (single_name == null) try requestedNameMap(allocator, names, 4, context) else .empty;
    const requested_count = if (single_name != null) 1 else requested.count();
    const records = allocator.alloc(IndexRecord, requested_count) catch return context.allocation(.records, @as(u64, requested_count) * @sizeOf(IndexRecord));
    const name_slices = allocator.alloc([]const u8, requested_count) catch return context.allocation(.name_slices, @as(u64, requested_count) * @sizeOf([]const u8));
    var side_words: std.ArrayList(u64) = .empty;
    var matched_count: usize = 0;
    var corrupt_metadata = false;
    var previous_side_end = records_end;
    var name_buffer: [64 * 1024]u8 = undefined;
    var side_buffer: [64 * 1024]u8 = undefined;
    var name_window = ZfiReadWindow{ .file = zfi_file, .io = io, .buffer = buffered_body orelse &name_buffer, .start = if (buffered_body != null) @sizeOf(ZfiHeader) else 0, .len = if (buffered_body) |body| body.len else 0 };
    var side_window = ZfiReadWindow{ .file = zfi_file, .io = io, .buffer = buffered_body orelse &side_buffer, .start = if (buffered_body != null) @sizeOf(ZfiHeader) else 0, .len = if (buffered_body) |body| body.len else 0 };
    for (0..header.record_count) |i| {
        context.at = .{ .source = .zfi, .record = i + 1, .location = .{ .offset = ZFI_HEADER_BYTES + i * ZFI_INDEX_RECORD_BYTES } };
        const bytes = reader.interface.takeArray(ZFI_INDEX_RECORD_BYTES) catch |err| switch (err) {
            error.ReadFailed => return context.ioFailure(.read, reader.err orelse err),
            error.EndOfStream => return context.shortRead(reader.logicalPos(), ZFI_INDEX_RECORD_BYTES, reader.interface.buffered().len),
        };
        var rec = std.mem.bytesToValue(IndexRecord, bytes);
        if (!rec.nameInZfi()) {
            context.at.location.offset.? += 10;
            return context.fail(.{ .@"index.name_not_embedded" = .{ .flags = rec._pad[0], .required_flag = NAME_IN_ZFI_FLAG } });
        }
        if (stale_identity or corrupt_metadata) continue;
        context.deferred = true;
        checkNameExtent(rec, blob_len, context) catch {
            corrupt_metadata = true;
            context.deferred = false;
            continue;
        };
        const name = try name_window.take(blob_start + rec.name_offset, rec.name_len, context);
        const matched = if (single_name == null) requested.getEntry(name) else null;
        const retain = if (single_name) |wanted| matched_count == 0 and std.mem.eql(u8, name, wanted) else if (matched) |entry| entry.value_ptr.* == UNMATCHED_RECORD_INDEX else false;
        if (rec.isUniformWidth()) {
            checkUniformGeometry(rec, fasta_stat.size, context) catch {
                corrupt_metadata = true;
                context.deferred = false;
                continue;
            };
        } else {
            checkSequenceStart(rec, fasta_stat.size, context) catch {
                corrupt_metadata = true;
                context.deferred = false;
                continue;
            };
            const wire_offset = rec.sideTableOffset();
            const retained_offset = side_words.items.len * @sizeOf(u64);
            const end = loadZfiSideTable(allocator, &side_window, rec, records_end, blob_start, fasta_stat.size, if (retain and wire_offset >= previous_side_end) &side_words else null, context) catch |err| switch (err) {
                error.InvalidIndex => {
                    corrupt_metadata = true;
                    context.deferred = false;
                    continue;
                },
                else => return err,
            };
            if (wire_offset < previous_side_end) {
                context.at.location.offset = wire_offset;
                context.deferFailure(.{ .@"index.side_table_order_invalid" = .{ .table_start = wire_offset, .minimum_start = previous_side_end } });
                corrupt_metadata = true;
                context.deferred = false;
                continue;
            }
            previous_side_end = end;
            context.deferred = false;
            if (retain) rec.markNonUniform(retained_offset) catch return context.fail(.{ .@"index.side_table_offset_limit_exceeded" = .{ .offset = retained_offset, .max_offset = MAX_SIDE_TABLE_OFFSET, .offset_role = .retained } });
        }
        context.deferred = false;
        if (!retain) continue;
        name_slices[matched_count] = allocator.dupe(u8, name) catch return context.allocation(.retained_names, name.len);
        records[matched_count] = rec;
        if (matched) |entry| {
            entry.key_ptr.* = name_slices[matched_count];
            entry.value_ptr.* = @intCast(matched_count);
        }
        matched_count += 1;
    }
    try context.failure.finishDeferred();

    if (single_name == null) discardUnmatchedNames(&requested, name_slices, matched_count);
    const side_data = std.mem.sliceAsBytes(side_words.items);
    return .{
        .io = io,
        .fasta_file = fasta_file,
        .records = records[0..matched_count],
        .name_map = if (single_name != null) null else requested,
        .name_slices = name_slices[0..matched_count],
        .fasta_size = fasta_stat.size,
        .zfi_data = side_data,
        .zfi_side_end = side_data.len,
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
    context: *LoadDiagnostics,
) diagnostics.Error!?LoadedIndex {
    const fai_file = std.Io.Dir.cwd().openFile(io, fai_path, .{}) catch |err| switch (err) {
        error.FileNotFound => return null,
        else => return context.ioFailure(.open, err),
    };
    defer fai_file.close(io);

    const fai_stat = fai_file.stat(io) catch |err| return context.ioFailure(.stat, err);

    // `.fai` identity is mtime-only: the text format has no source size or mtime field.
    // Same-size FASTA replacement is detected only when the FASTA mtime moves forward.
    if (fai_stat.mtime.nanoseconds < fasta_stat.mtime.nanoseconds) {
        return context.fail(.{ .@"index.older_than_source" = .{ .index_mtime_ns = fai_stat.mtime.nanoseconds, .source_mtime_ns = fasta_stat.mtime.nanoseconds } });
    }

    if (fai_stat.size == 0) {
        return context.fail(.{ .@"index.no_records" = .{ .format = .fai } });
    }
    if (mode == .matched and fai_stat.kind == .directory) return context.mapping(fai_stat.size, null);

    switch (mode) {
        .stats_scan, .positional, .matched => return @as(?LoadedIndex, try loadFaiStreamed(backing_allocator, io, fai_file, fasta_file, fasta_stat, mode, context)),
        .lookup_full_map => {},
    }

    const fai_size = std.math.cast(usize, fai_stat.size) orelse return context.capacity(.map_bytes, std.math.maxInt(usize), fai_stat.size);
    var fai_map = platform.mapFileReadOnly(io, fai_file, fai_size) catch |err| return context.mapping(fai_stat.size, err);

    errdefer fai_map.destroy(io);

    var arena = std.heap.ArenaAllocator.init(backing_allocator);
    errdefer arena.deinit();
    const allocator = arena.allocator();

    const fai_data: platform.MappedBytes = fai_map.memory;

    // Newline count undercounts by one when the final `.fai` line has no trailing `\n`.
    var approx_records = std.mem.count(u8, fai_data, "\n");
    if (fai_data.len > 0 and fai_data[fai_data.len - 1] != '\n') {
        approx_records = std.math.add(usize, approx_records, 1) catch return context.overflow(.record_count, .add, approx_records, 1, std.math.maxInt(usize));
    }
    // Reject sizes that would overflow the records allocation.
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
        // Advance past `\n` only when present; a final line may omit it.
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
    const name_map = try buildNameMapFast(allocator, loaded_records, fai_data, context);

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

const FaiNameRef = struct {
    offset: u64,
    len: usize,
    at: diagnostics.Reference,
};

// Keep one scan body for positional, stats, and file requests.
noinline fn loadFaiStreamed(
    backing_allocator: std.mem.Allocator,
    io: std.Io,
    fai_file: std.Io.File,
    fasta_file: std.Io.File,
    fasta_stat: std.Io.File.Stat,
    mode: LoadMode,
    context: *LoadDiagnostics,
) diagnostics.Error!LoadedIndex {
    var arena = std.heap.ArenaAllocator.init(backing_allocator);
    errdefer arena.deinit();
    const allocator = arena.allocator();

    var records_list: std.ArrayList(IndexRecord) = .empty;
    var slices_list: std.ArrayList([]const u8) = .empty;
    const single_name = if (mode == .matched and mode.matched.len == 1) mode.matched[0] else null;
    var requested: std.StringHashMapUnmanaged(u32) = .empty;
    switch (mode) {
        .positional, .matched => |names| {
            if (single_name == null) requested = try requestedNameMap(allocator, names, if (mode == .matched) 4 else 1, context);
            const count = if (single_name != null) 1 else requested.count();
            records_list.ensureTotalCapacityPrecise(allocator, count) catch return context.allocation(.records, @as(u64, count) * @sizeOf(IndexRecord));
            slices_list.ensureTotalCapacityPrecise(allocator, count) catch return context.allocation(.name_slices, @as(u64, count) * @sizeOf([]const u8));
        },
        .stats_scan => {},
        .lookup_full_map => unreachable,
    }

    var io_buf: [FAI_READER_BUFFER_BYTES]u8 = undefined;
    var file_reader = fai_file.reader(io, &io_buf);
    var saw_record = false;
    var long_line_buffer: [FAI_READER_BUFFER_BYTES]u8 = undefined;
    var shortest_ref: FaiNameRef = undefined;
    var longest_ref: FaiNameRef = undefined;
    var shortest_index: usize = 0;
    var longest_index: usize = 0;
    var shortest_len: u64 = 0;
    var longest_len: u64 = 0;

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
        const maybe_line = if (std.mem.findScalar(u8, available, '\n')) |end| buffered: {
            file_reader.interface.toss(end + 1);
            break :buffered available[0..end];
        } else file_reader.interface.takeDelimiter('\n') catch |err| switch (err) {
            error.ReadFailed => return context.ioFailure(.read, file_reader.err orelse err),
            error.StreamTooLong => if (mode == .matched) long: {
                context.at.record = try context.nextOrdinal(record_number, .records);
                const parsed = try takeLongFaiLine(&file_reader, &long_line_buffer, context);
                long_fields = parsed.fields;
                break :long parsed.name;
            } else return context.lineBuffer(file_reader.interface.buffered().len),
        };
        const line = maybe_line orelse break;
        row_number += 1;
        if (line.len == 0 and long_fields == null) continue;
        record_number = try context.nextOrdinal(record_number, .records);
        context.at.record = record_number;

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

        saw_record = true;
        const matched = if (mode != .stats_scan and single_name == null) requested.getEntry(line[0..fields.name_end]) else null;
        const retain = switch (mode) {
            .stats_scan => true,
            .positional, .matched => if (single_name) |wanted|
                records_list.items.len == 0 and std.mem.eql(u8, line[0..fields.name_end], wanted)
            else if (matched) |entry|
                entry.value_ptr.* == UNMATCHED_RECORD_INDEX
            else
                false,
            .lookup_full_map => unreachable,
        };
        if (!retain) continue;

        const rec_index = records_list.items.len;
        switch (mode) {
            .positional, .matched => {
                const name = allocator.dupe(u8, line[0..fields.name_end]) catch return context.allocation(.retained_names, fields.name_end);
                slices_list.append(allocator, name) catch return context.allocation(.name_slices, null);
                if (matched) |entry| {
                    entry.key_ptr.* = name;
                    entry.value_ptr.* = @intCast(rec_index);
                }
            },
            .stats_scan => {
                const name_ref = FaiNameRef{ .offset = current_offset, .len = fields.name_end, .at = context.at };
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
            .lookup_full_map => unreachable,
        }
        records_list.append(allocator, rec) catch return context.allocation(.records, null);
    }

    if (!saw_record) return context.fail(.{ .@"index.no_records" = .{ .format = .fai } });

    if (mode != .stats_scan and single_name == null) discardUnmatchedNames(&requested, slices_list.allocatedSlice(), records_list.items.len);

    const owned_records = records_list.toOwnedSlice(allocator) catch return context.allocation(.records, null);
    const owned_slices = switch (mode) {
        .positional, .matched => slices_list.toOwnedSlice(allocator) catch return context.allocation(.name_slices, null),
        .stats_scan => blk: {
            const names = allocator.alloc([]const u8, 2) catch return context.allocation(.name_slices, 2 * @sizeOf([]const u8));
            names[0] = try readFaiName(allocator, io, fai_file, shortest_ref, context);
            names[1] = if (shortest_index == longest_index)
                names[0]
            else
                try readFaiName(allocator, io, fai_file, longest_ref, context);
            break :blk names;
        },
        .lookup_full_map => unreachable,
    };

    const name_map: ?std.StringHashMapUnmanaged(u32) = switch (mode) {
        .positional, .matched => if (single_name != null) null else requested,
        .stats_scan => null,
        .lookup_full_map => unreachable,
    };

    return LoadedIndex{
        .io = io,
        .fasta_file = fasta_file,
        .records = owned_records,
        .name_map = name_map,
        .name_slices = owned_slices,
        .stats_name_indices = switch (mode) {
            .positional, .matched => null,
            .stats_scan => .{ shortest_index, longest_index },
            .lookup_full_map => unreachable,
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
    context: *LoadDiagnostics,
) diagnostics.Error![]const u8 {
    var name_context = context.*;
    name_context.at = name_ref.at;
    const name = allocator.alloc(u8, name_ref.len) catch return name_context.allocation(.retained_names, name_ref.len);
    if (name.len == 0) return name;
    const got = std.Io.File.readPositionalAll(fai_file, io, name, name_ref.offset) catch |err| return name_context.readFailure(name_ref.offset, name.len, err);
    if (got != name.len) return name_context.shortRead(name_ref.offset, name.len, got);
    return name;
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

fn faiFieldContext(context: *LoadDiagnostics, start: u64) LoadDiagnostics {
    var field_context = context.*;
    field_context.at.location.offset = context.at.location.offset.? + start;
    field_context.at.location.column = start + 1;
    return field_context;
}

fn missingFaiField(context: *LoadDiagnostics, field: diagnostics.FaiField) diagnostics.Error {
    return context.fail(.{ .@"index.fai.field_missing" = .{ .field = field, .column_number = @intCast(@intFromEnum(field) + 1) } });
}

fn checkFaiInteger(value: u64, lexical: bool, overflow: bool, bytes: u64, prefix: []const u8, field: diagnostics.FaiField, max: u64, context: *LoadDiagnostics) diagnostics.Error!u64 {
    if (bytes == 0) return missingFaiField(context, field);
    if (lexical) return context.fail(.{ .@"index.fai.invalid_integer" = .{ .field = field, .value = context.failure.retainPrefix(prefix, .{ .exact = bytes }, 32) } });
    if (overflow or value > max) return context.fail(.{ .@"index.fai.integer_out_of_range" = .{ .field = field, .max_value = max, .value = context.failure.retainPrefix(prefix, .{ .exact = bytes }, 32) } });
    return value;
}

// Keep successful decimal parsing independent of diagnostics. Only rejected
// fields need a complete lexical scan, value preview, and original location.
fn parseFaiAsciiU64(text: []const u8) ?u64 {
    if (text.len == 0) return null;
    var value: u64 = 0;
    for (text) |byte| {
        if (byte < '0' or byte > '9') return null;
        value = appendFaiDigit(value, byte) orelse return null;
    }
    return value;
}

fn faiIntegerFailure(text: []const u8, field: diagnostics.FaiField, max: u64, context: *LoadDiagnostics) diagnostics.Error {
    @branchHint(.cold);
    if (text.len == 0) return missingFaiField(context, field);
    const preview = context.failure.retain(text, 32);
    for (text) |byte| {
        if (byte < '0' or byte > '9') return context.fail(.{ .@"index.fai.invalid_integer" = .{ .field = field, .value = preview } });
    }
    return context.fail(.{ .@"index.fai.integer_out_of_range" = .{ .field = field, .max_value = max, .value = preview } });
}

fn parseFaiIndexLine(line: []const u8, context: *LoadDiagnostics) diagnostics.Error!FaiLineFields {
    const name_end = std.mem.findScalar(u8, line, '\t') orelse {
        var field_context = faiFieldContext(context, line.len);
        return missingFaiField(&field_context, .length);
    };
    const name_len = std.math.cast(u16, name_end) orelse return context.fail(.{ .@"index.identifier_too_long" = .{ .observed_bytes = .{ .exact = name_end }, .limit_bytes = std.math.maxInt(u16), .format = .fai } });
    var start = name_end;
    var values: [4]u64 = undefined;
    inline for (FAI_NUMERIC_FIELDS, 0..) |field, i| {
        if (start < line.len and line[start] == '\t') start += 1;
        const end = std.mem.findScalarPos(u8, line, start, '\t') orelse line.len;
        const text = line[start..end];
        const max = if (i < 2) std.math.maxInt(u64) else std.math.maxInt(u32);
        values[i] = parseFaiAsciiU64(text) orelse {
            var field_context = faiFieldContext(context, start);
            return faiIntegerFailure(text, field, max, &field_context);
        };
        if (values[i] > max) {
            var field_context = faiFieldContext(context, start);
            return faiIntegerFailure(text, field, max, &field_context);
        }
        start = end;
    }
    return .{ .name_end = name_end, .name_len = name_len, .seq_len = values[0], .seq_offset = values[1], .line_bases = @intCast(values[2]), .line_bytes = @intCast(values[3]) };
}

// Parse long rows from their original bytes with bounded field previews. Names
// remain in the caller's fixed buffer; numeric fields are not reconstructed.
fn takeLongFaiLine(file_reader: *std.Io.File.Reader, buffer: []u8, context: *LoadDiagnostics) diagnostics.Error!struct { name: []const u8, fields: FaiLineFields } {
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

fn validateZfiRecords(records: []const IndexRecord, fasta_len: u64, data: platform.MappedBytes, blob: []const u8, side_start: usize, side_end: usize, context: *LoadDiagnostics) diagnostics.Error!void {
    var previous_end = side_start;
    for (records, 0..) |rec, record_index| {
        context.at = .{ .source = .zfi, .record = record_index + 1, .location = .{ .offset = ZFI_HEADER_BYTES + record_index * ZFI_INDEX_RECORD_BYTES } };
        try checkNameExtent(rec, blob.len, context);
        if (rec.isUniformWidth()) {
            try checkUniformGeometry(rec, fasta_len, context);
            continue;
        }
        try checkSequenceStart(rec, fasta_len, context);
        const count_end = try sideTableExtent(rec, side_start, side_end, null, context);
        const offset: usize = @intCast(rec.sideTableOffset());
        context.at.location.offset = offset;
        const count = std.mem.readInt(u64, data[offset..][0..8], .little);
        const table_end = try sideTableExtent(rec, side_start, side_end, count, context);
        var validation = SideTableValidation{};
        for (0..@intCast(count)) |i| {
            context.at.location.offset = count_end + i * @sizeOf(SideTableLine);
            const line = std.mem.bytesToValue(SideTableLine, data[@intCast(context.at.location.offset.?)..][0..@sizeOf(SideTableLine)]);
            try validation.accept(line, rec, fasta_len, i + 1, context);
        }
        if (validation.expected_base_start != rec.seq_len) return context.fail(.{ .@"index.side_sequence_length_mismatch" = .{ .indexed_length = rec.seq_len, .side_table_bases = validation.expected_base_start } });
        context.at.location.offset = offset;
        if (offset < previous_end) return context.fail(.{ .@"index.side_table_order_invalid" = .{ .table_start = offset, .minimum_start = previous_end } });
        previous_end = table_end;
    }
}
