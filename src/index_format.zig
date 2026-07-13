const std = @import("std");
const posix = std.posix;

// ============================================================================
// Types - shared by indexer, getter, stats
// ============================================================================

/// ZFI binary format magic + header (`ZFI\x01` since first release; no users on a prior on-disk layout).
pub const ZFI_MAGIC: [4]u8 = .{ 'Z', 'F', 'I', 0x01 };

const NON_UNIFORM_WIDTH_FLAG: u8 = 1;
const NAME_IN_ZFI_FLAG: u8 = 2;
pub const name_in_zfi_flag = NAME_IN_ZFI_FLAG;
const SIDE_TABLE_OFFSET_BYTES = 5;
const MAX_SIDE_TABLE_OFFSET = (1 << (SIDE_TABLE_OFFSET_BYTES * 8)) - 1;

pub const ZFI_NAME_FOOTER_MAGIC: [4]u8 = .{ 'Z', 'F', 'N', 'M' };

/// On-disk name-table footer size (4-byte magic + little-endian u64 length).
pub const zfi_name_footer_bytes: usize = 12;

pub fn encodeZfiNameFooter(name_blob_len: u64) [zfi_name_footer_bytes]u8 {
    var out: [zfi_name_footer_bytes]u8 = undefined;
    @memcpy(out[0..4], &ZFI_NAME_FOOTER_MAGIC);
    std.mem.writeInt(u64, out[4..12], name_blob_len, .little);
    return out;
}

fn decodeZfiNameFooterBytes(footer_bytes: []const u8) ?u64 {
    if (footer_bytes.len < zfi_name_footer_bytes) return null;
    if (!std.mem.eql(u8, footer_bytes[0..4], &ZFI_NAME_FOOTER_MAGIC)) return null;
    return std.mem.readInt(u64, footer_bytes[4..12], .little);
}

/// Legacy on-disk footer written before tight layout (16 bytes: magic + 4 pad + u64).
const zfi_name_footer_legacy_bytes: usize = 16;

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

/// Index record for ZFI output (40 bytes padded).
///
/// For `.zfi` indexes, `name_offset` / `name_len` point into the mmap'd FASTA (byte after `>`)
/// unless `nameInZfi()` is set, in which case they point into the `.zfi` name blob.
/// For `.fai` fallback loads, `name_offset` / `name_len` point into the mmap'd `.fai`
/// line (name field before the first tab). Do not call `getName` with `fasta_data` on
/// `.fai` records; pass `LoadedIndex.fai_data` or use `getRecordName`.
///
/// v0.3 stores non-uniform-width metadata in `_pad` so uniform records remain
/// byte-identical to v0.2 records. `_pad[0] & 1 == 0` means the classic O(1)
/// line formula applies. Non-uniform records store a 40-bit little-endian
/// absolute `.zfi` side-table offset in `_pad[1..6]`.
pub const IndexRecord = extern struct {
    name_offset: u64,
    name_len: u16,
    _pad: [6]u8 = .{0} ** 6,
    seq_offset: u64,
    seq_len: u64,
    line_bases: u32,
    line_bytes: u32,

    pub fn getName(self: IndexRecord, data: []const u8) []const u8 {
        // `.zfi`: slice into FASTA (byte after `>`). `.fai`: slice into `.fai` mmap line.
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

/// Result of loading an index (from .zfi or .fai)
pub const LoadMode = enum {
    records_only,
    lookup_full_map,
    /// Stats composition: `.fai` rows streamed without retaining every name in RAM.
    stats_scan,
};

/// When a catalog exceeds this record count, multi-region GET (2..15) uses a by-record scan
/// instead of building a name hash map at query time.
pub const large_catalog_record_threshold: u32 = 4096;

/// Loaded FASTA + index state (from `.zfi` or samtools-compatible `.fai`).
///
/// `source` tells which index format was loaded. Names resolve through embedded
/// `name_offset` / `name_len` (`.zfi` → FASTA mmap, `.fai` → `.fai` mmap),
/// `name_slices` (arena-owned copies when a name hash map is built), or
/// `name_map` when loaded with `.lookup_full_map`.
pub const LoadedIndex = struct {
    records: []const IndexRecord,
    name_map: std.StringHashMap(usize),
    has_name_map: bool,
    /// Arena-owned name per record index; populated when `has_name_map`.
    name_slices: []const []const u8 = &.{},
    /// Populated when `.zfi` stores an embedded name blob (see `ZfiNameFooter`).
    name_blob: ?[]const u8 = null,
    fai_data: ?[]align(4096) const u8 = null,
    fasta_data: []align(4096) const u8,
    fasta_size: u64,
    zfi_data: ?[]align(4096) const u8,
    source: IndexSource,
    /// Optional `.fai` path for on-demand name reads (`stats_scan` without name table).
    sidecar_path: ?[]const u8 = null,
    /// Byte offsets into `sidecar_path` for each record line (`stats_scan` only).
    fai_line_offsets: []const u64 = &.{},
    arena: std.heap.ArenaAllocator,

    pub const IndexSource = enum { zfi, fai };

    fn nameBase(self: *const LoadedIndex) []const u8 {
        return if (self.source == .zfi) self.fasta_data else self.fai_data.?;
    }

    fn recordNameData(self: *const LoadedIndex, rec: IndexRecord) []const u8 {
        if (rec.nameInZfi()) return self.name_blob orelse self.zfi_data.?;
        return self.nameBase();
    }

    pub fn getRecordName(self: *const LoadedIndex, rec_idx: usize) []const u8 {
        if (self.name_slices.len > 0) {
            return self.name_slices[rec_idx];
        }
        const rec = self.records[rec_idx];
        if (rec.name_len > 0) {
            return rec.getName(self.recordNameData(rec));
        }
        var it = self.name_map.iterator();
        while (it.next()) |entry| {
            if (entry.value_ptr.* == rec_idx) {
                return entry.key_ptr.*;
            }
        }
        return "?";
    }

    pub fn getRecordNameWithIo(self: *LoadedIndex, io: std.Io, rec_idx: usize) []const u8 {
        if (self.sidecar_path) |path| {
            if (self.fai_line_offsets.len > rec_idx) {
                return readFaiNameAtOffset(io, self.arena.allocator(), path, self.fai_line_offsets[rec_idx]) catch "?";
            }
            return readFaiNameLine(io, self.arena.allocator(), path, rec_idx) catch "?";
        }
        return self.getRecordName(rec_idx);
    }

    pub fn deinit(self: *LoadedIndex) void {
        posix.munmap(@alignCast(@constCast(self.fasta_data)));
        if (self.zfi_data) |zd| {
            posix.munmap(@alignCast(@constCast(zd)));
        }
        if (self.fai_data) |fd| {
            posix.munmap(@alignCast(@constCast(fd)));
        }
        // Name-map keys live in `arena` when `name_slices` is populated.
        // Embedded `.zfi` names keep keys in the mmap'd index; skip key frees.
        if (self.name_slices.len == 0 and self.source == .zfi and self.name_blob == null) {
            self.name_map.deinit();
        }
        self.arena.deinit();
    }

    pub fn lookupName(self: *const LoadedIndex, name: []const u8) ?usize {
        if (self.has_name_map) {
            return self.name_map.get(name);
        }

        if (self.name_slices.len > 0) {
            var found: ?usize = null;
            for (self.name_slices, 0..) |rec_name, i| {
                if (std.mem.eql(u8, rec_name, name)) {
                    found = i;
                }
            }
            return found;
        }

        var found: ?usize = null;
        for (self.records, 0..) |rec, i| {
            if (rec.name_len == 0) continue;
            if (std.mem.eql(u8, rec.getName(self.recordNameData(rec)), name)) {
                found = i;
            }
        }
        return found;
    }

    pub fn sideTableLines(self: *const LoadedIndex, rec: IndexRecord) []const SideTableLine {
        if (rec.isUniformWidth()) return &.{};

        const zfi = self.zfi_data.?;
        const offset: usize = @intCast(rec.sideTableOffset());
        const count: *const u64 = @ptrCast(@alignCast(zfi[offset..].ptr));
        const line_bytes = zfi[offset + @sizeOf(u64) ..];
        return @as(
            [*]const SideTableLine,
            @ptrCast(@alignCast(line_bytes.ptr)),
        )[0..count.*];
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

const LoadAttempt = union(enum) {
    loaded: LoadedIndex,
    not_found,
    stale,
    corrupt,
};

// ============================================================================
// Error helper (shared by main, getter, indexer)
// ============================================================================

pub fn printErrorAndExit(comptime fmt: []const u8, args: anytype) noreturn {
    std.debug.print(fmt, args);
    std.process.exit(1);
}

// ============================================================================
// .zfi writer (used by indexer)
// ============================================================================

/// Writes the .zfi binary index file.
pub fn writeZfi(
    io: std.Io,
    path: []const u8,
    records: []const IndexRecord,
    source_size: u64,
) !void {
    const file = try std.Io.Dir.cwd().createFile(io, path, .{});
    errdefer std.Io.Dir.cwd().deleteFile(io, path) catch {};
    defer file.close(io);

    const header = ZfiHeader{
        .magic = ZFI_MAGIC,
        .record_count = @intCast(records.len),
        .source_size = source_size,
    };
    try std.Io.File.writeStreamingAll(file, io, std.mem.asBytes(&header));
    try std.Io.File.writeStreamingAll(file, io, std.mem.sliceAsBytes(records));
}

// ============================================================================
// Shared index loader - loads .zfi or falls back to .fai
// ============================================================================

/// Load the index for a FASTA file. Tries .zfi first, then .fai fallback.
/// The caller must call deinit() on the returned LoadedIndex.
pub fn loadIndex(io: std.Io, fasta_path: []const u8) LoadedIndex {
    return loadIndexWithMode(io, fasta_path, .lookup_full_map);
}

pub fn loadIndexWithMode(io: std.Io, fasta_path: []const u8, mode: LoadMode) LoadedIndex {
    return loadIndexCheckedWithMode(io, fasta_path, mode) catch |err| switch (err) {
        error.FileNotFound => printErrorAndExit("error: file not found: {s}\n", .{fasta_path}),
        error.AccessDenied => printErrorAndExit("error: access denied: {s}\n", .{fasta_path}),
        error.EmptyFile => printErrorAndExit("error: file is empty: {s}\n", .{fasta_path}),
        error.PathTooLong => printErrorAndExit("error: path too long\n", .{}),
        error.MmapFailed => printErrorAndExit("error: failed to mmap file: {s}\n", .{fasta_path}),
        error.CorruptIndex => printErrorAndExit("error: corrupt index file for: {s}\n", .{fasta_path}),
        error.StaleIndex => printErrorAndExit("error: index is stale (FASTA is newer than index). Re-run 'z-fasta index <path>'.\n", .{}),
        error.NoIndexFound => printErrorAndExit("error: no index found for {s}. Run 'z-fasta index {s}' first.\n", .{ fasta_path, fasta_path }),
        error.OutOfMemory => printErrorAndExit("error: out of memory loading index\n", .{}),
        error.Io => printErrorAndExit("error: failed to load index for: {s}\n", .{fasta_path}),
    };
}

/// Load the index for a FASTA file with typed errors for testing and fallback control.
pub fn loadIndexChecked(io: std.Io, fasta_path: []const u8) LoadIndexError!LoadedIndex {
    return loadIndexCheckedWithMode(io, fasta_path, .lookup_full_map);
}

pub fn loadIndexCheckedWithMode(io: std.Io, fasta_path: []const u8, mode: LoadMode) LoadIndexError!LoadedIndex {
    // Open FASTA file
    const fasta_file = std.Io.Dir.cwd().openFile(io, fasta_path, .{}) catch |err| switch (err) {
        error.FileNotFound => return error.FileNotFound,
        error.AccessDenied => return error.AccessDenied,
        else => return error.Io,
    };
    defer fasta_file.close(io);

    const fasta_stat = fasta_file.stat(io) catch return error.Io;

    if (fasta_stat.size == 0) {
        return error.EmptyFile;
    }

    // mmap FASTA
    const fasta_data = posix.mmap(
        null,
        fasta_stat.size,
        .{ .READ = true },
        .{ .TYPE = .PRIVATE },
        fasta_file.handle,
        0,
    ) catch return error.MmapFailed;

    // Try .zfi first
    var zfi_path_buf: [4096]u8 = undefined;
    const zfi_path = std.fmt.bufPrint(&zfi_path_buf, "{s}.zfi", .{fasta_path}) catch return error.PathTooLong;

    var zfi_failure: ?LoadIndexError = null;
    switch (tryLoadZfi(io, zfi_path, fasta_data, fasta_stat, mode) catch |err| {
        posix.munmap(@alignCast(@constCast(fasta_data)));
        return err;
    }) {
        .loaded => |result| return result,
        .stale => zfi_failure = error.StaleIndex,
        .corrupt => zfi_failure = error.CorruptIndex,
        .not_found => {},
    }

    // Try .fai fallback
    var fai_path_buf: [4096]u8 = undefined;
    const fai_path = std.fmt.bufPrint(&fai_path_buf, "{s}.fai", .{fasta_path}) catch {
        posix.munmap(@alignCast(@constCast(fasta_data)));
        return error.PathTooLong;
    };

    switch (tryLoadFai(io, fai_path, fasta_data, fasta_stat, mode) catch |err| {
        posix.munmap(@alignCast(@constCast(fasta_data)));
        return err;
    }) {
        .loaded => |result| return result,
        .stale => {
            posix.munmap(@alignCast(@constCast(fasta_data)));
            return error.StaleIndex;
        },
        .corrupt => {
            posix.munmap(@alignCast(@constCast(fasta_data)));
            return error.CorruptIndex;
        },
        .not_found => {
            posix.munmap(@alignCast(@constCast(fasta_data)));
            return zfi_failure orelse error.NoIndexFound;
        },
    }
}

fn parseZfiNameBlob(zfi_data: []const u8, records: []const IndexRecord) ?[]const u8 {
    if (records.len == 0 or !records[0].nameInZfi()) return null;

    const parsed = parseZfiNameFooter(zfi_data) orelse return null;
    const blob_len: usize = @intCast(parsed.name_blob_len);
    const footer_start = zfi_data.len - parsed.footer_bytes;
    if (footer_start < blob_len) return null;
    return zfi_data[footer_start - blob_len .. footer_start];
}

fn parseZfiNameFooter(zfi_data: []const u8) ?struct {
    name_blob_len: u64,
    footer_bytes: usize,
} {
    if (zfi_data.len >= zfi_name_footer_bytes) {
        if (decodeZfiNameFooterBytes(zfi_data[zfi_data.len - zfi_name_footer_bytes ..])) |name_blob_len| {
            return .{ .name_blob_len = name_blob_len, .footer_bytes = zfi_name_footer_bytes };
        }
    }

    if (zfi_data.len >= zfi_name_footer_legacy_bytes) {
        const tail = zfi_data[zfi_data.len - zfi_name_footer_legacy_bytes ..];
        if (std.mem.eql(u8, tail[0..4], &ZFI_NAME_FOOTER_MAGIC)) {
            const name_blob_len: u64 = @bitCast(tail[8..16].*);
            return .{ .name_blob_len = name_blob_len, .footer_bytes = zfi_name_footer_legacy_bytes };
        }
    }

    return null;
}

fn buildNameMapFast(
    records: []const IndexRecord,
    name_data: []const u8,
) LoadIndexError!struct {
    name_map: std.StringHashMap(usize),
    name_slices: []const []const u8,
} {
    var name_map = std.StringHashMap(usize).init(std.heap.page_allocator);
    name_map.ensureTotalCapacity(@intCast(records.len)) catch return error.OutOfMemory;
    for (records, 0..) |rec, i| {
        const name = rec.getName(name_data);
        name_map.putAssumeCapacity(name, i);
    }
    return .{ .name_map = name_map, .name_slices = &.{} };
}

fn buildNameTable(
    allocator: std.mem.Allocator,
    records: []const IndexRecord,
    name_data: []const u8,
    fasta_for_drop: ?[]align(4096) const u8,
) LoadIndexError!struct {
    name_map: std.StringHashMap(usize),
    name_slices: []const []const u8,
} {
    var name_map = std.StringHashMap(usize).init(allocator);
    const name_slices = try allocator.alloc([]const u8, records.len);

    name_map.ensureTotalCapacity(@intCast(records.len)) catch return error.OutOfMemory;
    const drop_stride: usize = 8 * 1024 * 1024;
    var next_drop_end: usize = drop_stride;
    for (records, 0..) |rec, i| {
        const name_src = rec.getName(name_data);
        const name_copy = try allocator.dupe(u8, name_src);
        name_slices[i] = name_copy;
        name_map.putAssumeCapacity(name_copy, i);
        if (fasta_for_drop) |fasta| {
            if (rec.name_offset > 0) {
                const header_start = @as(usize, @intCast(rec.name_offset - 1));
                if (header_start >= next_drop_end or i + 1 == records.len) {
                    dropFastaPrefix(fasta, header_start);
                    next_drop_end = header_start + drop_stride;
                }
            }
        }
    }

    if (fasta_for_drop) |fasta| {
        dropFastaCache(fasta);
    }

    return .{ .name_map = name_map, .name_slices = name_slices };
}

pub fn dropFastaPrefix(fasta_data: []const u8, end_exclusive: usize) void {
    dropFastaSpan(fasta_data, 0, end_exclusive);
}

pub fn dropFastaSpan(fasta_data: []const u8, start: usize, end_exclusive: usize) void {
    if (start >= end_exclusive or fasta_data.len == 0) return;
    const page = std.heap.page_size_min;
    const drop_start = std.mem.alignBackward(usize, start, page);
    const drop_end = @min(fasta_data.len, std.mem.alignForward(usize, end_exclusive, page));
    if (drop_end <= drop_start) return;
    posix.madvise(@alignCast(@constCast(fasta_data.ptr + drop_start)), drop_end - drop_start, posix.MADV.DONTNEED) catch {};
}

fn dropFastaCache(fasta_data: []align(4096) const u8) void {
    if (fasta_data.len == 0) return;
    const len = std.mem.alignBackward(usize, fasta_data.len, std.heap.page_size_min);
    if (len == 0) return;
    posix.madvise(@alignCast(@constCast(fasta_data.ptr)), len, posix.MADV.DONTNEED) catch {};
}

fn tryLoadZfi(
    io: std.Io,
    zfi_path: []const u8,
    fasta_data: []align(4096) const u8,
    fasta_stat: std.Io.File.Stat,
    mode: LoadMode,
) LoadIndexError!LoadAttempt {
    const zfi_file = std.Io.Dir.cwd().openFile(io, zfi_path, .{}) catch |err| switch (err) {
        error.FileNotFound => return .not_found,
        error.AccessDenied => return error.AccessDenied,
        else => return error.Io,
    };
    defer zfi_file.close(io);

    const zfi_stat = zfi_file.stat(io) catch return error.Io;

    // Check for 0-byte file before mmap (POSIX returns EINVAL)
    if (zfi_stat.size == 0) {
        return .corrupt;
    }

    // Staleness check: mtime
    if (zfi_stat.mtime.nanoseconds < fasta_stat.mtime.nanoseconds) {
        return .stale;
    }

    // mmap the .zfi
    const zfi_data = posix.mmap(
        null,
        zfi_stat.size,
        .{ .READ = true },
        .{ .TYPE = .PRIVATE },
        zfi_file.handle,
        0,
    ) catch return error.MmapFailed;
    errdefer posix.munmap(@alignCast(@constCast(zfi_data)));

    // Validate minimum size for header
    if (zfi_data.len < @sizeOf(ZfiHeader)) {
        return .corrupt;
    }

    // Validate magic
    const header: *const ZfiHeader = @ptrCast(@alignCast(zfi_data.ptr));
    if (!std.mem.eql(u8, &header.magic, &ZFI_MAGIC)) {
        return .corrupt;
    }

    // Validate source file size
    if (header.source_size != fasta_stat.size) {
        return .stale;
    }

    // Validate that the file has enough bytes for all records
    const expected_size = @sizeOf(ZfiHeader) + @as(usize, header.record_count) * @sizeOf(IndexRecord);
    if (zfi_data.len < expected_size) {
        return .corrupt;
    }

    // Cast record array from mmap bytes
    const record_bytes = zfi_data[@sizeOf(ZfiHeader)..];
    const records: []const IndexRecord = @as(
        [*]const IndexRecord,
        @ptrCast(@alignCast(record_bytes.ptr)),
    )[0..header.record_count];

    const name_blob = parseZfiNameBlob(zfi_data, records);

    for (records) |rec| {
        if (!isValidZfiRecordMetadata(rec, fasta_data.len, zfi_data, name_blob)) {
            return .corrupt;
        }
    }

    var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    errdefer arena.deinit();
    const allocator = arena.allocator();

    const build_name_map = mode == .lookup_full_map;

    var name_map = std.StringHashMap(usize).init(allocator);
    var name_slices: []const []const u8 = &.{};
    var has_name_map = false;
    if (build_name_map) {
        if (name_blob) |blob| {
            const built = buildNameMapFast(records, blob) catch return error.OutOfMemory;
            name_map = built.name_map;
            name_slices = built.name_slices;
        } else {
            const built = buildNameTable(allocator, records, fasta_data, fasta_data) catch return error.OutOfMemory;
            name_map = built.name_map;
            name_slices = built.name_slices;
        }
        has_name_map = true;
    }

    return .{ .loaded = LoadedIndex{
        .records = records,
        .name_map = name_map,
        .has_name_map = has_name_map,
        .name_slices = name_slices,
        .name_blob = name_blob,
        .fasta_data = fasta_data,
        .fasta_size = fasta_stat.size,
        .zfi_data = zfi_data,
        .source = .zfi,
        .arena = arena,
    } };
}

fn tryLoadFai(
    io: std.Io,
    fai_path: []const u8,
    fasta_data: []align(4096) const u8,
    fasta_stat: std.Io.File.Stat,
    mode: LoadMode,
) LoadIndexError!LoadAttempt {
    const fai_file = std.Io.Dir.cwd().openFile(io, fai_path, .{}) catch |err| switch (err) {
        error.FileNotFound => return .not_found,
        error.AccessDenied => return error.AccessDenied,
        else => return error.Io,
    };
    defer fai_file.close(io);

    const fai_stat = fai_file.stat(io) catch return error.Io;

    // Staleness check: mtime
    if (fai_stat.mtime.nanoseconds < fasta_stat.mtime.nanoseconds) {
        return .stale;
    }

    if (fai_stat.size == 0) {
        return .corrupt;
    }

    if (mode == .records_only) {
        return loadFaiRecordsOnly(io, fai_file, fasta_data, fasta_stat, true, null);
    }
    if (mode == .stats_scan) {
        return loadFaiStatsScan(io, fai_file, fai_path, fasta_data, fasta_stat);
    }

    const fai_data = posix.mmap(
        null,
        fai_stat.size,
        .{ .READ = true },
        .{ .TYPE = .PRIVATE },
        fai_file.handle,
        0,
    ) catch return error.MmapFailed;
    errdefer posix.munmap(@alignCast(@constCast(fai_data)));

    var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    errdefer arena.deinit();
    const allocator = arena.allocator();

    const approx_records = std.mem.count(u8, fai_data, "\n");
    const records = allocator.alloc(IndexRecord, approx_records) catch return error.OutOfMemory;

    var record_count: usize = 0;
    var pos: usize = 0;
    while (pos < fai_data.len) {
        const line_start = pos;
        const rel_eol = std.mem.indexOfScalar(u8, fai_data[pos..], '\n') orelse fai_data.len - pos;
        const line_len = rel_eol;
        pos += rel_eol + 1;
        if (line_len == 0) continue;

        const line = fai_data[line_start..][0..line_len];
        const name_end = std.mem.indexOfScalar(u8, line, '\t') orelse return .corrupt;
        var field_start: usize = name_end + 1;

        const seq_len = parseFaiFieldU64(line, &field_start) catch return .corrupt;
        const seq_offset = parseFaiFieldU64(line, &field_start) catch return .corrupt;
        const line_bases = parseFaiFieldU32(line, &field_start) catch return .corrupt;
        const line_bytes = parseFaiFieldU32(line, &field_start) catch return .corrupt;

        const rec = IndexRecord{
            .name_offset = line_start,
            .name_len = @intCast(name_end),
            .seq_offset = seq_offset,
            .seq_len = seq_len,
            .line_bases = line_bases,
            .line_bytes = line_bytes,
        };

        if (record_count >= records.len) return .corrupt;
        records[record_count] = rec;
        record_count += 1;
    }

    if (record_count == 0) {
        return .corrupt;
    }

    const loaded_records = records[0..record_count];
    const built = buildNameTable(allocator, loaded_records, fai_data, null) catch return error.OutOfMemory;

    return .{ .loaded = LoadedIndex{
        .records = loaded_records,
        .name_map = built.name_map,
        .has_name_map = true,
        .name_slices = built.name_slices,
        .fai_data = fai_data,
        .fasta_data = fasta_data,
        .fasta_size = fasta_stat.size,
        .zfi_data = null,
        .source = .fai,
        .arena = arena,
    } };
}

fn loadFaiRecordsOnly(
    io: std.Io,
    fai_file: std.Io.File,
    fasta_data: []align(4096) const u8,
    fasta_stat: std.Io.File.Stat,
    store_names: bool,
    sidecar_path: ?[]const u8,
) LoadIndexError!LoadAttempt {
    var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    errdefer arena.deinit();
    const allocator = arena.allocator();

    var records_list: std.ArrayListUnmanaged(IndexRecord) = .empty;
    errdefer records_list.deinit(allocator);
    var slices_list: std.ArrayListUnmanaged([]const u8) = .empty;
    errdefer slices_list.deinit(allocator);
    var offsets_list: std.ArrayListUnmanaged(u64) = .empty;
    errdefer offsets_list.deinit(allocator);
    const track_line_offsets = sidecar_path != null and !store_names;

    var io_buf: [65536]u8 = undefined;
    var file_reader = fai_file.reader(io, &io_buf);
    var file_offset: u64 = 0;

    while (true) {
        const maybe_line = file_reader.interface.takeDelimiter('\n') catch return error.Io;
        const line = maybe_line orelse break;
        const line_start_offset = file_offset;
        file_offset += @as(u64, @intCast(line.len)) + 1;
        if (line.len == 0) continue;

        const name_end = std.mem.indexOfScalar(u8, line, '\t') orelse return .corrupt;
        var field_start: usize = name_end + 1;

        const seq_len = parseFaiFieldU64(line, &field_start) catch return .corrupt;
        const seq_offset = parseFaiFieldU64(line, &field_start) catch return .corrupt;
        const line_bases = parseFaiFieldU32(line, &field_start) catch return .corrupt;
        const line_bytes = parseFaiFieldU32(line, &field_start) catch return .corrupt;

        const name = if (store_names)
            allocator.dupe(u8, line[0..name_end]) catch return error.OutOfMemory
        else
            @as([]const u8, &.{});
        if (store_names) try slices_list.append(allocator, name);
        if (track_line_offsets) try offsets_list.append(allocator, line_start_offset);
        try records_list.append(allocator, .{
            .name_offset = 0,
            .name_len = if (store_names) @intCast(name_end) else 0,
            .seq_offset = seq_offset,
            .seq_len = seq_len,
            .line_bases = line_bases,
            .line_bytes = line_bytes,
        });
    }

    if (records_list.items.len == 0) return .corrupt;

    const path_copy = if (sidecar_path) |path|
        allocator.dupe(u8, path) catch return error.OutOfMemory
    else
        null;

    return .{ .loaded = LoadedIndex{
        .records = records_list.toOwnedSlice(allocator) catch return error.OutOfMemory,
        .name_map = std.StringHashMap(usize).init(allocator),
        .has_name_map = false,
        .name_slices = if (store_names) slices_list.toOwnedSlice(allocator) catch return error.OutOfMemory else &.{},
        .fai_data = null,
        .fasta_data = fasta_data,
        .fasta_size = fasta_stat.size,
        .zfi_data = null,
        .source = .fai,
        .sidecar_path = path_copy,
        .fai_line_offsets = if (track_line_offsets) offsets_list.toOwnedSlice(allocator) catch return error.OutOfMemory else &.{},
        .arena = arena,
    } };
}

fn readFaiNameAtOffset(io: std.Io, allocator: std.mem.Allocator, fai_path: []const u8, offset: u64) LoadIndexError![]const u8 {
    const fai_file = std.Io.Dir.cwd().openFile(io, fai_path, .{}) catch return error.Io;
    defer fai_file.close(io);

    var io_buf: [65536]u8 = undefined;
    var file_reader = fai_file.reader(io, &io_buf);
    file_reader.seekTo(offset) catch return error.Io;
    const maybe_line = file_reader.interface.takeDelimiter('\n') catch return error.Io;
    const line = maybe_line orelse return error.CorruptIndex;
    if (line.len == 0) return error.CorruptIndex;
    const name_end = std.mem.indexOfScalar(u8, line, '\t') orelse return error.CorruptIndex;
    return allocator.dupe(u8, line[0..name_end]);
}

fn loadFaiStatsScan(
    io: std.Io,
    fai_file: std.Io.File,
    fai_path: []const u8,
    fasta_data: []align(4096) const u8,
    fasta_stat: std.Io.File.Stat,
) LoadIndexError!LoadAttempt {
    return loadFaiRecordsOnly(io, fai_file, fasta_data, fasta_stat, false, fai_path);
}

fn readFaiNameLine(io: std.Io, allocator: std.mem.Allocator, fai_path: []const u8, rec_idx: usize) LoadIndexError![]const u8 {
    const fai_file = std.Io.Dir.cwd().openFile(io, fai_path, .{}) catch return error.Io;
    defer fai_file.close(io);

    var io_buf: [65536]u8 = undefined;
    var file_reader = fai_file.reader(io, &io_buf);

    var line_no: usize = 0;
    while (true) {
        const maybe_line = file_reader.interface.takeDelimiter('\n') catch return error.Io;
        const line = maybe_line orelse return error.CorruptIndex;
        if (line.len == 0) continue;
        if (line_no == rec_idx) {
            const name_end = std.mem.indexOfScalar(u8, line, '\t') orelse return error.CorruptIndex;
            return allocator.dupe(u8, line[0..name_end]);
        }
        line_no += 1;
    }
}

fn parseFaiFieldU64(line: []const u8, field_start: *usize) LoadIndexError!u64 {
    if (field_start.* < line.len and line[field_start.*] == '\t') {
        field_start.* += 1;
    }
    if (field_start.* >= line.len) return error.CorruptIndex;

    const field_end = std.mem.indexOfScalarPos(u8, line, field_start.*, '\t') orelse line.len;
    const value = parseFaiAsciiU64(line[field_start.*..field_end]) orelse return error.CorruptIndex;
    field_start.* = field_end;
    return value;
}

fn parseFaiFieldU32(line: []const u8, field_start: *usize) LoadIndexError!u32 {
    if (field_start.* < line.len and line[field_start.*] == '\t') {
        field_start.* += 1;
    }
    if (field_start.* >= line.len) return error.CorruptIndex;

    const field_end = std.mem.indexOfScalarPos(u8, line, field_start.*, '\t') orelse line.len;
    const value = parseFaiAsciiU32(line[field_start.*..field_end]) orelse return error.CorruptIndex;
    field_start.* = field_end;
    return value;
}

fn parseFaiAsciiU64(text: []const u8) ?u64 {
    if (text.len == 0) return null;
    var value: u64 = 0;
    for (text) |byte| {
        if (byte < '0' or byte > '9') return null;
        value = value * 10 + (byte - '0');
    }
    return value;
}

fn parseFaiAsciiU32(text: []const u8) ?u32 {
    const wide = parseFaiAsciiU64(text) orelse return null;
    if (wide > std.math.maxInt(u32)) return null;
    return @intCast(wide);
}

fn isValidZfiRecordMetadata(
    rec: IndexRecord,
    fasta_len: usize,
    zfi_data: []align(4096) const u8,
    name_blob: ?[]const u8,
) bool {
    if (rec.name_len == 0) return false;
    if (rec.nameInZfi()) {
        const blob = name_blob orelse return false;
        if (rec.name_offset + rec.name_len > blob.len) return false;
    } else {
        if (rec.name_offset == 0 or rec.name_offset > fasta_len) return false;
        if (rec.name_offset + rec.name_len > fasta_len) return false;
    }

    if (rec.seq_len == 0) return false;
    if (rec.seq_offset >= fasta_len) return false;
    if (rec.line_bases == 0 or rec.line_bytes == 0) return false;
    if (rec.line_bytes < rec.line_bases) return false;

    if (!rec.isUniformWidth()) {
        return isValidSideTable(rec, fasta_len, zfi_data);
    }

    const full_lines: u128 = rec.seq_len / rec.line_bases;
    const remainder: u128 = rec.seq_len % rec.line_bases;
    var region_end: u128 = rec.seq_offset;
    if (remainder > 0) {
        region_end += full_lines * rec.line_bytes + remainder;
    } else {
        region_end += full_lines * rec.line_bytes;
    }

    return region_end <= fasta_len;
}

fn isValidZfiRecord(
    rec: IndexRecord,
    fasta_data: []align(4096) const u8,
    zfi_data: []align(4096) const u8,
) bool {
    if (!isValidZfiRecordMetadata(rec, fasta_data.len, zfi_data, null)) return false;
    if (rec.nameInZfi()) return true;
    if (fasta_data[rec.name_offset - 1] != '>') return false;
    return true;
}

fn isValidSideTable(rec: IndexRecord, fasta_len: usize, zfi_data: []align(4096) const u8) bool {
    const offset: usize = @intCast(rec.sideTableOffset());
    if (offset < @sizeOf(ZfiHeader)) return false;
    if (offset + @sizeOf(u64) > zfi_data.len) return false;

    const count_ptr: *const u64 = @ptrCast(@alignCast(zfi_data[offset..].ptr));
    const line_count = count_ptr.*;
    if (line_count == 0) return false;
    if (line_count > std.math.maxInt(usize) / @sizeOf(SideTableLine)) return false;

    const table_bytes = @as(usize, @intCast(line_count)) * @sizeOf(SideTableLine);
    const lines_offset = offset + @sizeOf(u64);
    if (lines_offset + table_bytes > zfi_data.len) return false;

    const lines = @as(
        [*]const SideTableLine,
        @ptrCast(@alignCast(zfi_data[lines_offset..].ptr)),
    )[0..@intCast(line_count)];

    var expected_base_start: u64 = 0;
    var previous_byte_offset: u64 = 0;
    for (lines, 0..) |line, i| {
        if (line.base_start != expected_base_start) return false;
        if (line.line_bases == 0 or line.line_bytes == 0) return false;
        if (line.line_bytes < line.line_bases) return false;
        if (line.byte_offset >= fasta_len) return false;
        if (line.byte_offset + line.line_bytes > fasta_len) return false;
        if (i > 0 and line.byte_offset <= previous_byte_offset) return false;

        expected_base_start += line.line_bases;
        previous_byte_offset = line.byte_offset;
    }

    return expected_base_start == rec.seq_len;
}
