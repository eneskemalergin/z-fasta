const std = @import("std");
const posix = std.posix;

// ============================================================================
// Types — shared by indexer, getter, stats
// ============================================================================

/// ZFI binary format magic + header
pub const ZFI_MAGIC: [4]u8 = .{ 'Z', 'F', 'I', 0x01 };

pub const ZfiHeader = extern struct {
    magic: [4]u8,
    record_count: u32,
    source_size: u64,
};

/// Index record for ZFI output (40 bytes padded)
pub const IndexRecord = extern struct {
    name_offset: u64,
    name_len: u16,
    _pad: [6]u8 = .{0} ** 6,
    seq_offset: u64,
    seq_len: u64,
    line_bases: u32,
    line_bytes: u32,

    pub fn getName(self: IndexRecord, data: []const u8) []const u8 {
        return data[self.name_offset..][0..self.name_len];
    }
};

/// Result of loading an index (from .zfi or .fai)
pub const LoadedIndex = struct {
    records: []const IndexRecord,
    name_map: std.StringHashMap(usize),
    fasta_data: []align(4096) const u8,
    fasta_size: u64,
    zfi_data: ?[]align(4096) const u8,
    source: IndexSource,
    arena: std.heap.ArenaAllocator,

    pub const IndexSource = enum { zfi, fai };

    pub fn deinit(self: *LoadedIndex) void {
        // munmap FASTA
        posix.munmap(@constCast(@alignCast(self.fasta_data)));
        // munmap .zfi if loaded
        if (self.zfi_data) |zd| {
            posix.munmap(@constCast(@alignCast(zd)));
        }
        self.name_map.deinit();
        self.arena.deinit();
    }

    pub fn lookupName(self: *const LoadedIndex, name: []const u8) ?usize {
        return self.name_map.get(name);
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
// Error helper (same pattern as main.zig)
// ============================================================================

pub fn printErrorAndExit(comptime fmt: []const u8, args: anytype) noreturn {
    std.io.getStdErr().writer().print(fmt, args) catch {};
    std.process.exit(1);
}

// ============================================================================
// .zfi writer (used by indexer)
// ============================================================================

/// Writes the .zfi binary index file.
pub fn writeZfi(
    path: []const u8,
    records: []const IndexRecord,
    source_size: u64,
) !void {
    const file = try std.fs.cwd().createFile(path, .{});
    errdefer std.fs.cwd().deleteFile(path) catch {};
    defer file.close();

    var writer = file.writer();

    const header = ZfiHeader{
        .magic = ZFI_MAGIC,
        .record_count = @intCast(records.len),
        .source_size = source_size,
    };
    try writer.writeAll(std.mem.asBytes(&header));
    try writer.writeAll(std.mem.sliceAsBytes(records));
}

// ============================================================================
// Shared index loader — loads .zfi or falls back to .fai
// ============================================================================

/// Load the index for a FASTA file. Tries .zfi first, then .fai fallback.
/// The caller must call deinit() on the returned LoadedIndex.
pub fn loadIndex(fasta_path: []const u8) LoadedIndex {
    return loadIndexChecked(fasta_path) catch |err| switch (err) {
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
pub fn loadIndexChecked(fasta_path: []const u8) LoadIndexError!LoadedIndex {
    // Open FASTA file
    const fasta_file = std.fs.cwd().openFile(fasta_path, .{}) catch |err| switch (err) {
        error.FileNotFound => return error.FileNotFound,
        error.AccessDenied => return error.AccessDenied,
        else => return error.Io,
    };
    defer fasta_file.close();

    const fasta_stat = fasta_file.stat() catch return error.Io;

    if (fasta_stat.size == 0) {
        return error.EmptyFile;
    }

    // mmap FASTA
    const fasta_data = posix.mmap(
        null,
        fasta_stat.size,
        posix.PROT.READ,
        .{ .TYPE = .PRIVATE },
        fasta_file.handle,
        0,
    ) catch return error.MmapFailed;

    // Try .zfi first
    var zfi_path_buf: [4096]u8 = undefined;
    const zfi_path = std.fmt.bufPrint(&zfi_path_buf, "{s}.zfi", .{fasta_path}) catch return error.PathTooLong;

    var zfi_failure: ?LoadIndexError = null;
    switch (tryLoadZfi(zfi_path, fasta_data, fasta_stat) catch |err| {
        posix.munmap(@constCast(@alignCast(fasta_data)));
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
        posix.munmap(@constCast(@alignCast(fasta_data)));
        return error.PathTooLong;
    };

    switch (tryLoadFai(fai_path, fasta_data, fasta_stat) catch |err| {
        posix.munmap(@constCast(@alignCast(fasta_data)));
        return err;
    }) {
        .loaded => |result| return result,
        .stale => {
            posix.munmap(@constCast(@alignCast(fasta_data)));
            return error.StaleIndex;
        },
        .corrupt => {
            posix.munmap(@constCast(@alignCast(fasta_data)));
            return error.CorruptIndex;
        },
        .not_found => {
            posix.munmap(@constCast(@alignCast(fasta_data)));
            return zfi_failure orelse error.NoIndexFound;
        },
    }
}

fn tryLoadZfi(
    zfi_path: []const u8,
    fasta_data: []align(4096) const u8,
    fasta_stat: std.fs.File.Stat,
) LoadIndexError!LoadAttempt {
    const zfi_file = std.fs.cwd().openFile(zfi_path, .{}) catch |err| switch (err) {
        error.FileNotFound => return .not_found,
        error.AccessDenied => return error.AccessDenied,
        else => return error.Io,
    };
    defer zfi_file.close();

    const zfi_stat = zfi_file.stat() catch return error.Io;

    // Check for 0-byte file before mmap (POSIX returns EINVAL)
    if (zfi_stat.size == 0) {
        return .corrupt;
    }

    // Staleness check: mtime
    if (zfi_stat.mtime < fasta_stat.mtime) {
        return .stale;
    }

    // mmap the .zfi
    const zfi_data = posix.mmap(
        null,
        zfi_stat.size,
        posix.PROT.READ,
        .{ .TYPE = .PRIVATE },
        zfi_file.handle,
        0,
    ) catch return error.MmapFailed;
    errdefer posix.munmap(@constCast(@alignCast(zfi_data)));

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

    for (records) |rec| {
        if (!isValidZfiRecord(rec, fasta_data)) {
            return .corrupt;
        }
    }

    // Build name map from FASTA mmap slices
    var name_map = std.StringHashMap(usize).init(std.heap.page_allocator);
    errdefer name_map.deinit();
    for (records, 0..) |rec, i| {
        const name = fasta_data[rec.name_offset..][0..rec.name_len];
        name_map.put(name, i) catch return error.OutOfMemory;
    }

    return .{ .loaded = LoadedIndex{
        .records = records,
        .name_map = name_map,
        .fasta_data = fasta_data,
        .fasta_size = fasta_stat.size,
        .zfi_data = zfi_data,
        .source = .zfi,
        .arena = std.heap.ArenaAllocator.init(std.heap.page_allocator),
    } };
}

fn tryLoadFai(
    fai_path: []const u8,
    fasta_data: []align(4096) const u8,
    fasta_stat: std.fs.File.Stat,
) LoadIndexError!LoadAttempt {
    const fai_file = std.fs.cwd().openFile(fai_path, .{}) catch |err| switch (err) {
        error.FileNotFound => return .not_found,
        error.AccessDenied => return error.AccessDenied,
        else => return error.Io,
    };
    defer fai_file.close();

    const fai_stat = fai_file.stat() catch return error.Io;

    // Staleness check: mtime
    if (fai_stat.mtime < fasta_stat.mtime) {
        return .stale;
    }

    if (fai_stat.size == 0) {
        return .corrupt;
    }

    // Read .fai into memory (they're small)
    var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    errdefer arena.deinit();
    const allocator = arena.allocator();

    const fai_contents = allocator.alloc(u8, fai_stat.size) catch return error.OutOfMemory;
    const bytes_read = fai_file.readAll(fai_contents) catch return error.Io;
    const fai_data = fai_contents[0..bytes_read];

    // Parse .fai lines: NAME\tLENGTH\tOFFSET\tLINE_BASES\tLINE_BYTES[\tQUAL_OFFSET]
    var records_list = std.ArrayList(IndexRecord).init(allocator);
    var name_map = std.StringHashMap(usize).init(allocator);

    var line_iter = std.mem.splitScalar(u8, fai_data, '\n');
    while (line_iter.next()) |line| {
        if (line.len == 0) continue;

        var fields = std.mem.splitScalar(u8, line, '\t');
        const name = fields.next() orelse continue;
        const len_str = fields.next() orelse {
            return .corrupt;
        };
        const offset_str = fields.next() orelse {
            return .corrupt;
        };
        const lb_str = fields.next() orelse {
            return .corrupt;
        };
        const lbytes_str = fields.next() orelse {
            return .corrupt;
        };
        // Ignore optional 6th field (qual_offset for FASTQ)

        const seq_len = std.fmt.parseInt(u64, len_str, 10) catch return .corrupt;
        const seq_offset = std.fmt.parseInt(u64, offset_str, 10) catch return .corrupt;
        const line_bases = std.fmt.parseInt(u32, lb_str, 10) catch return .corrupt;
        const line_bytes = std.fmt.parseInt(u32, lbytes_str, 10) catch return .corrupt;

        // For .fai records, name_offset/name_len refer to the name in fai_data
        // which is arena-owned. We DON'T scan the FASTA to recover name offsets.
        // Instead we store 0/0 and rely on the name_map for lookups.
        const rec = IndexRecord{
            .name_offset = 0,
            .name_len = 0,
            .seq_offset = seq_offset,
            .seq_len = seq_len,
            .line_bases = line_bases,
            .line_bytes = line_bytes,
        };

        const idx = records_list.items.len;
        records_list.append(rec) catch return error.OutOfMemory;

        // Store arena-owned name slice in the map
        const name_owned = allocator.dupe(u8, name) catch return error.OutOfMemory;
        name_map.put(name_owned, idx) catch return error.OutOfMemory;
    }

    if (records_list.items.len == 0) {
        return .corrupt;
    }

    return .{ .loaded = LoadedIndex{
        .records = records_list.items,
        .name_map = name_map,
        .fasta_data = fasta_data,
        .fasta_size = fasta_stat.size,
        .zfi_data = null,
        .source = .fai,
        .arena = arena,
    } };
}

fn isValidZfiRecord(rec: IndexRecord, fasta_data: []align(4096) const u8) bool {
    const fasta_len = fasta_data.len;

    if (rec.name_offset == 0 or rec.name_offset > fasta_len) return false;
    if (rec.name_len == 0) return false;
    if (rec.name_offset + rec.name_len > fasta_len) return false;
    if (fasta_data[rec.name_offset - 1] != '>') return false;

    if (rec.seq_len == 0 or rec.seq_offset > fasta_len) return false;
    if (rec.line_bases == 0 or rec.line_bytes == 0) return false;
    if (rec.line_bytes < rec.line_bases) return false;

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
