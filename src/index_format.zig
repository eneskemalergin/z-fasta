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
    // Open FASTA file
    const fasta_file = std.fs.cwd().openFile(fasta_path, .{}) catch |err| {
        switch (err) {
            error.FileNotFound => printErrorAndExit("error: file not found: {s}\n", .{fasta_path}),
            error.AccessDenied => printErrorAndExit("error: access denied: {s}\n", .{fasta_path}),
            else => printErrorAndExit("error: failed to open file: {s}\n", .{fasta_path}),
        }
    };
    defer fasta_file.close();

    const fasta_stat = fasta_file.stat() catch {
        printErrorAndExit("error: failed to stat file: {s}\n", .{fasta_path});
    };

    if (fasta_stat.size == 0) {
        printErrorAndExit("error: file is empty: {s}\n", .{fasta_path});
    }

    // mmap FASTA
    const fasta_data = posix.mmap(
        null,
        fasta_stat.size,
        posix.PROT.READ,
        .{ .TYPE = .PRIVATE },
        fasta_file.handle,
        0,
    ) catch {
        printErrorAndExit("error: failed to mmap file: {s}\n", .{fasta_path});
    };

    // Try .zfi first
    var zfi_path_buf: [4096]u8 = undefined;
    const zfi_path = std.fmt.bufPrint(&zfi_path_buf, "{s}.zfi", .{fasta_path}) catch {
        printErrorAndExit("error: path too long\n", .{});
    };

    if (tryLoadZfi(zfi_path, fasta_data, fasta_stat)) |result| {
        return result;
    }

    // Try .fai fallback
    var fai_path_buf: [4096]u8 = undefined;
    const fai_path = std.fmt.bufPrint(&fai_path_buf, "{s}.fai", .{fasta_path}) catch {
        printErrorAndExit("error: path too long\n", .{});
    };

    if (tryLoadFai(fai_path, fasta_data, fasta_stat)) |result| {
        return result;
    }

    printErrorAndExit("error: no index found for {s}. Run 'z-fasta index {s}' first.\n", .{ fasta_path, fasta_path });
}

fn tryLoadZfi(
    zfi_path: []const u8,
    fasta_data: []align(4096) const u8,
    fasta_stat: std.fs.File.Stat,
) ?LoadedIndex {
    const zfi_file = std.fs.cwd().openFile(zfi_path, .{}) catch {
        return null; // .zfi not found, not an error
    };
    defer zfi_file.close();

    const zfi_stat = zfi_file.stat() catch {
        printErrorAndExit("error: failed to stat index: {s}\n", .{zfi_path});
    };

    // Check for 0-byte file before mmap (POSIX returns EINVAL)
    if (zfi_stat.size == 0) {
        printErrorAndExit("error: corrupt index file: {s}\n", .{zfi_path});
    }

    // Staleness check: mtime
    if (zfi_stat.mtime < fasta_stat.mtime) {
        printErrorAndExit("error: index is stale (FASTA is newer than index). Re-run 'z-fasta index <path>'.\n", .{});
    }

    // mmap the .zfi
    const zfi_data = posix.mmap(
        null,
        zfi_stat.size,
        posix.PROT.READ,
        .{ .TYPE = .PRIVATE },
        zfi_file.handle,
        0,
    ) catch {
        printErrorAndExit("error: failed to mmap index: {s}\n", .{zfi_path});
    };

    // Validate minimum size for header
    if (zfi_data.len < @sizeOf(ZfiHeader)) {
        printErrorAndExit("error: corrupt index file: {s}\n", .{zfi_path});
    }

    // Validate magic
    const header: *const ZfiHeader = @ptrCast(@alignCast(zfi_data.ptr));
    if (!std.mem.eql(u8, &header.magic, &ZFI_MAGIC)) {
        printErrorAndExit("error: corrupt index file: {s}\n", .{zfi_path});
    }

    // Validate source file size
    if (header.source_size != fasta_stat.size) {
        printErrorAndExit("error: index is stale (file size changed). Re-run 'z-fasta index <path>'.\n", .{});
    }

    // Validate that the file has enough bytes for all records
    const expected_size = @sizeOf(ZfiHeader) + @as(usize, header.record_count) * @sizeOf(IndexRecord);
    if (zfi_data.len < expected_size) {
        printErrorAndExit("error: corrupt index file: {s}\n", .{zfi_path});
    }

    // Cast record array from mmap bytes
    const record_bytes = zfi_data[@sizeOf(ZfiHeader)..];
    const records: []const IndexRecord = @as(
        [*]const IndexRecord,
        @ptrCast(@alignCast(record_bytes.ptr)),
    )[0..header.record_count];

    // Build name map from FASTA mmap slices
    var name_map = std.StringHashMap(usize).init(std.heap.page_allocator);
    for (records, 0..) |rec, i| {
        const name = fasta_data[rec.name_offset..][0..rec.name_len];
        name_map.put(name, i) catch {
            printErrorAndExit("error: failed to build name lookup\n", .{});
        };
    }

    return LoadedIndex{
        .records = records,
        .name_map = name_map,
        .fasta_data = fasta_data,
        .fasta_size = fasta_stat.size,
        .zfi_data = zfi_data,
        .source = .zfi,
        .arena = std.heap.ArenaAllocator.init(std.heap.page_allocator),
    };
}

fn tryLoadFai(
    fai_path: []const u8,
    fasta_data: []align(4096) const u8,
    fasta_stat: std.fs.File.Stat,
) ?LoadedIndex {
    const fai_file = std.fs.cwd().openFile(fai_path, .{}) catch {
        return null; // .fai not found, not an error
    };
    defer fai_file.close();

    const fai_stat = fai_file.stat() catch {
        printErrorAndExit("error: failed to stat index: {s}\n", .{fai_path});
    };

    // Staleness check: mtime
    if (fai_stat.mtime < fasta_stat.mtime) {
        printErrorAndExit("error: index is stale (FASTA is newer than index). Re-run 'z-fasta index <path>'.\n", .{});
    }

    if (fai_stat.size == 0) {
        printErrorAndExit("error: corrupt index file: {s}\n", .{fai_path});
    }

    // Read .fai into memory (they're small)
    var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    errdefer arena.deinit();
    const allocator = arena.allocator();

    const fai_contents = allocator.alloc(u8, fai_stat.size) catch {
        printErrorAndExit("error: out of memory loading index\n", .{});
    };
    const bytes_read = fai_file.readAll(fai_contents) catch {
        printErrorAndExit("error: failed to read index: {s}\n", .{fai_path});
    };
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
            printErrorAndExit("error: corrupt index file: {s}\n", .{fai_path});
        };
        const offset_str = fields.next() orelse {
            printErrorAndExit("error: corrupt index file: {s}\n", .{fai_path});
        };
        const lb_str = fields.next() orelse {
            printErrorAndExit("error: corrupt index file: {s}\n", .{fai_path});
        };
        const lbytes_str = fields.next() orelse {
            printErrorAndExit("error: corrupt index file: {s}\n", .{fai_path});
        };
        // Ignore optional 6th field (qual_offset for FASTQ)

        const seq_len = std.fmt.parseInt(u64, len_str, 10) catch {
            printErrorAndExit("error: corrupt index file: {s}\n", .{fai_path});
        };
        const seq_offset = std.fmt.parseInt(u64, offset_str, 10) catch {
            printErrorAndExit("error: corrupt index file: {s}\n", .{fai_path});
        };
        const line_bases = std.fmt.parseInt(u32, lb_str, 10) catch {
            printErrorAndExit("error: corrupt index file: {s}\n", .{fai_path});
        };
        const line_bytes = std.fmt.parseInt(u32, lbytes_str, 10) catch {
            printErrorAndExit("error: corrupt index file: {s}\n", .{fai_path});
        };

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
        records_list.append(rec) catch {
            printErrorAndExit("error: out of memory loading index\n", .{});
        };

        // Store arena-owned name slice in the map
        const name_owned = allocator.dupe(u8, name) catch {
            printErrorAndExit("error: out of memory loading index\n", .{});
        };
        name_map.put(name_owned, idx) catch {
            printErrorAndExit("error: out of memory loading index\n", .{});
        };
    }

    if (records_list.items.len == 0) {
        printErrorAndExit("error: corrupt index file: {s}\n", .{fai_path});
    }

    return LoadedIndex{
        .records = records_list.items,
        .name_map = name_map,
        .fasta_data = fasta_data,
        .fasta_size = fasta_stat.size,
        .zfi_data = null,
        .source = .fai,
        .arena = arena,
    };
}
