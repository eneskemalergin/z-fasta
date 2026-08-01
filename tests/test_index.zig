//! Index unit and integration tests: scan parity, `.zfi`/`.fai` load validation, fixtures.
//!
//! Covers mmap vs streaming identity, malformed indexes, gates fixtures under
//! `tests/data/gates/`, and optional large-data / generated edge-case walks.

const std = @import("std");
const builtin = @import("builtin");
const main = @import("main");
const validateFasta = main.validateFasta;
const scanHeaders = main.scanHeaders;
const scanChunkedData = main.indexer.scanChunkedData;
const low_mem_chunk_size = main.indexer.low_mem_chunk_size;
const loadIndexChecked = main.index_format.loadIndexChecked;
const IndexRecord = main.IndexRecord;
const ZfiHeader = main.ZfiHeader;
const ZFI_MAGIC = main.ZFI_MAGIC;
const c = std.c;
const io = std.testing.io;

fn writeZfiFromRecords(path: []const u8, records: []const IndexRecord, source_size: u64, source_mtime_ns: u64, allocator: std.mem.Allocator) !void {
    var index = try main.indexer.zfiIndexFromRecords(records, allocator);
    defer index.deinit(allocator);
    try main.indexer.writeZfiIndexFile(io, path, &index, source_size, source_mtime_ns);
}

fn statMtimeNs(path: []const u8) !u64 {
    const file = try std.Io.Dir.cwd().openFile(io, path, .{});
    defer file.close(io);
    const st = try file.stat(io);
    return main.index_format.timestampToNs(st.mtime);
}

/// Write FASTA + raw `.zfi` bytes, patching a `ZFID` trailer (before the name footer)
/// to the FASTA mtime so load tests exercise structural checks instead of false staleness.
fn writeFastaAndRawZfi(allocator: std.mem.Allocator, stem: []const u8, fasta: []const u8, zfi_bytes: []const u8) !struct { fasta_path: []u8, zfi_path: []u8 } {
    const fasta_path = try uniqueArtifactPath(allocator, stem, "fa");
    const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{fasta_path});

    {
        const fasta_file = try std.Io.Dir.cwd().createFile(io, fasta_path, .{});
        defer fasta_file.close(io);
        try std.Io.File.writeStreamingAll(fasta_file, io, fasta);
    }
    // Re-open for mtime: Windows `stat` on the create handle can return AccessDenied.
    const fasta_mtime = try statMtimeNs(fasta_path);

    var zfi_copy = try allocator.dupe(u8, zfi_bytes);
    defer allocator.free(zfi_copy);
    const footer_len = main.index_format.zfi_name_footer_bytes;
    const id_len = main.index_format.zfi_source_id_bytes;
    // Production layout: `[name blob][ZFID][ZFNM footer]`.
    if (zfi_copy.len >= footer_len + id_len) {
        const id_off = zfi_copy.len - footer_len - id_len;
        if (std.mem.eql(u8, zfi_copy[id_off..][0..4], &main.index_format.ZFI_SOURCE_ID_MAGIC)) {
            std.mem.writeInt(u64, zfi_copy[id_off + 4 ..][0..8], fasta_mtime, .little);
        }
    }

    const zfi_file = try std.Io.Dir.cwd().createFile(io, zfi_path, .{ .truncate = true });
    defer zfi_file.close(io);
    try std.Io.File.writeStreamingAll(zfi_file, io, zfi_copy);
    return .{ .fasta_path = fasta_path, .zfi_path = zfi_path };
}

fn supportsPosixFutimens() bool {
    return switch (builtin.os.tag) {
        .linux, .dragonfly, .freebsd, .netbsd, .openbsd, .illumos, .macos, .ios, .tvos, .watchos, .visionos, .driverkit, .maccatalyst => true,
        else => false,
    };
}

fn uniqueArtifactPath(allocator: std.mem.Allocator, stem: []const u8, ext: []const u8) ![]u8 {
    try std.Io.Dir.cwd().createDirPath(io, "zig-cache/test-artifacts");
    const now = std.Io.Clock.Timestamp.now(io, .awake);
    const nanos: u64 = @intCast(now.raw.toNanoseconds());
    return std.fmt.allocPrint(allocator, "zig-cache/test-artifacts/{s}-{d}.{s}", .{
        stem,
        nanos,
        ext,
    });
}

fn markFileStaleOneHourAgo(file: std.Io.File) !void {
    if (!supportsPosixFutimens()) return error.SkipZigTest;

    switch (builtin.os.tag) {
        .linux => {
            var now: std.os.linux.timespec = undefined;
            if (std.os.linux.errno(std.os.linux.clock_gettime(.REALTIME, &now)) != .SUCCESS) {
                return error.Unexpected;
            }
            const stale: std.os.linux.timespec = .{ .sec = now.sec - 3600, .nsec = now.nsec };
            if (std.os.linux.errno(std.os.linux.futimens(file.handle, &.{ stale, stale })) != .SUCCESS) {
                return error.Unexpected;
            }
        },
        .dragonfly,
        .freebsd,
        .netbsd,
        .openbsd,
        .illumos,
        .macos,
        .ios,
        .tvos,
        .watchos,
        .visionos,
        .driverkit,
        .maccatalyst,
        => {
            var now: c.timespec = undefined;
            if (c.clock_gettime(c.CLOCK.REALTIME, &now) != 0) return error.Unexpected;
            const stale: c.timespec = .{ .sec = now.sec - 3600, .nsec = now.nsec };
            if (c.futimens(file.handle, &.{ stale, stale }) != 0) return error.Unexpected;
        },
        else => return error.SkipZigTest,
    }
}

// ============================================================================
// Test helper: read file into memory
// ============================================================================

fn readTestFile(allocator: std.mem.Allocator, path: []const u8) ![]u8 {
    const file = std.Io.Dir.cwd().openFile(io, path, .{}) catch |err| switch (err) {
        error.FileNotFound => {
            std.debug.print("required test fixture missing: {s}\n", .{path});
            return error.RequiredFixtureMissing;
        },
        else => |e| return e,
    };
    defer file.close(io);
    const stat = try file.stat(io);
    const data = try allocator.alloc(u8, @intCast(stat.size));
    const bytes_read = try std.Io.File.readPositionalAll(file, io, data, 0);
    return data[0..bytes_read];
}

/// In-repo fixtures required for CI. Missing path fails with a clear message.
fn requireFixturePath(path: []const u8) !void {
    std.Io.Dir.cwd().access(io, path, .{}) catch |err| switch (err) {
        error.FileNotFound => {
            std.debug.print("required test fixture missing: {s}\n", .{path});
            return error.RequiredFixtureMissing;
        },
        else => |e| return e,
    };
}

/// Optional large or generated fixtures. Prints why and skips the test.
fn skipOptionalFixture(path: []const u8, how_to_obtain: []const u8) error{SkipZigTest} {
    std.debug.print("optional fixture absent, skipping: {s}\n  obtain via: {s}\n", .{ path, how_to_obtain });
    return error.SkipZigTest;
}

fn expectZfiStreamingMatchesMmap(allocator: std.mem.Allocator, data: []const u8, label: []const u8) !void {
    var mmap_index = try main.indexer.scanZfiIndex(data, true, allocator);
    defer mmap_index.deinit(allocator);
    var stream_index = try main.indexer.scanZfiIndexStreamingData(data, true, allocator);
    defer stream_index.deinit(allocator);

    const mmap_bytes = try main.indexer.zfiIndexToBytes(&mmap_index, data.len, 0, allocator);
    defer allocator.free(mmap_bytes);
    const stream_bytes = try main.indexer.zfiIndexToBytes(&stream_index, data.len, 0, allocator);
    defer allocator.free(stream_bytes);

    std.testing.expectEqualSlices(u8, mmap_bytes, stream_bytes) catch {
        std.debug.print("ZFI mmap/stream mismatch: {s}\n", .{label});
        return error.TestExpectedEqual;
    };
}

fn walkFastaDirZfiParity(allocator: std.mem.Allocator, dir_path: []const u8, strict: bool) !void {
    var dir = try std.Io.Dir.cwd().openDir(io, dir_path, .{ .iterate = true });
    defer dir.close(io);
    var iter = dir.iterate();
    var checked: usize = 0;
    while (try iter.next(io)) |entry| {
        if (entry.kind != .file) continue;
        if (!std.mem.endsWith(u8, entry.name, ".fasta") and !std.mem.endsWith(u8, entry.name, ".fa")) continue;

        const path = try std.fmt.allocPrint(allocator, "{s}/{s}", .{ dir_path, entry.name });
        defer allocator.free(path);

        const data = try readTestFile(allocator, path);
        defer allocator.free(data);
        // Non-FASTA bytes under a .fasta name (e.g. empty file) are not parity targets.
        if (!validateFasta(data)) continue;

        var probe_arena = std.heap.ArenaAllocator.init(allocator);
        defer probe_arena.deinit();
        var probe_index = main.indexer.scanZfiIndex(data, true, probe_arena.allocator()) catch |err| {
            if (strict) {
                std.debug.print("required fixture failed to index: {s}\n", .{path});
                return err;
            }
            // Optional generated edge cases intentionally include unindexable inputs.
            continue;
        };
        probe_index.deinit(probe_arena.allocator());

        try expectZfiStreamingMatchesMmap(allocator, data, path);
        checked += 1;
    }
    if (strict and checked == 0) {
        std.debug.print("required fixture dir produced no parity checks: {s}\n", .{dir_path});
        return error.RequiredFixtureMissing;
    }
}

// ============================================================================
// validateFasta tests
// ============================================================================

test "validateFasta returns true for valid FASTA" {
    try std.testing.expect(validateFasta(">seq1\nACGT\n"));
}

test "validateFasta returns true for just header" {
    try std.testing.expect(validateFasta(">"));
}

test "validateFasta returns false for empty data" {
    try std.testing.expect(!validateFasta(""));
}

test "validateFasta returns false for non-FASTA" {
    try std.testing.expect(!validateFasta("not fasta"));
    try std.testing.expect(!validateFasta("ACGT"));
    try std.testing.expect(!validateFasta("\n>seq"));
    try std.testing.expect(!validateFasta(" >seq"));
}

test "validateFasta on real file - valid" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const data = try readTestFile(arena.allocator(), "tests/data/simple.fasta");
    try std.testing.expect(validateFasta(data));
}

test "validateFasta on real file - invalid" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const data = try readTestFile(arena.allocator(), "tests/data/not_fasta.txt");
    try std.testing.expect(!validateFasta(data));
}

// ============================================================================
// scanHeaders tests
// ============================================================================

test "scanHeaders basic" {
    const data = ">seq1\nACGT\n";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();

    const records = try scanHeaders(data, arena.allocator());
    try std.testing.expectEqual(@as(usize, 1), records.items.len);

    const rec = records.items[0];
    try std.testing.expectEqualStrings("seq1", rec.getName(data));
    try std.testing.expectEqual(@as(u64, 4), rec.seq_len);
    try std.testing.expectEqual(@as(u64, 6), rec.seq_offset);
    try std.testing.expectEqual(@as(u32, 4), rec.line_bases);
    try std.testing.expectEqual(@as(u32, 5), rec.line_bytes);
}

test "scanHeadersAt stores file-global offsets" {
    const prefix = "padding";
    const fasta = ">seq1\nACGT\n";
    const file_base: u64 = @intCast(prefix.len);

    var combined: std.ArrayList(u8) = .empty;
    defer combined.deinit(std.testing.allocator);
    try combined.appendSlice(std.testing.allocator, prefix);
    try combined.appendSlice(std.testing.allocator, fasta);

    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();

    const records = try main.indexer.scanHeadersAt(fasta, file_base, arena.allocator());
    try std.testing.expectEqual(@as(usize, 1), records.items.len);

    const rec = records.items[0];
    try std.testing.expectEqual(@as(u64, file_base + 1), rec.name_offset);
    try std.testing.expectEqual(@as(u64, file_base + 6), rec.seq_offset);
    try std.testing.expectEqualStrings("seq1", rec.getName(combined.items));
}

test "scanHeaders multiple sequences" {
    const data = ">seq1\nACGT\n>seq2\nGGGG\n";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();

    const records = try scanHeaders(data, arena.allocator());
    try std.testing.expectEqual(@as(usize, 2), records.items.len);

    try std.testing.expectEqualStrings("seq1", records.items[0].getName(data));
    try std.testing.expectEqualStrings("seq2", records.items[1].getName(data));
}

test "scanHeaders with header description" {
    const data = ">myseq some description here\nACGT\n";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();

    const records = try scanHeaders(data, arena.allocator());
    try std.testing.expectEqual(@as(usize, 1), records.items.len);
    try std.testing.expectEqualStrings("myseq", records.items[0].getName(data));
}

test "scanHeaders with pipe-delimited header (proteome style)" {
    const data = ">sp|P12345|PROT_NAME OS=Human\nACGT\n";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();

    const records = try scanHeaders(data, arena.allocator());
    try std.testing.expectEqual(@as(usize, 1), records.items.len);
    try std.testing.expectEqualStrings("sp|P12345|PROT_NAME", records.items[0].getName(data));
}

test "scanHeaders multiline sequence" {
    const data = ">seq1\nAAAA\nBBBB\nCCCC\n";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();

    const records = try scanHeaders(data, arena.allocator());
    try std.testing.expectEqual(@as(usize, 1), records.items.len);

    const rec = records.items[0];
    try std.testing.expectEqual(@as(u64, 12), rec.seq_len);
    try std.testing.expectEqual(@as(u32, 4), rec.line_bases);
    try std.testing.expectEqual(@as(u32, 5), rec.line_bytes);
}

test "scanHeaders counts wrapped final short line with trailing newline correctly" {
    const data = ">chrSynthetic\nAAAAAAAA\nCCCCCCCC\nGGGG\n";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();

    const records = try scanHeaders(data, arena.allocator());
    try std.testing.expectEqual(@as(usize, 1), records.items.len);

    const rec = records.items[0];
    try std.testing.expectEqual(@as(u64, 20), rec.seq_len);
    try std.testing.expectEqual(@as(u32, 8), rec.line_bases);
    try std.testing.expectEqual(@as(u32, 9), rec.line_bytes);
}

test "scanHeaders counts trailing whitespace lines via base scan fallback" {
    const data = ">seq1\nACGTACGT  \nACGT\n";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();

    const records = try scanHeaders(data, arena.allocator());
    try std.testing.expectEqual(@as(usize, 1), records.items.len);
    try std.testing.expectEqual(@as(u64, 12), records.items[0].seq_len);
}

test "scanHeaders counts record with trailing newline before next header" {
    const data = ">a\nAAAA\n>next\n";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();

    const records = try scanHeaders(data, arena.allocator());
    try std.testing.expectEqual(@as(usize, 1), records.items.len);
    try std.testing.expectEqual(@as(u64, 4), records.items[0].seq_len);
}

test "scanHeaders cross-checks non-dense first line against countBases" {
    const data = ">seq\nACGT    \nACGTACGT\n";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();

    const records = try scanHeaders(data, arena.allocator());
    try std.testing.expectEqual(@as(usize, 1), records.items.len);
    try std.testing.expectEqual(@as(u64, 12), records.items[0].seq_len);
}

test "scanHeaders handles CRLF line endings" {
    const data = ">seq1\r\nACGT\r\n";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();

    const records = try scanHeaders(data, arena.allocator());
    try std.testing.expectEqual(@as(usize, 1), records.items.len);

    const rec = records.items[0];
    try std.testing.expectEqual(@as(u64, 4), rec.seq_len);
    try std.testing.expectEqual(@as(u32, 4), rec.line_bases);
    try std.testing.expectEqual(@as(u32, 6), rec.line_bytes);
}

test "scanHeaders skips empty sequences (matches samtools)" {
    const data = ">empty\n>next\nACGT\n";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();

    const records = try scanHeaders(data, arena.allocator());
    try std.testing.expectEqual(@as(usize, 1), records.items.len);
    try std.testing.expectEqualStrings("next", records.items[0].getName(data));
}

test "scanHeaders skips duplicate names (matches samtools)" {
    const data = ">dup\nAAAA\n>other\nCCCC\n>dup\nGGGG\n";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();

    const records = try scanHeaders(data, arena.allocator());
    try std.testing.expectEqual(@as(usize, 2), records.items.len);
    try std.testing.expectEqualStrings("dup", records.items[0].getName(data));
    try std.testing.expectEqualStrings("other", records.items[1].getName(data));
}

test "NameDedup identity is string equality" {
    var seen = main.indexer.NameDedup.init(std.testing.allocator, true);
    defer seen.deinit();

    try std.testing.expect(!(try seen.observe("alpha")));
    try std.testing.expect(try seen.observe("alpha"));
    try std.testing.expect(!(try seen.observe("beta")));
    try std.testing.expect(try seen.observe("beta"));
}

test "mmap and streaming dedup agree and keep distinct names" {
    // First-wins on exact duplicate `dup`; both paths must keep `other` and must not
    // drop a distinct name solely because a hash collided (StringHashMap compares bytes).
    const data = ">dup\nAAAA\n>other\nCCCC\n>dup\nGGGG\n>also\nTTTT\n";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var mmap_fai = std.Io.Writer.Allocating.init(allocator);
    const mmap_count = try main.indexer.streamingScan(data, &mmap_fai.writer, .fai, true, allocator);

    var stream_fai = std.Io.Writer.Allocating.init(allocator);
    const stream_count = try scanChunkedData(data, &stream_fai.writer, true, allocator);

    try std.testing.expectEqual(@as(u32, 3), mmap_count);
    try std.testing.expectEqual(mmap_count, stream_count);
    try std.testing.expectEqualStrings(mmap_fai.written(), stream_fai.written());
    try std.testing.expect(std.mem.indexOf(u8, mmap_fai.written(), "other\t") != null);
    try std.testing.expect(std.mem.indexOf(u8, mmap_fai.written(), "also\t") != null);
}

test "scanZfiIndexStreaming matches mmap with dedup enabled" {
    const data = @embedFile("data/edge_cases.fasta");
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var mmap_index = try main.indexer.scanZfiIndex(data, true, allocator);
    defer mmap_index.deinit(allocator);
    var stream_index = try main.indexer.scanZfiIndexStreamingData(data, true, allocator);
    defer stream_index.deinit(allocator);

    const mmap_bytes = try main.indexer.zfiIndexToBytes(&mmap_index, data.len, 0, allocator);
    const stream_bytes = try main.indexer.zfiIndexToBytes(&stream_index, data.len, 0, allocator);
    try std.testing.expectEqualSlices(u8, mmap_bytes, stream_bytes);
}

test "scanHeaders on real simple.fasta file" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const data = try readTestFile(arena.allocator(), "tests/data/simple.fasta");

    const records = try scanHeaders(data, arena.allocator());
    try std.testing.expectEqual(@as(usize, 2), records.items.len);

    try std.testing.expectEqualStrings("seq1", records.items[0].getName(data));
    try std.testing.expectEqual(@as(u64, 24), records.items[0].seq_len);
}

test "scanHeaders on real proteome.fasta file" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const data = try readTestFile(arena.allocator(), "tests/data/proteome.fasta");

    const records = try scanHeaders(data, arena.allocator());
    try std.testing.expect(records.items.len > 0);

    // Proteome headers have pipes
    const name = records.items[0].getName(data);
    try std.testing.expect(std.mem.indexOf(u8, name, "|") != null);
}

// ============================================================================
// writeZfiIndex tests
// ============================================================================

test "writeZfiIndexFile creates valid production layout" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const records = [_]IndexRecord{
        .{ .name_offset = 1, .name_len = 4, .seq_offset = 10, .seq_len = 100, .line_bases = 80, .line_bytes = 81 },
    };

    const path = "tests/data/test_write.zfi";
    defer std.Io.Dir.cwd().deleteFile(io, path) catch {};

    try writeZfiFromRecords(path, &records, 1000, 0, allocator);

    const file = try std.Io.Dir.cwd().openFile(io, path, .{});
    defer file.close(io);
    const stat = try file.stat(io);

    var header_bytes: [@sizeOf(ZfiHeader)]u8 = undefined;
    _ = try std.Io.File.readPositionalAll(file, io, &header_bytes, 0);
    const header: ZfiHeader = @bitCast(header_bytes);

    try std.testing.expectEqualSlices(u8, &ZFI_MAGIC, &header.magic);
    try std.testing.expectEqual(@as(u32, 1), header.record_count);
    try std.testing.expectEqual(@as(u64, 1000), header.source_size);
    // header + one record + empty side/name regions + 12-byte footer + 12-byte source id
    try std.testing.expectEqual(
        @as(u64, 16 + 40 + main.index_format.zfi_name_footer_bytes + main.index_format.zfi_source_id_bytes),
        stat.size,
    );
}

test "ZfiHeader has correct size" {
    try std.testing.expectEqual(@as(usize, 16), @sizeOf(ZfiHeader));
}

test "IndexRecord has consistent size" {
    // name_offset(8) + name_len(2) + padding(6) + seq_offset(8) + seq_len(8) + line_bases(4) + line_bytes(4) = 40
    try std.testing.expectEqual(@as(usize, 40), @sizeOf(IndexRecord));
}

test "SideTableLine has explicit v0.3 size" {
    try std.testing.expectEqual(@as(usize, 32), @sizeOf(main.index_format.SideTableLine));
}

test "zfi wire encoders match asBytes for frozen little-endian layout" {
    const header = ZfiHeader{
        .magic = ZFI_MAGIC,
        .record_count = 0x01020304,
        .source_size = 0x0807060504030201,
    };
    try std.testing.expectEqualSlices(u8, std.mem.asBytes(&header), &main.index_format.encodeZfiHeader(header));

    var rec = IndexRecord{
        .name_offset = 0x1112131415161718,
        .name_len = 0x2122,
        .seq_offset = 0x3132333435363738,
        .seq_len = 0x4142434445464748,
        .line_bases = 0x51525354,
        .line_bytes = 0x61626364,
    };
    try rec.markNonUniform(0x000000AABBCCDDEE);
    try std.testing.expectEqualSlices(u8, std.mem.asBytes(&rec), &main.index_format.encodeIndexRecord(rec));
    try std.testing.expectEqual(@as(u8, 0xEE), std.mem.asBytes(&rec)[11]);
    try std.testing.expectEqual(@as(u64, 0x000000AABBCCDDEE), rec.sideTableOffset());

    const line = main.index_format.SideTableLine{
        .base_start = 1,
        .byte_offset = 0x90A0B0C0D0E0F000,
        .line_bytes = 10,
        .line_bases = 8,
    };
    try std.testing.expectEqualSlices(u8, std.mem.asBytes(&line), &main.index_format.encodeSideTableLine(line));

    const footer = main.index_format.encodeZfiNameFooter(0x0102030405060708);
    try std.testing.expectEqualSlices(u8, &main.index_format.ZFI_NAME_FOOTER_MAGIC, footer[0..4]);
    try std.testing.expectEqual(@as(u64, 0x0102030405060708), std.mem.readInt(u64, footer[4..12], .little));
}

test "zfiIndexToBytes matches encodeZfiHeader prefix and LE side-table count" {
    const data = ">seq\nAAA\nCCCC\nGG\n";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var zfi_index = try main.indexer.scanZfiIndex(data, true, allocator);
    defer zfi_index.deinit(allocator);
    const bytes = try main.indexer.zfiIndexToBytes(&zfi_index, data.len, 0, allocator);
    defer allocator.free(bytes);

    const expected_header = main.index_format.encodeZfiHeader(.{
        .magic = ZFI_MAGIC,
        .record_count = 1,
        .source_size = data.len,
    });
    try std.testing.expectEqualSlices(u8, &expected_header, bytes[0..16]);

    const rec = zfi_index.records.items[0];
    try std.testing.expectEqualSlices(
        u8,
        &main.index_format.encodeIndexRecord(rec),
        bytes[16 .. 16 + 40],
    );

    const side_off: usize = @intCast(rec.sideTableOffset());
    try std.testing.expectEqual(
        @as(u64, 3),
        std.mem.readInt(u64, bytes[side_off..][0..8], .little),
    );
}

test "scanZfiIndex marks non-uniform records and appends side table" {
    const data = ">seq\nAAA\nCCCC\nGG\n";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var zfi_index = try main.indexer.scanZfiIndex(data, true, allocator);
    defer zfi_index.deinit(allocator);

    try std.testing.expectEqual(@as(usize, 1), zfi_index.records.items.len);
    const rec = zfi_index.records.items[0];
    try std.testing.expect(!rec.isUniformWidth());
    try std.testing.expectEqual(
        @as(u64, @sizeOf(ZfiHeader) + @sizeOf(IndexRecord)),
        rec.sideTableOffset(),
    );

    try std.testing.expectEqual(
        @as(u64, 3),
        std.mem.readInt(u64, zfi_index.side_tables.items[0..@sizeOf(u64)], .little),
    );
}

test "scanZfiIndexStreaming zfi loads record names" {
    const data = @embedFile("data/simple.fasta");
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var stream_index = try main.indexer.scanZfiIndexStreamingData(data, true, allocator);
    defer stream_index.deinit(allocator);
    const zfi_bytes = try main.indexer.zfiIndexToBytes(&stream_index, data.len, 0, allocator);
    defer allocator.free(zfi_bytes);

    const paths = try writeFastaAndRawZfi(allocator, "stream-zfi", data, zfi_bytes);
    defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, paths.zfi_path) catch {};

    var idx = try loadIndexChecked(io, paths.fasta_path);
    defer idx.deinit(io);
    try std.testing.expectEqual(@as(?usize, 0), idx.lookupName("seq1"));
    try std.testing.expectEqual(@as(?usize, 1), idx.lookupName("seq2"));
}

test "scanZfiIndexStreaming matches mmap on tests/data fixtures" {
    try requireFixturePath("tests/data");
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    try walkFastaDirZfiParity(arena.allocator(), "tests/data", true);
}

test "scanZfiIndexStreaming matches mmap on edge cases (optional generated fixtures)" {
    // Generated by `bash bench/index/run.sh`; gitignored. Not required for CI.
    std.Io.Dir.cwd().access(io, "bench/index/edge_cases", .{}) catch |err| switch (err) {
        error.FileNotFound => return skipOptionalFixture(
            "bench/index/edge_cases",
            "bash bench/index/run.sh (generates edge_cases/)",
        ),
        else => |e| return e,
    };

    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    try walkFastaDirZfiParity(arena.allocator(), "bench/index/edge_cases", false);
}

test "scanZfiIndexStreaming matches mmap on messy variants" {
    try requireFixturePath("bench/shared/cache/messy_fixtures");
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    try walkFastaDirZfiParity(arena.allocator(), "bench/shared/cache/messy_fixtures", true);
}

test "scanZfiIndexStreaming matches mmap on REAL references (optional large data)" {
    const refs = [_][]const u8{
        "bench/shared/data/REAL_Genome.fa",
        "bench/shared/data/REAL_Transcriptome.fa",
        "bench/shared/data/REAL_Proteome.fasta",
    };
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var ran: usize = 0;
    for (refs) |path| {
        std.Io.Dir.cwd().access(io, path, .{}) catch |err| switch (err) {
            error.FileNotFound => {
                std.debug.print(
                    "optional large fixture absent, skipping file: {s}\n  obtain via: bash bench/shared/download_data.sh\n",
                    .{path},
                );
                continue;
            },
            else => |e| return e,
        };
        const data = try readTestFile(allocator, path);
        defer allocator.free(data);
        try expectZfiStreamingMatchesMmap(allocator, data, path);
        ran += 1;
    }
    if (ran == 0) {
        return skipOptionalFixture(
            "bench/shared/data/REAL_*.fa(sta)",
            "bash bench/shared/download_data.sh",
        );
    }
}

test "scanZfiIndexStreaming matches mmap across chunk boundaries" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var fasta: std.ArrayList(u8) = .empty;
    defer fasta.deinit(allocator);
    try fasta.appendSlice(allocator, ">seq1\n");
    try fasta.appendNTimes(allocator, 'A', low_mem_chunk_size + 10);
    try fasta.appendSlice(allocator, "\n>seq2\n");
    try fasta.appendNTimes(allocator, 'C', low_mem_chunk_size);
    try fasta.appendSlice(allocator, "GGTT\n");

    try expectZfiStreamingMatchesMmap(allocator, fasta.items, "synthetic-chunk-boundary");
}

test "scanZfiIndexStreaming matches mmap with dedup disabled" {
    const data = @embedFile("data/edge_cases.fasta");
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var mmap_index = try main.indexer.scanZfiIndex(data, false, allocator);
    defer mmap_index.deinit(allocator);
    var stream_index = try main.indexer.scanZfiIndexStreamingData(data, false, allocator);
    defer stream_index.deinit(allocator);

    const mmap_bytes = try main.indexer.zfiIndexToBytes(&mmap_index, data.len, 0, allocator);
    const stream_bytes = try main.indexer.zfiIndexToBytes(&stream_index, data.len, 0, allocator);
    try std.testing.expectEqualSlices(u8, mmap_bytes, stream_bytes);
}

test "low-mem indexing matches default across chunk boundary" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var fasta: std.ArrayList(u8) = .empty;
    defer fasta.deinit(allocator);
    try fasta.appendSlice(allocator, ">seq1\n");
    try fasta.appendNTimes(allocator, 'A', low_mem_chunk_size + 10);
    try fasta.appendSlice(allocator, "\n>seq2\nACGT\n");

    var expected = std.Io.Writer.Allocating.init(allocator);
    const expected_count = try main.indexer.streamingScan(fasta.items, &expected.writer, .fai, true, allocator);

    var actual = std.Io.Writer.Allocating.init(allocator);
    const actual_count = try scanChunkedData(fasta.items, &actual.writer, true, allocator);

    try std.testing.expectEqual(expected_count, actual_count);
    try std.testing.expectEqualStrings(expected.written(), actual.written());
}

test "low-mem indexing with header at chunk boundary" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const chunk = low_mem_chunk_size;
    var fasta: std.ArrayList(u8) = .empty;
    defer fasta.deinit(allocator);
    try fasta.appendSlice(allocator, ">seq1\n");
    try fasta.appendNTimes(allocator, 'A', chunk - 8);
    try fasta.appendSlice(allocator, "\n>seq2\nACGT\n");

    var expected = std.Io.Writer.Allocating.init(allocator);
    _ = try main.indexer.streamingScan(fasta.items, &expected.writer, .fai, true, allocator);

    var actual = std.Io.Writer.Allocating.init(allocator);
    _ = try scanChunkedData(fasta.items, &actual.writer, true, allocator);

    try std.testing.expectEqualStrings(expected.written(), actual.written());
}

test "low-mem indexing with sequence byte split across chunk boundary" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var fasta: std.ArrayList(u8) = .empty;
    defer fasta.deinit(allocator);
    try fasta.appendSlice(allocator, ">seq1\n");
    try fasta.appendNTimes(allocator, 'A', low_mem_chunk_size);
    try fasta.appendSlice(allocator, "ACGT\n");

    var expected = std.Io.Writer.Allocating.init(allocator);
    _ = try main.indexer.streamingScan(fasta.items, &expected.writer, .fai, true, allocator);

    var actual = std.Io.Writer.Allocating.init(allocator);
    _ = try scanChunkedData(fasta.items, &actual.writer, true, allocator);

    try std.testing.expectEqualStrings(expected.written(), actual.written());
}

test "streamingScan counts wrapped final short line with trailing newline correctly" {
    const data = ">chrSynthetic\nAAAAAAAA\nCCCCCCCC\nGGGG\n";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var out = std.Io.Writer.Allocating.init(allocator);
    const record_count = try main.indexer.streamingScan(data, &out.writer, .fai, true, allocator);

    try std.testing.expectEqual(@as(u32, 1), record_count);
    try std.testing.expectEqualStrings("chrSynthetic\t20\t14\t8\t9\n", out.written());
}

test "low-mem indexing rejects sequence names longer than u16" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var fasta: std.ArrayList(u8) = .empty;
    defer fasta.deinit(allocator);
    try fasta.append(allocator, '>');
    try fasta.appendNTimes(allocator, 'A', main.indexer.max_index_name_len + 1);
    try fasta.appendSlice(allocator, "\nACGT\n");

    var output = std.Io.Writer.Allocating.init(allocator);

    try std.testing.expectError(
        error.HeaderTooLong,
        scanChunkedData(fasta.items, &output.writer, true, allocator),
    );
}

test "mmap indexing rejects sequence names longer than u16" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var fasta: std.ArrayList(u8) = .empty;
    defer fasta.deinit(allocator);
    try fasta.append(allocator, '>');
    try fasta.appendNTimes(allocator, 'A', main.indexer.max_index_name_len + 1);
    try fasta.appendSlice(allocator, "\nACGT\n");

    var output = std.Io.Writer.Allocating.init(allocator);
    try std.testing.expectError(
        error.HeaderTooLong,
        main.indexer.streamingScan(fasta.items, &output.writer, .fai, true, allocator),
    );
}

test "mmap and streaming FAI agree with leading blank line after header" {
    const data = ">seq\n\nACGT\n";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var expected = std.Io.Writer.Allocating.init(allocator);
    defer expected.deinit();
    var actual = std.Io.Writer.Allocating.init(allocator);
    defer actual.deinit();

    _ = try main.indexer.streamingScan(data, &expected.writer, .fai, true, allocator);
    _ = try scanChunkedData(data, &actual.writer, true, allocator);
    try std.testing.expectEqualStrings(expected.written(), actual.written());

    var records = try scanHeaders(data, allocator);
    defer records.deinit(allocator);
    try std.testing.expectEqual(@as(usize, 1), records.items.len);
    // First base-bearing line starts after the blank line (`>seq\n\n` = 6 bytes).
    try std.testing.expectEqual(@as(u64, 6), records.items[0].seq_offset);
    try std.testing.expectEqual(@as(u32, 4), records.items[0].line_bases);
    try std.testing.expectEqual(@as(u32, 5), records.items[0].line_bytes);
}

test "mmap and streaming FAI agree when final line omits terminal newline" {
    const data = ">seq\nAAAA\nBB";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var expected = std.Io.Writer.Allocating.init(allocator);
    defer expected.deinit();
    var actual = std.Io.Writer.Allocating.init(allocator);
    defer actual.deinit();

    _ = try main.indexer.streamingScan(data, &expected.writer, .fai, true, allocator);
    _ = try scanChunkedData(data, &actual.writer, true, allocator);
    try std.testing.expectEqualStrings(expected.written(), actual.written());
}

test "mmap and streaming ZFI agree when blank lines sit between sequence lines" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const data = ">seq\nAAAA\n\nBBBB\n";
    try expectZfiStreamingMatchesMmap(allocator, data, "internal-blank-lines");

    var index = try main.indexer.scanZfiIndex(data, true, allocator);
    defer index.deinit(allocator);
    try std.testing.expectEqual(@as(usize, 1), index.records.items.len);
    try std.testing.expect(!index.records.items[0].isUniformWidth());
}

test "mmap and streaming skip empty records equivalently" {
    const data = ">empty\n>next\nACGT\n";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var expected = std.Io.Writer.Allocating.init(allocator);
    defer expected.deinit();
    var actual = std.Io.Writer.Allocating.init(allocator);
    defer actual.deinit();

    _ = try main.indexer.streamingScan(data, &expected.writer, .fai, true, allocator);
    _ = try scanChunkedData(data, &actual.writer, true, allocator);
    try std.testing.expectEqualStrings(expected.written(), actual.written());
    try expectZfiStreamingMatchesMmap(allocator, data, "empty-record-skip");
}

test "low-mem FAI matches mmap for long sequence names" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var fasta: std.ArrayList(u8) = .empty;
    defer fasta.deinit(allocator);
    try fasta.append(allocator, '>');
    try fasta.appendNTimes(allocator, 'A', 10_000);
    try fasta.appendSlice(allocator, " description\nACGT\n");

    var expected = std.Io.Writer.Allocating.init(allocator);
    var actual = std.Io.Writer.Allocating.init(allocator);
    _ = try main.indexer.streamingScan(fasta.items, &expected.writer, .fai, true, allocator);
    _ = try scanChunkedData(fasta.items, &actual.writer, true, allocator);
    try std.testing.expectEqualStrings(expected.written(), actual.written());
}

test "loadIndexChecked rejects corrupt zfi records" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_path = try uniqueArtifactPath(allocator, "corrupt-zfi", "fa");
    const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{fasta_path});
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, zfi_path) catch {};

    const fasta_file = try std.Io.Dir.cwd().createFile(io, fasta_path, .{});
    defer fasta_file.close(io);
    try std.Io.File.writeStreamingAll(fasta_file, io, ">seq1\nACGT\n");

    const bad_records = [_]IndexRecord{
        .{ .name_offset = 999_999, .name_len = 4, .seq_offset = 6, .seq_len = 4, .line_bases = 4, .line_bytes = 5 },
    };
    try writeZfiFromRecords(zfi_path, &bad_records, 11, try statMtimeNs(fasta_path), allocator);

    try std.testing.expectError(error.CorruptIndex, loadIndexChecked(io, fasta_path));
}

test "loadIndexChecked rejects seq_offset at end of file" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_path = try uniqueArtifactPath(allocator, "past-end-zfi", "fa");
    const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{fasta_path});
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, zfi_path) catch {};

    const fasta_data = ">seq1\nACGT\n";
    const fasta_file = try std.Io.Dir.cwd().createFile(io, fasta_path, .{});
    defer fasta_file.close(io);
    try std.Io.File.writeStreamingAll(fasta_file, io, fasta_data);

    const bad_records = [_]IndexRecord{
        .{ .name_offset = 1, .name_len = 4, .seq_offset = fasta_data.len, .seq_len = 4, .line_bases = 4, .line_bytes = 5 },
    };
    try writeZfiFromRecords(zfi_path, &bad_records, fasta_data.len, try statMtimeNs(fasta_path), allocator);

    try std.testing.expectError(error.CorruptIndex, loadIndexChecked(io, fasta_path));
}

fn appendZfiHeader(buf: *std.ArrayListUnmanaged(u8), allocator: std.mem.Allocator, record_count: u32, source_size: u64) !void {
    const header = ZfiHeader{
        .magic = ZFI_MAGIC,
        .record_count = record_count,
        .source_size = source_size,
    };
    try buf.appendSlice(allocator, std.mem.asBytes(&header));
}

test "loadIndexChecked rejects zfi with zero record_count" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var zfi_buf: std.ArrayListUnmanaged(u8) = .empty;
    defer zfi_buf.deinit(allocator);
    try appendZfiHeader(&zfi_buf, allocator, 0, 11);

    const paths = try writeFastaAndRawZfi(allocator, "zfi-zero-count", ">seq1\nACGT\n", zfi_buf.items);
    defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, paths.zfi_path) catch {};

    try std.testing.expectError(error.CorruptIndex, loadIndexChecked(io, paths.fasta_path));
}

test "loadIndexChecked rejects truncated zfi header" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const paths = try writeFastaAndRawZfi(allocator, "zfi-trunc-header", ">seq1\nACGT\n", ZFI_MAGIC[0..3]);
    defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, paths.zfi_path) catch {};

    try std.testing.expectError(error.CorruptIndex, loadIndexChecked(io, paths.fasta_path));
}

test "loadIndexChecked rejects zfi record_count larger than file" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var zfi_buf: std.ArrayListUnmanaged(u8) = .empty;
    defer zfi_buf.deinit(allocator);
    try appendZfiHeader(&zfi_buf, allocator, 2, 11);
    // Only one record body follows the header.
    const rec = IndexRecord{ .name_offset = 1, .name_len = 4, .seq_offset = 6, .seq_len = 4, .line_bases = 4, .line_bytes = 5 };
    try zfi_buf.appendSlice(allocator, std.mem.asBytes(&rec));

    const paths = try writeFastaAndRawZfi(allocator, "zfi-count-overflow", ">seq1\nACGT\n", zfi_buf.items);
    defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, paths.zfi_path) catch {};

    try std.testing.expectError(error.CorruptIndex, loadIndexChecked(io, paths.fasta_path));
}

test "loadIndexChecked rejects zfi with bad magic" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var zfi_buf: std.ArrayListUnmanaged(u8) = .empty;
    defer zfi_buf.deinit(allocator);
    var header = ZfiHeader{ .magic = .{ 'X', 'F', 'I', 0x01 }, .record_count = 1, .source_size = 11 };
    try zfi_buf.appendSlice(allocator, std.mem.asBytes(&header));
    const rec = IndexRecord{ .name_offset = 1, .name_len = 4, .seq_offset = 6, .seq_len = 4, .line_bases = 4, .line_bytes = 5 };
    try zfi_buf.appendSlice(allocator, std.mem.asBytes(&rec));

    const paths = try writeFastaAndRawZfi(allocator, "zfi-bad-magic", ">seq1\nACGT\n", zfi_buf.items);
    defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, paths.zfi_path) catch {};

    try std.testing.expectError(error.CorruptIndex, loadIndexChecked(io, paths.fasta_path));
}

test "loadIndexChecked rejects zfi name blob overlapping records" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_data = ">seq1\nACGT\n";
    var index = try main.indexer.scanZfiIndex(fasta_data, true, allocator);
    defer index.deinit(allocator);
    var zfi_bytes = try main.indexer.zfiIndexToBytes(&index, fasta_data.len, 0, allocator);
    defer allocator.free(zfi_bytes);

    // Inflate footer blob length so blob_start falls inside the record array.
    // Footer is at EOF; identity sits immediately before it.
    const footer = main.index_format.encodeZfiNameFooter(zfi_bytes.len);
    @memcpy(zfi_bytes[zfi_bytes.len - footer.len ..], &footer);

    const paths = try writeFastaAndRawZfi(allocator, "zfi-blob-overlap", fasta_data, zfi_bytes);
    defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, paths.zfi_path) catch {};

    try std.testing.expectError(error.CorruptIndex, loadIndexChecked(io, paths.fasta_path));
}

test "loadIndexChecked rejects zfi side table offset into records" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_data = ">seq\nAAA\nCCCC\nGG\n";
    var index = try main.indexer.scanZfiIndex(fasta_data, true, allocator);
    defer index.deinit(allocator);
    try std.testing.expect(!index.records.items[0].isUniformWidth());

    // Point the side table into the header/records region.
    try index.records.items[0].markNonUniform(@sizeOf(ZfiHeader));

    const zfi_bytes = try main.indexer.zfiIndexToBytes(&index, fasta_data.len, 0, allocator);
    defer allocator.free(zfi_bytes);

    const paths = try writeFastaAndRawZfi(allocator, "zfi-side-in-records", fasta_data, zfi_bytes);
    defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, paths.zfi_path) catch {};

    try std.testing.expectError(error.CorruptIndex, loadIndexChecked(io, paths.fasta_path));
}

test "loadIndexChecked rejects zfi side table with empty line_count" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_data = ">seq\nAAA\nCCCC\nGG\n";
    var index = try main.indexer.scanZfiIndex(fasta_data, true, allocator);
    defer index.deinit(allocator);

    const offset: usize = @intCast(index.records.items[0].sideTableOffset());
    const zfi_bytes = try main.indexer.zfiIndexToBytes(&index, fasta_data.len, 0, allocator);
    defer allocator.free(zfi_bytes);
    // Zero the line_count header at the side-table offset.
    @memset(zfi_bytes[offset .. offset + @sizeOf(u64)], 0);

    const paths = try writeFastaAndRawZfi(allocator, "zfi-side-zero-lines", fasta_data, zfi_bytes);
    defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, paths.zfi_path) catch {};

    try std.testing.expectError(error.CorruptIndex, loadIndexChecked(io, paths.fasta_path));
}

test "loadIndexChecked rejects mixed nameInZfi flags" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_data = ">a\nAAAA\n>b\nBBBB\n";
    var index = try main.indexer.scanZfiIndex(fasta_data, true, allocator);
    defer index.deinit(allocator);
    try std.testing.expect(index.records.items.len >= 2);
    try std.testing.expect(index.records.items[0].nameInZfi());

    // Clear NAME_IN_ZFI on the second record only.
    index.records.items[1]._pad[0] &= ~@as(u8, 2);

    const zfi_bytes = try main.indexer.zfiIndexToBytes(&index, fasta_data.len, 0, allocator);
    defer allocator.free(zfi_bytes);

    const paths = try writeFastaAndRawZfi(allocator, "zfi-mixed-name-flags", fasta_data, zfi_bytes);
    defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, paths.zfi_path) catch {};

    try std.testing.expectError(error.CorruptIndex, loadIndexChecked(io, paths.fasta_path));
}

test "loadIndexChecked accepts production zfi with side tables and name blob" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_data = ">seq\nAAA\nCCCC\nGG\n";
    var index = try main.indexer.scanZfiIndex(fasta_data, true, allocator);
    defer index.deinit(allocator);
    const zfi_bytes = try main.indexer.zfiIndexToBytes(&index, fasta_data.len, 0, allocator);
    defer allocator.free(zfi_bytes);

    const paths = try writeFastaAndRawZfi(allocator, "zfi-production-ok", fasta_data, zfi_bytes);
    defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, paths.zfi_path) catch {};

    var idx = try loadIndexChecked(io, paths.fasta_path);
    defer idx.deinit(io);
    try std.testing.expectEqual(main.index_format.LoadedIndex.IndexSource.zfi, idx.source);
    try std.testing.expect(idx.name_blob != null);
    try std.testing.expectEqual(@as(usize, 1), idx.records.len);
    try std.testing.expectEqualStrings("seq", idx.getRecordName(0));
}

test "loadIndexCheckedWithMode preserves duplicate lookup semantics" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_path = try uniqueArtifactPath(allocator, "duplicate-zfi", "fa");
    const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{fasta_path});
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, zfi_path) catch {};

    const fasta_data = ">dup\nAAAA\n>dup\nCCCC\n";
    const fasta_file = try std.Io.Dir.cwd().createFile(io, fasta_path, .{});
    defer fasta_file.close(io);
    try std.Io.File.writeStreamingAll(fasta_file, io, fasta_data);

    const records = [_]IndexRecord{
        .{ .name_offset = 1, .name_len = 3, .seq_offset = 5, .seq_len = 4, .line_bases = 4, .line_bytes = 5 },
        .{ .name_offset = 11, .name_len = 3, .seq_offset = 15, .seq_len = 4, .line_bases = 4, .line_bytes = 5 },
    };
    try writeZfiFromRecords(zfi_path, &records, fasta_data.len, try statMtimeNs(fasta_path), allocator);

    var full_map_idx = try main.index_format.loadIndexCheckedWithMode(io, fasta_path, .lookup_full_map);
    defer full_map_idx.deinit(io);
    var records_only_idx = try main.index_format.loadIndexCheckedWithMode(io, fasta_path, .records_only);
    defer records_only_idx.deinit(io);

    try std.testing.expectEqual(full_map_idx.lookupName("dup"), records_only_idx.lookupName("dup"));
    try std.testing.expectEqual(@as(?usize, 1), records_only_idx.lookupName("dup"));
}

test "loadIndexChecked falls back to fai when zfi is stale" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_path = try uniqueArtifactPath(allocator, "stale-zfi", "fa");
    const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{fasta_path});
    const fai_path = try std.fmt.allocPrint(allocator, "{s}.fai", .{fasta_path});
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, zfi_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, fai_path) catch {};

    const fasta_data = ">seq1\nACGTACGT\n>seq2\nGGGG\n";

    const fasta_file = try std.Io.Dir.cwd().createFile(io, fasta_path, .{ .truncate = true });
    defer fasta_file.close(io);
    try std.Io.File.writeStreamingAll(fasta_file, io, fasta_data);

    var records = try scanHeaders(fasta_data, allocator);
    defer records.deinit(allocator);
    try writeZfiFromRecords(zfi_path, records.items, fasta_data.len, try statMtimeNs(fasta_path), allocator);

    const fai_file = try std.Io.Dir.cwd().createFile(io, fai_path, .{ .truncate = true });
    defer fai_file.close(io);
    var fai_buf: [65536]u8 = undefined;
    var fai_fw = fai_file.writer(io, &fai_buf);
    _ = try main.indexer.streamingScan(fasta_data, &fai_fw.interface, .fai, true, allocator);
    try fai_fw.flush();

    const zfi_file = try std.Io.Dir.cwd().openFile(io, zfi_path, .{});
    defer zfi_file.close(io);
    try markFileStaleOneHourAgo(zfi_file);

    var idx = try loadIndexChecked(io, fasta_path);
    defer idx.deinit(io);

    try std.testing.expectEqual(main.index_format.LoadedIndex.IndexSource.fai, idx.source);
    try std.testing.expectEqual(@as(?usize, 0), idx.lookupName("seq1"));
}

test "loadIndexCheckedWithMode preserves fai duplicate lookup semantics" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_path = try uniqueArtifactPath(allocator, "duplicate-fai", "fa");
    const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{fasta_path});
    const fai_path = try std.fmt.allocPrint(allocator, "{s}.fai", .{fasta_path});
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, zfi_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, fai_path) catch {};

    const fasta_data = ">dup\nAAAA\n>dup\nCCCC\n";
    const fasta_file = try std.Io.Dir.cwd().createFile(io, fasta_path, .{});
    defer fasta_file.close(io);
    try std.Io.File.writeStreamingAll(fasta_file, io, fasta_data);

    const fai_file = try std.Io.Dir.cwd().createFile(io, fai_path, .{ .truncate = true });
    defer fai_file.close(io);
    try std.Io.File.writeStreamingAll(fai_file, io, "dup\t4\t5\t4\t5\ndup\t4\t15\t4\t5\n");

    const zfi_file = try std.Io.Dir.cwd().createFile(io, zfi_path, .{ .truncate = true });
    defer zfi_file.close(io);
    try std.Io.File.writeStreamingAll(zfi_file, io, "stale");

    var full_map_idx = try main.index_format.loadIndexCheckedWithMode(io, fasta_path, .lookup_full_map);
    defer full_map_idx.deinit(io);
    var records_only_idx = try main.index_format.loadIndexCheckedWithMode(io, fasta_path, .records_only);
    defer records_only_idx.deinit(io);

    try std.testing.expectEqual(main.index_format.LoadedIndex.IndexSource.fai, records_only_idx.source);
    try std.testing.expect(!records_only_idx.has_name_map);
    try std.testing.expectEqual(@as(u16, 3), records_only_idx.records[0].name_len);
    try std.testing.expectEqualStrings("dup", records_only_idx.getRecordName(0));
    try std.testing.expectEqual(full_map_idx.lookupName("dup"), records_only_idx.lookupName("dup"));
    try std.testing.expectEqual(@as(?usize, 1), records_only_idx.lookupName("dup"));
}

test "loadIndexChecked rejects fai with zero line_bases" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_path = try uniqueArtifactPath(allocator, "bad-fai-zero-bases", "fa");
    const fai_path = try std.fmt.allocPrint(allocator, "{s}.fai", .{fasta_path});
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, fai_path) catch {};

    const fasta_file = try std.Io.Dir.cwd().createFile(io, fasta_path, .{});
    defer fasta_file.close(io);
    try std.Io.File.writeStreamingAll(fasta_file, io, ">seq1\nACGT\n");

    const fai_file = try std.Io.Dir.cwd().createFile(io, fai_path, .{ .truncate = true });
    defer fai_file.close(io);
    try std.Io.File.writeStreamingAll(fai_file, io, "seq1\t4\t6\t0\t5\n");

    try std.testing.expectError(error.CorruptIndex, loadIndexChecked(io, fasta_path));
}

test "loadIndexChecked rejects fai with line_bytes less than line_bases" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_path = try uniqueArtifactPath(allocator, "bad-fai-line-bytes", "fa");
    const fai_path = try std.fmt.allocPrint(allocator, "{s}.fai", .{fasta_path});
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, fai_path) catch {};

    const fasta_file = try std.Io.Dir.cwd().createFile(io, fasta_path, .{});
    defer fasta_file.close(io);
    try std.Io.File.writeStreamingAll(fasta_file, io, ">seq1\nACGT\n");

    const fai_file = try std.Io.Dir.cwd().createFile(io, fai_path, .{ .truncate = true });
    defer fai_file.close(io);
    try std.Io.File.writeStreamingAll(fai_file, io, "seq1\t4\t6\t5\t4\n");

    try std.testing.expectError(error.CorruptIndex, loadIndexChecked(io, fasta_path));
}

test "loadIndexChecked rejects fai with seq_offset past EOF" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_path = try uniqueArtifactPath(allocator, "bad-fai-offset", "fa");
    const fai_path = try std.fmt.allocPrint(allocator, "{s}.fai", .{fasta_path});
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, fai_path) catch {};

    const fasta_file = try std.Io.Dir.cwd().createFile(io, fasta_path, .{});
    defer fasta_file.close(io);
    try std.Io.File.writeStreamingAll(fasta_file, io, ">seq1\nACGT\n");

    const fai_file = try std.Io.Dir.cwd().createFile(io, fai_path, .{ .truncate = true });
    defer fai_file.close(io);
    try std.Io.File.writeStreamingAll(fai_file, io, "seq1\t4\t99\t4\t5\n");

    try std.testing.expectError(error.CorruptIndex, loadIndexChecked(io, fasta_path));
}

test "loadIndexChecked rejects fai whose sequence span exceeds FASTA" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_path = try uniqueArtifactPath(allocator, "bad-fai-span", "fa");
    const fai_path = try std.fmt.allocPrint(allocator, "{s}.fai", .{fasta_path});
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, fai_path) catch {};

    const fasta_file = try std.Io.Dir.cwd().createFile(io, fasta_path, .{});
    defer fasta_file.close(io);
    try std.Io.File.writeStreamingAll(fasta_file, io, ">seq1\nACGT\n");

    // seq_len 100 from offset 6 cannot fit in an 11-byte FASTA.
    const fai_file = try std.Io.Dir.cwd().createFile(io, fai_path, .{ .truncate = true });
    defer fai_file.close(io);
    try std.Io.File.writeStreamingAll(fai_file, io, "seq1\t100\t6\t4\t5\n");

    try std.testing.expectError(error.CorruptIndex, loadIndexChecked(io, fasta_path));
}

test "loadIndexCheckedWithMode records_only rejects impossible fai geometry" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_path = try uniqueArtifactPath(allocator, "bad-fai-records-only", "fa");
    const fai_path = try std.fmt.allocPrint(allocator, "{s}.fai", .{fasta_path});
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, fai_path) catch {};

    const fasta_file = try std.Io.Dir.cwd().createFile(io, fasta_path, .{});
    defer fasta_file.close(io);
    try std.Io.File.writeStreamingAll(fasta_file, io, ">seq1\nACGT\n");

    const fai_file = try std.Io.Dir.cwd().createFile(io, fai_path, .{ .truncate = true });
    defer fai_file.close(io);
    try std.Io.File.writeStreamingAll(fai_file, io, "seq1\t4\t6\t0\t5\n");

    try std.testing.expectError(
        error.CorruptIndex,
        main.index_format.loadIndexCheckedWithMode(io, fasta_path, .records_only),
    );
}

test "loadIndexCheckedWithMode stats_scan rejects impossible fai geometry" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_path = try uniqueArtifactPath(allocator, "bad-fai-stats-scan", "fa");
    const fai_path = try std.fmt.allocPrint(allocator, "{s}.fai", .{fasta_path});
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, fai_path) catch {};

    const fasta_file = try std.Io.Dir.cwd().createFile(io, fasta_path, .{});
    defer fasta_file.close(io);
    try std.Io.File.writeStreamingAll(fasta_file, io, ">seq1\nACGT\n");

    const fai_file = try std.Io.Dir.cwd().createFile(io, fai_path, .{ .truncate = true });
    defer fai_file.close(io);
    try std.Io.File.writeStreamingAll(fai_file, io, "seq1\t4\t6\t4\t3\n");

    try std.testing.expectError(
        error.CorruptIndex,
        main.index_format.loadIndexCheckedWithMode(io, fasta_path, .stats_scan),
    );
}

test "loadIndexChecked rejects fai decimal overflow in seq_len" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_path = try uniqueArtifactPath(allocator, "bad-fai-u64-overflow", "fa");
    const fai_path = try std.fmt.allocPrint(allocator, "{s}.fai", .{fasta_path});
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, fai_path) catch {};

    const fasta_file = try std.Io.Dir.cwd().createFile(io, fasta_path, .{});
    defer fasta_file.close(io);
    try std.Io.File.writeStreamingAll(fasta_file, io, ">seq1\nACGT\n");

    // 20 nines overflows u64 during decimal accumulation (max u64 is 20 digits starting with 1).
    const fai_file = try std.Io.Dir.cwd().createFile(io, fai_path, .{ .truncate = true });
    defer fai_file.close(io);
    try std.Io.File.writeStreamingAll(fai_file, io, "seq1\t99999999999999999999\t6\t4\t5\n");

    try std.testing.expectError(error.CorruptIndex, loadIndexChecked(io, fasta_path));
}

test "loadIndexChecked rejects fai line_bases above u32" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_path = try uniqueArtifactPath(allocator, "bad-fai-u32-overflow", "fa");
    const fai_path = try std.fmt.allocPrint(allocator, "{s}.fai", .{fasta_path});
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, fai_path) catch {};

    const fasta_file = try std.Io.Dir.cwd().createFile(io, fasta_path, .{});
    defer fasta_file.close(io);
    try std.Io.File.writeStreamingAll(fasta_file, io, ">seq1\nACGT\n");

    const fai_file = try std.Io.Dir.cwd().createFile(io, fai_path, .{ .truncate = true });
    defer fai_file.close(io);
    try std.Io.File.writeStreamingAll(fai_file, io, "seq1\t4\t6\t4294967296\t5\n");

    try std.testing.expectError(error.CorruptIndex, loadIndexChecked(io, fasta_path));
}

test "loadIndexChecked rejects fai name longer than u16" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_path = try uniqueArtifactPath(allocator, "bad-fai-long-name", "fa");
    const fai_path = try std.fmt.allocPrint(allocator, "{s}.fai", .{fasta_path});
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, fai_path) catch {};

    const fasta_file = try std.Io.Dir.cwd().createFile(io, fasta_path, .{});
    defer fasta_file.close(io);
    try std.Io.File.writeStreamingAll(fasta_file, io, ">seq1\nACGT\n");

    var fai_buf: std.ArrayListUnmanaged(u8) = .empty;
    defer fai_buf.deinit(allocator);
    try fai_buf.appendNTimes(allocator, 'N', 65536);
    try fai_buf.appendSlice(allocator, "\t4\t6\t4\t5\n");

    const fai_file = try std.Io.Dir.cwd().createFile(io, fai_path, .{ .truncate = true });
    defer fai_file.close(io);
    try std.Io.File.writeStreamingAll(fai_file, io, fai_buf.items);

    try std.testing.expectError(error.CorruptIndex, loadIndexChecked(io, fasta_path));
}

fn expectFaiLoaderModesAgree(fasta_path: []const u8) !void {
    var full = try main.index_format.loadIndexCheckedWithMode(io, fasta_path, .lookup_full_map);
    defer full.deinit(io);
    var records_only = try main.index_format.loadIndexCheckedWithMode(io, fasta_path, .records_only);
    defer records_only.deinit(io);
    var stats_scan = try main.index_format.loadIndexCheckedWithMode(io, fasta_path, .stats_scan);
    defer stats_scan.deinit(io);

    try std.testing.expectEqual(main.index_format.LoadedIndex.IndexSource.fai, full.source);
    try std.testing.expectEqual(full.source, records_only.source);
    try std.testing.expectEqual(full.source, stats_scan.source);
    try std.testing.expectEqual(full.records.len, records_only.records.len);
    try std.testing.expectEqual(full.records.len, stats_scan.records.len);
    try std.testing.expectEqual(full.records.len, stats_scan.fai_line_offsets.len);

    for (full.records, 0..) |rec, i| {
        const streamed = records_only.records[i];
        const scanned = stats_scan.records[i];

        // Geometry and length fields must match across mmap and streaming paths.
        try std.testing.expectEqual(rec.seq_offset, streamed.seq_offset);
        try std.testing.expectEqual(rec.seq_len, streamed.seq_len);
        try std.testing.expectEqual(rec.line_bases, streamed.line_bases);
        try std.testing.expectEqual(rec.line_bytes, streamed.line_bytes);

        try std.testing.expectEqual(rec.seq_offset, scanned.seq_offset);
        try std.testing.expectEqual(rec.seq_len, scanned.seq_len);
        try std.testing.expectEqual(rec.line_bases, scanned.line_bases);
        try std.testing.expectEqual(rec.line_bytes, scanned.line_bytes);

        // stats_scan sidecar offsets are the same byte positions mmap stores in name_offset.
        try std.testing.expectEqual(rec.name_offset, stats_scan.fai_line_offsets[i]);

        try std.testing.expectEqualStrings(full.getRecordName(i), records_only.getRecordName(i));
        try std.testing.expectEqualStrings(full.getRecordName(i), stats_scan.getRecordNameWithIo(io, i));
    }
}

test "fai loaders agree when index omits terminal newline" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_path = try uniqueArtifactPath(allocator, "fai-no-term-nl", "fa");
    const fai_path = try std.fmt.allocPrint(allocator, "{s}.fai", .{fasta_path});
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, fai_path) catch {};

    const fasta_file = try std.Io.Dir.cwd().createFile(io, fasta_path, .{});
    defer fasta_file.close(io);
    try std.Io.File.writeStreamingAll(fasta_file, io, ">s1\nACGT\n>s2\nGGGG\n");

    // Two records; final `.fai` line has no trailing `\n`.
    const fai_file = try std.Io.Dir.cwd().createFile(io, fai_path, .{ .truncate = true });
    defer fai_file.close(io);
    try std.Io.File.writeStreamingAll(fai_file, io, "s1\t4\t4\t4\t5\ns2\t4\t13\t4\t5");

    try expectFaiLoaderModesAgree(fasta_path);
}

test "fai loaders agree when FASTA omits terminal newline" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_path = try uniqueArtifactPath(allocator, "fasta-no-term-nl", "fa");
    const fai_path = try std.fmt.allocPrint(allocator, "{s}.fai", .{fasta_path});
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, fai_path) catch {};

    // Samtools-style line_bytes still counts the separator even when the final line has none.
    const fasta_file = try std.Io.Dir.cwd().createFile(io, fasta_path, .{});
    defer fasta_file.close(io);
    try std.Io.File.writeStreamingAll(fasta_file, io, ">seq1\nACGT");

    const fai_file = try std.Io.Dir.cwd().createFile(io, fai_path, .{ .truncate = true });
    defer fai_file.close(io);
    try std.Io.File.writeStreamingAll(fai_file, io, "seq1\t4\t6\t4\t5\n");

    try expectFaiLoaderModesAgree(fasta_path);
}

test "fai loaders agree when FASTA and index both omit terminal newline" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_path = try uniqueArtifactPath(allocator, "both-no-term-nl", "fa");
    const fai_path = try std.fmt.allocPrint(allocator, "{s}.fai", .{fasta_path});
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, fai_path) catch {};

    const fasta_file = try std.Io.Dir.cwd().createFile(io, fasta_path, .{});
    defer fasta_file.close(io);
    try std.Io.File.writeStreamingAll(fasta_file, io, ">seq1\nAAAA\nBBBB");

    const fai_file = try std.Io.Dir.cwd().createFile(io, fai_path, .{ .truncate = true });
    defer fai_file.close(io);
    try std.Io.File.writeStreamingAll(fai_file, io, "seq1\t8\t6\t4\t5");

    try expectFaiLoaderModesAgree(fasta_path);
}

test "loadIndexChecked accepts zfi when FASTA omits terminal newline" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_path = try uniqueArtifactPath(allocator, "zfi-fasta-no-term-nl", "fa");
    const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{fasta_path});
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, zfi_path) catch {};

    const fasta_data = ">seq1\nACGT";
    const fasta_file = try std.Io.Dir.cwd().createFile(io, fasta_path, .{});
    defer fasta_file.close(io);
    try std.Io.File.writeStreamingAll(fasta_file, io, fasta_data);

    var records = try scanHeaders(fasta_data, allocator);
    defer records.deinit(allocator);
    try writeZfiFromRecords(zfi_path, records.items, fasta_data.len, try statMtimeNs(fasta_path), allocator);

    var idx = try loadIndexChecked(io, fasta_path);
    defer idx.deinit(io);
    try std.testing.expectEqual(main.index_format.LoadedIndex.IndexSource.zfi, idx.source);
    try std.testing.expectEqual(@as(u64, 4), idx.records[0].seq_len);
}

test "emit-fai rejects non-uniform sequence layout" {
    const data =
        \\>mixed
        \\AAAA
        \\BBBBBB
        \\CC
        \\
    ;
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var out = std.Io.Writer.Allocating.init(allocator);
    defer out.deinit();
    try std.testing.expectError(
        error.NonUniformFai,
        main.indexer.streamingScan(data, &out.writer, .fai, true, allocator),
    );
    try std.testing.expectError(
        error.NonUniformFai,
        scanChunkedData(data, &out.writer, true, allocator),
    );
    try std.testing.expectEqual(@as(usize, 0), out.written().len);
}

test "loadIndexChecked rejects zfi after FASTA size change" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_path = try uniqueArtifactPath(allocator, "stale-size", "fa");
    const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{fasta_path});
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, zfi_path) catch {};

    const fasta_data = ">seq1\nACGT\n";
    {
        const fasta_file = try std.Io.Dir.cwd().createFile(io, fasta_path, .{});
        defer fasta_file.close(io);
        try std.Io.File.writeStreamingAll(fasta_file, io, fasta_data);
    }
    var index = try main.indexer.scanZfiIndex(fasta_data, true, allocator);
    defer index.deinit(allocator);
    try main.indexer.writeZfiIndexFile(io, zfi_path, &index, fasta_data.len, try statMtimeNs(fasta_path));

    {
        const fasta_file = try std.Io.Dir.cwd().createFile(io, fasta_path, .{ .truncate = true });
        defer fasta_file.close(io);
        try std.Io.File.writeStreamingAll(fasta_file, io, ">seq1\nACGTT\n");
    }

    try std.testing.expectError(error.StaleIndex, loadIndexChecked(io, fasta_path));
}

test "loadIndexChecked rejects zfi after same-size FASTA replacement" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_path = try uniqueArtifactPath(allocator, "stale-same-size", "fa");
    const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{fasta_path});
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, zfi_path) catch {};

    const fasta_a = ">seq1\nACGT\n";
    const fasta_b = ">seq1\nTTTT\n";
    try std.testing.expectEqual(fasta_a.len, fasta_b.len);

    {
        const fasta_file = try std.Io.Dir.cwd().createFile(io, fasta_path, .{});
        defer fasta_file.close(io);
        try std.Io.File.writeStreamingAll(fasta_file, io, fasta_a);
    }
    const mtime_a = try statMtimeNs(fasta_path);
    var index = try main.indexer.scanZfiIndex(fasta_a, true, allocator);
    defer index.deinit(allocator);
    try main.indexer.writeZfiIndexFile(io, zfi_path, &index, fasta_a.len, mtime_a);

    // Replacement updates content; force mtime forward so ZFID staleness is visible
    // even when the clock has not moved (Windows CI can keep the same stamp).
    {
        const fasta_file = try std.Io.Dir.cwd().createFile(io, fasta_path, .{ .truncate = true });
        defer fasta_file.close(io);
        try std.Io.File.writeStreamingAll(fasta_file, io, fasta_b);
        const bumped = std.Io.Timestamp.fromNanoseconds(@as(i96, @intCast(mtime_a + std.time.ns_per_s)));
        try fasta_file.setTimestamps(io, .{
            .modify_timestamp = .{ .new = bumped },
        });
    }

    try std.testing.expectError(error.StaleIndex, loadIndexChecked(io, fasta_path));
}

test "loadIndexChecked accepts legacy zfi without source-identity trailer" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_data = ">seq1\nACGT\n";
    var index = try main.indexer.scanZfiIndex(fasta_data, true, allocator);
    defer index.deinit(allocator);
    var zfi_bytes = try main.indexer.zfiIndexToBytes(&index, fasta_data.len, 0, allocator);
    defer allocator.free(zfi_bytes);

    // Strip production ZFID (between blob and footer): legacy weak identity.
    const footer_len = main.index_format.zfi_name_footer_bytes;
    const id_len = main.index_format.zfi_source_id_bytes;
    try std.testing.expect(zfi_bytes.len >= footer_len + id_len);
    const id_off = zfi_bytes.len - footer_len - id_len;
    try std.testing.expectEqualSlices(
        u8,
        &main.index_format.ZFI_SOURCE_ID_MAGIC,
        zfi_bytes[id_off..][0..4],
    );
    var legacy: std.ArrayList(u8) = .empty;
    defer legacy.deinit(allocator);
    try legacy.appendSlice(allocator, zfi_bytes[0..id_off]);
    try legacy.appendSlice(allocator, zfi_bytes[zfi_bytes.len - footer_len ..]);

    const paths = try writeFastaAndRawZfi(allocator, "zfi-legacy-no-id", fasta_data, legacy.items);
    defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, paths.zfi_path) catch {};

    var idx = try loadIndexChecked(io, paths.fasta_path);
    defer idx.deinit(io);
    try std.testing.expectEqual(main.index_format.LoadedIndex.IndexSource.zfi, idx.source);
}

test "loadIndexChecked rejects zfi with wrong embedded source mtime" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_data = ">seq1\nACGT\n";
    const fasta_path = try uniqueArtifactPath(allocator, "zfi-bad-mtime", "fa");
    const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{fasta_path});
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, zfi_path) catch {};

    {
        const fasta_file = try std.Io.Dir.cwd().createFile(io, fasta_path, .{});
        defer fasta_file.close(io);
        try std.Io.File.writeStreamingAll(fasta_file, io, fasta_data);
    }
    const real_mtime = try statMtimeNs(fasta_path);
    const wrong_mtime = real_mtime +% 1;

    var index = try main.indexer.scanZfiIndex(fasta_data, true, allocator);
    defer index.deinit(allocator);
    const zfi_bytes = try main.indexer.zfiIndexToBytes(&index, fasta_data.len, wrong_mtime, allocator);
    defer allocator.free(zfi_bytes);

    {
        const zfi_file = try std.Io.Dir.cwd().createFile(io, zfi_path, .{ .truncate = true });
        defer zfi_file.close(io);
        try std.Io.File.writeStreamingAll(zfi_file, io, zfi_bytes);
    }

    try std.testing.expectError(error.StaleIndex, loadIndexChecked(io, fasta_path));
}

test "loadIndexChecked rejects truncated zfi mid-record" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var zfi_buf: std.ArrayListUnmanaged(u8) = .empty;
    defer zfi_buf.deinit(allocator);
    try appendZfiHeader(&zfi_buf, allocator, 1, 11);
    try zfi_buf.appendSlice(allocator, &[_]u8{0} ** 20); // half an IndexRecord

    const paths = try writeFastaAndRawZfi(allocator, "zfi-trunc-record", ">seq1\nACGT\n", zfi_buf.items);
    defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, paths.zfi_path) catch {};

    try std.testing.expectError(error.CorruptIndex, loadIndexChecked(io, paths.fasta_path));
}

test "loadIndexChecked rejects zfi with zero line geometry" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_data = ">seq1\nACGT\n";
    const bad_records = [_]IndexRecord{
        .{ .name_offset = 1, .name_len = 4, .seq_offset = 6, .seq_len = 4, .line_bases = 0, .line_bytes = 5 },
    };
    const fasta_path = try uniqueArtifactPath(allocator, "zfi-zero-geom", "fa");
    const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{fasta_path});
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, zfi_path) catch {};

    {
        const fasta_file = try std.Io.Dir.cwd().createFile(io, fasta_path, .{});
        defer fasta_file.close(io);
        try std.Io.File.writeStreamingAll(fasta_file, io, fasta_data);
    }
    try writeZfiFromRecords(zfi_path, &bad_records, fasta_data.len, try statMtimeNs(fasta_path), allocator);
    try std.testing.expectError(error.CorruptIndex, loadIndexChecked(io, fasta_path));
}

test "low-mem zfi matches mmap when record ends with blank line" {
    const data = ">seq\nACGT\n\n";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    try expectZfiStreamingMatchesMmap(arena.allocator(), data, "trailing-blank");
}

test "low-mem zfi matches mmap when blank line precedes next header" {
    const data = ">a\nAAAA\n\n>b\nCCCC\n";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    try expectZfiStreamingMatchesMmap(arena.allocator(), data, "blank-before-next-header");
}

test "streaming scan keeps trailing-blank record uniform like mmap" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var stream_index = try main.indexer.scanZfiIndexStreamingData(">seq\nACGT\n\n", true, allocator);
    defer stream_index.deinit(allocator);
    try std.testing.expectEqual(@as(usize, 1), stream_index.records.items.len);
    try std.testing.expect(stream_index.records.items[0].isUniformWidth());
    try std.testing.expectEqual(@as(usize, 0), stream_index.side_tables.items.len);
}

test "streaming scan marks interior blank non-uniform like mmap" {
    const data = ">seq\nAAAA\n\nCCCC\n";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    try expectZfiStreamingMatchesMmap(allocator, data, "interior-blank");

    var stream_index = try main.indexer.scanZfiIndexStreamingData(data, true, allocator);
    defer stream_index.deinit(allocator);
    try std.testing.expect(!stream_index.records.items[0].isUniformWidth());
}

test "scanZfiIndexStreaming keeps uniform-prefix side rows after late width break" {
    // Two equal stride lines, then a wider line: stream must keep both prefix rows.
    const data = ">seq\nAAAA\nAAAA\nCCCCCC\nAAAA\n";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    try expectZfiStreamingMatchesMmap(arena.allocator(), data, "late-width-break");
}

test "streaming short final line stays uniform without side table (Track A)" {
    // Many full-width lines + short last wrap: mmap is formula-uniform; streaming must match
    // without materializing O(lines) side-table rows (Genome --low-mem RSS bug).
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var fasta: std.ArrayList(u8) = .empty;
    defer fasta.deinit(allocator);
    try fasta.appendSlice(allocator, ">chrom\n");
    var line: [61]u8 = undefined;
    @memset(line[0..60], 'A');
    line[60] = '\n';
    var i: usize = 0;
    while (i < 4000) : (i += 1) {
        try fasta.appendSlice(allocator, line[0..]);
    }
    try fasta.appendSlice(allocator, "ACGT\n");

    try expectZfiStreamingMatchesMmap(allocator, fasta.items, "short-final-wrap");

    var stream_index = try main.indexer.scanZfiIndexStreamingData(fasta.items, true, allocator);
    defer stream_index.deinit(allocator);
    try std.testing.expectEqual(@as(usize, 1), stream_index.records.items.len);
    try std.testing.expect(stream_index.records.items[0].isUniformWidth());
    try std.testing.expectEqual(@as(usize, 0), stream_index.side_tables.items.len);
    try std.testing.expectEqual(@as(u64, 4000 * 60 + 4), stream_index.records.items[0].seq_len);
}

test "streaming interior short line then resume matches mmap side table" {
    const data = ">seq\nAAAA\nAAAA\nAA\nAAAA\n";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    try expectZfiStreamingMatchesMmap(allocator, data, "interior-short-then-resume");

    var stream_index = try main.indexer.scanZfiIndexStreamingData(data, true, allocator);
    defer stream_index.deinit(allocator);
    try std.testing.expect(!stream_index.records.items[0].isUniformWidth());
    try std.testing.expect(stream_index.side_tables.items.len > 0);
}

test "streaming long final line is non-uniform like mmap" {
    const data = ">seq\nAAAA\nAAAA\nAAAAAA\n";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    try expectZfiStreamingMatchesMmap(allocator, data, "long-final-line");

    var stream_index = try main.indexer.scanZfiIndexStreamingData(data, true, allocator);
    defer stream_index.deinit(allocator);
    try std.testing.expect(!stream_index.records.items[0].isUniformWidth());
    try std.testing.expect(stream_index.side_tables.items.len > 0);
}

test "streaming short final line without terminal newline matches mmap" {
    const data = ">seq\nAAAA\nAAAA\nAC";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    try expectZfiStreamingMatchesMmap(allocator, data, "short-final-no-terminal-nl");

    var stream_index = try main.indexer.scanZfiIndexStreamingData(data, true, allocator);
    defer stream_index.deinit(allocator);
    try std.testing.expect(stream_index.records.items[0].isUniformWidth());
    try std.testing.expectEqual(@as(usize, 0), stream_index.side_tables.items.len);
}

test "streaming short final line across chunk boundary matches mmap" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var fasta: std.ArrayList(u8) = .empty;
    defer fasta.deinit(allocator);
    try fasta.appendSlice(allocator, ">seq\n");
    const full_line = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"; // 60 + nl
    while (fasta.items.len + full_line.len < low_mem_chunk_size + 200) {
        try fasta.appendSlice(allocator, full_line);
    }
    try fasta.appendSlice(allocator, "ACGT\n");

    try expectZfiStreamingMatchesMmap(allocator, fasta.items, "short-final-chunk-boundary");
    var stream_index = try main.indexer.scanZfiIndexStreamingData(fasta.items, true, allocator);
    defer stream_index.deinit(allocator);
    try std.testing.expect(stream_index.records.items[0].isUniformWidth());
    try std.testing.expectEqual(@as(usize, 0), stream_index.side_tables.items.len);
}

test "streaming CRLF uniform body with short final matches mmap (Track A2)" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var fasta: std.ArrayList(u8) = .empty;
    defer fasta.deinit(allocator);
    try fasta.appendSlice(allocator, ">chrom\r\n");
    var line: [62]u8 = undefined;
    @memset(line[0..60], 'A');
    line[60] = '\r';
    line[61] = '\n';
    var i: usize = 0;
    while (i < 200) : (i += 1) {
        try fasta.appendSlice(allocator, line[0..]);
    }
    try fasta.appendSlice(allocator, "ACGT\r\n");

    try expectZfiStreamingMatchesMmap(allocator, fasta.items, "crlf-short-final");
    var stream_index = try main.indexer.scanZfiIndexStreamingData(fasta.items, true, allocator);
    defer stream_index.deinit(allocator);
    try std.testing.expect(stream_index.records.items[0].isUniformWidth());
    try std.testing.expectEqual(@as(usize, 0), stream_index.side_tables.items.len);
    try std.testing.expectEqual(@as(u32, 62), stream_index.records.items[0].line_bytes);
}

test "streaming multi-record uniform bodies stay mmap-identical (Track A2)" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var fasta: std.ArrayList(u8) = .empty;
    defer fasta.deinit(allocator);
    const full = "AAAA\n";
    var r: usize = 0;
    while (r < 3) : (r += 1) {
        var hdr_buf: [16]u8 = undefined;
        const hdr = try std.fmt.bufPrint(&hdr_buf, ">r{d}\n", .{r});
        try fasta.appendSlice(allocator, hdr);
        var i: usize = 0;
        while (i < 50) : (i += 1) {
            try fasta.appendSlice(allocator, full);
        }
        try fasta.appendSlice(allocator, "AC\n");
    }

    try expectZfiStreamingMatchesMmap(allocator, fasta.items, "multi-record-stride");
    var stream_index = try main.indexer.scanZfiIndexStreamingData(fasta.items, true, allocator);
    defer stream_index.deinit(allocator);
    try std.testing.expectEqual(@as(usize, 3), stream_index.records.items.len);
    for (stream_index.records.items) |rec| {
        try std.testing.expect(rec.isUniformWidth());
    }
}

test "streaming interior blank after stride prefix matches mmap" {
    const data = ">seq\nAAAA\nAAAA\nAAAA\n\nAAAA\n";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    try expectZfiStreamingMatchesMmap(arena.allocator(), data, "interior-blank-after-stride");
}

test "streaming short final plus next header spanning stride width matches mmap (Track A2)" {
    // Genome chr10→chr11: `NN\\n` + `>11 ... REF\\n` is exactly line_bytes (61). Stride must
    // not treat that span as one wrap or the next record is swallowed into the previous.
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var fasta: std.ArrayList(u8) = .empty;
    defer fasta.deinit(allocator);
    try fasta.appendSlice(allocator, ">a\n");
    const full = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"; // 60+nl
    var i: usize = 0;
    while (i < 5) : (i += 1) {
        try fasta.appendSlice(allocator, full);
    }
    // short final (2) + nl + header line crafted so short+header == 61 bytes
    try fasta.appendSlice(allocator, "NN\n");
    // 58-byte header line (incl. `\n`) so `NN\n` + header == 61 == line_bytes.
    const hdr = ">b xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n";
    try std.testing.expectEqual(@as(usize, 58), hdr.len);
    try fasta.appendSlice(allocator, hdr);
    try fasta.appendSlice(allocator, full);
    try fasta.appendSlice(allocator, "AC\n");

    try expectZfiStreamingMatchesMmap(allocator, fasta.items, "short-plus-header-eq-stride");
    var stream_index = try main.indexer.scanZfiIndexStreamingData(fasta.items, true, allocator);
    defer stream_index.deinit(allocator);
    try std.testing.expectEqual(@as(usize, 2), stream_index.records.items.len);
    const n0 = stream_index.name_blob.items[stream_index.records.items[0].name_offset..][0..stream_index.records.items[0].name_len];
    const n1 = stream_index.name_blob.items[stream_index.records.items[1].name_offset..][0..stream_index.records.items[1].name_len];
    try std.testing.expectEqualStrings("a", n0);
    try std.testing.expectEqualStrings("b", n1);
}

// ============================================================================
// Stable gates fixtures (tests/data/gates): malformed, boundary, portability
// ============================================================================

test "gates fixture: duplicate names keep first-wins on mmap and streaming" {
    try requireFixturePath("tests/data/gates/duplicates.fasta");
    const data = try readTestFile(std.testing.allocator, "tests/data/gates/duplicates.fasta");
    defer std.testing.allocator.free(data);

    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var mmap_index = try main.indexer.scanZfiIndex(data, true, allocator);
    defer mmap_index.deinit(allocator);
    var stream_index = try main.indexer.scanZfiIndexStreamingData(data, true, allocator);
    defer stream_index.deinit(allocator);

    try std.testing.expectEqual(@as(usize, 2), mmap_index.records.items.len);
    try std.testing.expectEqual(mmap_index.records.items.len, stream_index.records.items.len);
    try expectZfiStreamingMatchesMmap(allocator, data, "gates/duplicates.fasta");

    const n0 = mmap_index.name_blob.items[mmap_index.records.items[0].name_offset..][0..mmap_index.records.items[0].name_len];
    const n1 = mmap_index.name_blob.items[mmap_index.records.items[1].name_offset..][0..mmap_index.records.items[1].name_len];
    try std.testing.expectEqualStrings("dup", n0);
    try std.testing.expectEqualStrings("other", n1);
    try std.testing.expectEqual(@as(u64, 4), mmap_index.records.items[0].seq_len);
    try std.testing.expectEqual(@as(u64, 4), mmap_index.records.items[1].seq_len);
}

test "gates fixture: long header indexes; over-u16 name rejects on both paths" {
    try requireFixturePath("tests/data/gates/long_header.fasta");
    try requireFixturePath("tests/data/gates/long_header_reject.fasta");

    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const ok_data = try readTestFile(allocator, "tests/data/gates/long_header.fasta");
    defer allocator.free(ok_data);
    var ok_index = try main.indexer.scanZfiIndex(ok_data, true, allocator);
    defer ok_index.deinit(allocator);
    try std.testing.expectEqual(@as(usize, 1), ok_index.records.items.len);
    try std.testing.expectEqual(@as(u16, 1025), ok_index.records.items[0].name_len);
    try expectZfiStreamingMatchesMmap(allocator, ok_data, "gates/long_header.fasta");

    const reject_data = try readTestFile(allocator, "tests/data/gates/long_header_reject.fasta");
    defer allocator.free(reject_data);
    try std.testing.expect(reject_data[0] == '>');
    // Fixture must stay LF-only; a Windows CRLF checkout makes indexOf('\\n')-1 count the '\\r'.
    try std.testing.expect(std.mem.indexOfScalar(u8, reject_data, '\r') == null);
    try std.testing.expectEqual(main.indexer.max_index_name_len + 1, std.mem.indexOfScalar(u8, reject_data, '\n').? - 1);
    try std.testing.expectError(error.HeaderTooLong, main.indexer.scanZfiIndex(reject_data, true, allocator));
    try std.testing.expectError(error.HeaderTooLong, main.indexer.scanZfiIndexStreamingData(reject_data, true, allocator));
}

test "gates fixture: invalid UTF-8 header bytes index and stay embeddable" {
    try requireFixturePath("tests/data/gates/invalid_utf8_header.fasta");
    const data = try readTestFile(std.testing.allocator, "tests/data/gates/invalid_utf8_header.fasta");
    defer std.testing.allocator.free(data);

    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    try expectZfiStreamingMatchesMmap(allocator, data, "gates/invalid_utf8_header.fasta");
    var index = try main.indexer.scanZfiIndex(data, true, allocator);
    defer index.deinit(allocator);
    // Duplicate headers: first-wins keeps one record with the raw name bytes.
    try std.testing.expectEqual(@as(usize, 1), index.records.items.len);
    const name0 = index.name_blob.items[index.records.items[0].name_offset..][0..index.records.items[0].name_len];
    try std.testing.expectEqualSlices(u8, "bad\xffname", name0);
}

test "gates fixture: zero geometry and large offset FAI rejected" {
    try requireFixturePath("tests/data/gates/seq1_zero_geometry.fasta");
    try requireFixturePath("tests/data/gates/seq1_zero_geometry.fasta.fai");
    try requireFixturePath("tests/data/gates/seq1_large_offset.fasta");
    try requireFixturePath("tests/data/gates/seq1_large_offset.fasta.fai");

    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    // Copy into temp paths so a stray checked-in `.zfi` cannot mask the `.fai` case.
    try expectCorruptFaiFixture(allocator, "tests/data/gates/seq1_zero_geometry");
    try expectCorruptFaiFixture(allocator, "tests/data/gates/seq1_large_offset");
}

fn expectCorruptFaiFixture(allocator: std.mem.Allocator, stem: []const u8) !void {
    const src_fa = try std.fmt.allocPrint(allocator, "{s}.fasta", .{stem});
    const src_fai = try std.fmt.allocPrint(allocator, "{s}.fasta.fai", .{stem});
    const fasta_data = try readTestFile(allocator, src_fa);
    defer allocator.free(fasta_data);
    const fai_data = try readTestFile(allocator, src_fai);
    defer allocator.free(fai_data);

    const fasta_path = try uniqueArtifactPath(allocator, "gates-fai", "fa");
    const fai_path = try std.fmt.allocPrint(allocator, "{s}.fai", .{fasta_path});
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, fai_path) catch {};

    {
        const fasta_file = try std.Io.Dir.cwd().createFile(io, fasta_path, .{});
        defer fasta_file.close(io);
        try std.Io.File.writeStreamingAll(fasta_file, io, fasta_data);
    }
    {
        const fai_file = try std.Io.Dir.cwd().createFile(io, fai_path, .{ .truncate = true });
        defer fai_file.close(io);
        try std.Io.File.writeStreamingAll(fai_file, io, fai_data);
    }
    try std.testing.expectError(error.CorruptIndex, loadIndexChecked(io, fasta_path));
}

test "gates fixture: zfi zero geometry from stable seq1 FASTA rejected" {
    try requireFixturePath("tests/data/gates/seq1_zero_geometry.fasta");
    const fasta_data = try readTestFile(std.testing.allocator, "tests/data/gates/seq1_zero_geometry.fasta");
    defer std.testing.allocator.free(fasta_data);

    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_path = try uniqueArtifactPath(allocator, "gates-zfi-zero", "fa");
    const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{fasta_path});
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, zfi_path) catch {};

    {
        const fasta_file = try std.Io.Dir.cwd().createFile(io, fasta_path, .{});
        defer fasta_file.close(io);
        try std.Io.File.writeStreamingAll(fasta_file, io, fasta_data);
    }
    const bad_records = [_]IndexRecord{
        .{ .name_offset = 1, .name_len = 4, .seq_offset = 6, .seq_len = 4, .line_bases = 0, .line_bytes = 5 },
    };
    try writeZfiFromRecords(zfi_path, &bad_records, fasta_data.len, try statMtimeNs(fasta_path), allocator);
    try std.testing.expectError(error.CorruptIndex, loadIndexChecked(io, fasta_path));
}

test "NameDedup identity holds under forced hash collisions (P0 injectable check)" {
    // Injectable hasher: every key lands in one bucket. NameDedupWith.observe must still
    // keep distinct names via Context.eql (same code path as production NameDedup).
    const CollisionCtx = struct {
        pub fn hash(_: @This(), _: []const u8) u64 {
            return 0xC0FFEE;
        }
        pub fn eql(_: @This(), a: []const u8, b: []const u8) bool {
            return std.mem.eql(u8, a, b);
        }
    };

    const CollidingDedup = main.indexer.NameDedupWith(CollisionCtx);
    var seen = CollidingDedup.init(std.testing.allocator, true);
    defer seen.deinit();

    const names = [_][]const u8{ "alpha", "beta", "alpha", "gamma", "beta" };
    var kept: usize = 0;
    for (names) |name| {
        if (!(try seen.observe(name))) kept += 1;
    }
    try std.testing.expectEqual(@as(usize, 3), kept);
    try std.testing.expectEqual(@as(usize, 3), seen.map.count());

    var prod = main.indexer.NameDedup.init(std.testing.allocator, true);
    defer prod.deinit();
    var prod_kept: usize = 0;
    for (names) |name| {
        if (!(try prod.observe(name))) prod_kept += 1;
    }
    try std.testing.expectEqual(kept, prod_kept);
}
