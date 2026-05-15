const std = @import("std");
const main = @import("main");
const validateFasta = main.validateFasta;
const scanHeaders = main.scanHeaders;
const scanChunkedData = main.indexer.scanChunkedData;
const loadIndexChecked = main.index_format.loadIndexChecked;
const IndexRecord = main.IndexRecord;
const writeZfi = main.writeZfi;
const ZfiHeader = main.ZfiHeader;
const ZFI_MAGIC = main.ZFI_MAGIC;
const io = std.Io.Threaded.global_single_threaded.io();

// ============================================================================
// Test helper: read file into memory
// ============================================================================

fn readTestFile(allocator: std.mem.Allocator, path: []const u8) ![]u8 {
    const file = try std.Io.Dir.cwd().openFile(io, path, .{});
    defer file.close(io);
    const stat = try file.stat(io);
    const data = try allocator.alloc(u8, @intCast(stat.size));
    const bytes_read = try std.Io.File.readPositionalAll(file, io, data, 0);
    return data[0..bytes_read];
}

fn uniqueArtifactPath(allocator: std.mem.Allocator, stem: []const u8, ext: []const u8) ![]u8 {
    try std.Io.Dir.cwd().createDirPath(io, "zig-cache/test-artifacts");
    var ts: std.os.linux.timespec = undefined;
    _ = std.os.linux.clock_gettime(.MONOTONIC, &ts);
    const nanos: u64 = @intCast(@as(i128, ts.sec) * std.time.ns_per_s + ts.nsec);
    return std.fmt.allocPrint(allocator, "zig-cache/test-artifacts/{s}-{d}.{s}", .{
        stem,
        nanos,
        ext,
    });
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
// writeZfi tests
// ============================================================================

test "writeZfi creates valid file" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();

    const records = [_]IndexRecord{
        .{ .name_offset = 1, .name_len = 4, .seq_offset = 10, .seq_len = 100, .line_bases = 80, .line_bytes = 81 },
    };

    const path = "tests/data/test_write.zfi";
    defer std.Io.Dir.cwd().deleteFile(io, path) catch {};

    try writeZfi(io, path, &records, 1000);

    // Read and verify
    const file = try std.Io.Dir.cwd().openFile(io, path, .{});
    defer file.close(io);

    var header_bytes: [@sizeOf(ZfiHeader)]u8 = undefined;
    _ = try std.Io.File.readPositionalAll(file, io, &header_bytes, 0);
    const header: ZfiHeader = @bitCast(header_bytes);

    try std.testing.expectEqualSlices(u8, &ZFI_MAGIC, &header.magic);
    try std.testing.expectEqual(@as(u32, 1), header.record_count);
    try std.testing.expectEqual(@as(u64, 1000), header.source_size);
}

test "ZfiHeader has correct size" {
    try std.testing.expectEqual(@as(usize, 16), @sizeOf(ZfiHeader));
}

test "IndexRecord has consistent size" {
    // name_offset(8) + name_len(2) + padding(6) + seq_offset(8) + seq_len(8) + line_bases(4) + line_bytes(4) = 40
    try std.testing.expectEqual(@as(usize, 40), @sizeOf(IndexRecord));
}

test "low-mem indexing matches default across chunk boundary" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var fasta: std.ArrayList(u8) = .empty;
    defer fasta.deinit(allocator);
    try fasta.appendSlice(allocator, ">seq1\n");
    try fasta.appendNTimes(allocator, 'A', 4 * 1024 * 1024 + 10);
    try fasta.appendSlice(allocator, "\n>seq2\nACGT\n");

    var expected = std.Io.Writer.Allocating.init(allocator);
    const expected_count = try main.indexer.streamingScan(fasta.items, &expected.writer, .fai, true, allocator);

    var actual = std.Io.Writer.Allocating.init(allocator);
    const actual_count = try scanChunkedData(fasta.items, &actual.writer, allocator);

    try std.testing.expectEqual(expected_count, actual_count);
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

test "low-mem indexing rejects overlong sequence names" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var fasta: std.ArrayList(u8) = .empty;
    defer fasta.deinit(allocator);
    try fasta.append(allocator, '>');
    try fasta.appendNTimes(allocator, 'A', 4097);
    try fasta.appendSlice(allocator, "\nACGT\n");

    var output = std.Io.Writer.Allocating.init(allocator);

    try std.testing.expectError(
        error.HeaderTooLong,
        scanChunkedData(fasta.items, &output.writer, allocator),
    );
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
    try writeZfi(io, zfi_path, &bad_records, 11);

    try std.testing.expectError(error.CorruptIndex, loadIndexChecked(io, fasta_path));
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
    try writeZfi(io, zfi_path, &records, fasta_data.len);

    var full_map_idx = try main.index_format.loadIndexCheckedWithMode(io, fasta_path, .lookup_full_map);
    defer full_map_idx.deinit();
    var records_only_idx = try main.index_format.loadIndexCheckedWithMode(io, fasta_path, .records_only);
    defer records_only_idx.deinit();

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
    try writeZfi(io, zfi_path, records.items, fasta_data.len);

    const fai_file = try std.Io.Dir.cwd().createFile(io, fai_path, .{ .truncate = true });
    defer fai_file.close(io);
    var fai_buf: [65536]u8 = undefined;
    var fai_fw = fai_file.writer(io, &fai_buf);
    _ = try main.indexer.streamingScan(fasta_data, &fai_fw.interface, .fai, true, allocator);
    try fai_fw.flush();

    var now: std.os.linux.timespec = undefined;
    _ = std.os.linux.clock_gettime(.REALTIME, &now);
    const stale_time: std.os.linux.timespec = .{ .sec = now.sec - 3600, .nsec = now.nsec };
    const zfi_file = try std.Io.Dir.cwd().openFile(io, zfi_path, .{});
    defer zfi_file.close(io);
    try std.testing.expectEqual(@as(usize, 0), std.os.linux.futimens(zfi_file.handle, &.{ stale_time, stale_time }));

    var idx = try loadIndexChecked(io, fasta_path);
    defer idx.deinit();

    try std.testing.expectEqual(main.index_format.LoadedIndex.IndexSource.fai, idx.source);
    try std.testing.expectEqual(@as(?usize, 0), idx.lookupName("seq1"));
}
