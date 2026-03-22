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

// ============================================================================
// Test helper: read file into memory
// ============================================================================

fn readTestFile(allocator: std.mem.Allocator, path: []const u8) ![]u8 {
    const file = try std.fs.cwd().openFile(path, .{});
    defer file.close();
    const stat = try file.stat();
    const data = try allocator.alloc(u8, stat.size);
    const bytes_read = try file.readAll(data);
    return data[0..bytes_read];
}

fn uniqueArtifactPath(allocator: std.mem.Allocator, stem: []const u8, ext: []const u8) ![]u8 {
    try std.fs.cwd().makePath("zig-cache/test-artifacts");
    return std.fmt.allocPrint(allocator, "zig-cache/test-artifacts/{s}-{d}.{s}", .{
        stem,
        std.time.nanoTimestamp(),
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
    defer std.fs.cwd().deleteFile(path) catch {};

    try writeZfi(path, &records, 1000);

    // Read and verify
    const file = try std.fs.cwd().openFile(path, .{});
    defer file.close();

    var header_bytes: [@sizeOf(ZfiHeader)]u8 = undefined;
    _ = try file.readAll(&header_bytes);
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

    var fasta = std.ArrayList(u8).init(allocator);
    defer fasta.deinit();
    try fasta.appendSlice(">seq1\n");
    try fasta.appendNTimes('A', 4 * 1024 * 1024 + 10);
    try fasta.appendSlice("\n>seq2\nACGT\n");

    var expected = std.ArrayList(u8).init(allocator);
    defer expected.deinit();
    const expected_writer = expected.writer();
    const expected_count = try main.indexer.streamingScan(fasta.items, expected_writer, .fai, true, allocator);

    var actual = std.ArrayList(u8).init(allocator);
    defer actual.deinit();
    const actual_writer = actual.writer();
    const actual_count = try scanChunkedData(fasta.items, actual_writer, allocator);

    try std.testing.expectEqual(expected_count, actual_count);
    try std.testing.expectEqualStrings(expected.items, actual.items);
}

test "low-mem indexing rejects overlong sequence names" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var fasta = std.ArrayList(u8).init(allocator);
    defer fasta.deinit();
    try fasta.append('>');
    try fasta.appendNTimes('A', 4097);
    try fasta.appendSlice("\nACGT\n");

    var output = std.ArrayList(u8).init(allocator);
    defer output.deinit();

    try std.testing.expectError(
        error.HeaderTooLong,
        scanChunkedData(fasta.items, output.writer(), allocator),
    );
}

test "loadIndexChecked rejects corrupt zfi records" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_path = try uniqueArtifactPath(allocator, "corrupt-zfi", "fa");
    const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{fasta_path});
    defer std.fs.cwd().deleteFile(fasta_path) catch {};
    defer std.fs.cwd().deleteFile(zfi_path) catch {};

    const fasta_file = try std.fs.cwd().createFile(fasta_path, .{});
    defer fasta_file.close();
    try fasta_file.writeAll(">seq1\nACGT\n");

    const bad_records = [_]IndexRecord{
        .{ .name_offset = 999_999, .name_len = 4, .seq_offset = 6, .seq_len = 4, .line_bases = 4, .line_bytes = 5 },
    };
    try writeZfi(zfi_path, &bad_records, 11);

    try std.testing.expectError(error.CorruptIndex, loadIndexChecked(fasta_path));
}

test "loadIndexChecked falls back to fai when zfi is stale" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_path = try uniqueArtifactPath(allocator, "stale-zfi", "fa");
    const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{fasta_path});
    const fai_path = try std.fmt.allocPrint(allocator, "{s}.fai", .{fasta_path});
    defer std.fs.cwd().deleteFile(fasta_path) catch {};
    defer std.fs.cwd().deleteFile(zfi_path) catch {};
    defer std.fs.cwd().deleteFile(fai_path) catch {};

    try std.fs.cwd().copyFile("tests/data/simple.fasta", std.fs.cwd(), fasta_path, .{});
    try std.fs.cwd().copyFile("tests/data/simple.fasta.zfi", std.fs.cwd(), zfi_path, .{});
    try std.fs.cwd().copyFile("tests/data/simple.fasta.fai", std.fs.cwd(), fai_path, .{});

    const stale_time = std.time.timestamp() - 3600;
    const zfi_file = try std.fs.cwd().openFile(zfi_path, .{});
    defer zfi_file.close();
    try zfi_file.updateTimes(stale_time, stale_time);

    var idx = try loadIndexChecked(fasta_path);
    defer idx.deinit();

    try std.testing.expectEqual(main.index_format.LoadedIndex.IndexSource.fai, idx.source);
    try std.testing.expectEqual(@as(?usize, 0), idx.lookupName("seq1"));
}
