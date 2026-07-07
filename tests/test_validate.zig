const std = @import("std");
const main = @import("main");

const io = std.Io.Threaded.global_single_threaded.io();
const validator = main.validator;

fn countKind(summary: *const validator.Summary, kind: validator.Kind) usize {
    var count: usize = 0;
    for (summary.events.items) |event| {
        if (event.kind == kind) count += 1;
    }
    return count;
}

fn expectFixIdempotent(allocator: std.mem.Allocator, broken: []const u8) !void {
    var before = try validator.validateData(allocator, broken, .{});
    defer before.deinit(allocator);
    try std.testing.expect(validator.fixRejection(&before, .{}) == null);

    const fixed = try validator.fixData(allocator, broken, before.record_widths.items);
    defer allocator.free(fixed);

    var after = try validator.validateData(allocator, fixed, .{});
    defer after.deinit(allocator);
    try std.testing.expectEqual(@as(usize, 0), after.events.items.len);
}

// --- unit tests (moved from src/validator.zig) ---

test "validateData reports duplicate and empty sequence" {
    var summary = try validator.validateData(std.testing.allocator, ">dup\nAAAA\n>dup\n", .{});
    defer summary.deinit(std.testing.allocator);

    try std.testing.expectEqual(@as(usize, 1), countKind(&summary, .duplicate_name));
    try std.testing.expectEqual(@as(usize, 1), countKind(&summary, .empty_sequence));
    try std.testing.expectEqual(@as(usize, 0), countKind(&summary, .no_sequences));
    try std.testing.expectEqual(@as(usize, 1), summary.sequence_count);
    try std.testing.expectEqual(@as(usize, 2), summary.header_count);
}

test "validateData reports missing terminal newline and invalid nucleotide character" {
    var summary = try validator.validateData(std.testing.allocator, ">seq\nACGTZ", .{});
    defer summary.deinit(std.testing.allocator);

    try std.testing.expectEqual(@as(usize, 1), countKind(&summary, .missing_terminal_newline));
    try std.testing.expectEqual(@as(usize, 1), countKind(&summary, .invalid_character));
}

test "validateData custom alphabet overrides built-in alphabet" {
    var summary = try validator.validateData(std.testing.allocator, ">seq\nACGTZ\n", .{
        .custom_alphabet = "ACGTZ",
    });
    defer summary.deinit(std.testing.allocator);

    try std.testing.expectEqual(@as(usize, 0), countKind(&summary, .invalid_character));
}

test "validateData tracks line-width and trailing whitespace warnings" {
    var summary = try validator.validateData(std.testing.allocator, ">seq\nAAAA\nCC \nGGGG\nTT\n", .{});
    defer summary.deinit(std.testing.allocator);

    try std.testing.expectEqual(@as(usize, 1), countKind(&summary, .trailing_whitespace));
    try std.testing.expectEqual(@as(usize, 1), countKind(&summary, .inconsistent_line_widths));
}

test "validateData checks schemas" {
    var uniprot = try validator.validateData(std.testing.allocator, ">sp|P12345|PROT_HUMAN\nMAV\n", .{
        .schema = .uniprot,
    });
    defer uniprot.deinit(std.testing.allocator);

    var refseq = try validator.validateData(std.testing.allocator, ">bad\nACGT\n", .{
        .schema = .refseq,
    });
    defer refseq.deinit(std.testing.allocator);

    try std.testing.expectEqual(@as(usize, 0), countKind(&uniprot, .schema_violation));
    try std.testing.expectEqual(@as(usize, 1), countKind(&refseq, .schema_violation));
}

test "exitCode implements strict warning promotion" {
    try std.testing.expectEqual(@as(u8, 0), validator.exitCode(0, 0, false));
    try std.testing.expectEqual(@as(u8, 2), validator.exitCode(0, 1, false));
    try std.testing.expectEqual(@as(u8, 1), validator.exitCode(0, 1, true));
    try std.testing.expectEqual(@as(u8, 1), validator.exitCode(1, 0, false));
}

test "exitCodeForOptions ignores format warnings fixed by rewrite" {
    var summary = try validator.validateData(std.testing.allocator, ">seq\nACGT \n", .{});
    defer summary.deinit(std.testing.allocator);

    try std.testing.expectEqual(@as(usize, 1), countKind(&summary, .trailing_whitespace));
    try std.testing.expectEqual(@as(u8, 2), validator.exitCodeForOptions(&summary, .{}));
    try std.testing.expectEqual(@as(u8, 0), validator.exitCodeForOptions(&summary, .{ .fix = true }));
}

test "fixData idempotent: BOM, CRLF, trailing whitespace, width, missing newline" {
    const broken = "\xEF\xBB\xBF>seq\r\nAAAA\r\nCC \r\nGGGG\r\nTT";
    try expectFixIdempotent(std.testing.allocator, broken);
}

test "fixData idempotent: mixed line endings" {
    const broken = ">seq\r\nAAAA\nCC\n";
    try expectFixIdempotent(std.testing.allocator, broken);
}

test "fixData idempotent: already clean FASTA" {
    const clean = ">seq\nACGT\n";
    try expectFixIdempotent(std.testing.allocator, clean);

    var summary = try validator.validateData(std.testing.allocator, clean, .{});
    defer summary.deinit(std.testing.allocator);
    const fixed = try validator.fixData(std.testing.allocator, clean, summary.record_widths.items);
    defer std.testing.allocator.free(fixed);
    try std.testing.expectEqualStrings(clean, fixed);
}

test "fixRejection blocks unfixable errors" {
    var dup = try validator.validateData(std.testing.allocator, ">dup\nACGT\n>dup\nGCTA\n", .{});
    defer dup.deinit(std.testing.allocator);
    try std.testing.expectEqual(.duplicate_name, validator.fixRejection(&dup, .{}).?);

    var none = try validator.validateData(std.testing.allocator, "not fasta\n", .{});
    defer none.deinit(std.testing.allocator);
    try std.testing.expectEqual(.no_sequences, validator.fixRejection(&none, .{}).?);

    var bad = try validator.validateData(std.testing.allocator, ">seq\nACGTZ\n", .{});
    defer bad.deinit(std.testing.allocator);
    try std.testing.expectEqual(.invalid_character, validator.fixRejection(&bad, .{}).?);
    try std.testing.expect(validator.fixRejection(&bad, .{ .fix_format_only = true }) == null);

    var nulls = try validator.validateData(std.testing.allocator, ">seq\nACGT\x00TT\n", .{});
    defer nulls.deinit(std.testing.allocator);
    try std.testing.expectEqual(.null_byte, validator.fixRejection(&nulls, .{}).?);
}

test "fixData with fix-format-only preserves invalid characters" {
    const broken = ">seq\nACGTZ\n";
    var summary = try validator.validateData(std.testing.allocator, broken, .{});
    defer summary.deinit(std.testing.allocator);
    try std.testing.expect(validator.fixRejection(&summary, .{ .fix_format_only = true }) == null);

    const fixed = try validator.fixData(std.testing.allocator, broken, summary.record_widths.items);
    defer std.testing.allocator.free(fixed);
    try std.testing.expect(std.mem.indexOf(u8, fixed, "ACGTZ") != null);

    var after = try validator.validateData(std.testing.allocator, fixed, .{});
    defer after.deinit(std.testing.allocator);
    try std.testing.expectEqual(@as(usize, 1), countKind(&after, .invalid_character));
}

test "fixData does not remove empty sequence warnings" {
    const broken = ">empty\n";
    var summary = try validator.validateData(std.testing.allocator, broken, .{});
    defer summary.deinit(std.testing.allocator);
    try std.testing.expectEqual(@as(usize, 1), countKind(&summary, .empty_sequence));

    const fixed = try validator.fixData(std.testing.allocator, broken, summary.record_widths.items);
    defer std.testing.allocator.free(fixed);

    var after = try validator.validateData(std.testing.allocator, fixed, .{});
    defer after.deinit(std.testing.allocator);
    try std.testing.expectEqual(@as(usize, 1), countKind(&after, .empty_sequence));
}

// --- integration tests ---

fn readTestFile(allocator: std.mem.Allocator, path: []const u8) ![]u8 {
    const file = try std.Io.Dir.cwd().openFile(io, path, .{});
    defer file.close(io);
    const stat = try file.stat(io);
    const data = try allocator.alloc(u8, stat.size);
    _ = try file.readPositionalAll(io, data, 0);
    return data;
}

fn uniqueArtifactPath(allocator: std.mem.Allocator, stem: []const u8) ![]u8 {
    try std.Io.Dir.cwd().createDirPath(io, "zig-cache/test-artifacts");
    const now = std.Io.Clock.Timestamp.now(io, .awake);
    const nanos: u64 = @intCast(now.raw.toNanoseconds());
    return std.fmt.allocPrint(allocator, "zig-cache/test-artifacts/{s}-{d}.fa", .{ stem, nanos });
}

fn writeFastaArtifact(allocator: std.mem.Allocator, stem: []const u8, data: []const u8) ![]const u8 {
    const path = try uniqueArtifactPath(allocator, stem);
    const file = try std.Io.Dir.cwd().createFile(io, path, .{});
    defer file.close(io);
    try std.Io.File.writeStreamingAll(file, io, data);
    return path;
}

fn scanZfi(allocator: std.mem.Allocator, data: []const u8) !main.indexer.ZfiIndex {
    return main.indexer.scanZfiIndex(data, true, allocator);
}

fn expectHasSideTable(index: *const main.indexer.ZfiIndex) !void {
    var any_non_uniform = false;
    for (index.records.items) |rec| {
        if (!rec.isUniformWidth()) any_non_uniform = true;
    }
    try std.testing.expect(any_non_uniform);
    try std.testing.expect(index.side_tables.items.len > 0);
}

fn expectAllUniformNoSideTable(index: *const main.indexer.ZfiIndex) !void {
    for (index.records.items) |rec| {
        try std.testing.expect(rec.isUniformWidth());
    }
    try std.testing.expectEqual(@as(usize, 0), index.side_tables.items.len);
}

fn writeZfiForData(allocator: std.mem.Allocator, fasta_path: []const u8, data: []const u8) !void {
    var arena = std.heap.ArenaAllocator.init(allocator);
    defer arena.deinit();

    var index = try scanZfi(arena.allocator(), data);
    defer index.deinit(arena.allocator());

    var zfi_path_buf: [4096]u8 = undefined;
    const zfi_path = try std.fmt.bufPrint(&zfi_path_buf, "{s}.zfi", .{fasta_path});
    try main.indexer.writeZfiIndexFile(io, zfi_path, &index, data.len);
}

fn captureExtractRegion(allocator: std.mem.Allocator, fasta_path: []const u8, region: []const u8) ![]u8 {
    var idx = try main.index_format.loadIndexChecked(io, fasta_path);
    defer idx.deinit();

    var out = std.Io.Writer.Allocating.init(allocator);
    main.getter.extractRegion(&idx, region, &out.writer);
    return out.toOwnedSlice();
}

fn fixFormatOnly(allocator: std.mem.Allocator, broken: []const u8) ![]u8 {
    var summary = try validator.validateData(allocator, broken, .{});
    defer summary.deinit(allocator);
    try std.testing.expect(validator.fixRejection(&summary, .{}) == null);
    return validator.fixData(allocator, broken, summary.record_widths.items);
}

fn expectValidateClean(allocator: std.mem.Allocator, data: []const u8) !void {
    var summary = try validator.validateData(allocator, data, .{});
    defer summary.deinit(allocator);
    try std.testing.expectEqual(@as(usize, 0), summary.events.items.len);
}

test "validate --fix then index uses uniform O(1) path on crafted messy FASTA" {
    const broken =
        \\>messy_seq widths and trailing ws
        \\AAAA    
        \\CCCC
        \\GGGGTT
        \\
    ;

    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var broken_index = try scanZfi(allocator, broken);
    defer broken_index.deinit(allocator);
    try expectHasSideTable(&broken_index);

    const fixed = try fixFormatOnly(allocator, broken);
    try expectValidateClean(allocator, fixed);

    var fixed_index = try scanZfi(allocator, fixed);
    defer fixed_index.deinit(allocator);
    try expectAllUniformNoSideTable(&fixed_index);

    const path = try writeFastaArtifact(allocator, "validate-fix-crafted", fixed);
    const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{path});
    defer std.Io.Dir.cwd().deleteFile(io, path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, zfi_path) catch {};

    try writeZfiForData(allocator, path, fixed);

    const got = try captureExtractRegion(allocator, path, "messy_seq:1-14");
    const expected =
        \\>messy_seq:1-14
        \\AAAACCCCGGGGTT
        \\
    ;
    try std.testing.expectEqualStrings(expected, got);
}

test "validate --fix then index then get on mixed_line_widths fixture" {
    const broken = try readTestFile(
        std.testing.allocator,
        "bench/index/messy_variants/mixed_line_widths.fasta",
    );
    defer std.testing.allocator.free(broken);

    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var broken_index = try scanZfi(allocator, broken);
    defer broken_index.deinit(allocator);
    try expectHasSideTable(&broken_index);

    const fixed = try fixFormatOnly(allocator, broken);
    try expectValidateClean(allocator, fixed);

    var fixed_index = try scanZfi(allocator, fixed);
    defer fixed_index.deinit(allocator);
    try expectAllUniformNoSideTable(&fixed_index);

    const messy_path = try writeFastaArtifact(allocator, "validate-fix-messy", broken);
    const fixed_path = try writeFastaArtifact(allocator, "validate-fix-fixed", fixed);
    const messy_zfi = try std.fmt.allocPrint(allocator, "{s}.zfi", .{messy_path});
    const fixed_zfi = try std.fmt.allocPrint(allocator, "{s}.zfi", .{fixed_path});
    defer std.Io.Dir.cwd().deleteFile(io, messy_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, fixed_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, messy_zfi) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, fixed_zfi) catch {};

    try writeZfiForData(allocator, messy_path, broken);
    try writeZfiForData(allocator, fixed_path, fixed);

    const region = "mixed_line_widths:3-24";
    const messy_out = try captureExtractRegion(allocator, messy_path, region);
    const fixed_out = try captureExtractRegion(allocator, fixed_path, region);
    try std.testing.expectEqualStrings(messy_out, fixed_out);

    const expected =
        \\>mixed_line_widths:3-24
        \\AACCCCGGGGTTTTAAAACCCC
        \\
    ;
    try std.testing.expectEqualStrings(expected, fixed_out);
}

test "validate --fix then index then get on trailing_whitespace fixture" {
    const broken = try readTestFile(
        std.testing.allocator,
        "bench/index/messy_variants/trailing_whitespace.fasta",
    );
    defer std.testing.allocator.free(broken);

    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var broken_index = try scanZfi(allocator, broken);
    defer broken_index.deinit(allocator);
    try expectHasSideTable(&broken_index);

    const fixed = try fixFormatOnly(allocator, broken);
    try expectValidateClean(allocator, fixed);

    var fixed_index = try scanZfi(allocator, fixed);
    defer fixed_index.deinit(allocator);
    try expectAllUniformNoSideTable(&fixed_index);

    const fixed_path = try writeFastaArtifact(allocator, "validate-fix-trail", fixed);
    const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{fixed_path});
    defer std.Io.Dir.cwd().deleteFile(io, fixed_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, zfi_path) catch {};

    try writeZfiForData(allocator, fixed_path, fixed);

    const got = try captureExtractRegion(allocator, fixed_path, "trailing_whitespace:1-16");
    const expected =
        \\>trailing_whitespace:1-16
        \\AAAACCCCGGGGTTTT
        \\
    ;
    try std.testing.expectEqualStrings(expected, got);
}
