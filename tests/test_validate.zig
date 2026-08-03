//! Validate unit and CLI tests: issue kinds, JSON rendering, fix flow, and gates fixtures.
//!
//! Also covers validator/indexer agreement on `tests/data/validator_indexer_agreement.fasta`.

const std = @import("std");
const builtin = @import("builtin");
const main = @import("main");

const io = std.testing.io;
const validator = main.validator;

const ZFASTA_BIN = if (builtin.os.tag == .windows) "zig-out\\bin\\z-fasta.exe" else "zig-out/bin/z-fasta";

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

test "validateData caps retained events and sets truncated" {
    const allocator = std.testing.allocator;

    var at_cap: std.ArrayList(u8) = .empty;
    defer at_cap.deinit(allocator);
    try at_cap.appendSlice(allocator, ">seq\n");
    var i: usize = 0;
    while (i < validator.max_validate_events) : (i += 1) {
        try at_cap.appendSlice(allocator, "ACGT \n");
    }

    var capped = try validator.validateData(allocator, at_cap.items, .{});
    defer capped.deinit(allocator);
    try std.testing.expect(!capped.truncated);
    try std.testing.expectEqual(validator.max_validate_events, capped.events.items.len);

    try at_cap.appendSlice(allocator, "ACGT \n");
    var over = try validator.validateData(allocator, at_cap.items, .{});
    defer over.deinit(allocator);
    try std.testing.expect(over.truncated);
    try std.testing.expectEqual(validator.max_validate_events, over.events.items.len);
}

test "renderJsonEvent keeps long names without 256-byte truncation" {
    const allocator = std.testing.allocator;
    const long_name = "n" ** 400;

    const json = try validator.renderJsonEvent(allocator, .{
        .level = .warning,
        .kind = .empty_sequence,
        .line = 1,
        .name = long_name,
    });
    defer allocator.free(json);

    try std.testing.expect(std.unicode.utf8ValidateSlice(json));
    const parsed = try std.json.parseFromSlice(std.json.Value, allocator, json, .{});
    defer parsed.deinit();
    const message = parsed.value.object.get("message").?.string;
    try std.testing.expect(std.mem.indexOf(u8, message, long_name) != null);
    try std.testing.expectEqualStrings(long_name, parsed.value.object.get("name").?.string);
}

test "renderJsonEvent preserves valid non-ASCII UTF-8 names" {
    const allocator = std.testing.allocator;
    const name = "seq_\u{20ac}_α";

    const json = try validator.renderJsonEvent(allocator, .{
        .level = .warning,
        .kind = .empty_sequence,
        .line = 1,
        .name = name,
    });
    defer allocator.free(json);

    try std.testing.expect(std.unicode.utf8ValidateSlice(json));
    const parsed = try std.json.parseFromSlice(std.json.Value, allocator, json, .{});
    defer parsed.deinit();
    try std.testing.expectEqualStrings(name, parsed.value.object.get("name").?.string);
}

test "renderJsonEvent escapes quotes and backslashes in names" {
    const allocator = std.testing.allocator;
    const name = "seq\"a\\b";

    const json = try validator.renderJsonEvent(allocator, .{
        .level = .warning,
        .kind = .empty_sequence,
        .line = 1,
        .name = name,
    });
    defer allocator.free(json);

    try std.testing.expect(std.unicode.utf8ValidateSlice(json));
    try std.testing.expect(std.mem.indexOf(u8, json, "\\\"") != null);
    try std.testing.expect(std.mem.indexOf(u8, json, "\\\\") != null);
    const parsed = try std.json.parseFromSlice(std.json.Value, allocator, json, .{});
    defer parsed.deinit();
    try std.testing.expectEqualStrings(name, parsed.value.object.get("name").?.string);
}

test "renderJsonEvent escapes invalid UTF-8 name bytes" {
    const allocator = std.testing.allocator;
    const bad_name = "seq\xffname";

    const json = try validator.renderJsonEvent(allocator, .{
        .level = .error_level,
        .kind = .duplicate_name,
        .line = 3,
        .name = bad_name,
        .first_line = 1,
    });
    defer allocator.free(json);

    try std.testing.expect(std.unicode.utf8ValidateSlice(json));
    try std.testing.expect(std.mem.indexOf(u8, json, "\\u00ff") != null);

    const parsed = try std.json.parseFromSlice(std.json.Value, allocator, json, .{});
    defer parsed.deinit();
    // JSON \u00ff decodes to Unicode U+00FF (ÿ), not the raw invalid FASTA byte.
    try std.testing.expectEqualStrings("seq\u{00ff}name", parsed.value.object.get("name").?.string);
}

test "validateData to JSON keeps invalid UTF-8 header bytes escapable" {
    const allocator = std.testing.allocator;
    const fasta = ">bad\xffname\nACGT\n>bad\xffname\nTTTT\n";

    var summary = try validator.validateData(allocator, fasta, .{});
    defer summary.deinit(allocator);
    try std.testing.expectEqual(@as(usize, 1), countKind(&summary, .duplicate_name));

    const event = summary.events.items[0];
    try std.testing.expectEqual(validator.Kind.duplicate_name, event.kind);
    try std.testing.expectEqualSlices(u8, "bad\xffname", event.name);

    const json = try validator.renderJsonEvent(allocator, event);
    defer allocator.free(json);
    try std.testing.expect(std.unicode.utf8ValidateSlice(json));
    const parsed = try std.json.parseFromSlice(std.json.Value, allocator, json, .{});
    defer parsed.deinit();
    try std.testing.expectEqualStrings("bad\u{00ff}name", parsed.value.object.get("name").?.string);
}

test "gates fixture: invalid UTF-8 and long header are visible to validate" {
    const allocator = std.testing.allocator;

    const utf8_data = try readTestFile(allocator, "tests/data/gates/invalid_utf8_header.fasta");
    defer allocator.free(utf8_data);
    var utf8_summary = try validator.validateData(allocator, utf8_data, .{});
    defer utf8_summary.deinit(allocator);
    try std.testing.expectEqual(@as(usize, 2), utf8_summary.header_count);
    try std.testing.expectEqual(@as(usize, 1), countKind(&utf8_summary, .duplicate_name));
    const event = utf8_summary.events.items[0];
    try std.testing.expectEqual(validator.Kind.duplicate_name, event.kind);
    try std.testing.expectEqualSlices(u8, "bad\xffname", event.name);
    const json = try validator.renderJsonEvent(allocator, event);
    defer allocator.free(json);
    try std.testing.expect(std.unicode.utf8ValidateSlice(json));
    try std.testing.expect(std.mem.indexOf(u8, json, "\\u00ff") != null);

    const long_data = try readTestFile(allocator, "tests/data/gates/long_header.fasta");
    defer allocator.free(long_data);
    var long_summary = try validator.validateData(allocator, long_data, .{});
    defer long_summary.deinit(allocator);
    try std.testing.expect(countKind(&long_summary, .long_header) >= 1);
}

test "validateData records type sample metadata" {
    var summary = try validator.validateData(std.testing.allocator, ">seq\nACGT\n", .{});
    defer summary.deinit(std.testing.allocator);

    try std.testing.expectEqual(main.stats.SequenceType.nucleotide, summary.sequence_type);
    try std.testing.expectEqual(@as(u64, 4), summary.type_bases_sampled);
}

test "validate --json --summary reports sequence type sample fields" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_path = try writeFastaArtifact(allocator, "validate-type-json", ">seq\nACGT\n");
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};

    var threaded = std.Io.Threaded.init(allocator, .{});
    defer threaded.deinit();
    const spawn_io = threaded.io();

    const result = try std.process.run(allocator, spawn_io, .{
        .argv = &.{ ZFASTA_BIN, "validate", "--json", "--summary", fasta_path },
        .stdout_limit = .limited(64 * 1024),
        .stderr_limit = .limited(64 * 1024),
    });
    defer allocator.free(result.stdout);
    defer allocator.free(result.stderr);

    switch (result.term) {
        .exited => |code| try std.testing.expectEqual(@as(u8, 0), code),
        else => return error.ChildProcessFailed,
    }

    const parsed = try std.json.parseFromSlice(std.json.Value, allocator, std.mem.trim(u8, result.stdout, " \t\r\n"), .{});
    defer parsed.deinit();
    try std.testing.expectEqualStrings("nucleotide", parsed.value.object.get("sequence_type").?.string);
    try std.testing.expectEqual(@as(i64, 4), parsed.value.object.get("type_bases_sampled").?.integer);
    try std.testing.expectEqual(
        @as(i64, @intCast(main.stats.validate_type_sample_bases)),
        parsed.value.object.get("type_sample_cap").?.integer,
    );
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
    const file = std.Io.Dir.cwd().openFile(io, path, .{}) catch |err| switch (err) {
        error.FileNotFound => {
            std.debug.print("required test fixture missing: {s}\n", .{path});
            return error.RequiredFixtureMissing;
        },
        else => |e| return e,
    };
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
    return main.indexer.scanZfiData(data, true, allocator);
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

    const fasta_file = try std.Io.Dir.cwd().openFile(io, fasta_path, .{});
    defer fasta_file.close(io);
    const mtime_ns = main.index_format.timestampToNs((try fasta_file.stat(io)).mtime);

    var zfi_path_buf: [4096]u8 = undefined;
    const zfi_path = try std.fmt.bufPrint(&zfi_path_buf, "{s}.zfi", .{fasta_path});
    try main.indexer.writeZfiIndexFile(io, zfi_path, &index, data.len, mtime_ns);
}

fn captureExtractRegion(allocator: std.mem.Allocator, fasta_path: []const u8, region: []const u8) ![]u8 {
    var idx = try main.index_format.loadIndexChecked(io, fasta_path);
    defer idx.deinit(io);

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

test "validate --fix -o matches fixData rewrite" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const broken =
        \\>seq
        \\AAAA  
        \\CCCC
        \\
    ;

    var summary = try validator.validateData(allocator, broken, .{});
    defer summary.deinit(allocator);
    const expected = try validator.fixData(allocator, broken, summary.record_widths.items);

    const in_path = try writeFastaArtifact(allocator, "validate-fix-cli-in", broken);
    const out_path = try uniqueArtifactPath(allocator, "validate-fix-cli-out");
    defer std.Io.Dir.cwd().deleteFile(io, in_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, out_path) catch {};

    var threaded = std.Io.Threaded.init(allocator, .{});
    defer threaded.deinit();
    const spawn_io = threaded.io();

    const result = try std.process.run(allocator, spawn_io, .{
        .argv = &.{ ZFASTA_BIN, "validate", "--fix", "-o", out_path, in_path },
        .stdout_limit = .limited(1024 * 1024),
        .stderr_limit = .limited(1024 * 1024),
    });
    defer allocator.free(result.stdout);
    defer allocator.free(result.stderr);

    switch (result.term) {
        .exited => |code| try std.testing.expectEqual(@as(u8, 0), code),
        else => return error.ChildProcessFailed,
    }

    const got = try readTestFile(allocator, out_path);
    try std.testing.expectEqualStrings(expected, got);
}

test "validate --fix succeeds after event list truncation" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var input: std.ArrayList(u8) = .empty;
    defer input.deinit(allocator);
    try input.appendSlice(allocator, ">seq\n");
    var i: usize = 0;
    while (i <= validator.max_validate_events) : (i += 1) {
        try input.appendSlice(allocator, "A \n");
    }

    var summary = try validator.validateData(allocator, input.items, .{});
    defer summary.deinit(allocator);
    try std.testing.expect(summary.truncated);

    const in_path = try writeFastaArtifact(allocator, "validate-fix-truncated-in", input.items);
    const out_path = try uniqueArtifactPath(allocator, "validate-fix-truncated-out");
    defer std.Io.Dir.cwd().deleteFile(io, in_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, out_path) catch {};

    var threaded = std.Io.Threaded.init(allocator, .{});
    defer threaded.deinit();
    const spawn_io = threaded.io();

    const result = try std.process.run(allocator, spawn_io, .{
        .argv = &.{ ZFASTA_BIN, "validate", "--fix", "-o", out_path, in_path },
        .stdout_limit = .limited(4 * 1024 * 1024),
        .stderr_limit = .limited(64 * 1024),
    });
    defer allocator.free(result.stdout);
    defer allocator.free(result.stderr);

    switch (result.term) {
        .exited => |code| try std.testing.expectEqual(@as(u8, 0), code),
        else => return error.ChildProcessFailed,
    }

    const fixed = try readTestFile(allocator, out_path);
    var fixed_summary = try validator.validateData(allocator, fixed, .{});
    defer fixed_summary.deinit(allocator);
    try std.testing.expectEqual(@as(usize, 0), fixed_summary.events.items.len);
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

test "validate --fix then index then get on mixed_widths fixture" {
    const broken = try readTestFile(
        std.testing.allocator,
        "bench/shared/cache/messy_fixtures/mixed_widths.fasta",
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

    const region = "mixed_widths:3-24";
    const messy_out = try captureExtractRegion(allocator, messy_path, region);
    const fixed_out = try captureExtractRegion(allocator, fixed_path, region);
    try std.testing.expectEqualStrings(messy_out, fixed_out);

    const expected =
        \\>mixed_widths:3-24
        \\AACCCCGGGGTTTTAAAACCCC
        \\
    ;
    try std.testing.expectEqualStrings(expected, fixed_out);
}

test "validate --fix then index then get on trailing_whitespace fixture" {
    const broken = try readTestFile(
        std.testing.allocator,
        "bench/shared/cache/messy_fixtures/trailing_whitespace.fasta",
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

fn zfiEmbeddedName(index: *const main.indexer.ZfiIndex, rec_idx: usize) []const u8 {
    const rec = index.records.items[rec_idx];
    return index.name_blob.items[rec.name_offset..][0..rec.name_len];
}

test "validator and indexer agree on tests/data/validator_indexer_agreement.fasta" {
    const data = try readTestFile(std.testing.allocator, "tests/data/validator_indexer_agreement.fasta");
    defer std.testing.allocator.free(data);

    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    // Fixture bytes must actually contain each layout feature (not only related warnings).
    try std.testing.expect(std.mem.indexOf(u8, data, "\r\n") != null);
    try std.testing.expect(std.mem.indexOf(u8, data, "\n\n") != null);
    try std.testing.expect(!std.mem.endsWith(u8, data, "\n"));
    try std.testing.expect(!std.mem.endsWith(u8, data, "\r\n"));
    try std.testing.expect(std.mem.indexOf(u8, data, ">empty_rec\n") != null);
    try std.testing.expect(std.mem.indexOf(u8, data, " \n") != null or std.mem.indexOf(u8, data, "\t\n") != null);

    // --- validator: every listed layout issue is visible ---
    var summary = try validator.validateData(allocator, data, .{});
    defer summary.deinit(allocator);

    try std.testing.expect(countKind(&summary, .empty_sequence) >= 1);
    try std.testing.expect(countKind(&summary, .inconsistent_line_widths) >= 1);
    try std.testing.expect(countKind(&summary, .trailing_whitespace) >= 1);
    try std.testing.expect(countKind(&summary, .mixed_line_endings) >= 1);
    try std.testing.expect(countKind(&summary, .missing_terminal_newline) >= 1);
    try std.testing.expect(countKind(&summary, .long_header) >= 1);

    try std.testing.expectEqual(@as(usize, 6), summary.header_count);
    try std.testing.expectEqual(@as(usize, 5), summary.sequence_count);

    // --- indexer: same non-empty catalog, side tables for messy layout ---
    var index = try scanZfi(allocator, data);
    defer index.deinit(allocator);

    try std.testing.expectEqual(summary.sequence_count, index.records.items.len);
    try expectHasSideTable(&index);

    const expected_names = [_][]const u8{
        "var_widths",
        "blank_lines",
        "trailing_ws",
        "mixed_crlf",
    };
    for (expected_names, 0..) |want, i| {
        try std.testing.expectEqualStrings(want, zfiEmbeddedName(&index, i));
    }
    // Empty source record is validated but not indexed.
    for (0..index.records.items.len) |i| {
        try std.testing.expect(!std.mem.eql(u8, zfiEmbeddedName(&index, i), "empty_rec"));
    }
    // Near validate default --max-header-len (1024): still indexable, triggers long_header.
    const long_name = zfiEmbeddedName(&index, 4);
    try std.testing.expect(long_name.len > 1024);
    try std.testing.expectEqual(@as(usize, 1025), long_name.len);
}
