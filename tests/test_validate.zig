//! Validate tests: issue kinds, JSON output, `--fix`, and agreement with the indexer.

const std = @import("std");
const main = @import("main");
const utility = @import("utility.zig");

const io = std.testing.io;
const validator = main.validator;

const ZFASTA_BIN = utility.ZFASTA_BIN;
const expectCliResult = utility.expectCliResult;
const readTestFile = utility.readRequiredFile;
const writeFastaArtifact = utility.writeFastaArtifact;
const writeZfi = utility.writeZfi;
const captureExtractRegion = utility.captureExtractRegion;

fn validateBytes(allocator: std.mem.Allocator, data: []const u8, options: validator.Options) main.index_format.diagnostics.Error!validator.Summary {
    var failure: main.index_format.diagnostics.Failure = undefined;
    return validator.validateData(allocator, data, options, &failure);
}

fn fixBytes(allocator: std.mem.Allocator, data: []const u8, record_widths: []const u32) main.index_format.diagnostics.Error![]u8 {
    var failure: main.index_format.diagnostics.Failure = undefined;
    return validator.fixData(allocator, data, record_widths, &failure);
}

fn renderEvent(allocator: std.mem.Allocator, event: validator.Event) main.index_format.diagnostics.Error![]u8 {
    var failure: main.index_format.diagnostics.Failure = undefined;
    return validator.renderJsonEvent(allocator, event, &failure);
}

/// The blocker `checkFixAllowed` records, or null when the fix may proceed.
fn fixBlocker(summary: *const validator.Summary, options: validator.Options) !?std.meta.Tag(main.index_format.diagnostics.Blocker) {
    var failure: main.index_format.diagnostics.Failure = undefined;
    failure.reset(.validate);
    validator.checkFixAllowed(summary, options, &failure) catch |err| {
        if (err != error.RepairRefused) return err;
        return std.meta.activeTag(failure.primary().?.details.@"repair.blocked".blocker);
    };
    return null;
}

fn countKind(summary: *const validator.Summary, kind: validator.Kind) usize {
    var count: usize = 0;
    for (summary.events.items) |event| {
        if (event.kind == kind) count += 1;
    }
    return count;
}

fn expectFixIdempotent(allocator: std.mem.Allocator, broken: []const u8) !void {
    var before = try validateBytes(allocator, broken, .{});
    defer before.deinit(allocator);
    try std.testing.expectEqual(null, try fixBlocker(&before, .{}));

    const fixed = try fixBytes(allocator, broken, before.record_widths.items);
    defer allocator.free(fixed);

    var after = try validateBytes(allocator, fixed, .{});
    defer after.deinit(allocator);
    try std.testing.expectEqual(@as(usize, 0), after.events.items.len);

    const fixed_again = try fixBytes(allocator, fixed, after.record_widths.items);
    defer allocator.free(fixed_again);
    try std.testing.expectEqualStrings(fixed, fixed_again);
}

fn validateAndDeinitForAllocationCheck(allocator: std.mem.Allocator, data: []const u8) !void {
    var summary = try validateBytes(allocator, data, .{});
    defer summary.deinit(allocator);
}

fn fixAndFreeForAllocationCheck(allocator: std.mem.Allocator, data: []const u8) !void {
    const widths = [_]u32{4};
    const fixed = try fixBytes(allocator, data, &widths);
    defer allocator.free(fixed);
}

fn renderAndFreeForAllocationCheck(allocator: std.mem.Allocator, name: []const u8) !void {
    const json = try renderEvent(allocator, .{
        .level = .error_level,
        .kind = .duplicate_name,
        .line = 2,
        .name = name,
        .first_line = 1,
    });
    defer allocator.free(json);
}

fn fuzzValidator(_: void, smith: *std.testing.Smith) !void {
    var storage: [512]u8 = undefined;
    try checkValidatorProperties(storage[0..smith.slice(&storage)]);
}

fn checkValidatorProperties(data: []const u8) !void {
    var summary = try validateBytes(std.testing.allocator, data, .{});
    defer summary.deinit(std.testing.allocator);

    try std.testing.expect(summary.events.items.len <= validator.MAX_VALIDATE_EVENTS);
    try std.testing.expect(summary.sequence_count <= summary.header_count);
    try std.testing.expectEqual(summary.header_count, summary.record_widths.items.len);
    var total_kind_count: usize = 0;
    for (summary.kind_counts) |kind_count| total_kind_count += kind_count;
    try std.testing.expectEqual(summary.error_count + summary.warning_count, total_kind_count);
    try std.testing.expect(summary.format_warning_count <= summary.warning_count);
    for (summary.record_widths.items) |width| try std.testing.expect(width > 0);

    var failure: main.index_format.diagnostics.Failure = undefined;
    failure.reset(.validate);
    validator.checkFixAllowed(&summary, .{}, &failure) catch |err| switch (err) {
        error.RepairRefused => return,
        else => return err,
    };
    const fixed = try fixBytes(std.testing.allocator, data, summary.record_widths.items);
    defer std.testing.allocator.free(fixed);

    var fixed_summary = try validateBytes(std.testing.allocator, fixed, .{});
    defer fixed_summary.deinit(std.testing.allocator);
    failure.reset(.validate);
    try validator.checkFixAllowed(&fixed_summary, .{}, &failure);
    const fixed_again = try fixBytes(
        std.testing.allocator,
        fixed,
        fixed_summary.record_widths.items,
    );
    defer std.testing.allocator.free(fixed_again);
    try std.testing.expectEqualStrings(fixed, fixed_again);
}

fn scanZfi(allocator: std.mem.Allocator, data: []const u8) !main.indexer.ZfiIndex {
    return main.indexer.scanZfiData(allocator, data, true);
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

fn fixFormatOnly(allocator: std.mem.Allocator, broken: []const u8) ![]u8 {
    var summary = try validateBytes(allocator, broken, .{});
    defer summary.deinit(allocator);
    try std.testing.expectEqual(null, try fixBlocker(&summary, .{}));
    return fixBytes(allocator, broken, summary.record_widths.items);
}

fn expectValidateClean(allocator: std.mem.Allocator, data: []const u8) !void {
    var summary = try validateBytes(allocator, data, .{});
    defer summary.deinit(allocator);
    try std.testing.expectEqual(@as(usize, 0), summary.events.items.len);
}

fn expectFormatFixPreservesRegion(
    allocator: std.mem.Allocator,
    stem: []const u8,
    broken: []const u8,
    region: []const u8,
    expected: []const u8,
) !void {
    var broken_index = try scanZfi(allocator, broken);
    defer broken_index.deinit(allocator);
    try expectHasSideTable(&broken_index);

    const fixed = try fixFormatOnly(allocator, broken);
    defer allocator.free(fixed);
    try expectValidateClean(allocator, fixed);

    var fixed_index = try scanZfi(allocator, fixed);
    defer fixed_index.deinit(allocator);
    try expectAllUniformNoSideTable(&fixed_index);

    const broken_stem = try std.fmt.allocPrint(allocator, "{s}-broken", .{stem});
    defer allocator.free(broken_stem);
    const fixed_stem = try std.fmt.allocPrint(allocator, "{s}-fixed", .{stem});
    defer allocator.free(fixed_stem);
    const broken_path = try writeFastaArtifact(allocator, broken_stem, broken);
    defer allocator.free(broken_path);
    defer std.Io.Dir.cwd().deleteFile(io, broken_path) catch {};
    const fixed_path = try writeFastaArtifact(allocator, fixed_stem, fixed);
    defer allocator.free(fixed_path);
    defer std.Io.Dir.cwd().deleteFile(io, fixed_path) catch {};
    const broken_zfi = try std.fmt.allocPrint(allocator, "{s}.zfi", .{broken_path});
    defer allocator.free(broken_zfi);
    defer std.Io.Dir.cwd().deleteFile(io, broken_zfi) catch {};
    const fixed_zfi = try std.fmt.allocPrint(allocator, "{s}.zfi", .{fixed_path});
    defer allocator.free(fixed_zfi);
    defer std.Io.Dir.cwd().deleteFile(io, fixed_zfi) catch {};

    try writeZfi(allocator, broken_path, broken, true);
    try writeZfi(allocator, fixed_path, fixed, true);

    const broken_output = try captureExtractRegion(allocator, broken_path, region);
    defer allocator.free(broken_output);
    const fixed_output = try captureExtractRegion(allocator, fixed_path, region);
    defer allocator.free(fixed_output);

    try std.testing.expectEqualStrings(broken_output, fixed_output);
    try std.testing.expectEqualStrings(expected, fixed_output);
}

fn zfiEmbeddedName(index: *const main.indexer.ZfiIndex, rec_idx: usize) []const u8 {
    const rec = index.records.items[rec_idx];
    return index.name_blob.items[rec.name_offset..][0..rec.name_len];
}

fn checkNativeValidatorAllocations(allocator: std.mem.Allocator, data: []const u8) !void {
    const diagnostics = main.index_format.diagnostics;
    var failure: diagnostics.Failure = undefined;
    var summary = validator.validateData(allocator, data, .{}, &failure) catch |err| {
        try std.testing.expectEqual(error.OutOfMemory, err);
        try std.testing.expectEqual(diagnostics.Code.@"resource.out_of_memory", failure.primary().?.code());
        try std.testing.expectEqual(error.OutOfMemory, failure.primary().?.cause.?.zig);
        return error.OutOfMemory;
    };
    defer summary.deinit(allocator);
    const fixed = validator.fixData(allocator, data, summary.record_widths.items, &failure) catch |err| {
        try std.testing.expectEqual(error.OutOfMemory, err);
        try std.testing.expectEqual(diagnostics.AllocationRole.rewrite_bytes, failure.primary().?.details.@"resource.out_of_memory".allocation_role);
        return error.OutOfMemory;
    };
    defer allocator.free(fixed);
    const json = validator.renderJsonEvent(allocator, .{ .level = .warning, .kind = .empty_sequence, .line = 1, .name = "seq" }, &failure) catch |err| {
        try std.testing.expectEqual(error.OutOfMemory, err);
        try std.testing.expectEqual(diagnostics.AllocationRole.serialized_event, failure.primary().?.details.@"resource.out_of_memory".allocation_role);
        return error.OutOfMemory;
    };
    defer allocator.free(json);
}

// --- Validation and report contracts ---
test "[cli] - [validate]: rejects unknown options regardless of position" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const fasta = "tests/data/simple.fasta";

    try expectCliResult(
        allocator,
        &.{ ZFASTA_BIN, "validate", "--not-a-flag", fasta },
        1,
        "",
        "error [usage.unknown_option] validate: Unknown option. option=\"--not-a-flag\"\nhint: Correct the command arguments; use --help for accepted options.\n",
    );
    try expectCliResult(
        allocator,
        &.{ ZFASTA_BIN, "validate", fasta, "--not-a-flag" },
        1,
        "",
        "error [usage.unknown_option] validate: Unknown option. option=\"--not-a-flag\"\nhint: Correct the command arguments; use --help for accepted options.\n",
    );
}

test "[cli] - [validate]: rejects invalid arguments with exact diagnostics" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const cases = [_]struct { argv: []const []const u8, stderr: []const u8 }{
        .{
            .argv = &.{ ZFASTA_BIN, "validate", "--summary", "tests/data/simple.fasta" },
            .stderr = "error [usage.option_requires] validate: An option prerequisite is missing. option=summary required=json\nhint: Correct the command arguments; use --help for accepted options.\n",
        },
        .{
            .argv = &.{ ZFASTA_BIN, "validate", "--fix-format-only", "tests/data/simple.fasta" },
            .stderr = "error [usage.option_requires] validate: An option prerequisite is missing. option=fix_format_only required=fix\nhint: Correct the command arguments; use --help for accepted options.\n",
        },
        .{
            .argv = &.{ ZFASTA_BIN, "validate", "-o", "zig-cache/test-artifacts/unused.fasta", "tests/data/simple.fasta" },
            .stderr = "error [usage.option_requires] validate: An option prerequisite is missing. option=output required=fix\nhint: Correct the command arguments; use --help for accepted options.\n",
        },
        .{
            .argv = &.{ ZFASTA_BIN, "validate" },
            .stderr = "error [usage.missing_argument] validate: A required argument is missing. argument=fasta_path\nhint: Correct the command arguments; use --help for accepted options.\n",
        },
        .{
            .argv = &.{ ZFASTA_BIN, "validate", "tests/data/definitely-missing-cli-failure.fasta" },
            .stderr = "error [io.file_not_found] validate: \"tests/data/definitely-mi\" (47 bytes): Required file was not found.; cause=\"FileNotFound\"\nhint: Check the named input/output and retained I/O cause.\n",
        },
        .{
            .argv = &.{ ZFASTA_BIN, "validate", "tests/data/simple.fasta/child" },
            .stderr = "error [io.open_failed] validate: \"tests/data/simple.fasta/\" (29 bytes): Could not open the required input.; cause=\"NotDir\"\nhint: Check the named input/output and retained I/O cause.\n",
        },
        .{
            .argv = &.{ ZFASTA_BIN, "validate", "tests/data/simple.fasta", "tests/data/single.fasta" },
            .stderr = "error [usage.extra_argument] validate: Too many arguments. argument=fasta_path allowed_count=1\nhint: Correct the command arguments; use --help for accepted options.\n",
        },
        .{
            .argv = &.{ ZFASTA_BIN, "validate", "-o" },
            .stderr = "error [usage.missing_option_value] validate: A required option value is missing. option=output expected=path\nhint: Correct the command arguments; use --help for accepted options.\n",
        },
        .{
            .argv = &.{ ZFASTA_BIN, "validate", "--schema" },
            .stderr = "error [usage.missing_option_value] validate: A required option value is missing. option=schema expected=schema\nhint: Correct the command arguments; use --help for accepted options.\n",
        },
        .{
            .argv = &.{ ZFASTA_BIN, "validate", "--schema", "ena", "tests/data/simple.fasta" },
            .stderr = "error [usage.invalid_option_value] validate: Invalid option value. option=schema expected=schema value=\"ena\"\nhint: Correct the command arguments; use --help for accepted options.\n",
        },
        .{
            .argv = &.{ ZFASTA_BIN, "validate", "--custom-alphabet" },
            .stderr = "error [usage.missing_option_value] validate: A required option value is missing. option=custom_alphabet expected=alphabet\nhint: Correct the command arguments; use --help for accepted options.\n",
        },
        .{
            .argv = &.{ ZFASTA_BIN, "validate", "--max-header-len" },
            .stderr = "error [usage.missing_option_value] validate: A required option value is missing. option=max_header_len expected=positive_usize\nhint: Correct the command arguments; use --help for accepted options.\n",
        },
        .{
            .argv = &.{ ZFASTA_BIN, "validate", "--max-header-len", "0", "tests/data/simple.fasta" },
            .stderr = "error [usage.invalid_option_value] validate: Invalid option value. option=max_header_len expected=positive_usize value=\"0\"\nhint: Correct the command arguments; use --help for accepted options.\n",
        },
        .{
            .argv = &.{ ZFASTA_BIN, "validate", "--fix", "tests/data/simple.fasta" },
            .stderr = "error [usage.option_requires] validate: \"tests/data/simple.fasta\": An option prerequisite is missing. option=fix required=output\nhint: Correct the command arguments; use --help for accepted options.\n",
        },
    };
    for (cases) |case| {
        try expectCliResult(allocator, case.argv, 1, "", case.stderr);
    }
}

test "[cli] - [validate report]: returns exact status and output for clean and warning reports" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const path = try writeFastaArtifact(allocator, "cli-validate-warn", ">empty_rec\n");
    defer allocator.free(path);
    defer std.Io.Dir.cwd().deleteFile(io, path) catch {};
    const invalid_path = try writeFastaArtifact(allocator, "cli-validate-error", "not FASTA\n");
    defer allocator.free(invalid_path);
    defer std.Io.Dir.cwd().deleteFile(io, invalid_path) catch {};

    try expectCliResult(
        allocator,
        &.{ ZFASTA_BIN, "validate", path },
        2,
        "WARNING: line 1: empty sequence 'empty_rec'\n",
        "",
    );
    try expectCliResult(
        allocator,
        &.{ ZFASTA_BIN, "validate", "--strict", path },
        1,
        "WARNING: line 1: empty sequence 'empty_rec'\n",
        "",
    );
    try expectCliResult(
        allocator,
        &.{ ZFASTA_BIN, "validate", "--json", path },
        2,
        "{\"schema_version\":\"v1\",\"level\":\"warning\",\"line\":1,\"kind\":\"empty_sequence\",\"message\":\"empty sequence 'empty_rec'\",\"name\":\"empty_rec\"}\n",
        "",
    );
    try expectCliResult(
        allocator,
        &.{ ZFASTA_BIN, "validate", "tests/data/simple.fasta" },
        0,
        "OK: no issues found\n",
        "",
    );
    try expectCliResult(
        allocator,
        &.{ ZFASTA_BIN, "validate", invalid_path },
        1,
        "ERROR: line 1: no sequences found\n",
        "",
    );
}

test "[cli] - [validate options]: applies schema, alphabet, and header limits together" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const path = try writeFastaArtifact(allocator, "cli-validate-options", ">bad\nACGTZ\n");
    defer std.Io.Dir.cwd().deleteFile(io, path) catch {};

    try expectCliResult(
        allocator,
        &.{
            ZFASTA_BIN,          "validate",
            "--schema",          "refseq",
            "--custom-alphabet", "ACGTZ",
            "--max-header-len",  "2",
            path,
        },
        2,
        "WARNING: line 1: header exceeds 2 bytes\n" ++
            "WARNING: line 1: header for 'bad' does not match schema\n",
        "",
    );
}

test "[unit] - [validate scan]: reports duplicate names and empty records" {
    var summary = try validateBytes(std.testing.allocator, ">dup\nAAAA\n>dup\n", .{});
    defer summary.deinit(std.testing.allocator);

    try std.testing.expectEqual(@as(usize, 1), countKind(&summary, .duplicate_name));
    try std.testing.expectEqual(@as(usize, 1), countKind(&summary, .empty_sequence));
    try std.testing.expectEqual(@as(usize, 0), countKind(&summary, .no_sequences));
    try std.testing.expectEqual(@as(usize, 1), summary.sequence_count);
    try std.testing.expectEqual(@as(usize, 2), summary.header_count);
}

test "[edge] - [validate scan]: reports the complete empty-file contract" {
    var summary = try validateBytes(std.testing.allocator, "", .{});
    defer summary.deinit(std.testing.allocator);

    try std.testing.expectEqual(@as(usize, 0), summary.header_count);
    try std.testing.expectEqual(@as(usize, 0), summary.sequence_count);
    try std.testing.expectEqual(@as(usize, 0), summary.record_widths.items.len);
    try std.testing.expectEqual(@as(usize, 1), summary.error_count);
    try std.testing.expectEqual(@as(usize, 1), summary.warning_count);
    try std.testing.expectEqual(@as(usize, 1), countKind(&summary, .no_sequences));
    try std.testing.expectEqual(@as(usize, 1), countKind(&summary, .missing_terminal_newline));
}

test "[failure] - [validate operations]: release partial allocations" {
    try std.testing.checkAllAllocationFailures(
        std.testing.allocator,
        validateAndDeinitForAllocationCheck,
        .{">dup\nAAAA \n>dup\n"},
    );
    try std.testing.checkAllAllocationFailures(
        std.testing.allocator,
        fixAndFreeForAllocationCheck,
        .{">seq\nAAAA \n"},
    );
    try std.testing.checkAllAllocationFailures(
        std.testing.allocator,
        renderAndFreeForAllocationCheck,
        .{"seq\"\\\xff"},
    );
}

test "[regression] - [validate scan]: repair properties accept safe preamble refusal" {
    try checkValidatorProperties("junk\n>seq1\nACGT\n");
}

test "[fuzz] - [validate scan]: returns bounded summaries and idempotent allowed fixes" {
    try std.testing.fuzz({}, fuzzValidator, .{ .corpus = &.{
        "",
        ">a\nA\n",
        ">a description\r\nAC\r\nGT",
        "junk\n>seq1\nACGT\n",
        ">a\nAAA \nCC\nGGGG\n",
        ">bad\xffname\nA\x00C\n",
    } });
}

test "[edge] - [validate events]: caps retained events without losing counts" {
    const allocator = std.testing.allocator;

    var at_cap: std.ArrayList(u8) = .empty;
    defer at_cap.deinit(allocator);
    try at_cap.appendSlice(allocator, ">seq\n");
    var i: usize = 0;
    while (i < validator.MAX_VALIDATE_EVENTS) : (i += 1) {
        try at_cap.appendSlice(allocator, "ACGT \n");
    }

    var capped = try validateBytes(allocator, at_cap.items, .{});
    defer capped.deinit(allocator);
    try std.testing.expect(!capped.truncated);
    try std.testing.expectEqual(validator.MAX_VALIDATE_EVENTS, capped.events.items.len);

    try at_cap.appendSlice(allocator, "ACGT \n");
    var over = try validateBytes(allocator, at_cap.items, .{});
    defer over.deinit(allocator);
    try std.testing.expect(over.truncated);
    try std.testing.expectEqual(validator.MAX_VALIDATE_EVENTS, over.events.items.len);
    try std.testing.expectEqual(
        validator.MAX_VALIDATE_EVENTS + 1,
        over.kind_counts[@intFromEnum(validator.Kind.trailing_whitespace)],
    );
}

test "[edge] - [validate events]: retains fix blockers after truncation" {
    const allocator = std.testing.allocator;
    var input: std.ArrayList(u8) = .empty;
    defer input.deinit(allocator);
    try input.appendSlice(allocator, ">dup\n");
    for (0..validator.MAX_VALIDATE_EVENTS) |_| {
        try input.appendSlice(allocator, "A \n");
    }
    try input.appendSlice(allocator, ">dup\nA\n");

    var summary = try validateBytes(allocator, input.items, .{});
    defer summary.deinit(allocator);

    try std.testing.expect(summary.truncated);
    try std.testing.expectEqual(@as(usize, 1), summary.error_count);
    try std.testing.expectEqual(.duplicate_identifier, (try fixBlocker(&summary, .{})).?);
    try std.testing.expectEqual(.duplicate_identifier, (try fixBlocker(&summary, .{ .fix_format_only = true })).?);
}

test "[edge] - [validate events]: retains warnings left after a truncated format fix" {
    const allocator = std.testing.allocator;
    var input: std.ArrayList(u8) = .empty;
    defer input.deinit(allocator);
    try input.appendSlice(allocator, ">seq\n");
    for (0..validator.MAX_VALIDATE_EVENTS) |_| {
        try input.appendSlice(allocator, "A \n");
    }
    try input.appendSlice(allocator, ">empty\n");

    var summary = try validateBytes(allocator, input.items, .{});
    defer summary.deinit(allocator);

    try std.testing.expect(summary.truncated);
    try std.testing.expectEqual(validator.MAX_VALIDATE_EVENTS + 1, summary.warning_count);
    try std.testing.expectEqual(@as(u8, 2), validator.exitCodeForOptions(&summary, .{ .fix = true }));
    try std.testing.expectEqual(@as(u8, 1), validator.exitCodeForOptions(&summary, .{ .fix = true, .strict = true }));
}

test "[edge] - [validate JSON]: preserves names longer than 256 bytes" {
    const allocator = std.testing.allocator;
    const long_name = "n" ** 400;

    const json = try renderEvent(allocator, .{
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

test "[unit] - [validate JSON]: preserves valid non-ASCII names" {
    const allocator = std.testing.allocator;
    const name = "seq_\u{20ac}_α";

    const json = try renderEvent(allocator, .{
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

test "[unit] - [validate JSON]: escapes JSON syntax in names" {
    const allocator = std.testing.allocator;
    const name = "seq\"a\\b";

    const json = try renderEvent(allocator, .{
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

test "[edge] - [validate JSON]: escapes invalid UTF-8 name bytes" {
    const allocator = std.testing.allocator;
    const bad_name = "seq\xffname";

    const json = try renderEvent(allocator, .{
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
    try std.testing.expectEqualStrings("seq\u{00ff}name", parsed.value.object.get("name").?.string);
}

test "[integration] - [validate JSON]: reports invalid UTF-8 fixture names as valid JSON" {
    const allocator = std.testing.allocator;

    const utf8_data = try readTestFile(allocator, "tests/data/gates/invalid_utf8_header.fasta");
    defer allocator.free(utf8_data);
    var utf8_summary = try validateBytes(allocator, utf8_data, .{});
    defer utf8_summary.deinit(allocator);
    try std.testing.expectEqual(@as(usize, 2), utf8_summary.header_count);
    try std.testing.expectEqual(@as(usize, 1), countKind(&utf8_summary, .duplicate_name));
    const event = utf8_summary.events.items[0];
    try std.testing.expectEqual(validator.Kind.duplicate_name, event.kind);
    try std.testing.expectEqualSlices(u8, "bad\xffname", event.name);
    const json = try renderEvent(allocator, event);
    defer allocator.free(json);
    try std.testing.expect(std.unicode.utf8ValidateSlice(json));
    try std.testing.expect(std.mem.indexOf(u8, json, "\\u00ff") != null);
}

test "[property] - [validate alphabets]: accepts representative nucleotide and protein symbols" {
    const cases = [_]struct {
        sequence: []const u8,
        sequence_type: main.fasta.SequenceType,
    }{
        .{ .sequence = "ACGTNRYWSMKHBVDacgtnrywsmkhbvd", .sequence_type = .nucleotide },
        .{ .sequence = "ACGUNRYWSMKHBVDacgunrywsmkhbvd", .sequence_type = .nucleotide },
        .{ .sequence = "ACGTUacgtu", .sequence_type = .nucleotide },
        .{ .sequence = "EFILPQBZJXUO*-efilpqbzjxuo", .sequence_type = .protein },
    };

    for (cases) |case| {
        const fasta = try std.fmt.allocPrint(std.testing.allocator, ">seq\n{s}\n", .{case.sequence});
        defer std.testing.allocator.free(fasta);
        var summary = try validateBytes(std.testing.allocator, fasta, .{});
        defer summary.deinit(std.testing.allocator);

        try std.testing.expectEqual(@as(usize, 0), summary.events.items.len);
        try std.testing.expectEqual(case.sequence_type, summary.sequence_type);
        try std.testing.expectEqual(@as(u64, @intCast(case.sequence.len)), summary.type_bases_sampled);
    }
}

test "[edge] - [validate headers]: warns only above the configured byte limit" {
    const cases = [_]struct { fasta: []const u8, warning_count: usize }{
        .{ .fasta = ">abc\nA\n", .warning_count = 0 },
        .{ .fasta = ">abcd\nA\n", .warning_count = 1 },
    };

    for (cases) |case| {
        var summary = try validateBytes(std.testing.allocator, case.fasta, .{ .max_header_len = 3 });
        defer summary.deinit(std.testing.allocator);
        try std.testing.expectEqual(case.warning_count, countKind(&summary, .long_header));
        if (case.warning_count == 1) {
            const event = summary.events.items[0];
            try std.testing.expectEqual(@as(usize, 1), event.line);
            try std.testing.expectEqual(@as(usize, 3), event.limit);
        }
    }
}

test "[cli] - [validate JSON summary]: returns the complete stable schema" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_path = try writeFastaArtifact(allocator, "validate-type-json", ">seq\nACGT\n");
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};

    const expected = try std.fmt.allocPrint(
        allocator,
        "{{\"schema_version\":\"v1\",\"truncated\":false,\"sequence_type\":\"nucleotide\",\"type_bases_sampled\":4,\"type_sample_cap\":{d},\"counts\":{{\"no_sequences\":0,\"duplicate_name\":0,\"invalid_character\":0,\"null_byte\":0,\"utf8_bom\":0,\"inconsistent_line_widths\":0,\"trailing_whitespace\":0,\"empty_sequence\":0,\"missing_terminal_newline\":0,\"mixed_line_endings\":0,\"long_header\":0,\"schema_violation\":0}},\"first_examples\":{{}}}}\n",
        .{main.validator.TYPE_SAMPLE_BASES},
    );
    try expectCliResult(
        allocator,
        &.{ ZFASTA_BIN, "validate", "--json", "--summary", fasta_path },
        0,
        expected,
        "",
    );
}

test "[edge] - [validate scan]: reports missing terminal newline and invalid nucleotide byte" {
    var summary = try validateBytes(std.testing.allocator, ">seq\nACGTACGTAC?", .{});
    defer summary.deinit(std.testing.allocator);

    try std.testing.expectEqual(main.fasta.SequenceType.nucleotide, summary.sequence_type);
    try std.testing.expectEqual(@as(usize, 1), countKind(&summary, .missing_terminal_newline));
    try std.testing.expectEqual(@as(usize, 1), countKind(&summary, .invalid_character));
}

test "[unit] - [validate alphabet]: custom symbols replace the detected alphabet" {
    var summary = try validateBytes(std.testing.allocator, ">seq\nACGTZ\n", .{
        .custom_alphabet = "ACGTZ",
    });
    defer summary.deinit(std.testing.allocator);

    try std.testing.expectEqual(@as(usize, 0), countKind(&summary, .invalid_character));
}

test "[unit] - [validate layout]: reports width and trailing-whitespace warnings" {
    var summary = try validateBytes(std.testing.allocator, ">seq\nAAAA\nCC \nGGGG\nTT\n", .{});
    defer summary.deinit(std.testing.allocator);

    try std.testing.expectEqual(@as(usize, 1), countKind(&summary, .trailing_whitespace));
    try std.testing.expectEqual(@as(usize, 1), countKind(&summary, .inconsistent_line_widths));
}

test "[edge] - [validate layout]: reports a final line wider than the established width" {
    var summary = try validateBytes(std.testing.allocator, ">seq\nAA\nAAAA\n", .{});
    defer summary.deinit(std.testing.allocator);

    try std.testing.expectEqual(@as(usize, 1), countKind(&summary, .inconsistent_line_widths));
    const event = summary.events.items[0];
    try std.testing.expectEqual(@as(u32, 2), event.expected_width);
    try std.testing.expectEqual(@as(u32, 4), event.actual_width);
    try std.testing.expectEqual(@as(usize, 3), event.line);
}

test "[unit] - [validate schemas]: distinguishes valid and invalid UniProt and RefSeq headers" {
    const cases = [_]struct {
        fasta: []const u8,
        schema: validator.Schema,
        violations: usize,
    }{
        .{ .fasta = ">sp|P12345|PROT_HUMAN\nMAV\n", .schema = .uniprot, .violations = 0 },
        .{ .fasta = ">P12345\nMAV\n", .schema = .uniprot, .violations = 1 },
        .{ .fasta = ">NC_000001.11 Homo sapiens\nACGT\n", .schema = .refseq, .violations = 0 },
        .{ .fasta = ">bad\nACGT\n", .schema = .refseq, .violations = 1 },
    };

    for (cases) |case| {
        var summary = try validateBytes(std.testing.allocator, case.fasta, .{ .schema = case.schema });
        defer summary.deinit(std.testing.allocator);
        try std.testing.expectEqual(case.violations, countKind(&summary, .schema_violation));
    }
}

test "[unit] - [validate status]: promotes warnings in strict mode" {
    try std.testing.expectEqual(@as(u8, 0), validator.exitCode(0, 0, false));
    try std.testing.expectEqual(@as(u8, 2), validator.exitCode(0, 1, false));
    try std.testing.expectEqual(@as(u8, 1), validator.exitCode(0, 1, true));
    try std.testing.expectEqual(@as(u8, 1), validator.exitCode(1, 0, false));
}

test "[unit] - [validate status]: ignores warnings repaired by format rewrite" {
    var summary = try validateBytes(std.testing.allocator, ">seq\nACGT \n", .{});
    defer summary.deinit(std.testing.allocator);

    try std.testing.expectEqual(@as(usize, 1), countKind(&summary, .trailing_whitespace));
    try std.testing.expectEqual(@as(u8, 2), validator.exitCodeForOptions(&summary, .{}));
    try std.testing.expectEqual(@as(u8, 0), validator.exitCodeForOptions(&summary, .{ .fix = true }));
}

test "[property] - [validate fix]: is idempotent across supported layout rewrites" {
    const cases = [_][]const u8{
        "\xEF\xBB\xBF>seq\r\nAAAA\r\nCC \r\nGGGG\r\nTT",
        ">seq\r\nAAAA\nCC\n",
        ">seq\nACGT\n",
    };
    for (cases) |input| try expectFixIdempotent(std.testing.allocator, input);

    const clean = cases[2];
    var summary = try validateBytes(std.testing.allocator, clean, .{});
    defer summary.deinit(std.testing.allocator);
    const fixed = try fixBytes(std.testing.allocator, clean, summary.record_widths.items);
    defer std.testing.allocator.free(fixed);
    try std.testing.expectEqualStrings(clean, fixed);
}

test "[failure] - [validate fix]: rejects errors outside the format-only contract" {
    var dup = try validateBytes(std.testing.allocator, ">dup\nACGT\n>dup\nGCTA\n", .{});
    defer dup.deinit(std.testing.allocator);
    try std.testing.expectEqual(.duplicate_identifier, (try fixBlocker(&dup, .{})).?);

    var preamble_only = try validateBytes(std.testing.allocator, "not fasta\n", .{});
    defer preamble_only.deinit(std.testing.allocator);
    try std.testing.expectEqual(.content_before_header, (try fixBlocker(&preamble_only, .{})).?);

    var none = try validateBytes(std.testing.allocator, "", .{});
    defer none.deinit(std.testing.allocator);
    try std.testing.expectEqual(.no_records, (try fixBlocker(&none, .{})).?);

    var bad = try validateBytes(std.testing.allocator, ">seq\nACGT?\n", .{});
    defer bad.deinit(std.testing.allocator);
    try std.testing.expectEqual(.invalid_symbol, (try fixBlocker(&bad, .{})).?);
    try std.testing.expectEqual(null, try fixBlocker(&bad, .{ .fix_format_only = true }));

    var nulls = try validateBytes(std.testing.allocator, ">seq\nACGT\x00TT\n", .{});
    defer nulls.deinit(std.testing.allocator);
    try std.testing.expectEqual(.sequence_control_byte, (try fixBlocker(&nulls, .{})).?);
}

test "[property] - [validate format-only fix]: preserves invalid sequence characters" {
    const broken = ">seq\nACGT?\n";
    var summary = try validateBytes(std.testing.allocator, broken, .{});
    defer summary.deinit(std.testing.allocator);
    try std.testing.expectEqual(null, try fixBlocker(&summary, .{ .fix_format_only = true }));

    const fixed = try fixBytes(std.testing.allocator, broken, summary.record_widths.items);
    defer std.testing.allocator.free(fixed);
    try std.testing.expect(std.mem.indexOf(u8, fixed, "ACGT?") != null);

    var after = try validateBytes(std.testing.allocator, fixed, .{});
    defer after.deinit(std.testing.allocator);
    try std.testing.expectEqual(@as(usize, 1), countKind(&after, .invalid_character));
}

test "[property] - [validate fix]: preserves empty-record warnings" {
    const cases = [_][]const u8{ ">empty\n", ">empty" };
    for (cases) |broken| {
        var summary = try validateBytes(std.testing.allocator, broken, .{});
        defer summary.deinit(std.testing.allocator);
        try std.testing.expectEqual(@as(usize, 1), countKind(&summary, .empty_sequence));

        const fixed = try fixBytes(std.testing.allocator, broken, summary.record_widths.items);
        defer std.testing.allocator.free(fixed);
        try std.testing.expectEqualStrings(">empty\n", fixed);

        var after = try validateBytes(std.testing.allocator, fixed, .{});
        defer after.deinit(std.testing.allocator);
        try std.testing.expectEqual(@as(usize, 1), countKind(&after, .empty_sequence));
        try std.testing.expectEqual(@as(usize, 0), countKind(&after, .missing_terminal_newline));
    }
}

// --- Fix and cross-module contracts ---
test "[cli] - [validate fix]: replaces the output with the exact library rewrite" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const broken =
        \\>seq
        \\AAAA  
        \\CCCC
        \\
    ;

    var summary = try validateBytes(allocator, broken, .{});
    defer summary.deinit(allocator);
    const expected = try fixBytes(allocator, broken, summary.record_widths.items);

    const in_path = try writeFastaArtifact(allocator, "validate-fix-cli-in", broken);
    defer std.Io.Dir.cwd().deleteFile(io, in_path) catch {};
    const out_path = try writeFastaArtifact(allocator, "validate-fix-cli-out", "previous output\n");
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
    try std.testing.expectEqualStrings(
        "WARNING: line 2: trailing whitespace on sequence line in 'seq'\n",
        result.stdout,
    );
    try std.testing.expectEqual(@as(usize, 0), result.stderr.len);

    const got = try readTestFile(allocator, out_path);
    try std.testing.expectEqualStrings(expected, got);
}

test "[cli] - [validate format-only fix]: preserves invalid bytes while normal fix refuses" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const in_path = try writeFastaArtifact(allocator, "validate-format-only-in", ">seq\nACGT? \n");
    defer std.Io.Dir.cwd().deleteFile(io, in_path) catch {};
    const rejected_path = try utility.uniqueArtifactPath(allocator, "validate-fix-rejected", "fa");
    defer std.Io.Dir.cwd().deleteFile(io, rejected_path) catch {};
    const fixed_path = try utility.uniqueArtifactPath(allocator, "validate-format-only-out", "fa");
    defer std.Io.Dir.cwd().deleteFile(io, fixed_path) catch {};

    const expected_delivery = try std.fmt.allocPrint(allocator, "error [repair.blocked] validate: \"zig-cache/test-artifacts\" ({d} bytes):2 column=5 offset=9 record=1: Rewrite is blocked by the reported data condition. reason_code=sequence.invalid_symbol byte=0x3f policy=alphabet=protein selection=automatic custom_symbols=unknown sample=symbols=5 rule=nucleotide_above_90_percent interval=unknown request=unknown repair_policy=legacy_safe\nhint: Resolve the reported rewrite blocker before publishing output.\n", .{in_path.len});
    defer allocator.free(expected_delivery);
    try expectCliResult(
        allocator,
        &.{ ZFASTA_BIN, "validate", "--fix", "-o", rejected_path, in_path },
        1,
        "",
        expected_delivery,
    );
    try std.testing.expectError(
        error.FileNotFound,
        std.Io.Dir.cwd().openFile(io, rejected_path, .{}),
    );

    try expectCliResult(
        allocator,
        &.{ ZFASTA_BIN, "validate", "--fix", "--fix-format-only", "-o", fixed_path, in_path },
        1,
        "WARNING: line 2: trailing whitespace on sequence line in 'seq'\n" ++
            "ERROR: line 2: invalid sequence character 0x3f\n",
        "",
    );
    const fixed = try readTestFile(allocator, fixed_path);
    try std.testing.expectEqualStrings(">seq\nACGT?\n", fixed);
}

test "[cli] - [validate fix]: rejects an aliased input path without modifying it" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const broken = ">seq\nAAAA  \n";
    const in_path = try writeFastaArtifact(allocator, "validate-fix-input-alias", broken);
    defer std.Io.Dir.cwd().deleteFile(io, in_path) catch {};
    const output_alias = try std.fmt.allocPrint(allocator, "./{s}", .{in_path});

    var threaded = std.Io.Threaded.init(allocator, .{});
    defer threaded.deinit();
    const spawn_io = threaded.io();

    const result = try std.process.run(allocator, spawn_io, .{
        .argv = &.{ ZFASTA_BIN, "validate", "--error-format=json", "--fix", "-o", output_alias, in_path },
        .stdout_limit = .limited(64 * 1024),
        .stderr_limit = .limited(64 * 1024),
    });
    defer allocator.free(result.stdout);
    defer allocator.free(result.stderr);

    switch (result.term) {
        .exited => |code| try std.testing.expectEqual(@as(u8, 1), code),
        else => return error.ChildProcessFailed,
    }
    try std.testing.expectEqual(@as(usize, 0), result.stdout.len);
    try utility.expectDiagnostic(allocator, result.stderr, .@"repair.output_is_input");

    const input_after = try readTestFile(allocator, in_path);
    try std.testing.expectEqualStrings(broken, input_after);
}

test "[cli] - [validate fix]: writes a complete rewrite after report truncation" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var input: std.ArrayList(u8) = .empty;
    defer input.deinit(allocator);
    try input.appendSlice(allocator, ">seq\n");
    var i: usize = 0;
    while (i <= validator.MAX_VALIDATE_EVENTS) : (i += 1) {
        try input.appendSlice(allocator, "A \n");
    }

    var summary = try validateBytes(allocator, input.items, .{});
    defer summary.deinit(allocator);
    try std.testing.expect(summary.truncated);

    const in_path = try writeFastaArtifact(allocator, "validate-fix-truncated-in", input.items);
    defer std.Io.Dir.cwd().deleteFile(io, in_path) catch {};
    const out_path = try utility.uniqueArtifactPath(allocator, "validate-fix-truncated-out", "fa");
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
    var expected_stdout: std.ArrayList(u8) = .empty;
    defer expected_stdout.deinit(allocator);
    for (0..validator.MAX_VALIDATE_EVENTS) |event_index| {
        const line = try std.fmt.allocPrint(
            allocator,
            "WARNING: line {d}: trailing whitespace on sequence line in 'seq'\n",
            .{event_index + 2},
        );
        try expected_stdout.appendSlice(allocator, line);
    }
    try std.testing.expectEqualStrings(expected_stdout.items, result.stdout);
    const expected_stderr = try std.fmt.allocPrint(
        allocator,
        "warning: validate event list truncated at {d}; --fix output was still written\n",
        .{validator.MAX_VALIDATE_EVENTS},
    );
    try std.testing.expectEqualStrings(expected_stderr, result.stderr);

    const fixed = try readTestFile(allocator, out_path);
    var fixed_summary = try validateBytes(allocator, fixed, .{});
    defer fixed_summary.deinit(allocator);
    try std.testing.expectEqual(@as(usize, 0), fixed_summary.events.items.len);
}

test "[integration] - [validate fix]: preserves indexed retrieval across layout rewrites" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const cases = [_]struct {
        stem: []const u8,
        broken: []const u8,
        region: []const u8,
        expected: []const u8,
    }{
        .{
            .stem = "validate-fix-crafted",
            .broken = ">messy_seq widths and trailing ws\nAAAA    \nCCCC\nGGGGTT\n",
            .region = "messy_seq:1-14",
            .expected = ">messy_seq:1-14\nAAAACCCCGGGGTT\n",
        },
        .{
            .stem = "validate-fix-mixed",
            .broken = ">mixed_widths internal line widths vary\nAAAACCCCGGGG\nTTTTAA\nAACCCCGGGGTT\nTT\n",
            .region = "mixed_widths:3-24",
            .expected = ">mixed_widths:3-24\nAACCCCGGGGTTTTAAAACCCC\n",
        },
        .{
            .stem = "validate-fix-trailing",
            .broken = ">trailing_whitespace spaces and tabs after sequence bytes\nAAAACCCC    \nGGGGTTTT\t\nCCCCAAAA\n",
            .region = "trailing_whitespace:1-16",
            .expected = ">trailing_whitespace:1-16\nAAAACCCCGGGGTTTT\n",
        },
    };

    for (cases) |case| {
        try expectFormatFixPreservesRegion(
            allocator,
            case.stem,
            case.broken,
            case.region,
            case.expected,
        );
    }
}

test "[integration] - [validator and indexer]: agree on messy record catalog and geometry" {
    const data = try readTestFile(std.testing.allocator, "tests/data/validator_indexer_agreement.fasta");
    defer std.testing.allocator.free(data);

    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    try std.testing.expect(std.mem.indexOf(u8, data, "\r\n") != null);
    try std.testing.expect(std.mem.indexOf(u8, data, "\n\n") != null);
    try std.testing.expect(!std.mem.endsWith(u8, data, "\n"));
    try std.testing.expect(!std.mem.endsWith(u8, data, "\r\n"));
    try std.testing.expect(std.mem.indexOf(u8, data, ">empty_rec\n") != null);
    try std.testing.expect(std.mem.indexOf(u8, data, " \n") != null or std.mem.indexOf(u8, data, "\t\n") != null);

    var summary = try validateBytes(allocator, data, .{});
    defer summary.deinit(allocator);

    try std.testing.expect(countKind(&summary, .empty_sequence) >= 1);
    try std.testing.expect(countKind(&summary, .inconsistent_line_widths) >= 1);
    try std.testing.expect(countKind(&summary, .trailing_whitespace) >= 1);
    try std.testing.expect(countKind(&summary, .mixed_line_endings) >= 1);
    try std.testing.expect(countKind(&summary, .missing_terminal_newline) >= 1);
    try std.testing.expect(countKind(&summary, .long_header) >= 1);

    try std.testing.expectEqual(@as(usize, 6), summary.header_count);
    try std.testing.expectEqual(@as(usize, 5), summary.sequence_count);

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
    for (0..index.records.items.len) |i| {
        try std.testing.expect(!std.mem.eql(u8, zfiEmbeddedName(&index, i), "empty_rec"));
    }
    const long_name = zfiEmbeddedName(&index, 4);
    try std.testing.expect(long_name.len > 1024);
    try std.testing.expectEqual(@as(usize, 1025), long_name.len);
}

test "[cli] - [validate safety]: refuses preamble loss without changing legacy counts or destination" {
    const allocator = std.testing.allocator;
    const data = "junk\n>seq1\nACGT\n";
    var summary = try validateBytes(allocator, data, .{});
    defer summary.deinit(allocator);
    try std.testing.expectEqual(@as(usize, 0), summary.events.items.len);
    try std.testing.expectError(error.RepairRefused, fixBytes(allocator, data, summary.record_widths.items));
    const input = try writeFastaArtifact(allocator, "validate-preamble-input", data);
    defer allocator.free(input);
    defer std.Io.Dir.cwd().deleteFile(io, input) catch {};
    const output = try writeFastaArtifact(allocator, "validate-preamble-output", "previous destination\n");
    defer allocator.free(output);
    defer std.Io.Dir.cwd().deleteFile(io, output) catch {};
    for ([_]bool{ false, true }) |format_only| {
        var threaded = std.Io.Threaded.init(allocator, .{});
        defer threaded.deinit();
        const argv: []const []const u8 = if (format_only)
            &.{ ZFASTA_BIN, "validate", "--error-format=json", "--fix", "--fix-format-only", "-o", output, input }
        else
            &.{ ZFASTA_BIN, "validate", "--error-format=json", "--fix", "-o", output, input };
        const result = try std.process.run(allocator, threaded.io(), .{ .argv = argv });
        defer allocator.free(result.stdout);
        defer allocator.free(result.stderr);
        try std.testing.expectEqual(@as(u8, 1), result.term.exited);
        try std.testing.expectEqualStrings("", result.stdout);
        try utility.expectDiagnostic(allocator, result.stderr, .@"repair.blocked");
        const after = try readTestFile(allocator, output);
        defer allocator.free(after);
        try std.testing.expectEqualStrings("previous destination\n", after);
    }
}

test "[cli] - [validate safety]: rejects a symlink alias and safely replaces a separate hard link" {
    const allocator = std.testing.allocator;
    const broken = ">seq1\nACGT";
    const input = try writeFastaArtifact(allocator, "validate-link-input", broken);
    defer allocator.free(input);
    defer std.Io.Dir.cwd().deleteFile(io, input) catch {};
    const hard = try utility.uniqueArtifactPath(allocator, "validate-hard-output", "fa");
    defer allocator.free(hard);
    defer std.Io.Dir.cwd().deleteFile(io, hard) catch {};
    const sym = try utility.uniqueArtifactPath(allocator, "validate-sym-output", "fa");
    defer allocator.free(sym);
    defer std.Io.Dir.cwd().deleteFile(io, sym) catch {};
    const absolute_input = try std.Io.Dir.cwd().realPathFileAlloc(io, input, allocator);
    defer allocator.free(absolute_input);
    try std.Io.Dir.cwd().symLink(io, absolute_input, sym, .{});
    try std.Io.Dir.cwd().hardLink(input, std.Io.Dir.cwd(), hard, io, .{});
    var threaded = std.Io.Threaded.init(allocator, .{});
    defer threaded.deinit();
    const rejected = try std.process.run(allocator, threaded.io(), .{ .argv = &.{ ZFASTA_BIN, "validate", "--fix", "-o", sym, input } });
    defer allocator.free(rejected.stdout);
    defer allocator.free(rejected.stderr);
    try std.testing.expectEqual(@as(u8, 1), rejected.term.exited);
    const accepted = try std.process.run(allocator, threaded.io(), .{ .argv = &.{ ZFASTA_BIN, "validate", "--fix", "-o", hard, input } });
    defer allocator.free(accepted.stdout);
    defer allocator.free(accepted.stderr);
    try std.testing.expectEqual(@as(u8, 0), accepted.term.exited);
    const input_after = try readTestFile(allocator, input);
    defer allocator.free(input_after);
    const output_after = try readTestFile(allocator, hard);
    defer allocator.free(output_after);
    try std.testing.expectEqualStrings(broken, input_after);
    try std.testing.expectEqualStrings(">seq1\nACGT\n", output_after);
}

test "[failure] - [validate safety]: refuses an output open failure before staging" {
    const diagnostics = main.index_format.diagnostics;
    const allocator = std.testing.allocator;
    const broken = ">seq1\nACGT";
    const input = try writeFastaArtifact(allocator, "validate-identity-input", broken);
    defer allocator.free(input);
    defer std.Io.Dir.cwd().deleteFile(io, input) catch {};
    const output = try utility.uniqueArtifactPath(allocator, "validate-identity-loop", "fa");
    defer allocator.free(output);
    defer std.Io.Dir.cwd().deleteFile(io, output) catch {};
    const directory = try std.Io.Dir.cwd().realPathFileAlloc(io, ".", allocator);
    defer allocator.free(directory);
    const absolute_output = try std.fs.path.join(allocator, &.{ directory, output });
    defer allocator.free(absolute_output);
    try std.Io.Dir.cwd().symLink(io, absolute_output, output, .{});

    var failure: diagnostics.Failure = undefined;
    try std.testing.expectError(error.RepairRefused, validator.runValidate(allocator, io, input, .{ .fix = true, .output_path = output }, &failure));
    const primary = failure.primary().?;
    try std.testing.expectEqual(diagnostics.Code.@"repair.output_identity_unknown", primary.code());
    try std.testing.expectEqual(diagnostics.Stage.check_identity, primary.stage);
    try std.testing.expectEqual(diagnostics.SourceId.output, primary.at.source.?);
    try std.testing.expectEqual(error.SymLinkLoop, primary.cause.?.zig);
    try std.testing.expectEqual(.not_started, failure.outputs.file.state);
    try std.testing.expectEqual(.not_started, failure.outputs.stdout.state);
    try std.testing.expect(failure.work_complete);
    try std.testing.expect(!failure.complete);
    const after = try readTestFile(allocator, input);
    defer allocator.free(after);
    try std.testing.expectEqualStrings(broken, after);
    try std.testing.expectError(error.SymLinkLoop, std.Io.Dir.cwd().openFile(io, output, .{}));
}

test "[failure] - [native validator blockers]: owns duplicate facts beyond retained events" {
    const diagnostics = main.index_format.diagnostics;
    const allocator = std.testing.allocator;
    var data: std.ArrayList(u8) = .empty;
    defer data.deinit(allocator);
    try data.appendSlice(allocator, ">dup\n");
    for (0..validator.MAX_VALIDATE_EVENTS) |_| try data.appendSlice(allocator, "A \n");
    const duplicate_offset = data.items.len;
    try data.appendSlice(allocator, ">dup\nA\n");
    var failure: diagnostics.Failure = undefined;
    var summary = try validator.validateData(allocator, data.items, .{}, &failure);
    try std.testing.expect(summary.truncated);
    try std.testing.expectEqual(@as(usize, 1), summary.error_count);
    try std.testing.expectError(error.RepairRefused, validator.checkFixAllowed(&summary, .{ .fix_format_only = true }, &failure));
    const primary = failure.primary().?;
    try std.testing.expectEqual(@as(?u64, duplicate_offset), primary.at.location.offset);
    try std.testing.expectEqual(@as(?u64, validator.MAX_VALIDATE_EVENTS + 2), primary.at.location.line);
    try std.testing.expectEqual(@as(?u64, 1), primary.at.location.column);
    try std.testing.expectEqual(@as(?u64, 2), primary.at.record);
    try std.testing.expectEqual(@as(?u64, 1), primary.related.?.location.line);
    const preview = primary.details.@"repair.blocked".blocker.duplicate_identifier.identifier.?;
    summary.deinit(allocator);
    @memset(data.items, '?');
    try std.testing.expectEqualStrings("dup", failure.preview(preview));
}

test "[failure] - [native validator blockers]: distinguishes header and sequence bytes and exact alphabet evidence" {
    const diagnostics = main.index_format.diagnostics;
    const cases = [_][]const u8{ ">na\x00me\nACGT\n", ">seq\nA\x00C \n" };
    for (cases, 0..) |data, i| {
        var failure: diagnostics.Failure = undefined;
        var summary = try validator.validateData(std.testing.allocator, data, .{}, &failure);
        defer summary.deinit(std.testing.allocator);
        try std.testing.expectError(error.RepairRefused, validator.checkFixAllowed(&summary, .{}, &failure));
        const primary = failure.primary().?;
        try std.testing.expectEqual(@as(?u64, if (i == 0) 3 else 6), primary.at.location.offset);
        try std.testing.expectEqual(@as(?u64, if (i == 0) 4 else 2), primary.at.location.column);
        try std.testing.expectEqual(@as(?u64, 1), primary.at.record);
        const blocker = primary.details.@"repair.blocked".blocker;
        try std.testing.expectEqual(if (i == 0) diagnostics.Code.@"fasta.header_control_byte" else diagnostics.Code.@"fasta.sequence_control_byte", blocker.code());
    }
    var failure: diagnostics.Failure = undefined;
    const options = validator.Options{ .custom_alphabet = "ACGT" };
    var summary = try validator.validateData(std.testing.allocator, ">seq\n?A#?\n", options, &failure);
    defer summary.deinit(std.testing.allocator);
    try std.testing.expectEqual(@as(usize, 2), summary.kind_counts[@intFromEnum(validator.Kind.invalid_character)]);
    try std.testing.expectError(error.RepairRefused, validator.checkFixAllowed(&summary, options, &failure));
    const primary = failure.primary().?;
    try std.testing.expectEqual(@as(?u64, 7), primary.at.location.offset);
    try std.testing.expectEqual(@as(?u64, 3), primary.at.location.column);
    const detail = primary.details.@"repair.blocked".blocker.invalid_symbol;
    try std.testing.expectEqual(@as(u8, '#'), detail.byte);
    try std.testing.expectEqual(diagnostics.Alphabet.custom, detail.policy.alphabet);
    try std.testing.expect(detail.policy.sample == null);
    try std.testing.expect(detail.policy.custom_symbols.?['A' / 8] & (@as(u8, 1) << ('A' % 8)) != 0);
    failure.reset(.validate);
    try validator.checkFixAllowed(&summary, .{ .fix_format_only = true, .custom_alphabet = "ACGT" }, &failure);
    try std.testing.expect(failure.primary() == null);
}

test "[failure] - [native validator allocations]: releases scanner, rewrite and event storage" {
    try std.testing.checkAllAllocationFailures(std.testing.allocator, checkNativeValidatorAllocations, .{">seq\nACG \nT\n"});
}

test "[edge] - [native validator BOM]: retains original first-line byte columns" {
    const cases = [_]struct { data: []const u8, offset: u64, column: u64 }{
        .{ .data = "\xef\xbb\xbf>na\x00me\nACGT\n", .offset = 6, .column = 7 },
        .{ .data = "\xef\xbb\xbf \tjunk\n>seq\nACGT\n", .offset = 5, .column = 6 },
    };
    for (cases) |case| {
        var failure: main.index_format.diagnostics.Failure = undefined;
        var summary = try validator.validateData(std.testing.allocator, case.data, .{}, &failure);
        defer summary.deinit(std.testing.allocator);
        try std.testing.expectError(error.RepairRefused, validator.checkFixAllowed(&summary, .{}, &failure));
        try std.testing.expectEqual(@as(?u64, 1), failure.primary().?.at.location.line);
        try std.testing.expectEqual(@as(?u64, case.offset), failure.primary().?.at.location.offset);
        try std.testing.expectEqual(@as(?u64, case.column), failure.primary().?.at.location.column);
        if (summary.preamble != null) {
            try std.testing.expectError(error.RepairRefused, validator.fixData(std.testing.allocator, case.data, summary.record_widths.items, &failure));
            try std.testing.expectEqual(@as(?u64, case.offset), failure.primary().?.at.location.offset);
            try std.testing.expectEqual(@as(?u64, case.column), failure.primary().?.at.location.column);
        }
    }
}

test "[cli] - [validate text]: escapes unsafe names and keeps v1 JSON escaping" {
    const allocator = std.testing.allocator;
    const input = try writeFastaArtifact(allocator, "validate-text-name", ">na\x1b\x7f\xffme\nA\n>na\x1b\x7f\xffme\nC\n");
    defer allocator.free(input);
    defer std.Io.Dir.cwd().deleteFile(io, input) catch {};
    try expectCliResult(allocator, &.{ ZFASTA_BIN, "validate", input }, 1, "ERROR: line 3: duplicate name 'na\\x1b\\x7f\\xffme' (first seen line 1)\n", "");
    const json = try renderEvent(allocator, .{ .level = .error_level, .kind = .duplicate_name, .line = 3, .name = "\xff", .first_line = 1 });
    defer allocator.free(json);
    try std.testing.expect(std.mem.indexOf(u8, json, "\\u00ff") != null);
}

test "[failure] - [native validator input]: retains stat and mapping causes before scanning" {
    const diagnostics = main.index_format.diagnostics;
    const Fault = struct {
        fn stat(_: ?*anyopaque, _: std.Io.File) std.Io.File.StatError!std.Io.File.Stat {
            return error.AccessDenied;
        }
        fn map(_: ?*anyopaque, _: std.Io.File, _: std.Io.File.MemoryMap.CreateOptions) std.Io.File.MemoryMap.CreateError!std.Io.File.MemoryMap {
            return error.AccessDenied;
        }
    };
    var vtable = io.vtable.*;
    const fault_io = std.Io{ .userdata = io.userdata, .vtable = &vtable };
    var failure: diagnostics.Failure = undefined;
    vtable.fileStat = Fault.stat;
    try std.testing.expectError(error.Io, validator.runValidate(std.testing.allocator, fault_io, "tests/data/simple.fasta", .{}, &failure));
    try std.testing.expectEqual(diagnostics.Code.@"io.stat_failed", failure.primary().?.code());
    try std.testing.expectEqual(error.AccessDenied, failure.primary().?.cause.?.zig);
    vtable.fileStat = io.vtable.fileStat;
    vtable.fileMemoryMapCreate = Fault.map;
    try std.testing.expectError(error.Io, validator.runValidate(std.testing.allocator, fault_io, "tests/data/simple.fasta", .{}, &failure));
    try std.testing.expectEqual(diagnostics.Code.@"io.map_failed", failure.primary().?.code());
    try std.testing.expectEqual(error.AccessDenied, failure.primary().?.cause.?.zig);
    try std.testing.expect(!failure.work_complete and !failure.complete);
    try std.testing.expectEqual(.not_started, failure.outputs.stdout.state);
}

test "[failure] - [native validator repair]: preserves published output after report flush failure" {
    const diagnostics = main.index_format.diagnostics;
    const Fault = struct {
        fn operate(userdata: ?*anyopaque, operation: std.Io.Operation) std.Io.Cancelable!std.Io.Operation.Result {
            return switch (operation) {
                .file_write_streaming => |write| if (write.file.handle == std.Io.File.stdout().handle)
                    .{ .file_write_streaming = error.BrokenPipe }
                else
                    io.vtable.operate(userdata, operation),
                else => io.vtable.operate(userdata, operation),
            };
        }
    };
    const allocator = std.testing.allocator;
    const input = try writeFastaArtifact(allocator, "validate-published-input", ">seq\nACGT");
    defer allocator.free(input);
    defer std.Io.Dir.cwd().deleteFile(io, input) catch {};
    const output = try writeFastaArtifact(allocator, "validate-published-output", "old destination\n");
    defer allocator.free(output);
    defer std.Io.Dir.cwd().deleteFile(io, output) catch {};
    var vtable = io.vtable.*;
    vtable.operate = Fault.operate;
    const fault_io = std.Io{ .userdata = io.userdata, .vtable = &vtable };
    var failure: diagnostics.Failure = undefined;
    try std.testing.expectError(error.Io, validator.runValidate(allocator, fault_io, input, .{ .fix = true, .output_path = output }, &failure));
    try std.testing.expectEqual(diagnostics.Code.@"io.flush_failed", failure.primary().?.code());
    try std.testing.expectEqual(diagnostics.Stage.flush, failure.primary().?.stage);
    try std.testing.expectEqual(diagnostics.SourceId.stdout, failure.primary().?.at.source);
    try std.testing.expectEqual(error.BrokenPipe, failure.primary().?.cause.?.zig);
    try std.testing.expect(failure.work_complete and !failure.complete);
    try std.testing.expectEqual(.published, failure.outputs.file.state);
    try std.testing.expectEqual(.succeeded, failure.outputs.file.final_flush);
    try std.testing.expectEqual(.failed, failure.outputs.stdout.final_flush);
    const written = try readTestFile(allocator, output);
    defer allocator.free(written);
    try std.testing.expectEqualStrings(">seq\nACGT\n", written);
    const original = try readTestFile(allocator, input);
    defer allocator.free(original);
    try std.testing.expectEqualStrings(">seq\nACGT", original);
}

test "[failure] - [native validator stdout]: retains final flush cause and resets for reuse" {
    const diagnostics = main.index_format.diagnostics;
    const Fault = struct {
        fn operate(userdata: ?*anyopaque, operation: std.Io.Operation) std.Io.Cancelable!std.Io.Operation.Result {
            return switch (operation) {
                .file_write_streaming => .{ .file_write_streaming = error.BrokenPipe },
                else => io.vtable.operate(userdata, operation),
            };
        }
        fn discard(userdata: ?*anyopaque, operation: std.Io.Operation) std.Io.Cancelable!std.Io.Operation.Result {
            return switch (operation) {
                .file_write_streaming => |write| blk: {
                    var size = write.header.len;
                    for (write.data, 0..) |data, i| size += data.len * (if (i + 1 == write.data.len) write.splat else 1);
                    break :blk .{ .file_write_streaming = size };
                },
                else => io.vtable.operate(userdata, operation),
            };
        }
    };
    var vtable = io.vtable.*;
    vtable.operate = Fault.operate;
    const fault_io = std.Io{ .userdata = io.userdata, .vtable = &vtable };
    var failure: diagnostics.Failure = undefined;
    try std.testing.expectError(error.Io, validator.runValidate(std.testing.allocator, fault_io, "tests/data/simple.fasta", .{ .output_mode = .json_summary }, &failure));
    try std.testing.expectEqual(diagnostics.Code.@"io.flush_failed", failure.primary().?.code());
    try std.testing.expectEqual(error.BrokenPipe, failure.primary().?.cause.?.zig);
    try std.testing.expect(failure.work_complete and !failure.complete);
    try std.testing.expectEqual(.failed, failure.outputs.stdout.final_flush);
    try std.testing.expect(failure.outputs.stdout.may_be_partial);
    try std.testing.expectError(error.Io, validator.runValidate(std.testing.allocator, fault_io, "tests/data/validate-file-does-not-exist.fa", .{}, &failure));
    try std.testing.expectEqual(diagnostics.Code.@"io.file_not_found", failure.primary().?.code());
    try std.testing.expect(!failure.work_complete and !failure.complete);
    try std.testing.expectEqual(.not_started, failure.outputs.stdout.state);
    vtable.operate = Fault.discard;
    const result = try validator.runValidate(std.testing.allocator, fault_io, "tests/data/simple.fasta", .{ .output_mode = .json_summary }, &failure);
    try std.testing.expectEqual(@as(u8, 0), result.status);
    try std.testing.expect(failure.work_complete and failure.complete);
    try std.testing.expect(failure.primary() == null);
    try std.testing.expectEqual(.succeeded, failure.outputs.stdout.final_flush);
}

test "[cli] - [validate cap]: keeps complete counts and explains retained examples" {
    const allocator = std.testing.allocator;
    var data: std.ArrayList(u8) = .empty;
    defer data.deinit(allocator);
    try data.appendSlice(allocator, ">seq\n");
    for (0..validator.MAX_VALIDATE_EVENTS + 1) |_| try data.appendSlice(allocator, "A \n");
    const input = try writeFastaArtifact(allocator, "validate-cap-counts", data.items);
    defer allocator.free(input);
    defer std.Io.Dir.cwd().deleteFile(io, input) catch {};
    const result = try std.process.run(allocator, io, .{ .argv = &.{ ZFASTA_BIN, "validate", "--json", "--summary", input } });
    defer allocator.free(result.stdout);
    defer allocator.free(result.stderr);
    switch (result.term) {
        .exited => |code| try std.testing.expectEqual(@as(u8, 1), code),
        else => return error.ChildProcessFailed,
    }
    try std.testing.expectEqualStrings("error: validate retained 10000 events; scanning and counting completed; fix reported issues and re-run\n", result.stderr);
    const parsed = try std.json.parseFromSlice(std.json.Value, allocator, result.stdout, .{});
    defer parsed.deinit();
    try std.testing.expect(parsed.value.object.get("truncated").?.bool);
    try std.testing.expectEqual(@as(i64, validator.MAX_VALIDATE_EVENTS + 1), parsed.value.object.get("counts").?.object.get("trailing_whitespace").?.integer);
    const machine = try std.process.run(allocator, io, .{ .argv = &.{ ZFASTA_BIN, "validate", "--json", "--summary", "--error-format=json", input } });
    defer allocator.free(machine.stdout);
    defer allocator.free(machine.stderr);
    try std.testing.expectEqual(@as(u8, 1), machine.term.exited);
    try std.testing.expectEqualStrings(result.stdout, machine.stdout);
    try std.testing.expectEqualStrings("", machine.stderr);

    const Fault = struct {
        fn operate(userdata: ?*anyopaque, operation: std.Io.Operation) std.Io.Cancelable!std.Io.Operation.Result {
            return switch (operation) {
                .file_write_streaming => .{ .file_write_streaming = error.BrokenPipe },
                else => io.vtable.operate(userdata, operation),
            };
        }
    };
    var vtable = io.vtable.*;
    vtable.operate = Fault.operate;
    var failure: main.index_format.diagnostics.Failure = undefined;
    try std.testing.expectError(error.Io, validator.runValidate(allocator, .{ .userdata = io.userdata, .vtable = &vtable }, input, .{}, &failure));
    try std.testing.expectEqual(main.index_format.diagnostics.Code.@"io.write_failed", failure.primary().?.code());
    try std.testing.expectEqual(error.BrokenPipe, failure.primary().?.cause.?.zig);
    try std.testing.expect(failure.work_complete and !failure.complete);
    try std.testing.expectEqual(.not_attempted, failure.outputs.stdout.final_flush);
}

test "[edge] - [validate fix widths]: equal frequencies select the first seen width" {
    const cases = [_]struct { input: []const u8, width: u32, fixed: []const u8 }{
        .{ .input = ">seq\nAAA\nCC\nGGGG\n", .width = 3, .fixed = ">seq\nAAA\nCCG\nGGG\n" },
        .{ .input = ">seq\nCC\nAAA\nGGGG\n", .width = 2, .fixed = ">seq\nCC\nAA\nAG\nGG\nG\n" },
    };
    for (cases) |case| {
        var summary = try validateBytes(std.testing.allocator, case.input, .{});
        defer summary.deinit(std.testing.allocator);
        try std.testing.expectEqual(case.width, summary.record_widths.items[0]);
        const fixed = try fixBytes(std.testing.allocator, case.input, summary.record_widths.items);
        defer std.testing.allocator.free(fixed);
        try std.testing.expectEqualStrings(case.fixed, fixed);
    }
}
