//! Diagnostic facts, storage bounds, source lifetimes, and error-union cleanup.

const std = @import("std");
const diagnostics = @import("diagnostics");
const caller = @import("diagnostics_caller.zig");

const testing = std.testing;

fn retainTemporaryLabel(failure: *diagnostics.Failure) void {
    var local_path = [_]u8{ 't', 'e', 'm', 'p', '/', 'a', 0xff };
    const retained = failure.retain(&local_path, 32);
    failure.sources.set(.staged_output, .{ .label = .{ .retained = retained }, .kind = .path });
    @memset(&local_path, 'x');
}

fn writeFailure(stage: diagnostics.Stage, cause: anyerror) diagnostics.Diagnostic {
    return .{
        .details = .{ .@"io.write_failed" = {} },
        .stage = stage,
        .at = .{ .source = .stdout },
        .cause = .{ .zig = cause },
    };
}

fn allocationCheck(allocator: std.mem.Allocator) !void {
    var failure: diagnostics.Failure = undefined;
    failure.reset(.stats);
    const first = allocator.alloc(u8, 8) catch |err| {
        try testing.expect(failure.fail(.{
            .details = .{ .@"resource.out_of_memory" = .{ .allocation_role = .names, .requested_bytes = 8 } },
            .stage = .allocate,
            .cause = .{ .zig = err },
        }) == error.OutOfMemory);
        try testing.expect(failure.primary().?.cause.?.zig == error.OutOfMemory);
        return err;
    };
    defer allocator.free(first);
    const second = allocator.alloc(u8, 16) catch |err| {
        try testing.expect(failure.fail(.{
            .details = .{ .@"resource.out_of_memory" = .{ .allocation_role = .records, .requested_bytes = 16 } },
            .stage = .allocate,
            .cause = .{ .zig = err },
        }) == error.OutOfMemory);
        try testing.expect(failure.primary().?.details.@"resource.out_of_memory".requested_bytes.? == 16);
        const count = std.math.maxInt(u64);
        try testing.expectError(error.OutOfMemory, diagnostics.oneBased(&failure, count, .{}));
        try testing.expect(failure.secondaryItems()[0].code() == .@"resource.count_overflow");
        return err;
    };
    defer allocator.free(second);
    failure.succeed();
    try testing.expect(failure.primary() == null and failure.complete);
}

test "[unit] - [diagnostics]: copied source labels survive cleanup and failure moves" {
    var failure: diagnostics.Failure = undefined;
    failure.reset(.index);
    retainTemporaryLabel(&failure);
    try testing.expectEqualSlices(u8, "temp/a\xff", failure.sourceLabel(.staged_output).?);
    try testing.expect(failure.sourceLabel(.name_spool) == null);
    var moved = failure;
    failure.reset(.get);
    @memset(&failure.previews, 'z');
    try testing.expectEqualSlices(u8, "temp/a\xff", moved.sourceLabel(.staged_output).?);
    const label = moved.sources.get(.staged_output).?.label.retained;
    try testing.expect(label.length_kind == .exact and label.original_bytes == 7);
}
test "[edge] - [diagnostics]: streamed prefixes preserve original lengths and unknown suffixes" {
    var failure: diagnostics.Failure = undefined;
    failure.reset(.index);
    const maximum = std.math.maxInt(u64);
    const exact = failure.retainPrefix("000x", .{ .exact = maximum }, 4);
    const lower_bound = failure.retainPrefix("12", .{ .at_least = 2 }, 2);
    const disabled = failure.retainPrefix("", .{ .exact = maximum }, 0);
    try testing.expectEqual(maximum, exact.original_bytes);
    try testing.expect(exact.length_kind == .exact and exact.truncated);
    try testing.expectEqualStrings("000x", failure.preview(exact));
    try testing.expect(lower_bound.length_kind == .at_least and lower_bound.truncated);
    try testing.expect(disabled.length_kind == .exact and disabled.truncated);
    try testing.expect(disabled.len == 0 and disabled.original_bytes == maximum);
}

test "[unit] - [diagnostics]: observed source changes and cancellation do not invent causes" {
    var failure: diagnostics.Failure = undefined;
    failure.reset(.index);
    try testing.expect(failure.fail(.{
        .details = .{ .@"input.changed" = .{
            .changed_fields = .{ .path = true },
            .before = .{ .bytes = 19, .mtime_ns = -1 },
            .after = null,
            .disappeared = true,
        } },
        .stage = .compare_source,
    }) == error.InputChanged);
    const change = failure.primary().?.details.@"input.changed";
    try testing.expect(change.after == null and change.disappeared);
    try testing.expect(change.changed_fields.path and !change.changed_fields.size);
    try testing.expectEqual(@as(i128, -1), change.before.mtime_ns.?);
    try testing.expect(failure.primary().?.cause == null);

    failure.reset(.get);
    try testing.expect(failure.fail(.{
        .details = .{ .@"io.read_failed" = .{} },
        .stage = .read,
        .cause = .{ .zig = error.Canceled },
    }) == error.Io);
    try testing.expect(failure.primary().?.code() == .@"io.read_failed");
    try testing.expect(failure.primary().?.cause.?.zig == error.Canceled);
}
test "[failure] - [diagnostics]: deferred failure survives success and yields to a nested blocker" {
    var failure: diagnostics.Failure = undefined;
    failure.reset(.index);
    failure.deferFailure(writeFailure(.write, error.NoSpaceLeft));
    try testing.expectEqual(@as(u64, 1), try diagnostics.oneBased(&failure, 0, .{}));
    try testing.expect(failure.primary_deferred);
    try testing.expectError(error.Io, failure.finishDeferred());
    try testing.expect(!failure.primary_deferred);
    try testing.expect(failure.primary().?.cause.?.zig == error.NoSpaceLeft);

    failure.reset(.index);
    failure.deferFailure(writeFailure(.write, error.NoSpaceLeft));
    failure.deferFailure(writeFailure(.flush, error.WriteFailed));
    try testing.expectError(error.CountOverflow, diagnostics.oneBased(&failure, std.math.maxInt(u64), .{}));
    try testing.expect(!failure.primary_deferred);
    try testing.expect(failure.primary().?.code() == .@"resource.count_overflow");
    try testing.expect(failure.secondaryItems()[0].cause.?.zig == error.NoSpaceLeft);
    try testing.expect(failure.secondaryItems()[1].stage == .flush);
}

test "[unit] - [diagnostics]: findings and unfinished work retain separate meanings" {
    const finding = diagnostics.Diagnostic{
        .details = .{ .@"index.source_size_mismatch" = .{ .stored_bytes = 35, .actual_bytes = 36 } },
        .stage = .compare_source,
        .at = .{ .source = .zfi, .location = .{ .offset = 4 } },
        .related = .{ .source = .fasta },
    };
    const operation_failure = writeFailure(.write, error.WriteFailed);
    const data = diagnostics.metadata(finding.code());
    const io = diagnostics.metadata(operation_failure.code());
    try testing.expect(data.kind == .finding and data.unit == .file);
    try testing.expect(io.kind == .failure and io.unit == .none);
    try testing.expect(data.severity == .@"error" and io.severity == .@"error");
    try testing.expect(data.returned_error == error.InvalidIndex);
    try testing.expect(io.returned_error == error.Io);
    try testing.expectEqualStrings("index.source_size_mismatch", data.name);
    try testing.expect(finding.related.?.location.offset == null);
}

test "[edge] - [diagnostics]: zero and maximum offsets differ from unknown" {
    var failure: diagnostics.Failure = undefined;
    failure.reset(.get);
    const maximum = std.math.maxInt(u64);
    try testing.expectEqual(@as(u64, 1), try diagnostics.oneBased(&failure, 0, .{}));
    try testing.expectEqual(maximum, try diagnostics.oneBased(&failure, maximum - 1, .{}));
    try testing.expect(failure.primary() == null);

    const diagnostic = diagnostics.Diagnostic{
        .details = .{ .@"io.unexpected_eof" = .{ .exact = .{ .requested_bytes = 1, .received_bytes = 0 } } },
        .stage = .read,
        .at = .{ .source = .fasta, .location = .{ .offset = 0 } },
        .related = .{ .source = .zfi, .location = .{ .offset = maximum }, .record = maximum },
    };
    try testing.expect(failure.fail(diagnostic) == error.Io);
    try testing.expectEqual(@as(u64, 0), failure.primary().?.at.location.offset.?);
    try testing.expect(failure.primary().?.at.location.line == null);
    try testing.expect(failure.primary().?.at.location.column == null);
    try testing.expect(failure.primary().?.at.record == null);
    try testing.expect(failure.primary().?.at.request == null);
    try testing.expectEqual(maximum, failure.primary().?.related.?.location.offset.?);
    try testing.expectEqual(maximum, failure.primary().?.related.?.record.?);
}

test "[failure] - [diagnostics]: ordinal conversion cannot wrap to zero" {
    var failure: diagnostics.Failure = undefined;
    failure.reset(.get);
    const maximum = std.math.maxInt(u64);
    try testing.expectError(error.CountOverflow, diagnostics.oneBased(&failure, maximum, .{ .source = .zfi }));
    const facts = failure.primary().?.details.@"resource.count_overflow";
    try testing.expect(facts.counter == .report_ordinal);
    try testing.expectEqual(@as(u128, maximum), facts.max_value);
    try testing.expectEqual(@as(u128, maximum), facts.current.?);
    try testing.expectEqual(@as(u128, 1), facts.increment.?);
    try testing.expect(failure.primary().?.at.record == null);
}

test "[unit] - [diagnostics]: caller sources and static suffixes survive helper return" {
    const caller_path = [_]u8{ 'a', 0xff, '\n', '.', 'f', 'a' };
    var failure: diagnostics.Failure = undefined;
    failure.reset(.index);
    failure.sources.set(.fasta, .{ .label = .{ .borrowed = &caller_path }, .kind = .path });
    failure.sources.set(.zfi, .{ .label = .{ .borrowed = &caller_path }, .kind = .path, .suffix = .zfi });
    try testing.expect(failure.sources.get(.bed) == null);
    try testing.expectEqualSlices(u8, &caller_path, failure.sources.get(.zfi).?.label.borrowed);
    try testing.expect(failure.sources.get(.zfi).?.label.borrowed.ptr == caller_path[0..].ptr);
    try testing.expect(failure.sources.get(.zfi).?.suffix == .zfi);
    failure.reset(.get);
    try testing.expect(failure.sources.get(.zfi) == null);
    failure.sources.set(.fasta, .{ .label = .{ .borrowed = "next.fa" }, .kind = .path });
    try testing.expectEqualStrings("next.fa", failure.sources.get(.fasta).?.label.borrowed);
}

test "[edge] - [diagnostics]: previews own bytes within one shared budget" {
    var failure: diagnostics.Failure = undefined;
    failure.reset(.get);
    var bytes: [300]u8 = undefined;
    @memset(&bytes, 0xff);
    bytes[0] = 0;
    bytes[1] = '\n';
    const first = failure.retain(&bytes, 200);
    var more_bytes: [100]u8 = undefined;
    @memset(&more_bytes, 'a');
    const second = failure.retain(&more_bytes, 100);
    const exhausted = failure.retain("known", 5);
    @memset(&bytes, 'x');
    try testing.expectEqual(@as(usize, 200), failure.preview(first).len);
    try testing.expectEqual(@as(u8, 0), failure.preview(first)[0]);
    try testing.expectEqual(@as(u8, '\n'), failure.preview(first)[1]);
    try testing.expectEqual(@as(u8, 0xff), failure.preview(first)[2]);
    try testing.expectEqual(@as(u64, 300), first.original_bytes);
    try testing.expect(first.truncated and second.truncated and exhausted.truncated);
    try testing.expectEqual(@as(usize, 56), failure.preview(second).len);
    try testing.expectEqual(@as(usize, 0), failure.preview(exhausted).len);
    try testing.expectEqual(@as(u64, 5), exhausted.original_bytes);
    try testing.expectEqual(@as(u16, 256), failure.preview_len);

    var copied = failure;
    @memset(&failure.previews, 'z');
    try testing.expectEqual(@as(u8, 0xff), copied.preview(first)[2]);
    const empty = copied.retain("", 0);
    try testing.expect(!empty.truncated and empty.original_bytes == 0);
}

test "[unit] - [diagnostics]: operation entry clears live state without wiping arrays" {
    var failure: diagnostics.Failure = undefined;
    failure.reset(.index);
    @memset(&failure.previews, 0xa5);
    _ = failure.retain("old", 3);
    try testing.expect(failure.fail(writeFailure(.write, error.WriteFailed)) == error.Io);
    failure.addSecondary(writeFailure(.flush, error.FlushFailed));
    failure.secondary_truncated = true;
    failure.work_complete = true;
    failure.outputs.stdout = .{ .state = .streaming, .logical_bytes = 10, .may_be_partial = true };
    failure.outputs.file.state = .published;
    failure.sources.set(.fasta, .{ .label = .{ .borrowed = "old.fa" }, .kind = .path });
    failure.reset(.stats);
    try testing.expect(failure.operation == .stats and failure.primary() == null);
    try testing.expect(failure.secondaryItems().len == 0 and !failure.secondary_truncated);
    try testing.expect(failure.preview_len == 0 and failure.sources.registered == 0);
    try testing.expect(!failure.complete and !failure.work_complete);
    try testing.expect(failure.outputs.stdout.logical_bytes == null);
    try testing.expect(failure.outputs.file.state == .not_started);
    try testing.expectEqual(@as(u8, 0xa5), failure.previews[255]);
    try testing.expectEqual(@as(u8, 'o'), failure.previews[0]);
    try testing.expect(failure.primary_storage.stage == .write);
    try testing.expect(failure.primary_storage.cause.?.zig == error.WriteFailed);
    failure.succeed();
    try testing.expect(failure.complete and failure.work_complete);
}

test "[failure] - [diagnostics]: cleanup and reporting preserve primary and output" {
    var failure: diagnostics.Failure = undefined;
    failure.reset(.validate);
    failure.work_complete = true;
    failure.outputs.stdout = .{ .state = .streaming, .logical_bytes = 17, .records = 1, .may_be_partial = true };
    failure.outputs.file.state = .published;
    try testing.expect(failure.fail(.{
        .details = .{ .@"io.flush_failed" = {} },
        .stage = .flush,
        .at = .{ .source = .stderr },
        .cause = .{ .zig = error.WriteFailed },
    }) == error.Io);
    failure.addSecondary(.{
        .details = .{ .@"io.cleanup_failed" = .{ .resource_role = .name_spool } },
        .stage = .cleanup,
        .cause = .{ .zig = error.AccessDenied },
    });
    failure.addSecondary(writeFailure(.report, error.BrokenPipe));
    failure.addSecondary(writeFailure(.report, error.WriteFailed));
    try testing.expect(failure.primary().?.code() == .@"io.flush_failed");
    try testing.expect(failure.work_complete and !failure.complete);
    try testing.expect(failure.outputs.file.state == .published);
    try testing.expect(failure.outputs.stdout.may_be_partial);
    try testing.expectEqual(@as(usize, 2), failure.secondaryItems().len);
    try testing.expect(failure.secondary_truncated);
    try testing.expect(failure.secondaryItems()[0].code() == .@"io.cleanup_failed");
    try testing.expect(failure.secondaryItems()[0].cause.?.zig == error.AccessDenied);
}

test "[failure] - [diagnostics]: explicit priority retains the first deferred cause" {
    var failure: diagnostics.Failure = undefined;
    failure.reset(.index);
    const token = failure.retain("bad", 3);
    try testing.expect(failure.fail(writeFailure(.write, error.NoSpaceLeft)) == error.Io);
    failure.addSecondary(writeFailure(.flush, error.WriteFailed));
    failure.addSecondary(writeFailure(.cleanup, error.AccessDenied));
    const returned = failure.promote(.{
        .details = .{ .@"index.fai.invalid_integer" = .{ .field = .length, .value = token } },
        .stage = .parse,
        .at = .{ .source = .fai, .location = .{ .line = 4, .column = 2 } },
    });
    try testing.expect(returned == error.InvalidIndex);
    try testing.expect(failure.primary().?.code() == .@"index.fai.invalid_integer");
    try testing.expectEqualStrings("bad", failure.preview(failure.primary().?.details.@"index.fai.invalid_integer".value.?));
    try testing.expect(failure.secondaryItems()[0].stage == .write);
    try testing.expect(failure.secondaryItems()[0].cause.?.zig == error.NoSpaceLeft);
    try testing.expect(failure.secondaryItems()[1].stage == .flush);
    try testing.expect(failure.secondary_truncated);
    try testing.expect(failure.fail(writeFailure(.cleanup, error.AccessDenied)) == error.InvalidIndex);
}

test "[edge] - [diagnostics]: exact EOF and missing delimiter keep different facts" {
    const exact = diagnostics.Diagnostic{
        .details = .{ .@"io.unexpected_eof" = .{ .exact = .{ .requested_bytes = 64, .received_bytes = 3 } } },
        .stage = .read,
        .at = .{ .location = .{ .offset = 128 } },
    };
    const delimiter = diagnostics.Diagnostic{
        .details = .{ .@"io.unexpected_eof" = .{ .delimiter = .{ .delimiter = '\t', .observed_bytes = .{ .at_least = 3 } } } },
        .stage = .read,
    };
    try testing.expectEqual(@as(u64, 3), exact.details.@"io.unexpected_eof".exact.received_bytes);
    try testing.expectEqual(@as(u128, 3), delimiter.details.@"io.unexpected_eof".delimiter.observed_bytes.?.at_least);
    try testing.expect(delimiter.at.location.offset == null and delimiter.cause == null);
}

test "[edge] - [diagnostics]: wide arithmetic and signed timestamps stay native" {
    var failure: diagnostics.Failure = undefined;
    failure.reset(.index);
    const required = @as(u128, std.math.maxInt(u64)) + 40;
    try testing.expect(failure.fail(.{
        .details = .{ .@"index.truncated" = .{ .section = .records, .required_end = required, .available_end = 16 } },
        .stage = .parse,
        .at = .{ .location = .{ .offset = 0 } },
    }) == error.InvalidIndex);
    try testing.expectEqual(required, failure.primary().?.details.@"index.truncated".required_end);
    failure.reset(.index);
    try testing.expect(failure.fail(.{
        .details = .{ .@"index.source_timestamp_out_of_range" = .{
            .timestamp_ns = -1,
            .min_ns = 0,
            .max_ns = std.math.maxInt(u64),
        } },
        .stage = .compare_source,
    }) == error.InvalidIndex);
    try testing.expectEqual(@as(i128, -1), failure.primary().?.details.@"index.source_timestamp_out_of_range".timestamp_ns);
}

test "[unit] - [diagnostics]: replay scope and cached sample origin remain explicit" {
    const missing = diagnostics.Diagnostic{
        .details = .{ .@"request.sequence_not_found" = .{ .lookup_scope = .retained_matches, .index_source = .zfi } },
        .stage = .resolve,
        .at = .{ .source = .bed, .location = .{ .line = 8193 }, .request = 8190 },
    };
    const transform = diagnostics.Diagnostic{
        .details = .{ .@"transform.complement_on_protein" = .{
            .transform = .{ .reverse = false, .complement = true },
            .classification = .protein,
            .sample = .{ .symbols = 19, .rule = .nucleotide_above_90_percent, .interval = .{ .start = 5, .end = 23 }, .request = 2 },
        } },
        .stage = .transform,
        .at = .{ .request = 3 },
    };
    try testing.expect(missing.details.@"request.sequence_not_found".lookup_scope == .retained_matches);
    try testing.expectEqual(@as(u64, 8193), missing.at.location.line.?);
    const facts = transform.details.@"transform.complement_on_protein";
    try testing.expectEqual(@as(u32, 19), facts.sample.symbols);
    try testing.expectEqual(@as(u64, 2), facts.sample.request.?);
    try testing.expectEqual(@as(u64, 3), transform.at.request.?);
    try testing.expect(!facts.transform.reverse and facts.transform.complement);
}

test "[unit] - [diagnostics]: rewrite refusal retains the precise blocking byte" {
    const diagnostic = diagnostics.Diagnostic{
        .details = .{ .@"repair.blocked" = .{
            .blocker = .{ .invalid_symbol = .{ .byte = '?', .policy = .{ .alphabet = .dna, .selection = .explicit } } },
            .repair_policy = .legacy_safe,
        } },
        .stage = .scan,
        .at = .{ .source = .fasta, .location = .{ .offset = 18, .line = 3, .column = 5 }, .record = 2 },
    };
    const blocker = &diagnostic.details.@"repair.blocked".blocker;
    try testing.expect(blocker.code() == .@"sequence.invalid_symbol");
    try testing.expectEqual(@as(u8, '?'), blocker.invalid_symbol.byte);
    try testing.expect(diagnostics.metadata(diagnostic.code()).kind == .failure);
    try testing.expect(diagnostics.metadata(blocker.code()).kind == .finding);
}

test "[edge] - [diagnostics]: temporary attempts and native causes are bounded" {
    var attempts: diagnostics.TempAttempts = .{};
    for (0..7) |index| attempts.append(.{
        .place = .system,
        .stage = .create,
        .reason = .creation,
        .cause = .{ .native = .{ .domain = .posix_errno, .number = @intCast(13 + index) } },
    });
    try testing.expectEqual(@as(usize, 4), attempts.items().len);
    try testing.expect(attempts.truncated);
    try testing.expectEqual(@as(u32, 13), attempts.items()[0].cause.?.native.number);
    try testing.expectEqual(@as(u32, 16), attempts.items()[3].cause.?.native.number);
}

test "[failure] - [diagnostics]: every partial allocation cleans up without diagnostic allocation" {
    try testing.checkAllAllocationFailures(testing.allocator, allocationCheck, .{});
    try caller.exercise(testing.allocator);
}

test "[unit] - [diagnostics]: maps open and read errors to their codes" {
    const cases = [_]struct { err: anyerror, code: diagnostics.Code }{
        .{ .err = error.FileNotFound, .code = .@"io.file_not_found" },
        .{ .err = error.AccessDenied, .code = .@"io.access_denied" },
        .{ .err = error.NameTooLong, .code = .@"io.invalid_path" },
        .{ .err = error.NotDir, .code = .@"io.open_failed" },
        .{ .err = error.SymLinkLoop, .code = .@"io.open_failed" },
    };
    for (cases) |case| {
        const diagnostic = diagnostics.Diagnostic{ .details = diagnostics.openFailure(case.err), .stage = .open, .at = .{} };
        try testing.expectEqual(case.code, diagnostic.code());
    }
    try testing.expect(diagnostics.openFailure(error.NameTooLong).@"io.invalid_path".reason == .too_long);

    try testing.expectEqual(@as(?u64, 10), diagnostics.readFailure(10, null).@"io.read_failed".requested_bytes);
    const eof = diagnostics.readFailure(10, 4).@"io.unexpected_eof".exact;
    try testing.expectEqual(@as(u64, 10), eof.requested_bytes);
    try testing.expectEqual(@as(u64, 4), eof.received_bytes);
}
