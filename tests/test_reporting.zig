//! Renderer output, byte identity, limits, and failure delivery checks.

const std = @import("std");

const diagnostics = @import("diagnostics");
const reporting = @import("reporting");

const testing = std.testing;

fn object(value: std.json.Value, name: []const u8) std.json.Value {
    return value.object.get(name).?;
}

const LONG_CAUSE = error.@"\x01\x01\x01\x01\x01\x01\x01\x01\x01\x01\x01\x01\x01\x01\x01\x01\x01\x01\x01\x01\x01\x01\x01\x01\x01\x01\x01\x01\x01\x01\x01\x01long_error_name";

fn largest(comptime T: type, preview: diagnostics.Preview) T {
    if (T == diagnostics.Preview) return preview;
    if (T == diagnostics.Cause) return .{ .zig = LONG_CAUSE };
    if (T == diagnostics.Arguments) return .initFull();
    if (T == diagnostics.Measurement) return .{ .at_least = std.math.maxInt(u128) };
    if (T == diagnostics.TempAttempts) {
        var attempts: T = .{};
        for (0..diagnostics.TEMP_ATTEMPT_LIMIT) |_| attempts.append(largest(diagnostics.TempAttempt, preview));
        attempts.truncated = true;
        return attempts;
    }
    if (T == diagnostics.Eof) return .{ .delimiter = .{ .delimiter = 255, .observed_bytes = .{ .at_least = std.math.maxInt(u128) } } };
    if (T == diagnostics.Blocker) return .{ .invalid_symbol = largest(diagnostics.InvalidSymbol, preview) };
    return switch (@typeInfo(T)) {
        .void => {},
        .bool => true,
        .int => |integer| if (integer.signedness == .signed) std.math.minInt(T) else std.math.maxInt(T),
        .optional => |optional| largest(optional.child, preview),
        .array => @splat(1),
        .@"enum" => |enumeration| blk: {
            comptime var longest = enumeration.fields[0];
            inline for (enumeration.fields) |member| {
                if (member.name.len > longest.name.len) longest = member;
            }
            break :blk @enumFromInt(longest.value);
        },
        .@"struct" => blk: {
            var value: T = undefined;
            inline for (std.meta.fields(T)) |member| @field(value, member.name) = largest(member.type, preview);
            break :blk value;
        },
        else => @compileError("add a bounded renderer fixture for " ++ @typeName(T)),
    };
}

const BrokenWriter = struct {
    writer: std.Io.Writer,
    calls: usize = 0,

    fn drain(writer: *std.Io.Writer, data: []const []const u8, splat: usize) std.Io.Writer.Error!usize {
        _ = data;
        _ = splat;
        const self: *BrokenWriter = @fieldParentPtr("writer", writer);
        self.calls += 1;
        return error.WriteFailed;
    }
};

test "[unit] - [reporting]: file output states distinguish abort from staging and publication" {
    inline for (.{ .not_started, .staged, .aborted, .published, .streaming, .unknown }, [_][]const u8{ "not_started", "staged", "aborted", "published", "streaming", "unknown" }) |state, expected| {
        var failure: diagnostics.Failure = undefined;
        failure.reset(.index);
        try testing.expectEqual(error.Io, failure.fail(.{ .details = .{ .@"io.write_failed" = {} }, .stage = .write }));
        failure.outputs.file.state = state;
        var buffer: [reporting.JSON_BYTES_LIMIT]u8 = undefined;
        var writer = std.Io.Writer.fixed(&buffer);
        try reporting.writeJson(&writer, &failure, .{});
        var parsed = try std.json.parseFromSlice(std.json.Value, testing.allocator, writer.buffered(), .{});
        defer parsed.deinit();
        try testing.expectEqualStrings(expected, object(object(parsed.value, "outputs").array.items[1], "state").string);
    }
}

test "[unit] - [reporting]: machine facts preserve wide values and publication state" {
    var failure: diagnostics.Failure = undefined;
    failure.reset(.get);
    failure.sources.set(.fasta, .{ .kind = .path, .label = .{ .borrowed = "sample.fa" } });
    failure.sources.set(.zfi, .{ .kind = .path, .label = .{ .borrowed = "sample.fa" }, .suffix = .zfi });
    failure.deferFailure(.{
        .details = .{ .@"index.source_size_mismatch" = .{ .stored_bytes = 9007199254740993, .actual_bytes = 0 } },
        .stage = .parse,
        .at = .{ .source = .zfi, .location = .{ .offset = 0 } },
        .related = .{ .source = .fasta },
    });
    failure.work_complete = true;
    failure.outputs.file.state = .published;
    failure.addSecondary(.{ .details = .{ .@"io.flush_failed" = {} }, .stage = .flush, .at = .{ .source = .stderr }, .cause = .{ .native = .{ .domain = .posix_errno, .number = 32 } } });
    var buffer: [reporting.JSON_BYTES_LIMIT]u8 = undefined;
    var writer = std.Io.Writer.fixed(&buffer);
    try reporting.writeJson(&writer, &failure, .{ .version = "0.3.5", .accepted = true });

    try testing.expectEqual(@as(u8, '\n'), writer.buffered()[writer.end - 1]);
    var parsed = try std.json.parseFromSlice(std.json.Value, testing.allocator, writer.buffered(), .{});
    defer parsed.deinit();
    const envelope = parsed.value;
    try testing.expectEqualStrings("v2", object(envelope, "schema_version").string);
    try testing.expectEqualStrings("0.4.0-1", object(envelope, "catalog_revision").string);
    try testing.expectEqualStrings("0.3.5", object(object(envelope, "producer"), "version").string);
    try testing.expect(!object(envelope, "accepted").bool);
    try testing.expect(object(object(envelope, "operation"), "work_complete").bool);
    try testing.expect(!object(object(envelope, "operation"), "complete").bool);
    try testing.expect(object(envelope, "validation") == .null);
    const primary = object(object(envelope, "failure"), "primary");
    try testing.expectEqualStrings("index.source_size_mismatch", object(primary, "code").string);
    try testing.expectEqualStrings("9007199254740993", object(object(primary, "data"), "stored_bytes").string);
    try testing.expectEqualStrings("0", object(object(primary, "data"), "actual_bytes").string);
    try testing.expectEqualStrings("0", object(object(primary, "location"), "offset").string);
    try testing.expect(object(object(primary, "location"), "line") == .null);
    try testing.expectEqualStrings("published", object(object(envelope, "outputs").array.items[1], "state").string);
    const sources = object(envelope, "sources").array.items;
    try testing.expectEqual(@as(usize, 3), sources.len);
    try testing.expectEqualStrings("sample.fa.zfi", object(object(sources[1], "label"), "value").string);
    try testing.expect(object(sources[2], "kind") == .null);
    try testing.expectEqualStrings("rebuild_index", object(object(envelope, "action"), "id").string);
}

test "[edge] - [reporting]: opaque bytes and lower bounds retain their identity" {
    var failure: diagnostics.Failure = undefined;
    failure.reset(.get);
    const preview = failure.retainPrefix("bad\xff\x00", .{ .at_least = std.math.maxInt(u64) }, 256);
    failure.sources.set(.names, .{ .kind = .path, .label = .{ .retained = preview }, .suffix = .fai });
    failure.deferFailure(.{ .details = .{ .@"request.sequence_not_found" = .{ .lookup_scope = .retained_matches, .index_source = .zfi, .identifier = preview } }, .stage = .resolve, .at = .{ .source = .names } });
    var buffer: [reporting.JSON_BYTES_LIMIT]u8 = undefined;
    var writer = std.Io.Writer.fixed(&buffer);
    try reporting.writeJson(&writer, &failure, .{});
    var parsed = try std.json.parseFromSlice(std.json.Value, testing.allocator, writer.buffered(), .{});
    defer parsed.deinit();
    const label = object(object(parsed.value, "sources").array.items[1], "label");
    try testing.expectEqualStrings("base64", object(label, "encoding").string);
    try testing.expectEqualStrings("YmFk/wA=", object(label, "value").string);
    try testing.expectEqualStrings("18446744073709551619", object(label, "observed_bytes_min").string);
    try testing.expect(label.object.get("original_bytes") == null);
}

test "[unit] - [reporting]: text escapes names and supplies observed size facts" {
    var failure: diagnostics.Failure = undefined;
    failure.reset(.get);
    failure.sources.set(.zfi, .{ .kind = .path, .label = .{ .borrowed = "a\n\xff\"\\" } });
    failure.deferFailure(.{ .details = .{ .@"index.source_size_mismatch" = .{ .stored_bytes = 35, .actual_bytes = 36 } }, .stage = .parse, .at = .{ .source = .zfi, .location = .{ .offset = 8 } } });
    var buffer: [reporting.TEXT_BYTES_LIMIT]u8 = undefined;
    var writer = std.Io.Writer.fixed(&buffer);
    try reporting.writeText(&writer, &failure);
    try testing.expectEqualStrings("error [index.source_size_mismatch] get: \"a\\x0a\\xff\\\"\\\\\" offset=8: Index source size differs from the FASTA size. stored_bytes=35 actual_bytes=36\nhint: Rebuild the index after confirming the intended FASTA input.\n", writer.buffered());
}

test "[unit] - [reporting]: requested metrics stay strings and acceptance is supplied by the caller" {
    var failure: diagnostics.Failure = undefined;
    failure.reset(.get);
    failure.succeed();
    var buffer: [reporting.JSON_BYTES_LIMIT]u8 = undefined;
    var writer = std.Io.Writer.fixed(&buffer);
    try reporting.writeJson(&writer, &failure, .{ .accepted = true, .metrics = .{ .region_count = 1, .total_bases = 0, .elapsed_ns = std.math.maxInt(u64) } });
    var parsed = try std.json.parseFromSlice(std.json.Value, testing.allocator, writer.buffered(), .{});
    defer parsed.deinit();
    try testing.expectEqualStrings("summary", object(parsed.value, "kind").string);
    try testing.expect(object(parsed.value, "accepted").bool);
    try testing.expect(object(parsed.value, "failure") == .null);
    const metrics = object(parsed.value, "metrics");
    try testing.expectEqualStrings("1", object(metrics, "region_count").string);
    try testing.expectEqualStrings("0", object(metrics, "total_bases").string);
    try testing.expectEqualStrings("18446744073709551615", object(metrics, "elapsed_ns").string);
    writer.end = 0;
    try reporting.writeJson(&writer, &failure, .{});
    var rejected = try std.json.parseFromSlice(std.json.Value, testing.allocator, writer.buffered(), .{});
    defer rejected.deinit();
    try testing.expect(!object(rejected.value, "accepted").bool);
}

test "[property] - [reporting]: every code fits both budgets with full sources and secondaries" {
    @setEvalBranchQuota(100000);
    var failure: diagnostics.Failure = undefined;
    var preview_bytes: [diagnostics.PREVIEW_BYTES]u8 = @splat(1);
    var json_buffer: [reporting.JSON_BYTES_LIMIT]u8 = undefined;
    var text_buffer: [reporting.TEXT_BYTES_LIMIT]u8 = undefined;
    inline for (std.meta.fields(diagnostics.Details)) |member| {
        failure.reset(.validate);
        const preview = failure.retainPrefix(&preview_bytes, .{ .exact = std.math.maxInt(u64) }, 256);
        for (std.enums.values(diagnostics.SourceId)) |id| {
            failure.sources.set(id, .{ .kind = .path, .label = .{ .retained = preview }, .suffix = .zfi });
        }
        const reference = largest(diagnostics.Reference, preview);
        const diagnostic: diagnostics.Diagnostic = .{
            .details = @unionInit(diagnostics.Details, member.name, largest(member.type, preview)),
            .stage = .compare_source,
            .at = reference,
            .related = reference,
            .cause = .{ .zig = LONG_CAUSE },
        };
        failure.deferFailure(diagnostic);
        failure.addSecondary(diagnostic);
        failure.addSecondary(diagnostic);
        failure.addSecondary(diagnostic);
        failure.outputs = largest(diagnostics.Outputs, preview);
        var json_writer = std.Io.Writer.fixed(&json_buffer);
        try reporting.writeJson(&json_writer, &failure, .{ .version = "12345678901234567890123456789012", .metrics = largest(reporting.GetMetrics, preview) });
        var parsed = try std.json.parseFromSlice(std.json.Value, testing.allocator, json_writer.buffered(), .{});
        defer parsed.deinit();
        const primary = object(object(parsed.value, "failure"), "primary");
        try testing.expectEqualStrings(member.name, object(primary, "code").string);
        try testing.expectEqual(@as(usize, 16), object(parsed.value, "sources").array.items.len);
        try testing.expect(object(object(parsed.value, "failure"), "secondary_truncated").bool);
        var text_writer = std.Io.Writer.fixed(&text_buffer);
        reporting.writeText(&text_writer, &failure) catch |err| {
            std.debug.print("text budget exceeded for {s}\n", .{member.name});
            return err;
        };
        for (text_writer.buffered()) |byte| try testing.expect(byte == '\n' or (byte >= 32 and byte < 127));
    }
}

test "[edge] - [reporting]: typed exceptional encodings retain read kind and blocker facts" {
    var failure: diagnostics.Failure = undefined;
    failure.reset(.validate);
    failure.deferFailure(.{ .details = .{ .@"io.unexpected_eof" = .{ .delimiter = .{ .delimiter = 9, .observed_bytes = .{ .at_least = 9007199254740993 } } } }, .stage = .read });
    failure.addSecondary(.{ .details = .{ .@"index.unsupported_version" = .{ .found_version = 2, .supported_version = 1 } }, .stage = .parse });
    failure.addSecondary(.{ .details = .{ .@"repair.blocked" = .{ .blocker = .{ .invalid_symbol = .{ .byte = 63, .policy = .{ .alphabet = .dna, .selection = .automatic, .sample = .{ .symbols = 256, .rule = .nucleotide_above_90_percent } } } }, .repair_policy = .legacy_safe } }, .stage = .scan });
    var buffer: [reporting.JSON_BYTES_LIMIT]u8 = undefined;
    var writer = std.Io.Writer.fixed(&buffer);
    try reporting.writeJson(&writer, &failure, .{});
    var parsed = try std.json.parseFromSlice(std.json.Value, testing.allocator, writer.buffered(), .{});
    defer parsed.deinit();
    const result = object(parsed.value, "failure");
    const eof = object(object(result, "primary"), "data");
    try testing.expectEqualStrings("delimiter", object(eof, "read_kind").string);
    try testing.expectEqual(@as(i64, 9), object(eof, "delimiter").integer);
    try testing.expectEqualStrings("9007199254740993", object(eof, "observed_bytes_min").string);
    const version = object(object(result, "secondary").array.items[0], "data");
    try testing.expectEqual(@as(i64, 1), object(version, "supported_versions").array.items[0].integer);
    const repair = object(object(result, "secondary").array.items[1], "data");
    try testing.expectEqualStrings("sequence.invalid_symbol", object(repair, "reason_code").string);
    const policy = object(object(repair, "blocker"), "policy");
    try testing.expectEqualStrings("256", object(object(policy, "sample"), "symbols").string);
}

test "[failure] - [reporting]: writes and late flush fail once without changing the original result" {
    var failure: diagnostics.Failure = undefined;
    failure.reset(.index);
    failure.deferFailure(.{ .details = .{ .@"resource.out_of_memory" = .{ .allocation_role = .records } }, .stage = .allocate });
    var small: [8]u8 = undefined;
    var broken: BrokenWriter = .{ .writer = .{ .buffer = &small, .vtable = &.{ .drain = BrokenWriter.drain } } };
    try testing.expectError(error.WriteFailed, reporting.writeJson(&broken.writer, &failure, .{}));
    try testing.expectEqual(@as(usize, 1), broken.calls);
    try testing.expect(broken.writer.end > 0);
    broken.calls = 0;
    broken.writer.end = 0;
    try testing.expectError(error.WriteFailed, reporting.writeText(&broken.writer, &failure));
    try testing.expectEqual(@as(usize, 1), broken.calls);

    var buffer: [reporting.JSON_BYTES_LIMIT]u8 = undefined;
    broken = .{ .writer = .{ .buffer = &buffer, .vtable = &.{ .drain = BrokenWriter.drain } } };
    try reporting.writeJson(&broken.writer, &failure, .{});
    try testing.expectEqual(@as(usize, 0), broken.calls);
    try testing.expectError(error.WriteFailed, broken.writer.flush());
    try testing.expectEqual(@as(usize, 1), broken.calls);
    try testing.expect(failure.primary().?.code() == .@"resource.out_of_memory");
    try testing.expect(failure.secondary_len == 0 and !failure.complete);
}

test "[edge] - [reporting]: bootstrap messages fit fixed storage and return broken delivery" {
    for (std.enums.values(reporting.Bootstrap)) |reason| {
        var buffer: [256]u8 = undefined;
        var writer = std.Io.Writer.fixed(&buffer);
        try reporting.writeBootstrap(&writer, reason);
        try testing.expect(writer.end > 0);
        var tiny: [1]u8 = undefined;
        var broken = std.Io.Writer.fixed(&tiny);
        try testing.expectError(error.WriteFailed, reporting.writeBootstrap(&broken, reason));
    }
}

test "[edge] - [reporting]: UTF-8 and invalid bytes remain distinct with disabled previews" {
    var failure: diagnostics.Failure = undefined;
    var buffer: [reporting.JSON_BYTES_LIMIT]u8 = undefined;
    for ([_][]const u8{ "bad\xc3\xbf\n", "bad\xff\n", "" }, [_][]const u8{ "utf8", "base64", "utf8" }, [_][]const u8{ "bad\xc3\xbf\n", "YmFk/wo=", "" }) |bytes, encoding, expected| {
        failure.reset(.validate);
        const preview = failure.retain(bytes, 256);
        failure.deferFailure(.{ .details = .{ .@"fasta.duplicate_identifier" = .{ .identifier = preview } }, .stage = .scan });
        var writer = std.Io.Writer.fixed(&buffer);
        try reporting.writeJson(&writer, &failure, .{});
        var parsed = try std.json.parseFromSlice(std.json.Value, testing.allocator, writer.buffered(), .{});
        defer parsed.deinit();
        const identifier = object(object(object(object(parsed.value, "failure"), "primary"), "data"), "identifier");
        try testing.expectEqualStrings(encoding, object(identifier, "encoding").string);
        try testing.expectEqualStrings(expected, object(identifier, "value").string);
        try testing.expect(!object(identifier, "truncated").bool);
    }
    failure.reset(.validate);
    const disabled = failure.retainPrefix("", .{ .exact = 9007199254740993 }, 0);
    failure.deferFailure(.{ .details = .{ .@"fasta.duplicate_identifier" = .{ .identifier = disabled } }, .stage = .scan });
    var writer = std.Io.Writer.fixed(&buffer);
    try reporting.writeJson(&writer, &failure, .{});
    var parsed = try std.json.parseFromSlice(std.json.Value, testing.allocator, writer.buffered(), .{});
    defer parsed.deinit();
    const identifier = object(object(object(object(parsed.value, "failure"), "primary"), "data"), "identifier");
    try testing.expectEqualStrings("", object(identifier, "value").string);
    try testing.expectEqualStrings("9007199254740993", object(identifier, "original_bytes").string);
    try testing.expect(object(identifier, "truncated").bool);
}
