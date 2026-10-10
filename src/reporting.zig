//! Bounded text and NDJSON rendering of caller-owned diagnostic facts.
//! Rendering borrows the result, never allocates or flushes, and returns a
//! failed write to the adapter without changing the original failure.

const std = @import("std");
const diagnostics = @import("diagnostics");

pub const TEXT_BYTES_LIMIT = 4096;
pub const JSON_BYTES_LIMIT = 16384;
const SOURCE_PREFIX_BYTES = 48;
const CAUSE_PREFIX_BYTES = 32;

pub const GetMetrics = struct { region_count: u64, total_bases: u64, elapsed_ns: u64 };

/// Version is trusted UTF-8 product metadata of at most 32 bytes. The adapter
/// supplies acceptance; incomplete work or a blocking failure overrides it.
pub const Options = struct {
    version: ?[]const u8 = null,
    accepted: bool = false,
    metrics: ?GetMetrics = null,
};

pub const Bootstrap = enum { argv_oom, missing_command, unknown_command, invalid_format };

/// Static reporting before a normal result exists. A failed write returns once;
/// the adapter must not attempt another report to the same broken sink.
pub fn writeBootstrap(writer: *std.Io.Writer, reason: Bootstrap) std.Io.Writer.Error!void {
    const code: diagnostics.Code = switch (reason) {
        .argv_oom => .@"resource.out_of_memory",
        .missing_command => .@"usage.command_required",
        .unknown_command => .@"usage.unknown_command",
        .invalid_format => .@"usage.invalid_option_value",
    };
    try writer.print("error [{s}]: {s}\n", .{ @tagName(code), switch (reason) {
        .argv_oom => "Out of memory while reading command arguments.",
        .missing_command => "A command is required. Use --help to list commands.",
        .unknown_command => "Unknown command. Use --help to list commands.",
        .invalid_format => "Use --error-format=text or --error-format=json.",
    } });
}

/// Writes one newline-terminated object. Does not certify delivery of its own
/// bytes; the adapter must check its writer's final flush and process status.
pub fn writeJson(writer: *std.Io.Writer, failure: *const diagnostics.Failure, options: Options) std.Io.Writer.Error!void {
    if (options.version) |version| {
        std.debug.assert(version.len <= 32 and std.unicode.utf8ValidateSlice(version));
    }
    var json: std.json.Stringify = .{ .writer = writer };
    try json.beginObject();
    try field(&json, "schema_version", "v2");
    try field(&json, "catalog_revision", "draft-2026-10-08");
    try json.objectField("producer");
    try json.beginObject();
    try field(&json, "name", "z-fasta");
    try field(&json, "version", options.version);
    try json.endObject();
    try field(&json, "kind", if (failure.has_primary) "failure" else "summary");
    try json.objectField("operation");
    try json.beginObject();
    try field(&json, "name", @tagName(failure.operation));
    try field(&json, "complete", failure.complete);
    try field(&json, "work_complete", failure.work_complete);
    try json.endObject();
    try field(&json, "accepted", options.accepted and failure.complete and !failure.has_primary);
    try field(&json, "validation", null);
    try json.objectField("sources");
    try json.beginArray();
    const sources = sourceMask(failure);
    for (std.enums.values(diagnostics.SourceId)) |id| {
        if (sources & (@as(u16, 1) << @intFromEnum(id)) == 0) continue;
        try json.beginObject();
        try field(&json, "id", @tagName(id));
        if (failure.sources.get(id)) |source| {
            try field(&json, "kind", @tagName(source.kind));
            try json.objectField("label");
            var label_buffer: [SOURCE_PREFIX_BYTES]u8 = undefined;
            const label = sourcePrefix(failure, id, &label_buffer);
            try writeOpaque(&json, label.bytes, label.original, label.lower_bound, label.truncated);
        } else {
            try field(&json, "kind", null);
            try field(&json, "label", null);
        }
        try json.endObject();
    }
    try json.endArray();
    try json.objectField("failure");
    if (failure.primary()) |primary| {
        try json.beginObject();
        try json.objectField("primary");
        try jsonDiagnostic(&json, failure, primary);
        try json.objectField("secondary");
        try json.beginArray();
        for (failure.secondaryItems()) |secondary| try jsonDiagnostic(&json, failure, &secondary);
        try json.endArray();
        try field(&json, "secondary_truncated", failure.secondary_truncated);
        try json.endObject();
    } else try json.write(null);
    try json.objectField("outputs");
    try json.beginArray();
    inline for (std.meta.fields(diagnostics.Outputs)) |output| {
        try json.beginObject();
        try field(&json, "destination", output.name);
        try jsonStructFields(&json, failure, @field(failure.outputs, output.name));
        try json.endObject();
    }
    try json.endArray();
    try json.objectField("metrics");
    if (options.metrics) |metrics| try jsonValue(&json, failure, metrics) else try json.write(null);
    try json.objectField("action");
    if (failure.primary()) |primary| {
        try json.beginObject();
        try field(&json, "id", @tagName(diagnostics.metadata(primary.code()).action));
        try jsonField(&json, failure, "source", actionSource(primary));
        try json.endObject();
    } else try json.write(null);
    try json.endObject();
    try writer.writeByte('\n');
}

/// Writes escaped primary/secondary facts and a concrete action. A success
/// result has no failure text. The caller owns all write/flush failure policy.
pub fn writeText(writer: *std.Io.Writer, failure: *const diagnostics.Failure) std.Io.Writer.Error!void {
    const primary = failure.primary() orelse return;
    try textDiagnostic(writer, failure, primary, "error");
    for (failure.secondaryItems()) |secondary| try textDiagnostic(writer, failure, &secondary, "secondary");
    if (failure.secondary_truncated) try writer.writeAll("Further secondary failures were omitted.\n");
    try writer.print("hint: {s}\n", .{actionText(diagnostics.metadata(primary.code()).action)});
}

fn field(json: *std.json.Stringify, name: []const u8, value: anytype) std.Io.Writer.Error!void {
    try json.objectField(name);
    try json.write(value);
}

fn jsonField(json: *std.json.Stringify, failure: *const diagnostics.Failure, comptime name: []const u8, value: anytype) std.Io.Writer.Error!void {
    if (@TypeOf(value) == ?diagnostics.Measurement) {
        if (value) |present| try jsonField(json, failure, name, present) else try field(json, name, null);
    } else if (@TypeOf(value) == u32 and std.mem.eql(u8, name, "symbols")) {
        try jsonField(json, failure, name, @as(u64, value));
    } else if (@TypeOf(value) == diagnostics.Measurement) {
        switch (value) {
            .exact => |number| {
                try json.objectField(name);
                try jsonValue(json, failure, number);
            },
            .at_least => |number| {
                try json.objectField(name ++ "_min");
                try jsonValue(json, failure, number);
            },
        }
    } else {
        try json.objectField(name);
        try jsonValue(json, failure, value);
    }
}

fn jsonStructFields(json: *std.json.Stringify, failure: *const diagnostics.Failure, value: anytype) std.Io.Writer.Error!void {
    inline for (std.meta.fields(@TypeOf(value))) |member| {
        if (!std.mem.eql(u8, member.name, "reserved")) {
            try jsonField(json, failure, member.name, @field(value, member.name));
        }
    }
}

fn jsonValue(json: *std.json.Stringify, failure: *const diagnostics.Failure, value: anytype) std.Io.Writer.Error!void {
    const T = @TypeOf(value);
    if (T == diagnostics.Preview) {
        try writeOpaque(json, failure.preview(value), value.original_bytes, value.length_kind == .at_least, value.truncated);
    } else if (T == diagnostics.Cause) {
        try json.beginObject();
        switch (value) {
            .zig => |err| {
                try field(json, "domain", "zig");
                try json.objectField("name");
                const name = @errorName(err);
                try writeOpaque(json, name[0..@min(name.len, CAUSE_PREFIX_BYTES)], name.len, false, name.len > CAUSE_PREFIX_BYTES);
            },
            .native => |native| {
                try field(json, "domain", @tagName(native.domain));
                try field(json, "number", native.number);
            },
        }
        try json.endObject();
    } else if (T == diagnostics.Arguments) {
        try json.beginArray();
        for (std.enums.values(diagnostics.Argument)) |argument| {
            if (value.contains(argument)) try json.write(@tagName(argument));
        }
        try json.endArray();
    } else if (T == diagnostics.TempAttempts) {
        try json.beginObject();
        try json.objectField("entries");
        try json.beginArray();
        for (value.items()) |attempt| try jsonValue(json, failure, attempt);
        try json.endArray();
        try field(json, "truncated", value.truncated);
        try json.endObject();
    } else if (T == diagnostics.Eof) {
        try json.beginObject();
        try field(json, "read_kind", @tagName(value));
        switch (value) {
            inline else => |facts| try jsonStructFields(json, failure, facts),
        }
        try json.endObject();
    } else if (T == diagnostics.Version) {
        try json.beginObject();
        try field(json, "found_version", value.found_version);
        try json.objectField("supported_versions");
        try json.beginArray();
        try json.write(value.supported_version);
        try json.endArray();
        try json.endObject();
    } else if (T == diagnostics.IndexDiscovery) {
        try json.beginObject();
        try field(json, "attempted_formats", [_][]const u8{ "zfi", "fai" });
        try json.objectField("attempted_sources");
        try json.write([_]diagnostics.SourceId{ value.zfi_source, value.fai_source });
        try json.endObject();
    } else if (T == diagnostics.RepairBlock) {
        try json.beginObject();
        try field(json, "reason_code", @tagName(value.blocker.code()));
        try json.objectField("blocker");
        switch (value.blocker) {
            inline else => |facts| try jsonValue(json, failure, facts),
        }
        try field(json, "repair_policy", @tagName(value.repair_policy));
        try json.endObject();
    } else if (T == diagnostics.Measurement) {
        try json.beginObject();
        switch (value) {
            .exact => |number| try jsonField(json, failure, "exact", number),
            .at_least => |number| try jsonField(json, failure, "at_least", number),
        }
        try json.endObject();
    } else switch (@typeInfo(T)) {
        .void => {
            try json.beginObject();
            try json.endObject();
        },
        .int => |integer| {
            if (integer.bits > 32) {
                var buffer: [40]u8 = undefined;
                try json.write(std.fmt.bufPrint(&buffer, "{d}", .{value}) catch unreachable);
            } else try json.write(value);
        },
        .optional => if (value) |present| try jsonValue(json, failure, present) else try json.write(null),
        .@"struct" => {
            try json.beginObject();
            try jsonStructFields(json, failure, value);
            try json.endObject();
        },
        .array => |array| {
            if (array.child != u8) @compileError("diagnostic arrays must declare their encoding");
            try writeOpaque(json, &value, value.len, false, false);
        },
        .@"enum", .bool => try json.write(value),
        else => @compileError("diagnostic type requires an explicit JSON encoding: " ++ @typeName(T)),
    }
}

fn writeOpaque(json: *std.json.Stringify, bytes: []const u8, original: u128, lower_bound: bool, truncated: bool) std.Io.Writer.Error!void {
    std.debug.assert(bytes.len <= diagnostics.PREVIEW_BYTES);
    const utf8 = std.unicode.utf8ValidateSlice(bytes);
    try json.beginObject();
    try field(json, "encoding", if (utf8) "utf8" else "base64");
    try json.objectField("value");
    if (utf8) {
        try json.write(bytes);
    } else {
        var buffer: [std.base64.standard.Encoder.calcSize(diagnostics.PREVIEW_BYTES)]u8 = undefined;
        try json.write(std.base64.standard.Encoder.encode(&buffer, bytes));
    }
    try field(json, "truncated", truncated);
    try json.objectField(if (lower_bound) "observed_bytes_min" else "original_bytes");
    var number: [39]u8 = undefined;
    try json.write(std.fmt.bufPrint(&number, "{d}", .{original}) catch unreachable);
    try json.endObject();
}

fn jsonReference(json: *std.json.Stringify, failure: *const diagnostics.Failure, reference: diagnostics.Reference) std.Io.Writer.Error!void {
    try json.beginObject();
    try jsonStructFields(json, failure, reference);
    try json.endObject();
}

fn jsonDiagnostic(json: *std.json.Stringify, failure: *const diagnostics.Failure, diagnostic: *const diagnostics.Diagnostic) std.Io.Writer.Error!void {
    const metadata = diagnostics.metadata(diagnostic.code());
    try json.beginObject();
    try field(json, "code", metadata.name);
    try field(json, "category", @tagName(metadata.category));
    try field(json, "severity", @tagName(metadata.severity));
    try field(json, "kind", @tagName(metadata.kind));
    try field(json, "unit", @tagName(metadata.unit));
    try field(json, "stage", @tagName(diagnostic.stage));
    try jsonStructFields(json, failure, diagnostic.at);
    try json.objectField("related");
    try json.beginArray();
    if (diagnostic.related) |related| try jsonReference(json, failure, related);
    try json.endArray();
    try jsonField(json, failure, "cause", diagnostic.cause);
    try json.objectField("data");
    switch (diagnostic.details) {
        inline else => |facts| try jsonValue(json, failure, facts),
    }
    try json.endObject();
}

const Label = struct { bytes: []const u8, original: u128, lower_bound: bool, truncated: bool };

fn sourcePrefix(failure: *const diagnostics.Failure, id: diagnostics.SourceId, buffer: *[SOURCE_PREFIX_BYTES]u8) Label {
    const source = failure.sources.get(id).?;
    const bytes = failure.sourceLabel(id).?;
    const original: u128 = switch (source.label) {
        .borrowed => bytes.len,
        .retained => |preview| preview.original_bytes,
    };
    const lower_bound = source.label == .retained and source.label.retained.length_kind == .at_least;
    const suffix = switch (source.suffix) {
        .none => "",
        .zfi => ".zfi",
        .fai => ".fai",
    };
    var count: usize = @min(bytes.len, buffer.len);
    @memcpy(buffer[0..count], bytes[0..count]);
    // A suffix on an incomplete base would no longer be a true source prefix.
    if (!lower_bound and count == original) {
        const appended = @min(suffix.len, buffer.len - count);
        @memcpy(buffer[count..][0..appended], suffix[0..appended]);
        count += appended;
    }
    const total = original + suffix.len;
    return .{ .bytes = buffer[0..count], .original = total, .lower_bound = lower_bound, .truncated = lower_bound or count != total };
}

fn referenceMask(reference: diagnostics.Reference) u16 {
    return if (reference.source) |id| @as(u16, 1) << @intFromEnum(id) else 0;
}

fn diagnosticMask(diagnostic: *const diagnostics.Diagnostic) u16 {
    var mask = referenceMask(diagnostic.at);
    if (diagnostic.related) |related| mask |= referenceMask(related);
    switch (diagnostic.details) {
        .@"index.not_found" => |facts| mask |= (@as(u16, 1) << @intFromEnum(facts.zfi_source)) | (@as(u16, 1) << @intFromEnum(facts.fai_source)),
        .@"request.sequence_not_found" => |facts| mask |= @as(u16, 1) << @intFromEnum(facts.index_source),
        .@"io.cleanup_failed" => |facts| mask |= @as(u16, 1) << @intFromEnum(facts.resource_role),
        .@"repair.output_is_input" => |facts| mask |= (@as(u16, 1) << @intFromEnum(facts.input)) | (@as(u16, 1) << @intFromEnum(facts.output)),
        .@"repair.output_identity_unknown" => |facts| mask |= (@as(u16, 1) << @intFromEnum(facts.input)) | (@as(u16, 1) << @intFromEnum(facts.output)),
        else => {},
    }
    if (actionSource(diagnostic)) |id| mask |= @as(u16, 1) << @intFromEnum(id);
    return mask;
}

fn sourceMask(failure: *const diagnostics.Failure) u16 {
    var mask = failure.sources.registered;
    if (failure.primary()) |primary| mask |= diagnosticMask(primary);
    for (failure.secondaryItems()) |secondary| mask |= diagnosticMask(&secondary);
    return mask;
}

fn actionSource(diagnostic: *const diagnostics.Diagnostic) ?diagnostics.SourceId {
    return switch (diagnostics.metadata(diagnostic.code()).action) {
        .rebuild_index => .fasta,
        else => diagnostic.at.source,
    };
}

/// Writes printable ASCII verbatim and escapes controls, opaque bytes, quotes
/// and backslashes without adding enclosing quotes. Does not allocate.
pub fn writeEscapedBytes(writer: *std.Io.Writer, bytes: []const u8) std.Io.Writer.Error!void {
    for (bytes) |byte| switch (byte) {
        '"', '\\' => {
            try writer.writeByte('\\');
            try writer.writeByte(byte);
        },
        0x20...0x21, 0x23...0x5b, 0x5d...0x7e => try writer.writeByte(byte),
        else => try writer.print("\\x{x:0>2}", .{byte}),
    };
}

fn escaped(writer: *std.Io.Writer, bytes: []const u8) std.Io.Writer.Error!void {
    try writer.writeByte('"');
    try writeEscapedBytes(writer, bytes);
    try writer.writeByte('"');
}

fn textOpaque(writer: *std.Io.Writer, bytes: []const u8, original: u128, lower_bound: bool, truncated: bool) std.Io.Writer.Error!void {
    try escaped(writer, bytes);
    if (truncated) try writer.print(" ({s}{d} bytes)", .{ if (lower_bound) "at least " else "", original });
}

fn causeTextPrefix(bytes: []const u8) []const u8 {
    var escaped_bytes: usize = 0;
    for (bytes[0..@min(bytes.len, CAUSE_PREFIX_BYTES)], 0..) |byte, index| {
        const cost: usize = switch (byte) {
            '"', '\\' => 2,
            0x20...0x21, 0x23...0x5b, 0x5d...0x7e => 1,
            else => 4,
        };
        if (escaped_bytes + cost > CAUSE_PREFIX_BYTES) return bytes[0..index];
        escaped_bytes += cost;
    }
    return bytes[0..@min(bytes.len, CAUSE_PREFIX_BYTES)];
}

fn textValue(writer: *std.Io.Writer, failure: *const diagnostics.Failure, value: anytype) std.Io.Writer.Error!void {
    const T = @TypeOf(value);
    if (T == diagnostics.Preview) {
        const bytes = failure.preview(value);
        const count = @min(bytes.len, 32);
        try textOpaque(writer, bytes[0..count], value.original_bytes, value.length_kind == .at_least, value.truncated or count != bytes.len);
    } else if (T == diagnostics.Cause) {
        switch (value) {
            .zig => |err| {
                const name = @errorName(err);
                const prefix = causeTextPrefix(name);
                try textOpaque(writer, prefix, name.len, false, prefix.len != name.len);
            },
            .native => |native| try writer.print("{s}:{d}", .{ @tagName(native.domain), native.number }),
        }
    } else if (T == diagnostics.Arguments) {
        var separator: []const u8 = "";
        for (std.enums.values(diagnostics.Argument)) |argument| {
            if (value.contains(argument)) {
                try writer.print("{s}{s}", .{ separator, @tagName(argument) });
                separator = ",";
            }
        }
    } else if (T == diagnostics.TempAttempts) {
        for (value.items()) |attempt| {
            try writer.writeByte('(');
            try textValue(writer, failure, attempt);
            try writer.writeByte(')');
        }
        if (value.truncated) try writer.writeAll(" (further attempts omitted)");
    } else if (T == diagnostics.Measurement) {
        switch (value) {
            .exact => |number| try writer.print("{d}", .{number}),
            .at_least => |number| try writer.print("at least {d}", .{number}),
        }
    } else if (T == diagnostics.RepairBlock) {
        try writer.print("reason_code={s} ", .{@tagName(value.blocker.code())});
        switch (value.blocker) {
            inline else => |facts| try textValue(writer, failure, facts),
        }
        try writer.print(" repair_policy={s}", .{@tagName(value.repair_policy)});
    } else switch (@typeInfo(T)) {
        .void => {},
        .int => try writer.print("{d}", .{value}),
        .bool => try writer.writeAll(if (value) "true" else "false"),
        .@"enum" => try writer.writeAll(@tagName(value)),
        .optional => if (value) |present| try textValue(writer, failure, present) else try writer.writeAll("unknown"),
        .@"union" => switch (value) {
            inline else => |facts, tag| {
                try writer.print("{s} ", .{@tagName(tag)});
                try textValue(writer, failure, facts);
            },
        },
        .@"struct" => {
            var separator: []const u8 = "";
            inline for (std.meta.fields(T)) |member| {
                if (!std.mem.eql(u8, member.name, "reserved")) {
                    try writer.print("{s}{s}=", .{ separator, member.name });
                    if (comptime std.mem.eql(u8, member.name, "byte") or std.mem.eql(u8, member.name, "delimiter") or std.mem.eql(u8, member.name, "flags") or std.mem.eql(u8, member.name, "required_flag")) {
                        try writer.print("0x{x:0>2}", .{@field(value, member.name)});
                    } else try textValue(writer, failure, @field(value, member.name));
                    separator = " ";
                }
            }
        },
        .array => try escaped(writer, &value),
        else => @compileError("diagnostic type requires an explicit text encoding"),
    }
}

fn textReference(writer: *std.Io.Writer, failure: *const diagnostics.Failure, reference: diagnostics.Reference) std.Io.Writer.Error!void {
    if (reference.source) |id| {
        if (failure.sources.get(id) != null) {
            var buffer: [SOURCE_PREFIX_BYTES]u8 = undefined;
            const label = sourcePrefix(failure, id, &buffer);
            const count = @min(label.bytes.len, 24);
            try textOpaque(writer, label.bytes[0..count], label.original, label.lower_bound, label.truncated or count != label.bytes.len);
        } else try writer.writeAll(@tagName(id));
    }
    if (reference.location.line) |line| try writer.print(":{d}", .{line});
    if (reference.location.column) |column| try writer.print(" column={d}", .{column});
    if (reference.location.offset) |offset| try writer.print(" offset={d}", .{offset});
    if (reference.record) |record| try writer.print(" record={d}", .{record});
    if (reference.request) |request| try writer.print(" request={d}", .{request});
}

fn textDiagnostic(writer: *std.Io.Writer, failure: *const diagnostics.Failure, diagnostic: *const diagnostics.Diagnostic, prefix: []const u8) std.Io.Writer.Error!void {
    try writer.print("{s} [{s}] {s}: ", .{ prefix, @tagName(diagnostic.code()), @tagName(failure.operation) });
    if (referenceMask(diagnostic.at) != 0 or diagnostic.at.location.offset != null or diagnostic.at.location.line != null or diagnostic.at.location.column != null or diagnostic.at.record != null or diagnostic.at.request != null) {
        try textReference(writer, failure, diagnostic.at);
        try writer.writeAll(": ");
    }
    try writer.writeAll(message(diagnostic.code()));
    switch (diagnostic.details) {
        inline else => |facts| {
            if (@TypeOf(facts) != void) {
                try writer.writeByte(' ');
                try textValue(writer, failure, facts);
            }
        },
    }
    if (diagnostic.related) |related| {
        try writer.writeAll("; related: ");
        try textReference(writer, failure, related);
    }
    if (diagnostic.cause) |cause| {
        try writer.writeAll("; cause=");
        try textValue(writer, failure, cause);
    }
    try writer.writeByte('\n');
}

fn actionText(action: diagnostics.Action) []const u8 {
    return switch (action) {
        .correct_arguments => "Correct the command arguments; use --help for accepted options.",
        .inspect_io => "Check the named input/output and retained I/O cause.",
        .check_memory => "Check available memory for the named allocation.",
        .check_capacity => "Reduce the requested size or use a host with sufficient capacity.",
        .check_counts => "Check the reported counter and its numeric limit.",
        .check_source => "Confirm the intended input has finished changing before retrying.",
        .correct_data => "Inspect the reported data and the selected rule.",
        .select_index_format => "Use ZFI for sequence layout that FAI cannot represent.",
        .rebuild_index => "Rebuild the index after confirming the intended FASTA input.",
        .correct_request => "Correct the reported request or select the intended index.",
        .select_transform => "Select a transform appropriate for the sampled sequence type.",
        .select_output => "Choose a checked output path distinct from the input.",
        .resolve_blocker => "Resolve the reported rewrite blocker before publishing output.",
        .report_defect => "Report the named internal check and diagnostic facts.",
    };
}

fn message(code: diagnostics.Code) []const u8 {
    return switch (code) {
        .@"usage.command_required" => "A command is required.",
        .@"usage.unknown_command" => "Unknown command.",
        .@"usage.unknown_option" => "Unknown option.",
        .@"usage.missing_argument" => "A required argument is missing.",
        .@"usage.extra_argument" => "Too many arguments.",
        .@"usage.missing_option_value" => "A required option value is missing.",
        .@"usage.repeated_option" => "This option cannot be repeated.",
        .@"usage.conflicting_options" => "Selected options conflict.",
        .@"usage.option_requires" => "An option prerequisite is missing.",
        .@"usage.invalid_option_value" => "Invalid option value.",
        .@"io.file_not_found" => "Required file was not found.",
        .@"io.access_denied" => "Access was denied.",
        .@"io.invalid_path" => "The selected path is invalid.",
        .@"io.open_failed" => "Could not open the required input.",
        .@"io.stat_failed" => "Could not read required file metadata.",
        .@"io.seek_failed" => "Could not seek to the requested offset.",
        .@"io.read_failed" => "Could not read the required input.",
        .@"io.unexpected_eof" => "Required input ended early.",
        .@"io.map_failed" => "Could not map the required input.",
        .@"io.create_failed" => "Could not create the output.",
        .@"io.write_failed" => "Could not write the output.",
        .@"io.flush_failed" => "Could not flush the output.",
        .@"io.publish_failed" => "Could not publish the completed output file.",
        .@"io.temp_unavailable" => "No usable temporary file is available.",
        .@"io.cleanup_failed" => "Could not remove the temporary resource.",
        .@"resource.out_of_memory" => "Out of memory while allocating required storage.",
        .@"resource.capacity_limit_exceeded" => "Requested storage exceeds supported capacity.",
        .@"resource.count_overflow" => "An exact counter exceeds its numeric limit.",
        .@"resource.line_buffer_exceeded" => "Input line exceeds the reader limit.",
        .@"input.changed" => "Input changed during the operation.",
        .@"fasta.no_records" => "No FASTA records were found.",
        .@"fasta.content_before_header" => "Content appears before the first FASTA header.",
        .@"fasta.header_control_byte" => "Header contains a prohibited control byte.",
        .@"fasta.sequence_control_byte" => "Sequence contains a prohibited control byte.",
        .@"fasta.duplicate_identifier" => "Identifier duplicates an earlier record.",
        .@"sequence.invalid_symbol" => "Sequence byte is outside the selected alphabet.",
        .@"index.initial_header_required" => "Indexing requires a FASTA header at byte zero.",
        .@"index.no_indexable_records" => "No nonempty records remain under the selected naming policy.",
        .@"index.identifier_too_long" => "Identifier exceeds the selected index limit.",
        .@"index.sequence_line_too_long" => "Sequence line exceeds its width limit.",
        .@"index.record_count_limit_exceeded" => "Record count exceeds the selected index limit.",
        .@"index.side_table_offset_limit_exceeded" => "Side-table offset exceeds its encoded limit.",
        .@"index.source_timestamp_out_of_range" => "Source timestamp is outside the ZFI range.",
        .@"index.nonuniform_fai" => "Record has a nonuniform layout that FAI cannot represent.",
        .@"index.not_found" => "No index was found for the FASTA input.",
        .@"index.no_records" => "Index contains no records.",
        .@"index.invalid_magic" => "Index signature does not match the expected magic.",
        .@"index.unsupported_version" => "ZFI version is unsupported.",
        .@"index.truncated" => "A required index section extends beyond available bytes.",
        .@"index.name_not_embedded" => "ZFI record does not declare an embedded name.",
        .@"index.name_footer_missing" => "Index has no accepted name-table footer.",
        .@"index.name_blob_out_of_bounds" => "Name blob extends outside its permitted region.",
        .@"index.name_out_of_bounds" => "Record name lies outside the index name blob.",
        .@"index.zero_sequence_length" => "Indexed record has zero sequence length.",
        .@"index.invalid_line_geometry" => "Uniform line geometry is invalid.",
        .@"index.sequence_offset_out_of_bounds" => "Sequence offset is outside the source byte bound.",
        .@"index.sequence_span_out_of_bounds" => "Indexed sequence extends beyond the source byte bound.",
        .@"index.arithmetic_overflow" => "Index arithmetic exceeds the wire numeric range.",
        .@"index.fai.field_missing" => "A required FAI field is missing or empty.",
        .@"index.fai.invalid_integer" => "FAI field must contain unsigned decimal digits.",
        .@"index.fai.integer_out_of_range" => "FAI integer exceeds its accepted range.",
        .@"index.side_table_offset_out_of_bounds" => "Side-table offset is outside its permitted region.",
        .@"index.side_table_misaligned" => "Side-table offset is misaligned.",
        .@"index.side_table_empty" => "Nonuniform record has an empty side table.",
        .@"index.side_table_order_invalid" => "Side table is out of order or overlaps its predecessor.",
        .@"index.side_base_start_mismatch" => "Side-table base prefix differs from the expected value.",
        .@"index.side_line_geometry_invalid" => "Side-table entry has invalid line geometry.",
        .@"index.side_first_offset_mismatch" => "First side entry differs from the sequence offset.",
        .@"index.side_line_order_invalid" => "Side-table line is out of order or overlaps its predecessor.",
        .@"index.side_sequence_length_mismatch" => "Side-table bases differ from the indexed length.",
        .@"index.older_than_source" => "Index modification time is earlier than the FASTA time.",
        .@"index.source_size_mismatch" => "Index source size differs from the FASTA size.",
        .@"index.source_mtime_mismatch" => "Index source timestamp differs from the FASTA time.",
        .@"pair.region_length_mismatch" => "Selected region yielded fewer symbols than requested.",
        .@"pair.total_length_mismatch" => "Observed composition total differs from the indexed total.",
        .@"request.no_regions" => "No executable requests were found.",
        .@"request.sequence_not_found" => "Sequence was not found in the selected lookup.",
        .@"request.start_before_one" => "Inclusive request start must be at least one.",
        .@"request.start_after_sequence" => "Request start exceeds the indexed sequence length.",
        .@"request.end_before_start" => "Request end is before its start.",
        .@"request.identifier_too_long" => "Request identifier exceeds its byte limit.",
        .@"request.bed.field_missing" => "A required BED field is missing or empty.",
        .@"request.bed.invalid_coordinate" => "BED coordinate must contain unsigned decimal digits.",
        .@"request.bed.coordinate_out_of_range" => "BED coordinate exceeds its accepted range.",
        .@"request.bed.invalid_interval" => "BED end must be greater than start.",
        .@"request.bed.invalid_strand" => "Selected BED strand is invalid.",
        .@"transform.complement_on_protein" => "Cannot complement the sampled protein sequence.",
        .@"repair.output_is_input" => "Rewrite output resolves to the input.",
        .@"repair.output_identity_unknown" => "Could not establish a distinct rewrite output.",
        .@"repair.blocked" => "Rewrite is blocked by the reported data condition.",
        .@"internal.consistency_failed" => "An internal consistency check failed.",
    };
}
