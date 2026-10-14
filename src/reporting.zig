//! Renders a `Failure` as text or as one NDJSON object, and a validation report inside that
//! object or as text. Never allocates or flushes; a failed write is returned and the original
//! failure is left unchanged.

const std = @import("std");
const diagnostics = @import("diagnostics");

pub const TEXT_BYTES_LIMIT = 4096;
pub const JSON_BYTES_LIMIT = 16384;

pub const GetMetrics = struct { region_count: u64, total_bases: u64, elapsed_ns: u64 };

/// `version` must be trusted UTF-8 of at most 32 bytes.
pub const Options = struct {
    version: ?[]const u8 = null,
    accepted: bool = false,
    metrics: ?GetMetrics = null,
    validation: ?*const Validation = null,
};

/// What `validate` found. Counts are exact; each code keeps only its first occurrence.
pub const Validation = struct {
    completion: Completion,
    strict: bool,
    alphabet: diagnostics.AlphabetPolicy,
    max_header_bytes: u64,
    convention: ?diagnostics.Convention,
    examined: Examined,
    records_with_errors: u64,
    records_with_warnings: u64,
    /// In report order; a code that did not occur has no occurrences and no first occurrence.
    findings: []const Finding,
    /// Occurrences of each byte value outside the selected alphabet.
    invalid_symbols: *const [256]u64,
    /// Set when the duplicate check reached its memory budget and stopped.
    name_budget_stop: ?BudgetStop = null,
    /// Set when the FASTA was compared with its index.
    index: ?IndexComparison = null,

    pub const Completion = enum { complete, incomplete, interrupted };
    pub const Outcome = enum { clean, warnings, invalid, incomplete };
    /// `failed` when the check found an error, `warnings` when it found only warnings.
    pub const CheckState = enum { passed, warnings, failed, incomplete };
    /// The four checks every report has, by the category of the codes each reports; `index_match`
    /// joins them when the FASTA was compared with its index.
    const CHECKS = [_]diagnostics.Category{ .structure, .identifiers, .alphabet, .format };
    /// `records` counts headers seen; `records_completed`, records whose end was reached.
    pub const Examined = struct { bytes: u64 = 0, lines: u64 = 0, records: u64 = 0, records_completed: u64 = 0, symbols: u64 = 0 };
    const Totals = struct { errors: u64 = 0, warnings: u64 = 0 };
    pub const BudgetStop = struct { budget_bytes: u64, stopped_at_record: u64 };
    /// `omitted_duplicates`: records the index leaves out because an earlier record has the name.
    pub const IndexComparison = struct { format: diagnostics.Format, entries: u64, matched: u64, omitted_duplicates: u64 };

    pub fn codes(self: *const Validation, severity: diagnostics.Severity) u64 {
        var count: u64 = 0;
        for (self.findings) |finding| {
            if (finding.occurrences != 0 and diagnostics.metadata(finding.code).severity == severity) count += 1;
        }
        return count;
    }

    /// Occurrences by severity, over every code or over one category.
    pub fn totals(self: *const Validation, category: ?diagnostics.Category) Totals {
        var result: Totals = .{};
        for (self.findings) |finding| {
            const metadata = diagnostics.metadata(finding.code);
            if (category != null and metadata.category != category.?) continue;
            switch (metadata.severity) {
                .@"error" => result.errors += finding.occurrences,
                .warning => result.warnings += finding.occurrences,
            }
        }
        return result;
    }

    pub fn outcome(self: *const Validation) Outcome {
        if (self.completion != .complete) return .incomplete;
        const all = self.totals(null);
        if (all.errors != 0) return .invalid;
        return if (all.warnings != 0) .warnings else .clean;
    }

    pub fn check(self: *const Validation, category: diagnostics.Category) CheckState {
        if (self.completion == .interrupted) return .incomplete;
        if (category == .identifiers and self.name_budget_stop != null) return .incomplete;
        const found = self.totals(category);
        if (found.errors != 0) return .failed;
        return if (found.warnings != 0) .warnings else .passed;
    }
};

pub const Finding = struct {
    code: diagnostics.Code,
    occurrences: u64 = 0,
    /// Distinct records with an occurrence; file-level and preamble occurrences belong to none.
    records: u64 = 0,
    first: ?Occurrence = null,
};

pub const Occurrence = struct {
    details: diagnostics.Details,
    at: diagnostics.Reference,
    related: ?diagnostics.Reference = null,
    identifier: ?Identifier = null,
};

/// A copied identifier prefix and the identifier's full length.
pub const Identifier = struct {
    bytes: [diagnostics.PREVIEW_BYTES]u8,
    len: u16,
    original_bytes: u64,

    pub fn init(identifier: []const u8) Identifier {
        var value: Identifier = .{ .bytes = undefined, .len = @intCast(@min(identifier.len, diagnostics.PREVIEW_BYTES)), .original_bytes = identifier.len };
        @memcpy(value.bytes[0..value.len], identifier[0..value.len]);
        return value;
    }

    pub fn slice(self: *const Identifier) []const u8 {
        return self.bytes[0..self.len];
    }
};

pub const Bootstrap = enum { argv_oom, missing_command, unknown_command, invalid_format };

/// For failures before a `Failure` exists; do not retry on the same broken sink.
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

pub fn writeJson(writer: *std.Io.Writer, failure: *const diagnostics.Failure, options: Options) std.Io.Writer.Error!void {
    if (options.version) |version| {
        std.debug.assert(version.len <= 32 and std.unicode.utf8ValidateSlice(version));
    }
    var json: std.json.Stringify = .{ .writer = writer };
    try json.beginObject();
    try field(&json, "catalog_revision", "0.4.0-1");
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
    try json.objectField("validation");
    if (options.validation) |validation| try jsonValidation(&json, failure, validation) else try json.write(null);
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

pub fn writeText(writer: *std.Io.Writer, failure: *const diagnostics.Failure) std.Io.Writer.Error!void {
    const primary = failure.primary() orelse return;
    try textDiagnostic(writer, failure, primary, "error");
    for (failure.secondaryItems()) |secondary| try textDiagnostic(writer, failure, &secondary, "secondary");
    if (failure.secondary_truncated) try writer.writeAll("Further secondary failures were omitted.\n");
    try writer.print("hint: {s}\n", .{actionText(diagnostics.metadata(primary.code()).action)});
}

pub fn writeValidationText(writer: *std.Io.Writer, failure: *const diagnostics.Failure, validation: *const Validation) std.Io.Writer.Error!void {
    const all = validation.totals(null);
    try writer.writeAll(switch (validation.outcome()) {
        .incomplete => if (validation.completion == .interrupted) "INTERRUPTED " else "INCOMPLETE ",
        .invalid => "INVALID ",
        .warnings => "WARNINGS ",
        .clean => "OK ",
    });
    try writeEscapedBytes(writer, failure.sourceLabel(.fasta) orelse "input");
    if (all.errors + all.warnings == 0) {
        try writer.writeAll(": no findings\n");
    } else {
        const error_codes = validation.codes(.@"error");
        const warning_codes = validation.codes(.warning);
        try writer.print(": {d} error{s} in {d} code{s}, {d} warning{s} in {d} code{s}; {d} of {d} record{s} with errors, {d} with warnings\n", .{ all.errors, plural(all.errors), error_codes, plural(error_codes), all.warnings, plural(all.warnings), warning_codes, plural(warning_codes), validation.records_with_errors, validation.examined.records, plural(validation.examined.records), validation.records_with_warnings });
    }
    for (validation.findings) |finding| {
        if (finding.occurrences != 0) try textFinding(writer, failure, finding);
    }
    const examined = validation.examined;
    try writer.print("{s}: examined {d} byte{s}, {d} line{s}, {d} record{s}", .{ @tagName(validation.completion), examined.bytes, plural(examined.bytes), examined.lines, plural(examined.lines), examined.records, plural(examined.records) });
    if (examined.records_completed != examined.records) try writer.print(" ({d} completed)", .{examined.records_completed});
    try writer.print(", {d} symbol{s}; ", .{ examined.symbols, plural(examined.symbols) });
    if (validation.alphabet.sample) |sample| {
        try writer.print("{s} alphabet chosen from {d} symbol{s}", .{ @tagName(validation.alphabet.alphabet), sample.symbols, plural(sample.symbols) });
    } else try writer.print("{s} alphabet", .{@tagName(validation.alphabet.alphabet)});
    if (validation.name_budget_stop) |stop| try writer.print("; duplicate check stopped at record {d} at its {d}-byte name budget", .{ stop.stopped_at_record, stop.budget_bytes });
    if (validation.index) |index| {
        try writer.print("; {s} index: {d} of {d} entr{s} matched", .{ @tagName(index.format), index.matched, index.entries, if (index.entries == 1) "y" else "ies" });
        if (index.omitted_duplicates != 0) try writer.print(", {d} later duplicate{s} not indexed", .{ index.omitted_duplicates, plural(index.omitted_duplicates) });
    }
    try writer.writeByte('\n');
}

/// Printable ASCII as is; controls, non-ASCII bytes, quotes, and backslashes escaped.
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

const SOURCE_PREFIX_BYTES = 48;
const CAUSE_PREFIX_BYTES = 32;

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
                try jsonNamedValue(json, failure, name, number);
            },
            .at_least => |number| {
                try jsonNamedValue(json, failure, name ++ "_min", number);
            },
        }
    } else {
        try jsonNamedValue(json, failure, name, value);
    }
}

fn jsonNamedValue(json: *std.json.Stringify, failure: *const diagnostics.Failure, name: []const u8, value: anytype) std.Io.Writer.Error!void {
    try json.objectField(name);
    try jsonValue(json, failure, value);
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

fn jsonValidation(json: *std.json.Stringify, failure: *const diagnostics.Failure, validation: *const Validation) std.Io.Writer.Error!void {
    try json.beginObject();
    try field(json, "outcome", @tagName(validation.outcome()));
    try field(json, "completion", @tagName(validation.completion));
    try json.objectField("policy");
    try json.beginObject();
    try field(json, "strict", validation.strict);
    try jsonNamedValue(json, failure, "alphabet", validation.alphabet);
    try jsonNamedValue(json, failure, "max_header_bytes", validation.max_header_bytes);
    try jsonNamedValue(json, failure, "header_convention", validation.convention);
    try json.endObject();
    try json.objectField("checks");
    try json.beginObject();
    inline for (Validation.CHECKS) |category| try field(json, @tagName(category), @tagName(validation.check(category)));
    if (validation.index != null) try field(json, "index_match", @tagName(validation.check(.index_match)));
    try json.endObject();
    try jsonNamedValue(json, failure, "examined", validation.examined);
    const all = validation.totals(null);
    try json.objectField("findings");
    try json.beginObject();
    try jsonNamedValue(json, failure, "errors", all.errors);
    try jsonNamedValue(json, failure, "warnings", all.warnings);
    try jsonNamedValue(json, failure, "error_codes", validation.codes(.@"error"));
    try jsonNamedValue(json, failure, "warning_codes", validation.codes(.warning));
    try json.endObject();
    try json.objectField("categories");
    try json.beginObject();
    inline for (Validation.CHECKS) |category| try jsonNamedValue(json, failure, @tagName(category), validation.totals(category));
    if (validation.index != null) try jsonNamedValue(json, failure, "index_match", validation.totals(.index_match));
    try json.endObject();
    try jsonNamedValue(json, failure, "name_budget", validation.name_budget_stop);
    try jsonNamedValue(json, failure, "index", validation.index);
    try jsonNamedValue(json, failure, "records_with_errors", validation.records_with_errors);
    try jsonNamedValue(json, failure, "records_with_warnings", validation.records_with_warnings);
    try json.objectField("codes");
    try json.beginArray();
    for (validation.findings) |finding| {
        if (finding.occurrences == 0) continue;
        const metadata = diagnostics.metadata(finding.code);
        try json.beginObject();
        try field(json, "code", metadata.name);
        try field(json, "category", @tagName(metadata.category));
        try field(json, "severity", @tagName(metadata.severity));
        try field(json, "unit", @tagName(metadata.unit));
        try jsonNamedValue(json, failure, "occurrences", finding.occurrences);
        try jsonNamedValue(json, failure, "records", finding.records);
        try json.objectField("first");
        if (finding.first) |first| {
            try json.beginObject();
            try jsonStructFields(json, failure, first.at.location);
            try jsonNamedValue(json, failure, "record", first.at.record);
            try json.objectField("identifier");
            if (first.identifier) |identifier| {
                try writeOpaque(json, identifier.slice(), identifier.original_bytes, false, identifier.len != identifier.original_bytes);
            } else try json.write(null);
            try json.objectField("related");
            if (first.related) |related| {
                try json.beginObject();
                try jsonStructFields(json, failure, related.location);
                try jsonNamedValue(json, failure, "record", related.record);
                try json.endObject();
            } else try json.write(null);
            try json.objectField("data");
            try jsonDetails(json, failure, first.details);
            try json.endObject();
        } else try json.write(null);
        try json.endObject();
    }
    try json.endArray();
    try json.objectField("invalid_symbols");
    try json.beginArray();
    for (validation.invalid_symbols, 0..) |count, byte| {
        if (count == 0) continue;
        try json.beginObject();
        try field(json, "byte", @as(u8, @intCast(byte)));
        try jsonNamedValue(json, failure, "occurrences", count);
        try json.endObject();
    }
    try json.endArray();
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
    try jsonDetails(json, failure, diagnostic.details);
    try json.endObject();
}

fn jsonDetails(json: *std.json.Stringify, failure: *const diagnostics.Failure, details: diagnostics.Details) std.Io.Writer.Error!void {
    switch (details) {
        inline else => |facts| try jsonValue(json, failure, facts),
    }
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
    try textDetails(writer, failure, diagnostic.details);
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

fn plural(count: u64) []const u8 {
    return if (count == 1) "" else "s";
}

fn textLocation(writer: *std.Io.Writer, reference: diagnostics.Reference) std.Io.Writer.Error!void {
    if (reference.location.line) |line| try writer.print("line {d}", .{line});
    if (reference.location.column) |column| try writer.print(" column {d}", .{column});
}

fn textFinding(writer: *std.Io.Writer, failure: *const diagnostics.Failure, finding: Finding) std.Io.Writer.Error!void {
    const metadata = diagnostics.metadata(finding.code);
    try writer.print("{s} [{s}]", .{ @tagName(metadata.severity), metadata.name });
    if (metadata.unit != .file) {
        try writer.print(" {d} {s}{s}", .{ finding.occurrences, @tagName(metadata.unit), plural(finding.occurrences) });
        if (metadata.unit != .record and finding.records != 0) try writer.print(" in {d} record{s}", .{ finding.records, plural(finding.records) });
    }
    const first = finding.first orelse {
        try writer.print(": {s}\n", .{message(finding.code)});
        return;
    };
    if (first.at.location.line != null) {
        try writer.writeAll(if (metadata.unit == .file) " at " else "; first at ");
        try textLocation(writer, first.at);
    }
    if (first.at.source == .zfi or first.at.source == .fai) {
        try writer.print("; first at {s} entry {d}", .{ @tagName(first.at.source.?), first.at.record.? });
    } else if (first.at.record) |record| try writer.print(", record {d}", .{record});
    if (first.identifier) |identifier| {
        try writer.writeByte(' ');
        const shown = identifier.slice()[0..@min(identifier.len, 32)];
        try textOpaque(writer, shown, identifier.original_bytes, false, shown.len != identifier.original_bytes);
    }
    try writer.print(": {s}", .{message(finding.code)});
    switch (first.details) {
        .@"fasta.duplicate_identifier" => {},
        .@"sequence.invalid_symbol" => |facts| try writer.print(" byte=0x{x:0>2}", .{facts.byte}),
        else => try textDetails(writer, failure, first.details),
    }
    if (first.related) |related| {
        try writer.writeAll(" First defined at ");
        try textLocation(writer, related);
        try writer.writeByte('.');
    }
    try writer.writeByte('\n');
}

fn textDetails(writer: *std.Io.Writer, failure: *const diagnostics.Failure, details: diagnostics.Details) std.Io.Writer.Error!void {
    switch (details) {
        inline else => |facts| {
            if (@TypeOf(facts) != void) {
                try writer.writeByte(' ');
                try textValue(writer, failure, facts);
            }
        },
    }
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
        .@"input.unsupported_format" => "Input uses a compressed format that is not accepted here.",
        .@"fasta.no_records" => "No FASTA records were found.",
        .@"fasta.content_before_header" => "Content appears before the first FASTA header.",
        .@"fasta.header_control_byte" => "Header contains a prohibited control byte.",
        .@"fasta.sequence_control_byte" => "Sequence contains a prohibited control byte.",
        .@"fasta.empty_record" => "FASTA record contains no sequence symbols.",
        .@"fasta.empty_identifier" => "Header has an empty identifier.",
        .@"fasta.duplicate_identifier" => "Identifier duplicates an earlier record.",
        .@"fasta.identifier_too_long" => "Identifier is longer than an index can hold.",
        .@"fasta.long_header" => "Header exceeds the selected length limit.",
        .@"header.convention_mismatch" => "Header does not match the selected convention.",
        .@"sequence.invalid_symbol" => "Sequence byte is outside the selected alphabet.",
        .@"format.sequence_whitespace" => "Sequence line contains whitespace.",
        .@"format.blank_sequence_line" => "Blank line occurs inside a FASTA record.",
        .@"format.inconsistent_wrapping" => "Sequence line width differs from the record's wrap width.",
        .@"format.leading_bom" => "File starts with a UTF-8 BOM.",
        .@"format.missing_terminal_newline" => "Nonempty file has no terminal newline.",
        .@"format.mixed_line_endings" => "File mixes LF and CRLF line endings.",
        .@"index.initial_header_required" => "Indexing requires a FASTA header at byte zero.",
        .@"index.no_indexable_records" => "Every record is empty.",
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
        .@"index.name_footer_missing" => "Index does not end with its source identity and name-table footer.",
        .@"index.name_table_invalid" => "Index name-lookup table does not match its records.",
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
        .@"pair.record_missing_from_index" => "FASTA record has no index entry.",
        .@"pair.record_missing_from_fasta" => "Index entry has no FASTA record.",
        .@"pair.sequence_length_mismatch" => "Indexed sequence length differs from the FASTA.",
        .@"pair.sequence_offset_mismatch" => "Indexed sequence offset differs from the FASTA.",
        .@"pair.line_geometry_mismatch" => "Indexed line layout differs from the FASTA.",
        .@"pair.side_table_mismatch" => "Indexed side table differs from the FASTA lines.",
        .@"request.no_regions" => "No executable requests were found.",
        .@"request.sequence_not_found" => "Sequence was not found in the selected lookup.",
        .@"request.start_before_one" => "Inclusive request start must be at least one.",
        .@"request.start_after_sequence" => "Request start exceeds the indexed sequence length.",
        .@"request.end_before_start" => "Request end is before its start.",
        .@"request.end_after_sequence" => "Request end exceeds the indexed sequence length.",
        .@"request.ordinal_out_of_range" => "Record ordinal is outside the open index.",
        .@"request.output_buffer_too_small" => "Output requires more bytes than the destination holds.",
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
