//! Typed diagnostic facts and bounded caller-owned failure storage.
//! Fallible operations return Error!T, reset Failure once at entry, and retain
//! facts before returning. This module does not allocate, print, log, or exit.

const std = @import("std");

pub const SECONDARY_LIMIT = 2;
pub const PREVIEW_BYTES = 256;
pub const SOURCE_LIMIT = 16;
pub const TEMP_ATTEMPT_LIMIT = 4;
pub const FAILURE_BYTES_LIMIT = 2048;

/// Error tags control Zig propagation and errdefer. Inspect Failure for the
/// precise code and facts; multiple codes intentionally share an error tag.
pub const Error = error{
    InvalidArgument,
    Io,
    OutOfMemory,
    CapacityExceeded,
    CountOverflow,
    InputChanged,
    InvalidFasta,
    InvalidIndex,
    InvalidRequest,
    UnsupportedTransform,
    RepairRefused,
    Internal,
};

pub const Operation = enum { cli, index, get, stats, validate };
pub const Kind = enum { finding, failure };
pub const Severity = enum { @"error", warning };
pub const Unit = enum { none, file, record, line, byte, field };
pub const Category = enum {
    usage,
    io,
    resources,
    source_stability,
    structure,
    identifiers,
    alphabet,
    index_compatibility,
    index_integrity,
    index_match,
    requests,
    transforms,
    repair,
    internal,
};

/// Suggested action, independent of exit status, severity, and retry policy.
pub const Action = enum {
    correct_arguments,
    inspect_io,
    check_memory,
    check_capacity,
    check_counts,
    check_source,
    correct_data,
    select_index_format,
    rebuild_index,
    correct_request,
    select_transform,
    select_output,
    resolve_blocker,
    report_defect,
};

/// Stable roles within one operation; their labels need no registration heap.
pub const SourceId = enum(u4) {
    fasta,
    zfi,
    fai,
    names,
    bed,
    stdin,
    stdout,
    stderr,
    output,
    report,
    staged_output,
    fai_spool,
    name_spool,
    side_spool,
    catalog_spool,
    request_snapshot,
};
pub const SourceKind = enum { path, stream, memory };
pub const Suffix = enum { none, zfi, fai };

/// A borrowed label lives through handling/rendering or Failure.reset. A
/// retained label uses the owning Failure's preview pool. Derived paths use
/// the original caller label plus a static suffix, without a stack-built path.
pub const Source = struct {
    label: union(enum) { borrowed: []const u8, retained: Preview },
    kind: SourceKind,
    suffix: Suffix = .none,
};

/// Fixed role table. Unregistered roles return null, not an invented path.
pub const Sources = struct {
    entries: [SOURCE_LIMIT]Source = undefined,
    registered: u16 = 0,

    /// Register each role once after operation reset, before fallible work.
    pub fn set(self: *Sources, id: SourceId, source: Source) void {
        const bit = @as(u16, 1) << @intFromEnum(id);
        std.debug.assert(self.registered & bit == 0);
        self.entries[@intFromEnum(id)] = source;
        self.registered |= bit;
    }

    /// Borrowed bytes expire with their caller. Retained references expire at
    /// the owning Failure's reset; resolve those through Failure.sourceLabel.
    pub fn get(self: *const Sources, id: SourceId) ?Source {
        const bit = @as(u16, 1) << @intFromEnum(id);
        if (self.registered & bit == 0) return null;
        return self.entries[@intFromEnum(id)];
    }
};

/// Byte offset is zero-based. Text line and byte column are one-based.
/// null means unobserved; zero is a valid known byte offset.
pub const Location = struct {
    offset: ?u64 = null,
    line: ?u64 = null,
    column: ?u64 = null,
};

/// Record/request numbers are original one-based ordinals, not retained slots.
/// Binary fields can have only an offset. Unknown values stay null.
pub const Reference = struct {
    source: ?SourceId = null,
    location: Location = .{},
    record: ?u64 = null,
    request: ?u64 = null,
};

/// A fixed-pool reference, not a slice into Failure. Copying Failure preserves
/// this reference. Views borrow that Failure until its next reset or overwrite.
pub const Preview = struct {
    original_bytes: u64,
    start: u16,
    len: u16,
    truncated: bool,
    length_kind: std.meta.Tag(ByteLength),
};

/// Byte length of the original value, independent of the retained prefix.
pub const ByteLength = union(enum(u8)) { exact: u64, at_least: u64 };

/// An observed bound must not be presented as an exact measurement.
pub const Measurement = union(enum) { exact: u128, at_least: u128 };
pub const NativeError = struct {
    domain: enum { posix_errno, win32, ntstatus },
    number: u32,
};

/// Retained backend evidence; the Zig interface tag may itself be coarse.
/// null at the owning Diagnostic means that no cause was retained.
pub const Cause = union(enum) { zig: anyerror, native: NativeError };

pub const Stage = enum {
    arguments,
    open,
    stat,
    seek,
    read,
    map,
    parse,
    allocate,
    scan,
    resolve,
    sample,
    transform,
    create,
    write,
    flush,
    publish,
    cleanup,
    report,
    compare_source,
    check_identity,
};
pub const OutputState = struct {
    state: enum { not_started, staged, published, streaming, unknown, aborted } = .not_started,
    /// Records completed by the writer, not a promise about the visible prefix.
    records: ?u64 = null,
    /// Bytes accepted by the writer; buffering does not prove delivery.
    logical_bytes: ?u64 = null,
    final_flush: enum { not_attempted, succeeded, failed } = .not_attempted,
    may_be_partial: bool = false,
};

/// stdout, a requested destination, and stderr can reach different states.
/// `file` also tracks caller-supplied writers whose path/kind is unobserved.
pub const Outputs = struct {
    stdout: OutputState = .{},
    file: OutputState = .{},
    stderr: OutputState = .{},
};

pub const Argument = enum {
    fasta_path,
    request_source,
    positional,
    emit_fai,
    no_dedup,
    bed,
    names,
    strand_aware,
    rc,
    complement_only,
    reverse_only,
    annotate_rc,
    summary,
    strict,
    json,
    fix,
    output,
    fix_format_only,
    schema,
    custom_alphabet,
    max_header_len,
    error_format,
    read_buffer,
};
pub const Arguments = std.EnumSet(Argument);
pub const Expected = enum { path, request_source, schema, alphabet, positive_usize, error_format, nonempty_buffer };
pub const Requirement = enum { bed, transform_or_strand, json, fix, output };
pub const AllocationRole = enum {
    argv,
    arena,
    records,
    names,
    name_slices,
    name_map,
    dedup_table,
    side_tables,
    retained_names,
    name_cache,
    catalog,
    compressed_names,
    spool_buffer,
    reader_buffer,
    request_buffer,
    requests,
    resolved_requests,
    request_names,
    request_collection,
    lengths,
    traversal_order,
    validation_identifiers,
    validation_events,
    line_widths,
    rewrite_bytes,
    serialized_index,
    serialized_event,
};
pub const CapacityRole = enum {
    input_bytes,
    map_bytes,
    records,
    record_bytes,
    name_blob,
    name_map,
    name_offset,
    name_cache_offset,
    side_offset,
    positional_requests,
    request_names,
    request_reservation,
    disk_span,
    reader_buffer,
};
pub const Counter = enum {
    source_bytes,
    sequence_symbols,
    total_symbols,
    aun_numerator,
    physical_lines,
    records,
    requests,
    events,
    composition_bin,
    lowercase_symbols,
    report_ordinal,
};
pub const PathReason = enum { rejected, too_long, directory };
pub const TemporaryPurpose = enum { fai, names, side_tables, catalog };
pub const TemporaryPlace = enum { output_directory, environment, system, current_directory };
pub const TemporaryReason = enum { randomness, directory_open, creation, collisions };
pub const TempAttempt = struct {
    place: TemporaryPlace,
    stage: Stage,
    reason: TemporaryReason,
    cause: ?Cause = null,
};
pub const TempAttempts = struct {
    entries: [TEMP_ATTEMPT_LIMIT]TempAttempt = undefined,
    len: u8 = 0,
    truncated: bool = false,

    /// Keep the first attempts without allocating. Later attempts set truncation.
    pub fn append(self: *TempAttempts, attempt: TempAttempt) void {
        if (self.len == TEMP_ATTEMPT_LIMIT) {
            self.truncated = true;
            return;
        }
        self.entries[self.len] = attempt;
        self.len += 1;
    }

    pub fn items(self: *const TempAttempts) []const TempAttempt {
        return self.entries[0..self.len];
    }
};

pub const FileSnapshot = struct {
    bytes: ?u64 = null,
    /// Signed nanoseconds, preserving pre-epoch observations.
    mtime_ns: ?i128 = null,
};
pub const ChangedFields = packed struct(u8) {
    path: bool = false,
    size: bool = false,
    mtime: bool = false,
    reserved: u5 = 0,
};
pub const SourceChange = struct {
    changed_fields: ChangedFields,
    before: FileSnapshot,
    after: ?FileSnapshot,
    /// A required comparison observed disappearance, rather than access denial.
    disappeared: bool = false,
};

pub const Format = enum { zfi, fai };
pub const NamingPolicy = enum { first_nonempty, all_nonempty };
pub const Grammar = enum { legacy_nul_only };
pub const Alphabet = enum { dna, rna, nucleotide, protein, custom };
pub const AlphabetSelection = enum { explicit, automatic, custom };
pub const ClassifierRule = enum { nucleotide_above_90_percent };
/// One-based inclusive interval. No conversion from BED is implicit.
pub const Interval = struct { start: u64, end: u64 };
/// Zero-based half-open byte extent in the source, not a retained buffer.
pub const ByteSpan = struct { start: u64, end: u64 };
pub const Sample = struct {
    /// Actual classified symbols, not the configured cap.
    symbols: u32,
    rule: ClassifierRule,
    interval: ?Interval = null,
    request: ?u64 = null,
};
pub const AlphabetPolicy = struct {
    alphabet: Alphabet,
    selection: AlphabetSelection,
    /// Exact custom byte set when selected; bit n allows byte n.
    custom_symbols: ?[32]u8 = null,
    sample: ?Sample = null,
};
pub const Transform = struct { reverse: bool, complement: bool };
pub const IndexField = enum {
    length,
    offset,
    line_bases,
    line_bytes,
    name_offset,
    name_bytes,
    name_blob_bytes,
    record_count,
    side_offset,
    side_count,
    side_byte_end,
    base_start,
    sequence_span,
    source_timestamp,
};
pub const FaiField = enum { name, length, offset, line_bases, line_bytes };
pub const BedField = enum { chrom, start, end };
pub const BedCoordinate = enum { start, end };
pub const CoordinateSystem = enum { one_based_inclusive, zero_based_half_open };
pub const Arithmetic = enum { add, subtract, multiply };
pub const Section = enum { header, records, side_table_count, side_table_entries, name_blob };
pub const BoundSource = enum { declared, actual };
pub const OffsetRole = enum { wire, retained, local };
pub const CountLimitRole = enum { zfi_header, scanner_count };
pub const GeometryConstraint = enum { bases_zero, bytes_zero, bytes_below_bases };
pub const Geometry = struct { line_bases: u64, line_bytes: u64 };
pub const NonuniformReason = enum {
    interior_short_line,
    wider_line,
    mixed_line_endings,
    interior_blank_line,
    sequence_whitespace,
    nonuniform_layout,
};
pub const FooterLayouts = packed struct(u8) { tight: bool, legacy: bool, reserved: u6 = 0 };
pub const LookupScope = enum { complete_catalog, retained_matches };
pub const RepairPolicy = enum { legacy_safe, legacy_format_only };
pub const IdentityComparison = enum { lexical_path, canonical_path };
pub const Component = enum { cli, indexer, index_loader, getter, bed_parser, stats, validator, reporting, platform };
pub const Check = enum {
    missing_side_table,
    retained_name,
    retained_shortest_name,
    retained_longest_name,
    record_reservation,
    record_ordinal,
    composition_accounting,
    reader_buffer,
    traversal_position,
    side_table_state,
    codec_state,
    numeric_buffer,
};

pub const MissingValue = struct { option: Argument, expected: Expected };
pub const InvalidValue = struct { option: Argument, expected: Expected, value: ?Preview = null };
pub const InvalidPath = struct { reason: PathReason, limit_bytes: ?u64 = null };
pub const Seek = struct { requested_offset: ?u64 = null };
pub const Read = struct { requested_bytes: ?u64 = null };
pub const Eof = union(enum) {
    exact: struct { requested_bytes: u64, received_bytes: u64 },
    delimiter: struct { delimiter: u8, observed_bytes: ?Measurement = null },
};
pub const Mapping = struct { map_bytes: u64 };
pub const TemporaryUnavailable = struct { purpose: TemporaryPurpose, attempts: TempAttempts };
pub const Cleanup = struct { resource_role: SourceId };
pub const Allocation = struct { allocation_role: AllocationRole, requested_bytes: ?u64 = null };
pub const Capacity = struct { capacity_role: CapacityRole, limit: u128, required: Measurement };
pub const CountOverflow = struct { counter: Counter, max_value: u128, current: ?u128 = null, increment: ?u128 = null };
pub const LineBuffer = struct { buffer_bytes: u64, observed_bytes: Measurement };
pub const Preamble = struct { content_bytes: Measurement };
pub const ControlByte = struct { byte: u8, grammar: Grammar };
pub const Duplicate = struct { identifier: ?Preview = null };
pub const InvalidSymbol = struct { byte: u8, policy: AlphabetPolicy };
pub const IdentifierLimit = struct { observed_bytes: Measurement, limit_bytes: u64, format: Format };
pub const LineWidthRole = enum { index_geometry, synthetic_separator, validator_width };
pub const LineLimit = struct { field: IndexField, observed: Measurement, limit: u64, width_role: LineWidthRole, format: ?Format = null };
pub const RecordLimit = struct { record_count: Measurement, max_records: u64, format: Format, limit_role: CountLimitRole };
pub const SideOffsetLimit = struct { offset: u64, max_offset: u64, offset_role: OffsetRole };
pub const TimestampRange = struct { timestamp_ns: i128, min_ns: i128, max_ns: i128 };
pub const Nonuniform = struct { reason: NonuniformReason, expected: ?Geometry = null, observed: ?Geometry = null };
/// Both permitted formats were attempted, in ZFI then FAI order.
pub const IndexDiscovery = struct { zfi_source: SourceId, fai_source: SourceId };
pub const Magic = struct { expected_magic: [4]u8, actual: ?Preview = null };
pub const Version = struct { found_version: u8, supported_version: u8 };
pub const Truncation = struct { section: Section, required_end: u128, available_end: u64 };
pub const NameFlags = struct { flags: u8, required_flag: u8 };
pub const NameBlob = struct { name_blob_bytes: u64, payload_end: u64, records_end: u64 };
pub const NameExtent = struct { name_offset: u64, name_bytes: u64, blob_bytes: u64 };
pub const InvalidGeometry = struct { geometry: Geometry, constraint: GeometryConstraint };
pub const SequenceOffset = struct { sequence_offset: u64, source_bytes: u64, bound_source: BoundSource };
pub const SequenceSpan = struct { required_end: u128, source_bytes: u64, bound_source: BoundSource, entry: ?u64 = null };
pub const ArithmeticOverflow = struct { field: IndexField, arithmetic: Arithmetic, lhs: u128, rhs: u128, max_value: u128 };
pub const MissingFaiField = struct { field: FaiField, column_number: u8 };
pub const FaiInteger = struct { field: FaiField, value: ?Preview = null };
pub const FaiIntegerRange = struct { field: FaiField, max_value: u64, value: ?Preview = null };
pub const SideExtent = struct { offset: u64, region_start: u64, region_end: u64 };
pub const Alignment = struct { offset: u64, alignment_bytes: u8, section: Section };
pub const SideOrder = struct { table_start: u64, minimum_start: u64 };
pub const SidePrefix = struct { entry: u64, base_start: u64, expected_base_start: u64 };
pub const SideGeometry = struct { entry: u64, geometry: Geometry, constraint: GeometryConstraint };
pub const FirstOffset = struct { first_byte_offset: u64, sequence_offset: u64 };
pub const LineOrder = struct { entry: u64, byte_offset: u64, minimum_offset: u64 };
pub const LengthMismatch = struct { indexed_length: u64, side_table_bases: u64 };
pub const IndexAge = struct { index_mtime_ns: i128, source_mtime_ns: i128 };
pub const SourceSize = struct { stored_bytes: u64, actual_bytes: u64 };
pub const SourceTime = struct { stored_mtime_ns: u64, actual_mtime_ns: i128 };
pub const RegionMismatch = struct { requested_symbols: u64, observed_symbols: u64, interval: Interval, examined_span: ByteSpan };
pub const TotalMismatch = struct { indexed_total: u64, observed_total: u64 };
pub const MissingSequence = struct { lookup_scope: LookupScope, index_source: SourceId, identifier: ?Preview = null };
pub const StartAfter = struct { start: u64, sequence_length: u64 };
pub const ReversedInterval = struct { start: u64, end: u64 };
pub const RequestIdentifier = struct { identifier_bytes: u64, limit_bytes: u64 };
pub const MissingBedField = struct { field: BedField, column_number: u8 };
pub const BedInteger = struct { field: BedCoordinate, coordinate_system: CoordinateSystem, value: ?Preview = null };
pub const BedIntegerRange = struct { field: BedCoordinate, max_value: u64, coordinate_system: CoordinateSystem, value: ?Preview = null };
pub const BedInterval = struct { start: u64, end: u64, coordinate_system: CoordinateSystem };
pub const BedStrand = struct { strand: ?Preview = null, accepted: enum { plus_minus_dot_or_absent } = .plus_minus_dot_or_absent };
pub const ProteinTransform = struct { transform: Transform, classification: Alphabet, sample: Sample, identifier: ?Preview = null };
pub const OutputAlias = struct { input: SourceId, output: SourceId, comparison: IdentityComparison };
pub const UnknownIdentity = struct { input: SourceId, output: SourceId };

/// A rewrite refusal retains the observed blocker, not a second finding count.
/// Its location is Diagnostic.at; a duplicate's first definition is related.
pub const Blocker = union(enum) {
    no_records: void,
    content_before_header: Preamble,
    header_control_byte: ControlByte,
    sequence_control_byte: ControlByte,
    duplicate_identifier: Duplicate,
    invalid_symbol: InvalidSymbol,

    pub fn code(self: *const Blocker) Code {
        return switch (self.*) {
            .no_records => .@"fasta.no_records",
            .content_before_header => .@"fasta.content_before_header",
            .header_control_byte => .@"fasta.header_control_byte",
            .sequence_control_byte => .@"fasta.sequence_control_byte",
            .duplicate_identifier => .@"fasta.duplicate_identifier",
            .invalid_symbol => .@"sequence.invalid_symbol",
        };
    }
};
pub const RepairBlock = struct { blocker: Blocker, repair_policy: RepairPolicy };
pub const Consistency = struct { component: Component, check: Check };

/// The code is the payload tag. A different code requires its own typed facts.
pub const Details = union(enum(u8)) {
    @"usage.command_required": void,
    @"usage.unknown_command": struct { command: ?Preview = null },
    @"usage.unknown_option": struct { option: ?Preview = null },
    @"usage.missing_argument": struct { argument: Argument },
    @"usage.extra_argument": struct { argument: Argument, allowed_count: u64 },
    @"usage.missing_option_value": MissingValue,
    @"usage.repeated_option": struct { option: Argument },
    @"usage.conflicting_options": struct { arguments: Arguments },
    @"usage.option_requires": struct { option: Argument, required: Requirement },
    @"usage.invalid_option_value": InvalidValue,
    @"io.file_not_found": void,
    @"io.access_denied": void,
    @"io.invalid_path": InvalidPath,
    @"io.open_failed": void,
    @"io.stat_failed": void,
    @"io.seek_failed": Seek,
    @"io.read_failed": Read,
    @"io.unexpected_eof": Eof,
    @"io.map_failed": Mapping,
    @"io.create_failed": void,
    @"io.write_failed": void,
    @"io.flush_failed": void,
    @"io.publish_failed": void,
    @"io.temp_unavailable": TemporaryUnavailable,
    @"io.cleanup_failed": Cleanup,
    @"resource.out_of_memory": Allocation,
    @"resource.capacity_limit_exceeded": Capacity,
    @"resource.count_overflow": CountOverflow,
    @"resource.line_buffer_exceeded": LineBuffer,
    @"input.changed": SourceChange,
    @"fasta.no_records": void,
    @"fasta.content_before_header": Preamble,
    @"fasta.header_control_byte": ControlByte,
    @"fasta.sequence_control_byte": ControlByte,
    @"fasta.duplicate_identifier": Duplicate,
    @"sequence.invalid_symbol": InvalidSymbol,
    @"index.initial_header_required": struct { byte: u8, required_offset: u64 = 0 },
    @"index.no_indexable_records": struct { naming_policy: NamingPolicy, seen_records: ?u64 = null, empty_records: ?u64 = null, duplicate_records: ?u64 = null },
    @"index.identifier_too_long": IdentifierLimit,
    @"index.sequence_line_too_long": LineLimit,
    @"index.record_count_limit_exceeded": RecordLimit,
    @"index.side_table_offset_limit_exceeded": SideOffsetLimit,
    @"index.source_timestamp_out_of_range": TimestampRange,
    @"index.nonuniform_fai": Nonuniform,
    @"index.not_found": IndexDiscovery,
    @"index.no_records": struct { format: Format },
    @"index.invalid_magic": Magic,
    @"index.unsupported_version": Version,
    @"index.truncated": Truncation,
    @"index.name_not_embedded": NameFlags,
    @"index.name_footer_missing": FooterLayouts,
    @"index.name_blob_out_of_bounds": NameBlob,
    @"index.name_out_of_bounds": NameExtent,
    @"index.zero_sequence_length": struct { sequence_length: u64 = 0 },
    @"index.invalid_line_geometry": InvalidGeometry,
    @"index.sequence_offset_out_of_bounds": SequenceOffset,
    @"index.sequence_span_out_of_bounds": SequenceSpan,
    @"index.arithmetic_overflow": ArithmeticOverflow,
    @"index.fai.field_missing": MissingFaiField,
    @"index.fai.invalid_integer": FaiInteger,
    @"index.fai.integer_out_of_range": FaiIntegerRange,
    @"index.side_table_offset_out_of_bounds": SideExtent,
    @"index.side_table_misaligned": Alignment,
    @"index.side_table_empty": struct { line_count: u64 = 0 },
    @"index.side_table_order_invalid": SideOrder,
    @"index.side_base_start_mismatch": SidePrefix,
    @"index.side_line_geometry_invalid": SideGeometry,
    @"index.side_first_offset_mismatch": FirstOffset,
    @"index.side_line_order_invalid": LineOrder,
    @"index.side_sequence_length_mismatch": LengthMismatch,
    @"index.older_than_source": IndexAge,
    @"index.source_size_mismatch": SourceSize,
    @"index.source_mtime_mismatch": SourceTime,
    @"pair.region_length_mismatch": RegionMismatch,
    @"pair.total_length_mismatch": TotalMismatch,
    @"request.no_regions": struct { request_source: enum { positional, names, bed } },
    @"request.sequence_not_found": MissingSequence,
    @"request.start_before_one": struct { start: u64, coordinate_system: CoordinateSystem },
    @"request.start_after_sequence": StartAfter,
    @"request.end_before_start": ReversedInterval,
    @"request.identifier_too_long": RequestIdentifier,
    @"request.bed.field_missing": MissingBedField,
    @"request.bed.invalid_coordinate": BedInteger,
    @"request.bed.coordinate_out_of_range": BedIntegerRange,
    @"request.bed.invalid_interval": BedInterval,
    @"request.bed.invalid_strand": BedStrand,
    @"transform.complement_on_protein": ProteinTransform,
    @"repair.output_is_input": OutputAlias,
    @"repair.output_identity_unknown": UnknownIdentity,
    @"repair.blocked": RepairBlock,
    @"internal.consistency_failed": Consistency,
};

pub const Code = std.meta.Tag(Details);

pub const Metadata = struct {
    code: Code,
    name: []const u8,
    category: Category,
    kind: Kind,
    severity: Severity,
    unit: Unit,
    action: Action,
    returned_error: Error,
};

/// One compiled registry. Tag names supply the stable strings; the exhaustive
/// switch requires metadata for every payload. Neither formatters nor callers
/// maintain a second list of code strings.
pub const REGISTRY = blk: {
    const codes = std.enums.values(Code);
    var entries: [codes.len]Metadata = undefined;
    for (codes, 0..) |code, index| entries[index] = describeCode(code);
    break :blk entries;
};

/// Static metadata; no formatting or allocation occurs at lookup.
pub fn metadata(code: Code) *const Metadata {
    return &REGISTRY[@intFromEnum(code)];
}

/// Source, location, original ordinals, and one related location are independent
/// of optional byte previews. stage names the observed work, not a guessed cause.
pub const Diagnostic = struct {
    details: Details,
    stage: Stage,
    at: Reference = .{},
    related: ?Reference = null,
    cause: ?Cause = null,

    pub fn code(self: *const Diagnostic) Code {
        return std.meta.activeTag(self.details);
    }
};

/// Reserve once in the caller. No helper returns this object by value.
/// Previews and retained labels own copied bytes; borrowed labels live with
/// their caller through handling.
/// Reset invalidates all preview views and registered source references.
/// This stores ordinary failures, not validator finding counts or event owners.
pub const Failure = struct {
    primary_storage: Diagnostic = undefined,
    has_primary: bool = false,
    primary_deferred: bool = false,
    secondary: [SECONDARY_LIMIT]Diagnostic = undefined,
    secondary_len: u8 = 0,
    secondary_truncated: bool = false,
    previews: [PREVIEW_BYTES]u8 = undefined,
    preview_len: u16 = 0,
    sources: Sources = .{},
    operation: Operation,
    /// Work can finish before publication, cleanup, or report delivery fails.
    work_complete: bool = false,
    /// Includes required side effects. Any blocking failure sets this false.
    complete: bool = false,
    outputs: Outputs = .{},

    /// Initialize undefined caller storage or reset it once at operation entry.
    /// Never call this from a propagating helper; existing facts must survive.
    /// Only live fields reset; unused preview/source/secondary arrays stay untouched.
    pub fn reset(self: *Failure, operation: Operation) void {
        self.has_primary = false;
        self.primary_deferred = false;
        self.secondary_len = 0;
        self.secondary_truncated = false;
        self.preview_len = 0;
        self.sources.registered = 0;
        self.operation = operation;
        self.work_complete = false;
        self.complete = false;
        self.outputs = .{};
    }

    /// Copy up to limit bytes and the remaining shared budget. An exhausted
    /// budget returns a truncated preview; it cannot prevent a failure return.
    pub fn retain(self: *Failure, bytes: []const u8, limit: u16) Preview {
        return self.retainPrefix(bytes, .{ .exact = bytes.len }, limit);
    }

    /// A streaming parser can keep a prefix without retaining the full field.
    /// original records its exact length or observed lower bound. A zero limit
    /// retains that scalar evidence even when byte previews are disabled.
    pub fn retainPrefix(self: *Failure, prefix: []const u8, original: ByteLength, limit: u16) Preview {
        const original_bytes = switch (original) {
            .exact, .at_least => |bytes| bytes,
        };
        std.debug.assert(prefix.len <= original_bytes);
        const count = @min(prefix.len, limit, PREVIEW_BYTES - self.preview_len);
        const start = self.preview_len;
        @memcpy(self.previews[start..][0..count], prefix[0..count]);
        self.preview_len += @intCast(count);
        return .{
            .original_bytes = original_bytes,
            .start = start,
            .len = @intCast(count),
            .truncated = count != original_bytes or original == .at_least,
            .length_kind = std.meta.activeTag(original),
        };
    }

    /// A view borrows this Failure. Do not retain it across reset/overwrite.
    pub fn preview(self: *const Failure, reference: Preview) []const u8 {
        const end = @as(usize, reference.start) + reference.len;
        std.debug.assert(end <= self.preview_len);
        return self.previews[reference.start..end];
    }

    /// Resolve base label bytes, excluding the static suffix. Borrowed labels
    /// live with their caller; retained labels borrow this Failure until reset.
    /// A retained label can be truncated; Source.label retains its length facts.
    pub fn sourceLabel(self: *const Failure, id: SourceId) ?[]const u8 {
        const source = self.sources.get(id) orelse return null;
        return switch (source.label) {
            .borrowed => |bytes| bytes,
            .retained => |reference| self.preview(reference),
        };
    }

    /// The view borrows this Failure until reset or later mutation. A separate
    /// presence flag avoids clearing the inactive payload at entry.
    pub fn primary(self: *const Failure) ?*const Diagnostic {
        return if (self.has_primary) &self.primary_storage else null;
    }

    /// Capture the first failure, or retain a later one as secondary. Returning
    /// the primary's projection keeps tag and facts consistent during cleanup.
    /// A deferred failure uses the same storage; promote selects a later primary
    /// only where the owning operation explicitly requires that priority.
    pub fn fail(self: *Failure, diagnostic: Diagnostic) Error {
        @branchHint(.cold);
        if (self.primary_deferred) return self.promote(diagnostic);
        self.complete = false;
        if (!self.has_primary) {
            self.primary_storage = diagnostic;
            self.has_primary = true;
        } else {
            self.addSecondary(diagnostic);
        }
        return metadata(self.primary_storage.code()).returned_error;
    }

    /// Retain a deliberately deferred failure while work continues. A later
    /// ordinary fail promotes its own facts and preserves this first cause.
    /// Subsequent deferred observations use the bounded secondary slots.
    pub fn deferFailure(self: *Failure, diagnostic: Diagnostic) void {
        @branchHint(.cold);
        std.debug.assert(!self.has_primary or self.primary_deferred);
        self.complete = false;
        if (!self.has_primary) {
            self.primary_storage = diagnostic;
            self.has_primary = true;
            self.primary_deferred = true;
        } else {
            self.addSecondary(diagnostic);
        }
    }

    /// After higher-priority work, return the retained failure without creating
    /// another diagnostic. This does not finish I/O or mark work complete.
    pub fn finishDeferred(self: *Failure) Error!void {
        if (self.has_primary) {
            @branchHint(.cold);
            self.primary_deferred = false;
            return metadata(self.primary_storage.code()).returned_error;
        }
    }

    /// Preserve the primary when checked cleanup/reporting also fails.
    pub fn addSecondary(self: *Failure, diagnostic: Diagnostic) void {
        @branchHint(.cold);
        std.debug.assert(self.has_primary);
        self.complete = false;
        if (self.secondary_len == SECONDARY_LIMIT) {
            self.secondary_truncated = true;
            return;
        }
        self.secondary[self.secondary_len] = diagnostic;
        self.secondary_len += 1;
    }

    /// Select an explicitly higher-priority failure. The displaced primary
    /// becomes the first secondary; retained causes/previews are unchanged.
    pub fn promote(self: *Failure, diagnostic: Diagnostic) Error {
        @branchHint(.cold);
        if (self.has_primary) {
            if (self.secondary_len == SECONDARY_LIMIT) {
                self.secondary_truncated = true;
            } else {
                self.secondary_len += 1;
            }
            var index = self.secondary_len - 1;
            while (index > 0) : (index -= 1) self.secondary[index] = self.secondary[index - 1];
            self.secondary[0] = self.primary_storage;
        }
        self.primary_storage = diagnostic;
        self.has_primary = true;
        self.primary_deferred = false;
        self.complete = false;
        return metadata(diagnostic.code()).returned_error;
    }

    /// Borrowed diagnostic views expire at reset or later mutation.
    pub fn secondaryItems(self: *const Failure) []const Diagnostic {
        return self.secondary[0..self.secondary_len];
    }

    /// Mark completion only after work and all required finish/cleanup succeed.
    pub fn succeed(self: *Failure) void {
        std.debug.assert(!self.has_primary);
        self.work_complete = true;
        self.complete = true;
    }
};

/// Convert an observed zero-based ordinal with checked addition. The maximum
/// cannot wrap to an invented zero; unknown callers do not call this function.
pub fn oneBased(failure: *Failure, zero_based: u64, at: Reference) Error!u64 {
    return std.math.add(u64, zero_based, 1) catch {
        @branchHint(.cold);
        return failure.fail(.{
            .details = .{ .@"resource.count_overflow" = .{
                .counter = .report_ordinal,
                .max_value = std.math.maxInt(u64),
                .current = zero_based,
                .increment = 1,
            } },
            .stage = .report,
            .at = at,
        });
    };
}

fn describe(code: Code, category: Category, kind: Kind, unit: Unit, action: Action, returned_error: Error) Metadata {
    return .{ .code = code, .name = @tagName(code), .category = category, .kind = kind, .severity = .@"error", .unit = unit, .action = action, .returned_error = returned_error };
}

fn describeCode(code: Code) Metadata {
    return switch (code) {
        .@"usage.command_required" => describe(code, .usage, .failure, .none, .correct_arguments, error.InvalidArgument),
        .@"usage.unknown_command" => describe(code, .usage, .failure, .none, .correct_arguments, error.InvalidArgument),
        .@"usage.unknown_option" => describe(code, .usage, .failure, .none, .correct_arguments, error.InvalidArgument),
        .@"usage.missing_argument" => describe(code, .usage, .failure, .none, .correct_arguments, error.InvalidArgument),
        .@"usage.extra_argument" => describe(code, .usage, .failure, .none, .correct_arguments, error.InvalidArgument),
        .@"usage.missing_option_value" => describe(code, .usage, .failure, .none, .correct_arguments, error.InvalidArgument),
        .@"usage.repeated_option" => describe(code, .usage, .failure, .none, .correct_arguments, error.InvalidArgument),
        .@"usage.conflicting_options" => describe(code, .usage, .failure, .none, .correct_arguments, error.InvalidArgument),
        .@"usage.option_requires" => describe(code, .usage, .failure, .none, .correct_arguments, error.InvalidArgument),
        .@"usage.invalid_option_value" => describe(code, .usage, .failure, .none, .correct_arguments, error.InvalidArgument),
        .@"io.file_not_found" => describe(code, .io, .failure, .none, .inspect_io, error.Io),
        .@"io.access_denied" => describe(code, .io, .failure, .none, .inspect_io, error.Io),
        .@"io.invalid_path" => describe(code, .io, .failure, .none, .inspect_io, error.Io),
        .@"io.open_failed" => describe(code, .io, .failure, .none, .inspect_io, error.Io),
        .@"io.stat_failed" => describe(code, .io, .failure, .none, .inspect_io, error.Io),
        .@"io.seek_failed" => describe(code, .io, .failure, .none, .inspect_io, error.Io),
        .@"io.read_failed" => describe(code, .io, .failure, .none, .inspect_io, error.Io),
        .@"io.unexpected_eof" => describe(code, .io, .failure, .none, .inspect_io, error.Io),
        .@"io.map_failed" => describe(code, .io, .failure, .none, .inspect_io, error.Io),
        .@"io.create_failed" => describe(code, .io, .failure, .none, .inspect_io, error.Io),
        .@"io.write_failed" => describe(code, .io, .failure, .none, .inspect_io, error.Io),
        .@"io.flush_failed" => describe(code, .io, .failure, .none, .inspect_io, error.Io),
        .@"io.publish_failed" => describe(code, .io, .failure, .none, .inspect_io, error.Io),
        .@"io.temp_unavailable" => describe(code, .io, .failure, .none, .inspect_io, error.Io),
        .@"io.cleanup_failed" => describe(code, .io, .failure, .none, .inspect_io, error.Io),
        .@"resource.out_of_memory" => describe(code, .resources, .failure, .none, .check_memory, error.OutOfMemory),
        .@"resource.capacity_limit_exceeded" => describe(code, .resources, .failure, .none, .check_capacity, error.CapacityExceeded),
        .@"resource.count_overflow" => describe(code, .resources, .failure, .none, .check_counts, error.CountOverflow),
        .@"resource.line_buffer_exceeded" => describe(code, .resources, .failure, .none, .check_capacity, error.CapacityExceeded),
        .@"input.changed" => describe(code, .source_stability, .failure, .none, .check_source, error.InputChanged),
        .@"fasta.no_records" => describe(code, .structure, .finding, .file, .correct_data, error.InvalidFasta),
        .@"fasta.content_before_header" => describe(code, .structure, .finding, .line, .correct_data, error.InvalidFasta),
        .@"fasta.header_control_byte" => describe(code, .structure, .finding, .byte, .correct_data, error.InvalidFasta),
        .@"fasta.sequence_control_byte" => describe(code, .structure, .finding, .byte, .correct_data, error.InvalidFasta),
        .@"fasta.duplicate_identifier" => describe(code, .identifiers, .finding, .record, .correct_data, error.InvalidFasta),
        .@"sequence.invalid_symbol" => describe(code, .alphabet, .finding, .byte, .correct_data, error.InvalidFasta),
        .@"index.initial_header_required" => describe(code, .index_compatibility, .finding, .file, .correct_data, error.InvalidIndex),
        .@"index.no_indexable_records" => describe(code, .index_compatibility, .finding, .file, .correct_data, error.InvalidIndex),
        .@"index.identifier_too_long" => describe(code, .index_compatibility, .finding, .record, .correct_data, error.InvalidIndex),
        .@"index.sequence_line_too_long" => describe(code, .index_compatibility, .finding, .line, .correct_data, error.InvalidIndex),
        .@"index.record_count_limit_exceeded" => describe(code, .index_compatibility, .finding, .file, .correct_data, error.InvalidIndex),
        .@"index.side_table_offset_limit_exceeded" => describe(code, .index_compatibility, .finding, .record, .correct_data, error.InvalidIndex),
        .@"index.source_timestamp_out_of_range" => describe(code, .index_compatibility, .finding, .file, .check_source, error.InvalidIndex),
        .@"index.nonuniform_fai" => describe(code, .index_compatibility, .finding, .record, .select_index_format, error.InvalidIndex),
        .@"index.not_found" => describe(code, .index_integrity, .failure, .none, .rebuild_index, error.InvalidIndex),
        .@"index.no_records" => describe(code, .index_integrity, .finding, .file, .rebuild_index, error.InvalidIndex),
        .@"index.invalid_magic" => describe(code, .index_integrity, .finding, .file, .rebuild_index, error.InvalidIndex),
        .@"index.unsupported_version" => describe(code, .index_integrity, .finding, .file, .rebuild_index, error.InvalidIndex),
        .@"index.truncated" => describe(code, .index_integrity, .finding, .file, .rebuild_index, error.InvalidIndex),
        .@"index.name_not_embedded" => describe(code, .index_integrity, .finding, .record, .rebuild_index, error.InvalidIndex),
        .@"index.name_footer_missing" => describe(code, .index_integrity, .finding, .file, .rebuild_index, error.InvalidIndex),
        .@"index.name_blob_out_of_bounds" => describe(code, .index_integrity, .finding, .file, .rebuild_index, error.InvalidIndex),
        .@"index.name_out_of_bounds" => describe(code, .index_integrity, .finding, .record, .rebuild_index, error.InvalidIndex),
        .@"index.zero_sequence_length" => describe(code, .index_integrity, .finding, .record, .rebuild_index, error.InvalidIndex),
        .@"index.invalid_line_geometry" => describe(code, .index_integrity, .finding, .record, .rebuild_index, error.InvalidIndex),
        .@"index.sequence_offset_out_of_bounds" => describe(code, .index_integrity, .finding, .record, .rebuild_index, error.InvalidIndex),
        .@"index.sequence_span_out_of_bounds" => describe(code, .index_integrity, .finding, .record, .rebuild_index, error.InvalidIndex),
        .@"index.arithmetic_overflow" => describe(code, .index_integrity, .finding, .field, .rebuild_index, error.InvalidIndex),
        .@"index.fai.field_missing" => describe(code, .index_integrity, .finding, .field, .rebuild_index, error.InvalidIndex),
        .@"index.fai.invalid_integer" => describe(code, .index_integrity, .finding, .field, .rebuild_index, error.InvalidIndex),
        .@"index.fai.integer_out_of_range" => describe(code, .index_integrity, .finding, .field, .rebuild_index, error.InvalidIndex),
        .@"index.side_table_offset_out_of_bounds" => describe(code, .index_integrity, .finding, .record, .rebuild_index, error.InvalidIndex),
        .@"index.side_table_misaligned" => describe(code, .index_integrity, .finding, .record, .rebuild_index, error.InvalidIndex),
        .@"index.side_table_empty" => describe(code, .index_integrity, .finding, .record, .rebuild_index, error.InvalidIndex),
        .@"index.side_table_order_invalid" => describe(code, .index_integrity, .finding, .record, .rebuild_index, error.InvalidIndex),
        .@"index.side_base_start_mismatch" => describe(code, .index_integrity, .finding, .record, .rebuild_index, error.InvalidIndex),
        .@"index.side_line_geometry_invalid" => describe(code, .index_integrity, .finding, .record, .rebuild_index, error.InvalidIndex),
        .@"index.side_first_offset_mismatch" => describe(code, .index_integrity, .finding, .record, .rebuild_index, error.InvalidIndex),
        .@"index.side_line_order_invalid" => describe(code, .index_integrity, .finding, .record, .rebuild_index, error.InvalidIndex),
        .@"index.side_sequence_length_mismatch" => describe(code, .index_integrity, .finding, .record, .rebuild_index, error.InvalidIndex),
        .@"index.older_than_source" => describe(code, .index_match, .finding, .file, .rebuild_index, error.InvalidIndex),
        .@"index.source_size_mismatch" => describe(code, .index_match, .finding, .file, .rebuild_index, error.InvalidIndex),
        .@"index.source_mtime_mismatch" => describe(code, .index_match, .finding, .file, .rebuild_index, error.InvalidIndex),
        .@"pair.region_length_mismatch" => describe(code, .index_match, .finding, .record, .rebuild_index, error.InvalidIndex),
        .@"pair.total_length_mismatch" => describe(code, .index_match, .finding, .file, .rebuild_index, error.InvalidIndex),
        .@"request.no_regions" => describe(code, .requests, .failure, .none, .correct_request, error.InvalidRequest),
        .@"request.sequence_not_found" => describe(code, .requests, .failure, .none, .correct_request, error.InvalidRequest),
        .@"request.start_before_one" => describe(code, .requests, .failure, .none, .correct_request, error.InvalidRequest),
        .@"request.start_after_sequence" => describe(code, .requests, .failure, .none, .correct_request, error.InvalidRequest),
        .@"request.end_before_start" => describe(code, .requests, .failure, .none, .correct_request, error.InvalidRequest),
        .@"request.identifier_too_long" => describe(code, .requests, .failure, .none, .correct_request, error.InvalidRequest),
        .@"request.bed.field_missing" => describe(code, .requests, .failure, .none, .correct_request, error.InvalidRequest),
        .@"request.bed.invalid_coordinate" => describe(code, .requests, .failure, .none, .correct_request, error.InvalidRequest),
        .@"request.bed.coordinate_out_of_range" => describe(code, .requests, .failure, .none, .correct_request, error.InvalidRequest),
        .@"request.bed.invalid_interval" => describe(code, .requests, .failure, .none, .correct_request, error.InvalidRequest),
        .@"request.bed.invalid_strand" => describe(code, .requests, .failure, .none, .correct_request, error.InvalidRequest),
        .@"transform.complement_on_protein" => describe(code, .transforms, .failure, .none, .select_transform, error.UnsupportedTransform),
        .@"repair.output_is_input" => describe(code, .repair, .failure, .none, .select_output, error.RepairRefused),
        .@"repair.output_identity_unknown" => describe(code, .repair, .failure, .none, .select_output, error.RepairRefused),
        .@"repair.blocked" => describe(code, .repair, .failure, .none, .resolve_blocker, error.RepairRefused),
        .@"internal.consistency_failed" => describe(code, .internal, .failure, .none, .report_defect, error.Internal),
    };
}

comptime {
    @setEvalBranchQuota(30000);
    for (REGISTRY, 0..) |entry, index| {
        if (@intFromEnum(entry.code) != index) @compileError("diagnostic registry tag order differs");
        for (REGISTRY[0..index]) |previous| {
            if (std.mem.eql(u8, entry.name, previous.name)) @compileError("duplicate diagnostic code");
        }
    }
    if (std.enums.values(SourceId).len != SOURCE_LIMIT) @compileError("diagnostic source role limit differs");
    if (@sizeOf(Failure) > FAILURE_BYTES_LIMIT) @compileError("diagnostic failure storage exceeds 2048 bytes");
}
