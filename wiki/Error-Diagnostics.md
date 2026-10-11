# Error Diagnostics

z-fasta reports command failures on stderr with stable symbolic codes and observed facts. Successful FASTA, FAI and statistics output stays on stdout. Validation findings keep their existing stdout format, including JSON v1.

## Choose a renderer

```bash
z-fasta get --error-format=json genome.fa chr1:1-100
z-fasta --error-format=text stats genome.fa
```

The default is `--error-format=text`. The exact alternatives are `--error-format=text` and `--error-format=json`; the space-separated form is not accepted. The option can precede the command or appear among its arguments. The last valid repetition wins. Values of `--names`, `--bed`, `-o`, `--schema`, `--custom-alphabet` and `--max-header-len` remain option values even if their bytes resemble a renderer option.

Help and version retain their position-dependent behavior. An invalid renderer option uses a short text fallback and exit 1. Failure while allocating command arguments also uses text. A broken stderr sink can prevent any report from arriving. A consumer must handle missing JSON as a failure to obtain a report.

## Text and machine examples

Running `z-fasta get` returns exit 1, empty stdout and:

```text
error [usage.missing_argument] get: A required argument is missing. argument=fasta_path
hint: Correct the command arguments; use --help for accepted options.
```

Running `z-fasta get --error-format=json` returns exit 1 and one newline-terminated JSON object on stderr. This excerpt comes from that report; the complete object also includes `validation`, `sources`, `outputs`, `metrics` and `action`:

```json
{
  "schema_version": "v2",
  "catalog_revision": "0.4.0-1",
  "producer": {
    "name": "z-fasta",
    "version": "0.3.5"
  },
  "kind": "failure",
  "operation": {
    "name": "get",
    "complete": false,
    "work_complete": false
  },
  "accepted": false,
  "failure": {
    "primary": {
      "code": "usage.missing_argument",
      "category": "usage",
      "severity": "error",
      "kind": "failure",
      "unit": "none",
      "stage": "arguments",
      "source": null,
      "location": {
        "offset": null,
        "line": null,
        "column": null
      },
      "record": null,
      "request": null,
      "related": [],
      "cause": null,
      "data": {
        "argument": "fasta_path"
      }
    },
    "secondary": [],
    "secondary_truncated": false
  }
}
```

Read `failure.primary.code`, not the English message or the numeric Zig error tag. The primary records the selected condition; `failure.secondary` retains up to two later failures, such as failed cleanup or reporting. `secondary_truncated` says whether more secondary failures were omitted.

## Interpret the report

`schema_version` identifies the envelope format. `catalog_revision` identifies the code definitions (`0.4.0-1` for this catalog), and `producer.version` identifies the executable. They are separate values. Preserve unknown codes and fields when forwarding a report; do not turn an unknown code into success or infer a retry policy from its category or suggested `action`.

`operation.work_complete` records whether the command's work finished. `operation.complete` also requires its required side effects to finish. `accepted` cannot be true for incomplete work or a blocking failure. A failed notice can therefore return exit 1 after data has already reached stdout.

`outputs` describes stdout, a requested file or caller writer, and stderr separately. States are `not_started`, `staged`, `published`, `streaming`, `unknown` or `aborted`. `records` counts records completed by the writer. `logical_bytes` counts bytes accepted by the writer, including buffered bytes. Neither proves the visible prefix is complete. `may_be_partial` marks a possible prefix, and `final_flush` is `not_attempted`, `succeeded` or `failed`. The serialized object cannot certify its own final stderr flush. Always check the actual process status and whether the full JSON object arrived.

Source IDs name roles within this command: FASTA, the selected sidecar, names or BED input, stdout, stderr, requested output, or a temporary spool. Their labels can be truncated; an unavailable label or location stays null. Derived sidecar labels use the original FASTA path and suffix. Offsets are zero-based byte offsets; lines, byte columns, records and requests are one-based original positions, not slots in a retained array. A code's unit describes the reported condition, not an occurrence count or a promise that every location is known.

Integers wider than 32 bits are decimal strings, including lengths, offsets, ordinals, timestamps and GET summary metrics. Small integers such as a byte value, format version or flags remain JSON numbers. Missing observations are null. An exact measurement uses its ordinary field; an observed lower bound uses a `_min` field. Do not read a lower bound as an exact count.

Opaque paths, identifiers and byte previews use an object with `encoding`, `value`, `truncated` and either `original_bytes` or `observed_bytes_min`. Valid UTF-8 uses `encoding: utf8`; other bytes use `encoding: base64`. Decode according to that field. A truncated preview is only the retained prefix, even after decoding. Native causes retain the observed Zig name or native domain and number; no missing errno is guessed.

## Keep status, data and reports together

This shell job retains failure evidence and publishes extracted data only after exit 0:

```bash
work=$(mktemp -d)
status=0
z-fasta get --error-format=json genome.fa chr1:1-100 > "$work/result.fa" 2> "$work/report.json" || status=$?
if [ "$status" -ne 0 ]; then
    printf 'GET failed with exit %s; files retained in %s\n' "$status" "$work" >&2
    exit "$status"
fi
mv "$work/result.fa" result.fa
rm -r "$work"
```

GET can write earlier valid records before a later request fails. Do not publish its temporary stdout file merely because it is non-empty. Ordinary successful machine-mode commands need not emit a stderr report. `get --summary --error-format=json` requests a summary envelope with decimal-string `region_count`, `total_bases` and `elapsed_ns`; its successful FASTA still goes to stdout.

Validate needs a separate policy: exit 0 means no remaining issues, exit 2 means warnings only, and exit 1 means errors, strict warning rejection, an incomplete retained report without fix, or an operational failure. With `--json --summary`, parse its v1 stdout report and the actual status. Allow exit 2 only if your job explicitly permits warnings. `--error-format=json` selects operational stderr diagnostics; it does not convert v1 findings into v2. Machine mode omits the prose event-cap notice, while v1 `truncated` and the existing status still report the cap. See [Validation and repair](Validation-and-Repair) for counting and repair limits.

A local consumer example in `tests/check_diagnostics.sh` exercises these policies, including missing and truncated JSON. It does not establish integration with a workflow scheduler, HTTP service or telemetry system.

## Native callers

Core command functions return a finite `diagnostics.Error!T` and fill an explicit caller-owned `diagnostics.Failure`. They return after releasing their working resources. They do not select a renderer or terminate the process. Multiple codes share a Zig error tag, so handle the tag for propagation and inspect the captured code and facts for the exact reason.

`tests/diagnostics_commands_caller.zig` is a compiled example using the existing command exports. It handles a missing-file failure from each of index, GET, stats and validate, then successfully calls each again without running the CLI dispatcher. Run it with `zig build test-diagnostics`. These are internal command APIs, not a newly published stable FASTA library.

Failure construction allocates no message or growing chain. It retains previews in a shared 256-byte pool, up to two secondary diagnostics, and 16 source roles. Borrowed labels must live through handling or the next reset. Preview views borrow the owning `Failure` until reset or overwrite. A validator `Summary` owns its event and width arrays but borrows event names from input bytes; retain the input until you finish using those names and call `Summary.deinit` with the original allocator. These storage limits do not bound validator memory use: its mapping, names and per-record rewrite data still grow with input.

See [For contributors](For-Contributors) for the source owners and checks.

## Active code catalog

The registry contains 95 codes. This table lists the native detail fields from `diagnostics.Details`; shared `stage`, `source`, `location`, `record`, `request`, `related` and `cause` sit outside those fields. The finding codes also describe retained repair blockers. They do not imply that every validator finding is emitted as a v2 diagnostic.

The machine renderer preserves the native fields with these explicit conversions:

- `io.unexpected_eof` uses `read_kind` plus the selected read's fields.
- `index.unsupported_version` uses `found_version` and `supported_versions`.
- `index.not_found` uses `attempted_formats` and `attempted_sources`.
- `repair.blocked` uses `reason_code`, the selected `blocker` fields and `repair_policy`.
- Argument sets become arrays of selected argument names. Temporary attempts become `entries` plus `truncated`. Byte previews and exact/lower-bound measurements follow the encodings above.

A field ending in `offset` or describing a byte extent refers to bytes in its named source. Lengths and start/end fields identify the sequence or coordinate system named in that code's data. Expected/found values retain the original check's values. Null fields mean unobserved facts, not zero. The typed definitions and renderer are the source of field types and nested shapes.

| Code | Kind | Unit | Detail fields |
| --- | --- | --- | --- |
| `usage.command_required` | failure | none | none |
| `usage.unknown_command` | failure | none | `command` |
| `usage.unknown_option` | failure | none | `option` |
| `usage.missing_argument` | failure | none | `argument` |
| `usage.extra_argument` | failure | none | `argument`, `allowed_count` |
| `usage.missing_option_value` | failure | none | `option`, `expected` |
| `usage.repeated_option` | failure | none | `option` |
| `usage.conflicting_options` | failure | none | `arguments` |
| `usage.option_requires` | failure | none | `option`, `required` |
| `usage.invalid_option_value` | failure | none | `option`, `expected`, `value` |
| `io.file_not_found` | failure | none | none |
| `io.access_denied` | failure | none | none |
| `io.invalid_path` | failure | none | `reason`, `limit_bytes` |
| `io.open_failed` | failure | none | none |
| `io.stat_failed` | failure | none | none |
| `io.seek_failed` | failure | none | `requested_offset` |
| `io.read_failed` | failure | none | `requested_bytes` |
| `io.unexpected_eof` | failure | none | `exact`, `delimiter` |
| `io.map_failed` | failure | none | `map_bytes` |
| `io.create_failed` | failure | none | none |
| `io.write_failed` | failure | none | none |
| `io.flush_failed` | failure | none | none |
| `io.publish_failed` | failure | none | none |
| `io.temp_unavailable` | failure | none | `purpose`, `attempts` |
| `io.cleanup_failed` | failure | none | `resource_role` |
| `resource.out_of_memory` | failure | none | `allocation_role`, `requested_bytes` |
| `resource.capacity_limit_exceeded` | failure | none | `capacity_role`, `limit`, `required` |
| `resource.count_overflow` | failure | none | `counter`, `max_value`, `current`, `increment` |
| `resource.line_buffer_exceeded` | failure | none | `buffer_bytes`, `observed_bytes` |
| `input.changed` | failure | none | `changed_fields`, `before`, `after`, `disappeared` |
| `input.unsupported_format` | failure | none | `detected` |
| `fasta.no_records` | finding | file | none |
| `fasta.content_before_header` | finding | line | `content_bytes` |
| `fasta.header_control_byte` | finding | byte | `byte`, `grammar` |
| `fasta.sequence_control_byte` | finding | byte | `byte`, `grammar` |
| `fasta.duplicate_identifier` | finding | record | `identifier` |
| `sequence.invalid_symbol` | finding | byte | `byte`, `policy` |
| `index.initial_header_required` | finding | file | `byte`, `required_offset` |
| `index.no_indexable_records` | finding | file | `naming_policy`, `seen_records`, `empty_records`, `duplicate_records` |
| `index.identifier_too_long` | finding | record | `observed_bytes`, `limit_bytes`, `format` |
| `index.sequence_line_too_long` | finding | line | `field`, `observed`, `limit`, `width_role`, `format` |
| `index.record_count_limit_exceeded` | finding | file | `record_count`, `max_records`, `format`, `limit_role` |
| `index.side_table_offset_limit_exceeded` | finding | record | `offset`, `max_offset`, `offset_role` |
| `index.source_timestamp_out_of_range` | finding | file | `timestamp_ns`, `min_ns`, `max_ns` |
| `index.nonuniform_fai` | finding | record | `reason`, `expected`, `observed` |
| `index.not_found` | failure | none | `zfi_source`, `fai_source` |
| `index.no_records` | finding | file | `format` |
| `index.invalid_magic` | finding | file | `expected_magic`, `actual` |
| `index.unsupported_version` | finding | file | `found_version`, `supported_version` |
| `index.truncated` | finding | file | `section`, `required_end`, `available_end` |
| `index.name_not_embedded` | finding | record | `flags`, `required_flag` |
| `index.name_footer_missing` | finding | file | `tight`, `legacy`, `reserved` |
| `index.name_blob_out_of_bounds` | finding | file | `name_blob_bytes`, `payload_end`, `records_end` |
| `index.name_out_of_bounds` | finding | record | `name_offset`, `name_bytes`, `blob_bytes` |
| `index.zero_sequence_length` | finding | record | `sequence_length` |
| `index.invalid_line_geometry` | finding | record | `geometry`, `constraint` |
| `index.sequence_offset_out_of_bounds` | finding | record | `sequence_offset`, `source_bytes`, `bound_source` |
| `index.sequence_span_out_of_bounds` | finding | record | `required_end`, `source_bytes`, `bound_source`, `entry` |
| `index.arithmetic_overflow` | finding | field | `field`, `arithmetic`, `lhs`, `rhs`, `max_value` |
| `index.fai.field_missing` | finding | field | `field`, `column_number` |
| `index.fai.invalid_integer` | finding | field | `field`, `value` |
| `index.fai.integer_out_of_range` | finding | field | `field`, `max_value`, `value` |
| `index.side_table_offset_out_of_bounds` | finding | record | `offset`, `region_start`, `region_end` |
| `index.side_table_misaligned` | finding | record | `offset`, `alignment_bytes`, `section` |
| `index.side_table_empty` | finding | record | `line_count` |
| `index.side_table_order_invalid` | finding | record | `table_start`, `minimum_start` |
| `index.side_base_start_mismatch` | finding | record | `entry`, `base_start`, `expected_base_start` |
| `index.side_line_geometry_invalid` | finding | record | `entry`, `geometry`, `constraint` |
| `index.side_first_offset_mismatch` | finding | record | `first_byte_offset`, `sequence_offset` |
| `index.side_line_order_invalid` | finding | record | `entry`, `byte_offset`, `minimum_offset` |
| `index.side_sequence_length_mismatch` | finding | record | `indexed_length`, `side_table_bases` |
| `index.older_than_source` | finding | file | `index_mtime_ns`, `source_mtime_ns` |
| `index.source_size_mismatch` | finding | file | `stored_bytes`, `actual_bytes` |
| `index.source_mtime_mismatch` | finding | file | `stored_mtime_ns`, `actual_mtime_ns` |
| `pair.region_length_mismatch` | finding | record | `requested_symbols`, `observed_symbols`, `interval`, `examined_span` |
| `pair.total_length_mismatch` | finding | file | `indexed_total`, `observed_total` |
| `request.no_regions` | failure | none | `request_source` |
| `request.sequence_not_found` | failure | none | `lookup_scope`, `index_source`, `identifier` |
| `request.start_before_one` | failure | none | `start`, `coordinate_system` |
| `request.start_after_sequence` | failure | none | `start`, `sequence_length` |
| `request.end_before_start` | failure | none | `start`, `end` |
| `request.end_after_sequence` | failure | none | `end`, `sequence_length` |
| `request.ordinal_out_of_range` | failure | none | `ordinal`, `record_count` |
| `request.output_buffer_too_small` | failure | none | `required_bytes`, `available_bytes` |
| `request.identifier_too_long` | failure | none | `identifier_bytes`, `limit_bytes` |
| `request.bed.field_missing` | failure | none | `field`, `column_number` |
| `request.bed.invalid_coordinate` | failure | none | `field`, `coordinate_system`, `value` |
| `request.bed.coordinate_out_of_range` | failure | none | `field`, `max_value`, `coordinate_system`, `value` |
| `request.bed.invalid_interval` | failure | none | `start`, `end`, `coordinate_system` |
| `request.bed.invalid_strand` | failure | none | `strand`, `accepted` |
| `transform.complement_on_protein` | failure | none | `transform`, `classification`, `sample`, `identifier` |
| `repair.output_is_input` | failure | none | `input`, `output`, `comparison` |
| `repair.output_identity_unknown` | failure | none | `input`, `output` |
| `repair.blocked` | failure | none | `blocker`, `repair_policy` |
| `internal.consistency_failed` | failure | none | `component`, `check` |
