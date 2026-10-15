# Zig Package

z-fasta is also a Zig package. It gives a Zig program the same FASTA rules as the CLI, without starting a process:

- `IndexedReader` opens a FASTA and its index once and fetches many intervals from it.
- `RecordReader` streams whole records from any `std.Io.Reader`, with no index.
- `validation` checks a FASTA stream in one pass and can write a repaired copy in the same pass, with no index.

Choose by the work: stream when you read every record once, and open an index when you fetch selected records or intervals repeatedly. A program can do both at once, and streaming never needs an index.

## Add the dependency

The package needs Zig 0.16.0 or later. Fetch a release and add the module in `build.zig`:

```bash
zig fetch --save git+https://github.com/eneskemalergin/z-fasta
```

```zig
const z_fasta = b.dependency("z_fasta", .{ .target = target, .optimize = optimize });
exe.root_module.addImport("z-fasta", z_fasta.module("z-fasta"));
```

Then `const zfasta = @import("z-fasta");`. Every optimize mode builds the package; only the z-fasta executable itself rejects ReleaseSmall.

## Failures

Every call that can fail takes a caller-owned `zfasta.diagnostics.Failure` and returns a `zfasta.diagnostics.Error`. Nothing prints or exits. Several codes share one error tag, so read the exact reason from the failure:

```zig
var failure: diagnostics.Failure = undefined;
try testing.expectError(error.InvalidRequest, reader.fetchInto(0, .{ .start = 3, .end = 2 }, &buf, &failure));
try testing.expectEqual(diagnostics.Code.@"request.end_before_start", failure.primary().?.code());
```

[Error diagnostics](Error-Diagnostics) lists every code and its fields.

## Coordinates

The package counts from zero; diagnostics and the CLI count from one.

- An ordinal is a record's position in the index, from 0, and it reaches every record, duplicates included.
- An interval is zero-based and half-open, `[start, end)`, and is never clipped: an end past the record fails with `request.end_after_sequence`. The CLI region `chr1:5-9` is the interval `.{ .start = 4, .end = 9 }`.
- In a failure, `record`, `line`, and `column` are one-based, as the CLI prints them: the record at ordinal 0 is reported as record 1.

## Fetch intervals with `IndexedReader`

Build the index once with the CLI (`z-fasta index ref.fa`), then:

```zig
var failure: diagnostics.Failure = undefined;
var reader = try zfasta.IndexedReader.open(allocator, io, fasta, .{}, &failure);
defer reader.deinit();
try testing.expectEqual(@as(usize, 2), reader.recordCount());
const a = (try reader.lookup("a", &failure)).?;
try testing.expectEqualStrings("a", a.name);
try testing.expectEqual(@as(u64, 12), a.length);

var buf: [16]u8 = undefined;
try testing.expectEqualStrings("ACGTA", try reader.fetchInto(a.ordinal, .{ .start = 4, .end = 9 }, &buf, &failure));
try testing.expect(try reader.lookup("missing", &failure) == null);
try testing.expect(try reader.record(2, &failure) == null);
```

- `open` takes `.{ .index = .automatic }` by default, which prefers `ref.fa.zfi` and then `ref.fa.fai`, as `get` does. `.{ .index = .{ .fai = path } }` or `.{ .zfi = path }` opens exactly that file and never falls back.
- `open` maps a `.zfi` and checks its header, source identity, and layout, so it costs the same for any index size; `record`, `lookup`, and `fetchInto` check each record, side table, table entry, and name the first time they read it, and fail with its index diagnostic. A `.fai` is mapped and its names read into a name map at `open`; each line is parsed and checked the first time it is used. Run `z-fasta validate --index` to check a whole index.
- `record(ordinal, &failure)` and `lookup(name, &failure)` read no sequence and return null past the last record or for an unknown name. The returned `name` borrows the reader until `deinit`. `lookup` returns the first record with that name; later duplicates are reachable by ordinal.
- `fetchInto` checks the ordinal, the interval, and that `dst` holds `end - start` bytes before it writes anything; a read failure after that can leave a prefix. It returns the slice of `dst` it filled.
- A reader keeps the files open until `deinit`. It is not safe to share between threads, and its `io` must outlive it. A gzip input fails with `input.unsupported_format`.

## Stream records with `RecordReader`

```zig
var source: std.Io.Reader = .fixed(input);
var storage: [512]u8 = undefined;
var failure: diagnostics.Failure = undefined;
var records = zfasta.RecordReader.init(&source, &storage, &failure);
const first = (try records.next()).?;
try testing.expectEqualStrings("sp|P99999|TEST_HUMAN", first.id);
try testing.expectEqualStrings("Test protein one OS=Homo sapiens OX=9606 GN=TST1 PE=1 SV=2", first.description);
try testing.expectEqualStrings(line1 ++ line2 ++ line3, first.sequence);
```

- `storage` holds one record: its header line and its bases. A record that does not fit fails with `resource.line_buffer_exceeded`, which reports how many bytes it needed at least; nothing is allocated.
- A record borrows `storage` until the next call. After a failure, every later call returns the same error.
- Records follow `index`'s rules: the identifier ends at the first space, tab, CR, or LF; whitespace and line ends are dropped from the sequence; a record without bases is skipped; duplicates are all returned. Bytes before the first header, a BOM included, fail with `fasta.content_before_header`.

## Validate and repair with `validation`

```zig
const report = try testing.allocator.create(zfasta.validation.Report);
defer testing.allocator.destroy(report);
var buffer: [16]u8 = undefined;
var failure: diagnostics.Failure = undefined;

var reader = std.Io.Reader.fixed(data);
try zfasta.validation.validate(testing.allocator, &reader, &buffer, null, .{}, report, &failure);
try testing.expectEqual(@as(u64, 1), report.finding(.@"fasta.duplicate_identifier").occurrences);
try testing.expectEqual(@as(u8, 1), report.status(false, false));

reader = std.Io.Reader.fixed(data);
var output: std.Io.Writer.Allocating = .init(testing.allocator);
defer output.deinit();
try zfasta.validation.repair(testing.allocator, &reader, &buffer, null, .{}, 4, &output.writer, report, &failure);
try testing.expectEqualStrings(">a desc\nACGT\nACG\n>a\nAC?\n", output.written());
```

- The rules, counts, and repair are the CLI's; see [Validation and repair](Validation-and-Repair). `validation.Rules` selects the alphabet, header schema and limit, and the duplicate check's memory budget; `.{}` is the CLI's default.
- `buffer` is the read window; the CLI uses `validation.READ_BUFFER_BYTES` (256 KiB). Results do not depend on its size.
- `Report` is a fixed-size value of about 16 KiB, so allocate it rather than placing it on a small stack. `report.finding(code)` gives a code's occurrences, records, and first occurrence; `report.status(strict, repaired)` gives the CLI's exit status.
- Pass a `validation.NameSource` (the open file and its size) when the input is a file: the duplicate check then keeps a fingerprint and an offset per identifier and re-reads names from the file. With `null`, it copies the identifiers.
- `repair` writes the copy to any writer, wrapped at the given width, in the same pass. On any error, including `repair.blocked`, discard what the writer received.
- The index comparison (`validate --index`) is CLI-only for now.

## Related pages

- [Coordinates and names](Coordinates-and-Sequence-Names)
- [Error diagnostics](Error-Diagnostics)
- [For contributors](For-Contributors)
