//! Sequence extraction: positional regions, BED/names batching, strand, and RC transforms.
//!
//! Uniform records use O(1) byte-offset math; non-uniform records use side-table lookup.
//! A present `.zfi` is authoritative; `.fai` is used only when `.zfi` is absent.

const std = @import("std");
const complement = @import("complement.zig");
const bed_parser = @import("bed_parser.zig");
const index_format = @import("index_format.zig");
const platform = @import("platform.zig");
const stats = @import("stats.zig");

const IndexRecord = index_format.IndexRecord;
const LoadedIndex = index_format.LoadedIndex;
const printErrorAndExit = index_format.printErrorAndExit;

// ============================================================================
// Region parsing
// ============================================================================

pub const Region = struct {
    name: []const u8,
    start: u64, // 1-based inclusive
    end: ?u64, // 1-based inclusive, null = to end of sequence
    is_full: bool, // true if no :START-END specified
};

pub const Orientation = struct {
    reverse: bool = false,
    complement: bool = false,

    pub fn compose(self: Orientation, other: Orientation) Orientation {
        return .{
            .reverse = self.reverse != other.reverse,
            .complement = self.complement != other.complement,
        };
    }

    pub fn isIdentity(self: Orientation) bool {
        return !self.reverse and !self.complement;
    }

    pub fn reverseComplement() Orientation {
        return .{ .reverse = true, .complement = true };
    }

    pub fn reverseOnly() Orientation {
        return .{ .reverse = true, .complement = false };
    }

    pub fn complementOnly() Orientation {
        return .{ .reverse = false, .complement = true };
    }
};

pub const GetOptions = struct {
    region_strs: []const []const u8,
    bed_path: ?[]const u8 = null,
    names_path: ?[]const u8 = null,
    honor_strand: bool = false,
    summary: bool = false,
    chunk_size: usize = 4_096,
    orientation: Orientation = .{},
    annotate_transform: bool = false,
};

pub const chunk_size_all = std.math.maxInt(usize);

/// Maximum size for inputs loaded fully into memory: `--names` (always) and
/// BED with `--chunk-size -1`. Default BED chunking is not size-capped here.
pub const max_input_file_mib: usize = 512;
pub const max_input_file_bytes: usize = max_input_file_mib * 1024 * 1024;

comptime {
    // Keep get help / README "512 MiB" wording synchronized with this constant.
    if (max_input_file_mib != 512) {
        @compileError("update get help and README for max_input_file_mib");
    }
}

const AllInMemoryInputKind = enum {
    names_file,
    bed_all_in_memory,
};

/// Per-line buffer for chunked BED streaming (`takeDelimiter` limit).
const bed_line_reader_buffer_bytes = 4096;

/// Per-region output cap for the multi-region sort buffer path (≥16 regions).
const max_sort_path_region_output_bytes: u64 = 64 * 1024 * 1024;

/// Total intermediate output cap for the multi-region sort buffer path.
const max_sort_path_total_output_bytes: u64 = 256 * 1024 * 1024;

const ParsedRequest = struct {
    region: Region,
    orientation: Orientation,
};

const BatchStats = struct {
    region_count: usize = 0,
    total_bases: u64 = 0,
};

/// Parse a region string. Handles the Ensembl colon trap by parsing from the right.
/// Accepted formats:
///   NAME              full sequence
///   NAME:START-END    1-based, inclusive
///   NAME:START-       from START to end
pub fn parseRegion(input: []const u8) Region {
    // Only the rightmost colon can start a coordinate suffix. Any earlier suffix
    // contains that colon, so its coordinate fields cannot both be integers.
    var colon_pos: ?usize = null;
    {
        var i: usize = input.len;
        while (i > 0) {
            i -= 1;
            if (input[i] == ':') {
                colon_pos = i;
                break;
            }
        }
    }

    if (colon_pos) |cp| {
        const suffix = input[cp + 1 ..];

        if (parseRangeSuffix(suffix)) |range| {
            return Region{
                .name = input[0..cp],
                .start = range.start,
                .end = range.end,
                .is_full = false,
            };
        }
    }

    // No valid region suffix found; treat entire input as sequence name
    return Region{
        .name = input,
        .start = 1,
        .end = null,
        .is_full = true,
    };
}

const RangeParsed = struct {
    start: u64,
    end: ?u64,
};

fn parseRangeSuffix(suffix: []const u8) ?RangeParsed {
    // Expect: START-END or START-
    // Find the '-' separator
    const dash_pos = std.mem.indexOfScalar(u8, suffix, '-') orelse return null;

    const start_str = suffix[0..dash_pos];
    const end_str = suffix[dash_pos + 1 ..];

    // START must be a valid positive integer
    if (start_str.len == 0) return null;
    const start = std.fmt.parseInt(u64, start_str, 10) catch return null;

    if (end_str.len == 0) {
        // NAME:START- form (to end of sequence)
        return RangeParsed{ .start = start, .end = null };
    }

    const end = std.fmt.parseInt(u64, end_str, 10) catch return null;
    return RangeParsed{ .start = start, .end = end };
}

// ============================================================================
// Region resolution
// ============================================================================

/// A fully resolved, validated extraction request.
/// Byte offset and length are pre-computed so extraction is a single mmap
/// slice walk; no further index lookups required.
pub const ResolvedRegion = struct {
    name: []const u8, // sequence name for FASTA header
    start: u64, // 1-based inclusive (validated)
    display_end: u64, // end value for header (pre-clamp, samtools convention)
    is_full: bool, // true -> emit ">NAME", false -> emit ">NAME:start-display_end"
    start_byte: u64, // absolute byte offset into fasta_data for first base
    seq_offset: u64,
    num_bases: u64, // number of bases to extract
    line_bases: u32,
    line_bytes: u32,
    side_table: []const index_format.SideTableLine,
    orientation: Orientation,
    annotate_transform: bool,
    original_index: usize, // position in CLI argument list (preserves output order)
};

/// Resolve one region string against a loaded index.
/// Validates coordinates and pre-computes the O(1) byte offset.
/// Calls printErrorAndExit on any error (sequence not found, bad coordinates).
pub fn resolveRegion(idx: *const LoadedIndex, region_str: []const u8, original_index: usize) ResolvedRegion {
    const region = parseRegion(region_str);

    const rec_idx = idx.lookupName(region.name) orelse {
        printErrorAndExit("error: sequence not found: {s}\n", .{region.name});
    };

    return resolveParsedRequest(
        idx,
        .{ .region = region, .orientation = .{} },
        rec_idx,
        original_index,
        false,
    );
}

fn resolveParsedRequest(
    idx: *const LoadedIndex,
    request: ParsedRequest,
    rec_idx: usize,
    original_index: usize,
    annotate_transform: bool,
) ResolvedRegion {
    const rec = idx.records[rec_idx];
    const region = request.region;

    var start = region.start;
    var end = region.end orelse rec.seq_len;
    const display_end = end; // capture before clamping (samtools keeps unclamped in header)

    if (region.is_full) {
        start = 1;
        end = rec.seq_len;
    }

    if (start < 1) {
        printErrorAndExit("error: start position must be >= 1\n", .{});
    }
    if (start > rec.seq_len) {
        printErrorAndExit("error: start position {d} exceeds sequence length {d}\n", .{ start, rec.seq_len });
    }
    if (end < start) {
        printErrorAndExit("error: end position must be >= start position\n", .{});
    }

    // Clamp end to seq_len silently (samtools behavior)
    if (end > rec.seq_len) {
        end = rec.seq_len;
    }

    const num_bases = end - start + 1;

    const side_table = idx.sideTableLines(rec);
    const start_byte = byteOffsetForBase(idx.fasta_data, rec, side_table, start - 1);

    return ResolvedRegion{
        .name = region.name,
        .start = start,
        .display_end = display_end,
        .is_full = region.is_full,
        .start_byte = start_byte,
        .seq_offset = rec.seq_offset,
        .num_bases = num_bases,
        .line_bases = rec.line_bases,
        .line_bytes = rec.line_bytes,
        .side_table = side_table,
        .orientation = request.orientation,
        .annotate_transform = annotate_transform,
        .original_index = original_index,
    };
}

fn byteOffsetForBase(
    fasta: []const u8,
    rec: IndexRecord,
    side_table: []const index_format.SideTableLine,
    base_index: u64,
) u64 {
    if (side_table.len == 0) {
        const line_number = base_index / rec.line_bases;
        const column = base_index % rec.line_bases;
        return rec.seq_offset + (line_number * rec.line_bytes) + column;
    }

    var lo: usize = 0;
    var hi: usize = side_table.len;
    while (lo < hi) {
        const mid = lo + (hi - lo) / 2;
        if (side_table[mid].base_start <= base_index) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    const line_idx = if (lo == 0) 0 else lo - 1;
    const line = side_table[line_idx];
    const column = base_index - line.base_start;

    var seen: u64 = 0;
    var pos: usize = @intCast(line.byte_offset);
    const end = @min(fasta.len, pos + @as(usize, @intCast(line.line_bytes)));
    while (pos < end) : (pos += 1) {
        if (fasta[pos] <= ' ') continue;
        if (seen == column) return @intCast(pos);
        seen += 1;
    }

    printErrorAndExit("error: corrupt non-uniform index side table\n", .{});
}

/// Classify a record as nucleotide or protein from a short sequence prefix.
/// Uses `stats.get_type_sample_bases` (not full-file): protein vs IUPAC nucleotide
/// is clear within a few hundred bases, and sampling up to 100k on every `--rc`
/// path was a large share of Genome RC overhead.
fn detectRecordType(rec: IndexRecord, fasta: []const u8) stats.SequenceType {
    var counts = [_]u64{0} ** 256;
    var total: u64 = 0;
    var pos: usize = @intCast(rec.seq_offset);
    const sample_limit: u64 = @min(rec.seq_len, stats.get_type_sample_bases);

    while (pos < fasta.len and total < sample_limit) : (pos += 1) {
        const byte = fasta[pos];
        if (byte <= ' ') continue;
        counts[byte] += 1;
        total += 1;
    }

    return stats.detectType(&counts, total);
}

fn ensureComplementAllowed(idx: *const LoadedIndex, requests: []const ParsedRequest) void {
    var last_name: ?[]const u8 = null;
    var last_rec_idx: usize = 0;
    var last_checked_rec_idx: ?usize = null;
    var last_checked_type: stats.SequenceType = undefined;

    for (requests) |request| {
        if (!request.orientation.complement) continue;

        const rec_idx = if (last_name) |name|
            if (std.mem.eql(u8, name, request.region.name))
                last_rec_idx
            else
                idx.lookupName(request.region.name) orelse {
                    printErrorAndExit("error: sequence not found: {s}\n", .{request.region.name});
                }
        else
            idx.lookupName(request.region.name) orelse {
                printErrorAndExit("error: sequence not found: {s}\n", .{request.region.name});
            };

        last_name = request.region.name;
        last_rec_idx = rec_idx;

        const rec_type = if (last_checked_rec_idx) |cached_rec_idx|
            if (cached_rec_idx == rec_idx)
                last_checked_type
            else blk: {
                const detected = detectRecordType(idx.records[rec_idx], idx.fasta_data);
                last_checked_rec_idx = rec_idx;
                last_checked_type = detected;
                break :blk detected;
            }
        else blk: {
            const detected = detectRecordType(idx.records[rec_idx], idx.fasta_data);
            last_checked_rec_idx = rec_idx;
            last_checked_type = detected;
            break :blk detected;
        };

        if (rec_type == .protein) {
            printErrorAndExit(
                "error: reverse complement is not defined for protein sequences: {s} (classified from up to {d} bases)\n",
                .{ request.region.name, stats.get_type_sample_bases },
            );
        }
    }
}

// ============================================================================
// Sequence emission
// ============================================================================

/// Writes one region's FASTA output and releases cached FASTA pages for the span read.
pub fn extractRegion(idx: *const LoadedIndex, region_str: []const u8, writer: anytype) void {
    const resolved = resolveRegion(idx, region_str, 0);
    emitRegionAndRelease(resolved, idx.fasta_data, writer);
}

fn emitRegionAndRelease(resolved: ResolvedRegion, fasta: []const u8, writer: anytype) void {
    emitRegion(resolved, fasta, writer);
    const span = fastaSpanForRegion(resolved, fasta);
    index_format.dropFastaSpan(fasta, span.start, span.end);
}

/// Batch sequential FASTA page releases to balance RSS vs `madvise` syscall cost.
const fasta_release_batch_bytes: usize = 8 * 1024 * 1024;

/// Release FASTA mmap pages during sequential extraction; batches `madvise` to limit syscall overhead.
/// Only safe when regions are visited in file order with small gaps (dense catalogs).
const FastaReleaseCursor = struct {
    released_end: usize = 0,
    /// Highest byte offset covered by a completed region emit (may lag `released_end` while batching).
    scan_end: usize = 0,

    fn beforeRegion(self: *FastaReleaseCursor, fasta: []const u8, start_byte: usize) void {
        if (start_byte <= self.released_end) return;
        const pending = start_byte - self.released_end;
        if (pending >= fasta_release_batch_bytes) {
            index_format.dropFastaSpan(fasta, self.released_end, start_byte);
            self.released_end = start_byte;
        }
    }

    fn afterRegion(self: *FastaReleaseCursor, span_end: usize) void {
        if (span_end > self.scan_end) self.scan_end = span_end;
    }

    fn flush(self: *FastaReleaseCursor, fasta: []const u8) void {
        const drop_end = @max(self.scan_end, self.released_end);
        if (drop_end > self.released_end) {
            index_format.dropFastaSpan(fasta, self.released_end, drop_end);
            self.released_end = drop_end;
        }
    }
};

fn emitRegionAndReleaseSequential(
    resolved: ResolvedRegion,
    fasta: []const u8,
    writer: anytype,
    release: *FastaReleaseCursor,
) void {
    const start_byte: usize = @intCast(resolved.start_byte);
    release.beforeRegion(fasta, start_byte);
    emitRegion(resolved, fasta, writer);
    const span = fastaSpanForRegion(resolved, fasta);
    release.afterRegion(@max(span.end, start_byte + 1));
}

fn fastaSpanForRegion(resolved: ResolvedRegion, fasta: []const u8) struct { start: usize, end: usize } {
    const rec = IndexRecord{
        .name_offset = 0,
        .name_len = 0,
        .seq_offset = resolved.seq_offset,
        .seq_len = resolved.start - 1 + resolved.num_bases,
        .line_bases = resolved.line_bases,
        .line_bytes = resolved.line_bytes,
    };
    const start: usize = @intCast(resolved.start_byte);
    const last_base_index = resolved.start - 1 + resolved.num_bases - 1;
    const last_byte = byteOffsetForBase(fasta, rec, resolved.side_table, last_base_index);
    const end: usize = @min(fasta.len, @as(usize, @intCast(last_byte)) + resolved.line_bytes);
    return .{ .start = start, .end = @max(start + 1, end) };
}

/// Write FASTA output for one resolved region to `writer`.
/// Output is wrapped at 60 bases per line (samtools default).
fn emitRegion(resolved: ResolvedRegion, fasta: []const u8, writer: anytype) void {
    const annotation = headerAnnotation(resolved.orientation, resolved.annotate_transform);
    if (resolved.is_full) {
        writer.print(">{s}{s}\n", .{ resolved.name, annotation }) catch {
            printErrorAndExit("error: write failed\n", .{});
        };
    } else {
        writer.print(">{s}:{d}-{d}{s}\n", .{ resolved.name, resolved.start, resolved.display_end, annotation }) catch {
            printErrorAndExit("error: write failed\n", .{});
        };
    }

    if (resolved.orientation.reverse) {
        if (resolved.side_table.len == 0 and resolved.line_bases > 0) {
            emitRegionBackwardUniform(resolved, fasta, writer);
        } else {
            emitRegionBackward(resolved, fasta, writer);
        }
        return;
    }

    emitRegionForward(resolved, fasta, writer);
}

/// Max on-disk span for sparse positional emit. Larger / reverse / side-table
/// regions fall back to mmap.
const max_pread_span_bytes: usize = 64 * 1024;

/// Borrowed loader handle for sparse large-FASTA reads via portable positional I/O.
const SparseFastaSource = struct {
    io: std.Io,
    file: std.Io.File,
};

fn canEmitRegionViaPread(resolved: ResolvedRegion, fasta_len: usize) bool {
    if (resolved.orientation.reverse) return false;
    if (resolved.side_table.len != 0 or resolved.line_bases == 0) return false;
    const start: usize = @intCast(resolved.start_byte);
    const last_base = resolved.start - 1 + resolved.num_bases - 1;
    const last_line = last_base / resolved.line_bases;
    const end_exclusive: usize = @min(
        fasta_len,
        @as(usize, @intCast(resolved.seq_offset + (last_line + 1) * resolved.line_bytes)),
    );
    return end_exclusive > start and (end_exclusive - start) <= max_pread_span_bytes;
}

// Positional reads keep sparse large-FASTA requests from faulting mapped pages.
fn emitRegionViaPread(
    resolved: ResolvedRegion,
    sparse: SparseFastaSource,
    fasta_len: usize,
    file_buf: []u8,
    writer: anytype,
) void {
    const annotation = headerAnnotation(resolved.orientation, resolved.annotate_transform);
    if (resolved.is_full) {
        writer.print(">{s}{s}\n", .{ resolved.name, annotation }) catch {
            printErrorAndExit("error: write failed\n", .{});
        };
    } else {
        writer.print(">{s}:{d}-{d}{s}\n", .{ resolved.name, resolved.start, resolved.display_end, annotation }) catch {
            printErrorAndExit("error: write failed\n", .{});
        };
    }

    const start: usize = @intCast(resolved.start_byte);
    const line_bases: u64 = resolved.line_bases;
    const line_bytes: u64 = resolved.line_bytes;
    const first_base = resolved.start - 1;
    const last_base = first_base + resolved.num_bases - 1;
    const first_line = first_base / line_bases;
    const first_column = first_base % line_bases;
    const last_line = last_base / line_bases;
    const end_exclusive: usize = @min(
        fasta_len,
        @as(usize, @intCast(resolved.seq_offset + (last_line + 1) * line_bytes)),
    );
    const span_len = end_exclusive - start;
    if (span_len > file_buf.len) {
        printErrorAndExit("error: internal: positional span exceeds buffer\n", .{});
    }

    const got = std.Io.File.readPositionalAll(sparse.file, sparse.io, file_buf[0..span_len], start) catch {
        printErrorAndExit("error: failed to read FASTA\n", .{});
    };
    if (got != span_len) {
        printErrorAndExit("error: read past end of FASTA\n", .{});
    }
    const slice = file_buf[0..span_len];

    const wrap_width: usize = 60;
    var bases_written: u64 = 0;
    var line_pos: usize = 0;
    var out_buf: [65536]u8 = undefined;
    var out_len: usize = 0;
    var base_index: u64 = first_base;

    while (bases_written < resolved.num_bases) {
        const line_number = base_index / line_bases;
        const column = base_index % line_bases;
        const src_start: usize = @intCast(
            (line_number - first_line) * line_bytes + column - first_column,
        );
        const available = line_bases - column;
        const take_u64 = @min(available, resolved.num_bases - bases_written);
        const take: usize = @intCast(take_u64);
        if (src_start + take > slice.len) {
            printErrorAndExit("error: read past end of FASTA\n", .{});
        }
        const src = slice[src_start .. src_start + take];

        var i: usize = 0;
        while (i < take) {
            const chunk = @min(take - i, wrap_width - line_pos);
            const needed = chunk + @as(usize, @intFromBool(line_pos + chunk == wrap_width));
            if (out_buf.len - out_len < needed) {
                writer.writeAll(out_buf[0..out_len]) catch {
                    printErrorAndExit("error: write failed\n", .{});
                };
                out_len = 0;
            }
            if (resolved.orientation.complement) {
                complement.complementInto(out_buf[out_len .. out_len + chunk], src[i .. i + chunk]);
            } else {
                @memcpy(out_buf[out_len .. out_len + chunk], src[i .. i + chunk]);
            }
            out_len += chunk;
            i += chunk;
            line_pos += chunk;
            if (line_pos >= wrap_width) {
                out_buf[out_len] = '\n';
                out_len += 1;
                line_pos = 0;
            }
        }

        bases_written += take_u64;
        base_index += take_u64;
    }

    if (line_pos > 0) {
        if (out_len == out_buf.len) {
            writer.writeAll(&out_buf) catch {
                printErrorAndExit("error: write failed\n", .{});
            };
            out_len = 0;
        }
        out_buf[out_len] = '\n';
        out_len += 1;
    }
    if (out_len > 0) {
        writer.writeAll(out_buf[0..out_len]) catch {
            printErrorAndExit("error: write failed\n", .{});
        };
    }
}

fn headerAnnotation(orientation: Orientation, annotate_transform: bool) []const u8 {
    if (!annotate_transform or orientation.isIdentity()) return "";
    if (orientation.reverse and orientation.complement) return " (reverse complement)";
    if (orientation.reverse) return " (reverse)";
    if (orientation.complement) return " (complement)";
    return "";
}

fn emitRegionForward(resolved: ResolvedRegion, fasta: []const u8, writer: anytype) void {
    // Uniform records: copy whole line runs instead of per-byte whitespace scans.
    // Messy (side-table) and reverse paths keep the byte walker.
    if (resolved.side_table.len == 0 and resolved.line_bases > 0) {
        emitRegionForwardUniform(resolved, fasta, writer);
        return;
    }
    emitRegionForwardScan(resolved, fasta, writer);
}

fn emitRegionForwardUniform(resolved: ResolvedRegion, fasta: []const u8, writer: anytype) void {
    const wrap_width: usize = 60;
    const line_bases: u64 = resolved.line_bases;
    const line_bytes: u64 = resolved.line_bytes;
    var bases_written: u64 = 0;
    var line_pos: usize = 0;
    var out_buf: [65536]u8 = undefined;
    var out_len: usize = 0;
    var base_index: u64 = resolved.start - 1;

    while (bases_written < resolved.num_bases) {
        const line_number = base_index / line_bases;
        const column = base_index % line_bases;
        const line_start: usize = @intCast(resolved.seq_offset + line_number * line_bytes);
        const available = line_bases - column;
        const take_u64 = @min(available, resolved.num_bases - bases_written);
        const take: usize = @intCast(take_u64);
        const src_start = line_start + @as(usize, @intCast(column));
        if (src_start + take > fasta.len) {
            printErrorAndExit("error: read past end of FASTA\n", .{});
        }
        const src = fasta[src_start .. src_start + take];

        var i: usize = 0;
        while (i < take) {
            const chunk = @min(take - i, wrap_width - line_pos);
            const needed = chunk + @as(usize, @intFromBool(line_pos + chunk == wrap_width));
            if (out_buf.len - out_len < needed) {
                writer.writeAll(out_buf[0..out_len]) catch {
                    printErrorAndExit("error: write failed\n", .{});
                };
                out_len = 0;
            }
            if (resolved.orientation.complement) {
                complement.complementInto(out_buf[out_len .. out_len + chunk], src[i .. i + chunk]);
            } else {
                @memcpy(out_buf[out_len .. out_len + chunk], src[i .. i + chunk]);
            }
            out_len += chunk;
            i += chunk;
            line_pos += chunk;
            if (line_pos >= wrap_width) {
                out_buf[out_len] = '\n';
                out_len += 1;
                line_pos = 0;
            }
        }

        bases_written += take_u64;
        base_index += take_u64;
    }

    if (line_pos > 0) {
        if (out_len == out_buf.len) {
            writer.writeAll(&out_buf) catch {
                printErrorAndExit("error: write failed\n", .{});
            };
            out_len = 0;
        }
        out_buf[out_len] = '\n';
        out_len += 1;
    }

    if (out_len > 0) {
        writer.writeAll(out_buf[0..out_len]) catch {
            printErrorAndExit("error: write failed\n", .{});
        };
    }
}

fn emitRegionForwardScan(resolved: ResolvedRegion, fasta: []const u8, writer: anytype) void {
    const wrap_width: usize = 60;
    var pos: usize = @intCast(resolved.start_byte);
    var bases_written: u64 = 0;
    var line_pos: usize = 0;
    var out_buf: [65536]u8 = undefined;
    var out_len: usize = 0;

    while (bases_written < resolved.num_bases and pos < fasta.len) {
        const byte = fasta[pos];
        pos += 1;

        if (byte <= ' ') continue;

        if (out_len + 2 > out_buf.len) {
            writer.writeAll(out_buf[0..out_len]) catch {
                printErrorAndExit("error: write failed\n", .{});
            };
            out_len = 0;
        }

        out_buf[out_len] = if (resolved.orientation.complement) complement.complement(byte) else byte;
        out_len += 1;
        bases_written += 1;
        line_pos += 1;

        if (line_pos >= wrap_width) {
            out_buf[out_len] = '\n';
            out_len += 1;
            line_pos = 0;
        }
    }

    if (line_pos > 0) {
        if (out_len == out_buf.len) {
            writer.writeAll(&out_buf) catch {
                printErrorAndExit("error: write failed\n", .{});
            };
            out_len = 0;
        }
        out_buf[out_len] = '\n';
        out_len += 1;
    }

    if (out_len > 0) {
        writer.writeAll(out_buf[0..out_len]) catch {
            printErrorAndExit("error: write failed\n", .{});
        };
    }
}

/// Uniform reverse / RC: walk whole line runs backward (same geometry as forward uniform).
fn emitRegionBackwardUniform(resolved: ResolvedRegion, fasta: []const u8, writer: anytype) void {
    const wrap_width: usize = 60;
    const line_bases: u64 = resolved.line_bases;
    const line_bytes: u64 = resolved.line_bytes;
    var bases_written: u64 = 0;
    var line_pos: usize = 0;
    var out_buf: [65536]u8 = undefined;
    var out_len: usize = 0;
    // Emit from last base toward first (reverse order).
    var base_index: u64 = resolved.start - 1 + resolved.num_bases - 1;

    while (bases_written < resolved.num_bases) {
        const line_number = base_index / line_bases;
        const column = base_index % line_bases;
        const line_start: usize = @intCast(resolved.seq_offset + line_number * line_bytes);
        // How many bases we can take walking left on this line (inclusive of column).
        const available = column + 1;
        const take_u64 = @min(available, resolved.num_bases - bases_written);
        const take: usize = @intCast(take_u64);
        const src_start = line_start + @as(usize, @intCast(column + 1 - take_u64));
        if (src_start + take > fasta.len) {
            printErrorAndExit("error: read past end of FASTA\n", .{});
        }
        const src = fasta[src_start .. src_start + take];

        var i: usize = 0;
        while (i < take) {
            const chunk = @min(take - i, wrap_width - line_pos);
            const needed = chunk + @as(usize, @intFromBool(line_pos + chunk == wrap_width));
            if (out_buf.len - out_len < needed) {
                writer.writeAll(out_buf[0..out_len]) catch {
                    printErrorAndExit("error: write failed\n", .{});
                };
                out_len = 0;
            }
            // Walk src from the right end of this take window.
            const src_end = take - i;
            if (resolved.orientation.complement) {
                var j: usize = 0;
                while (j < chunk) : (j += 1) {
                    out_buf[out_len + j] = complement.complement(src[src_end - 1 - j]);
                }
            } else {
                var j: usize = 0;
                while (j < chunk) : (j += 1) {
                    out_buf[out_len + j] = src[src_end - 1 - j];
                }
            }
            out_len += chunk;
            i += chunk;
            line_pos += chunk;
            if (line_pos >= wrap_width) {
                out_buf[out_len] = '\n';
                out_len += 1;
                line_pos = 0;
            }
        }

        bases_written += take_u64;
        if (bases_written >= resolved.num_bases) break;
        base_index -= take_u64;
    }

    if (line_pos > 0) {
        if (out_len == out_buf.len) {
            writer.writeAll(&out_buf) catch {
                printErrorAndExit("error: write failed\n", .{});
            };
            out_len = 0;
        }
        out_buf[out_len] = '\n';
        out_len += 1;
    }
    if (out_len > 0) {
        writer.writeAll(out_buf[0..out_len]) catch {
            printErrorAndExit("error: write failed\n", .{});
        };
    }
}

fn emitRegionBackward(resolved: ResolvedRegion, fasta: []const u8, writer: anytype) void {
    const wrap_width: usize = 60;
    var bases_remaining = resolved.num_bases;
    var line_pos: usize = 0;
    var out_buf: [65536]u8 = undefined;
    var out_len: usize = 0;

    const last_base_index = resolved.start - 1 + resolved.num_bases - 1;
    const rec = IndexRecord{
        .name_offset = 1,
        .name_len = 1,
        .seq_offset = resolved.seq_offset,
        .seq_len = resolved.start - 1 + resolved.num_bases,
        .line_bases = resolved.line_bases,
        .line_bytes = resolved.line_bytes,
    };
    var pos: usize = @intCast(byteOffsetForBase(fasta, rec, resolved.side_table, last_base_index));

    while (bases_remaining > 0) {
        const byte = if (resolved.orientation.complement) complement.complement(fasta[pos]) else fasta[pos];

        if (out_len + 2 > out_buf.len) {
            writer.writeAll(out_buf[0..out_len]) catch {
                printErrorAndExit("error: write failed\n", .{});
            };
            out_len = 0;
        }

        out_buf[out_len] = byte;
        out_len += 1;
        bases_remaining -= 1;
        line_pos += 1;

        if (line_pos >= wrap_width) {
            out_buf[out_len] = '\n';
            out_len += 1;
            line_pos = 0;
        }

        if (bases_remaining > 0) {
            while (pos > 0) {
                pos -= 1;
                const prev = fasta[pos];
                if (prev > ' ') break;
            }
        }
    }

    if (line_pos > 0) {
        if (out_len == out_buf.len) {
            writer.writeAll(&out_buf) catch {
                printErrorAndExit("error: write failed\n", .{});
            };
            out_len = 0;
        }
        out_buf[out_len] = '\n';
        out_len += 1;
    }

    if (out_len > 0) {
        writer.writeAll(out_buf[0..out_len]) catch {
            printErrorAndExit("error: write failed\n", .{});
        };
    }
}

fn monotonicNs(io: std.Io) u64 {
    const now = std.Io.Clock.Timestamp.now(io, .awake);
    return @intCast(now.raw.toNanoseconds());
}

fn estimateRegionOutputBytes(resolved: ResolvedRegion) u64 {
    const wrap_lines = resolved.num_bases / 60 + @as(u64, @intFromBool(resolved.num_bases % 60 != 0));
    const header_len: u64 = if (resolved.is_full)
        @intCast(resolved.name.len + 1)
    else
        @intCast(resolved.name.len + 32);
    const annotation_len: u64 = if (resolved.annotate_transform) 24 else 0;
    return resolved.num_bases + wrap_lines + header_len + annotation_len;
}

/// Drop FASTA cache after a sparse batch on huge files to cap peak RSS.
const sparse_large_fasta_bytes: u64 = 256 * 1024 * 1024;
/// Few contigs (e.g. GRCh38): BED rows are sparse; file-order sort walks huge gaps.
const sparse_catalog_record_threshold: usize = 512;
/// Median byte gap above this → emit in request order with direct seeks.
const sparse_median_gap_bytes: u64 = 512 * 1024;
/// File-order sort only helps on FASTA large enough that random BED seeks lose.
const sort_by_offset_min_fasta_bytes: u64 = 64 * 1024 * 1024;

fn medianStartByteGap(allocator: std.mem.Allocator, resolved: []const ResolvedRegion) !u64 {
    if (resolved.len < 2) return 0;
    const starts = try allocator.alloc(u64, resolved.len);
    defer allocator.free(starts);
    for (resolved, 0..) |r, i| starts[i] = r.start_byte;
    std.mem.sort(u64, starts, {}, std.sort.asc(u64));
    const gap_count = starts.len - 1;
    var gaps = try allocator.alloc(u64, gap_count);
    defer allocator.free(gaps);
    for (0..gap_count) |i| gaps[i] = starts[i + 1] - starts[i];
    std.mem.sort(u64, gaps, {}, std.sort.asc(u64));
    return gaps[gap_count / 2];
}

fn shouldSortByFileOffset(
    idx: *const index_format.LoadedIndex,
    allocator: std.mem.Allocator,
    resolved: []const ResolvedRegion,
) !bool {
    return shouldSortByFileOffsetForBatch(
        idx.records.len,
        idx.fasta_size,
        allocator,
        resolved,
    );
}

fn shouldSortByFileOffsetForBatch(
    record_count: usize,
    fasta_size: u64,
    allocator: std.mem.Allocator,
    resolved: []const ResolvedRegion,
) !bool {
    if (resolved.len < 16) return false;
    if (record_count <= sparse_catalog_record_threshold) return false;
    if (fasta_size < sort_by_offset_min_fasta_bytes) return false;
    const gap = try medianStartByteGap(allocator, resolved);
    return gap <= sparse_median_gap_bytes;
}

fn shouldUseSortBuffers(resolved: []const ResolvedRegion) bool {
    var total: u64 = 0;
    for (resolved) |r| {
        const est = estimateRegionOutputBytes(r);
        if (est > max_sort_path_region_output_bytes) return false;
        total += est;
        if (total > max_sort_path_total_output_bytes) return false;
    }
    return true;
}

fn batchHasReverseReads(resolved: []const ResolvedRegion) bool {
    for (resolved) |r| {
        if (r.orientation.reverse) return true;
    }
    return false;
}

fn shouldAdviseSequentialMmap(
    requests_len: usize,
    sequential_scan: bool,
    has_reverse_reads: bool,
    allow_sequential_madvise: bool,
) bool {
    if (!allow_sequential_madvise) return false;
    if (requests_len < 16) return false;
    if (has_reverse_reads) return false;
    return sequential_scan;
}

fn readAllInput(allocator: std.mem.Allocator, io: std.Io, path: []const u8, kind: AllInMemoryInputKind) []u8 {
    if (std.mem.eql(u8, path, "-")) {
        var stdin_buf: [4096]u8 = undefined;
        var reader = std.Io.File.stdin().reader(io, &stdin_buf);
        return reader.interface.allocRemaining(allocator, .limited(max_input_file_bytes)) catch |err| switch (err) {
            error.StreamTooLong => switch (kind) {
                .names_file => printErrorAndExit(
                    "error: --names stdin exceeds {d} MiB limit (--chunk-size does not stream --names)\n",
                    .{max_input_file_mib},
                ),
                .bed_all_in_memory => printErrorAndExit(
                    "error: BED stdin exceeds {d} MiB limit for --chunk-size -1; use default --chunk-size\n",
                    .{max_input_file_mib},
                ),
            },
            else => printErrorAndExit("error: failed to read stdin\n", .{}),
        };
    }

    const file = std.Io.Dir.cwd().openFile(io, path, .{}) catch |err| switch (err) {
        error.FileNotFound => printErrorAndExit("error: file not found: {s}\n", .{path}),
        error.AccessDenied => printErrorAndExit("error: access denied: {s}\n", .{path}),
        else => printErrorAndExit("error: failed to open file: {s}\n", .{path}),
    };
    defer file.close(io);

    const stat = file.stat(io) catch {
        printErrorAndExit("error: failed to stat file: {s}\n", .{path});
    };

    if (stat.size > max_input_file_bytes) {
        switch (kind) {
            .names_file => printErrorAndExit(
                "error: names file exceeds {d} MiB limit: {s} (--chunk-size does not stream --names)\n",
                .{ max_input_file_mib, path },
            ),
            .bed_all_in_memory => printErrorAndExit(
                "error: BED file exceeds {d} MiB limit for --chunk-size -1: {s}; use default --chunk-size\n",
                .{ max_input_file_mib, path },
            ),
        }
    }

    const bytes = allocator.alloc(u8, stat.size) catch {
        printErrorAndExit("error: out of memory\n", .{});
    };

    var file_buf: [4096]u8 = undefined;
    var reader = file.reader(io, &file_buf);
    reader.interface.readSliceAll(bytes) catch {
        printErrorAndExit("error: failed to read file: {s}\n", .{path});
    };
    return bytes;
}

fn appendBedRegionRequest(
    requests: *std.ArrayList(ParsedRequest),
    region: bed_parser.BedRegion,
    honor_strand: bool,
    global_orientation: Orientation,
    allocator: std.mem.Allocator,
    name_allocator: std.mem.Allocator,
    duplicate_name: bool,
    last_duplicated_name: ?*?[]const u8,
) void {
    if (honor_strand and region.strand == .invalid) {
        printErrorAndExit("error: invalid BED line {d}: invalid strand\n", .{region.line_number});
    }

    const name = if (duplicate_name) blk: {
        if (last_duplicated_name) |cached_name| {
            if (cached_name.*) |existing| {
                if (std.mem.eql(u8, existing, region.chrom)) break :blk existing;
            }
        }

        const duplicated = name_allocator.dupe(u8, region.chrom) catch {
            printErrorAndExit("error: out of memory\n", .{});
        };
        if (last_duplicated_name) |cached_name| cached_name.* = duplicated;
        break :blk duplicated;
    } else region.chrom;

    requests.append(allocator, .{
        .region = .{
            .name = name,
            .start = region.start1Based(),
            .end = region.end1BasedInclusive(),
            .is_full = false,
        },
        .orientation = (if (honor_strand and region.strand == .minus) Orientation.reverseComplement() else Orientation{}).compose(global_orientation),
    }) catch {
        printErrorAndExit("error: out of memory\n", .{});
    };
}

fn appendBedLineRequest(
    requests: *std.ArrayList(ParsedRequest),
    line: []const u8,
    line_number: usize,
    honor_strand: bool,
    global_orientation: Orientation,
    allocator: std.mem.Allocator,
    name_allocator: std.mem.Allocator,
    duplicate_name: bool,
    last_duplicated_name: ?*?[]const u8,
) void {
    const parsed = bed_parser.parseBedLine(line, line_number) catch |err| switch (err) {
        error.MissingChrom => printErrorAndExit("error: invalid BED line {d}: missing chrom\n", .{line_number}),
        error.MissingStart => printErrorAndExit("error: invalid BED line {d}: missing start\n", .{line_number}),
        error.MissingEnd => printErrorAndExit("error: invalid BED line {d}: missing end\n", .{line_number}),
        error.InvalidStart => printErrorAndExit("error: invalid BED line {d}: invalid start\n", .{line_number}),
        error.InvalidEnd => printErrorAndExit("error: invalid BED line {d}: invalid end\n", .{line_number}),
        error.EmptyInterval => printErrorAndExit("error: invalid BED line {d}: end must be greater than start\n", .{line_number}),
    };

    switch (parsed) {
        .skip => {},
        .region => |region| appendBedRegionRequest(requests, region, honor_strand, global_orientation, allocator, name_allocator, duplicate_name, last_duplicated_name),
    }
}

fn appendBedRequests(requests: *std.ArrayList(ParsedRequest), bed_data: []const u8, honor_strand: bool, global_orientation: Orientation, allocator: std.mem.Allocator) void {
    var lines = std.mem.splitScalar(u8, bed_data, '\n');
    var line_number: usize = 0;

    while (lines.next()) |line| {
        line_number += 1;
        appendBedLineRequest(requests, line, line_number, honor_strand, global_orientation, allocator, allocator, false, null);
    }
}

fn processBedData(
    idx: *LoadedIndex,
    allocator: std.mem.Allocator,
    bed_data: []const u8,
    honor_strand: bool,
    global_orientation: Orientation,
    annotate_transform: bool,
    writer: anytype,
) BatchStats {
    var requests = std.ArrayList(ParsedRequest).empty;
    defer requests.deinit(allocator);

    appendBedRequests(&requests, bed_data, honor_strand, global_orientation, allocator);
    return processParsedRequests(idx, allocator, requests.items, annotate_transform, writer, false, null);
}

fn appendNamesRequests(requests: *std.ArrayList(ParsedRequest), names_data: []const u8, orientation: Orientation, allocator: std.mem.Allocator) void {
    var lines = std.mem.splitScalar(u8, names_data, '\n');
    while (lines.next()) |line| {
        const trimmed = if (line.len > 0 and line[line.len - 1] == '\r') line[0 .. line.len - 1] else line;
        if (trimmed.len == 0 or trimmed[0] == '#') continue;

        requests.append(allocator, .{
            .region = .{
                .name = trimmed,
                .start = 1,
                .end = null,
                .is_full = true,
            },
            .orientation = orientation,
        }) catch {
            printErrorAndExit("error: out of memory\n", .{});
        };
    }
}

fn appendCliRequests(requests: *std.ArrayList(ParsedRequest), region_strs: []const []const u8, orientation: Orientation, allocator: std.mem.Allocator) void {
    for (region_strs) |region_str| {
        requests.append(allocator, .{
            .region = parseRegion(region_str),
            .orientation = orientation,
        }) catch {
            printErrorAndExit("error: out of memory\n", .{});
        };
    }
}

fn writeSummary(writer: *std.Io.Writer, region_count: usize, total_bases: u64, elapsed_ns: u64) !void {
    const seconds = if (elapsed_ns == 0) 0.0 else @as(f64, @floatFromInt(elapsed_ns)) / @as(f64, @floatFromInt(std.time.ns_per_s));
    const regions_per_second = if (seconds == 0.0) 0.0 else @as(f64, @floatFromInt(region_count)) / seconds;

    try writer.print("summary: regions={d} total_bases={d} elapsed_s={d:.6} regions_per_s={d:.1}\n", .{ region_count, total_bases, seconds, regions_per_second });
}

fn processParsedRequests(
    idx: *index_format.LoadedIndex,
    allocator: std.mem.Allocator,
    requests: []const ParsedRequest,
    annotate_transform: bool,
    writer: anytype,
    allow_sequential_madvise: bool,
    sparse: ?SparseFastaSource,
) BatchStats {
    if (requests.len == 0) return .{};

    ensureComplementAllowed(idx, requests);

    const resolved = allocator.alloc(ResolvedRegion, requests.len) catch {
        printErrorAndExit("error: out of memory\n", .{});
    };
    defer allocator.free(resolved);
    var already_in_offset_order = requests.len >= 16;
    var prev_start_byte: u64 = 0;

    var last_name: ?[]const u8 = null;
    var last_rec_idx: usize = 0;
    for (requests, 0..) |request, i| {
        const rec_idx = if (last_name) |name|
            if (std.mem.eql(u8, name, request.region.name))
                last_rec_idx
            else
                idx.lookupName(request.region.name) orelse {
                    printErrorAndExit("error: sequence not found: {s}\n", .{request.region.name});
                }
        else
            idx.lookupName(request.region.name) orelse {
                printErrorAndExit("error: sequence not found: {s}\n", .{request.region.name});
            };

        last_name = request.region.name;
        last_rec_idx = rec_idx;

        resolved[i] = resolveParsedRequest(idx, request, rec_idx, i, annotate_transform);
        if (already_in_offset_order) {
            if (i > 0 and resolved[i].start_byte < prev_start_byte) {
                already_in_offset_order = false;
            }
            prev_start_byte = resolved[i].start_byte;
        }
    }

    const dense_sort = if (requests.len >= 16)
        shouldSortByFileOffset(idx, allocator, resolved) catch {
            printErrorAndExit("error: out of memory\n", .{});
        }
    else
        false;
    // Sorting widely spaced requests adds buffering without reducing RSS.
    // Eligible requests with small median gaps sort for sequential page release.
    const release_sorted = dense_sort;

    const use_sort_buffers = release_sorted and !already_in_offset_order and shouldUseSortBuffers(resolved);
    const sequential_scan = dense_sort and (already_in_offset_order or use_sort_buffers);
    const sequential_mmap = shouldAdviseSequentialMmap(
        requests.len,
        sequential_scan,
        batchHasReverseReads(resolved),
        allow_sequential_madvise,
    );
    // Sparse large-FASTA uses positional reads below; MADV on a 3 GiB map is wasted.
    const use_positional = !release_sorted and requests.len >= 16 and
        idx.fasta_size >= sort_by_offset_min_fasta_bytes and sparse != null;
    if (!use_positional) {
        const mmap_advice: platform.Advice = if (sequential_mmap) .sequential else .random;
        platform.advise(idx.fasta_data, mmap_advice);
    }

    var total_bases: u64 = 0;

    if (requests.len < 16) {
        for (resolved) |r| {
            total_bases += r.num_bases;
            emitRegionAndRelease(r, idx.fasta_data, writer);
        }
    } else if (!release_sorted) {
        // Eligible scattered requests use positional reads; other requests keep the mapped path.
        var pread_buf: [max_pread_span_bytes]u8 = undefined;
        for (resolved) |r| {
            total_bases += r.num_bases;
            if (use_positional) {
                if (canEmitRegionViaPread(r, idx.fasta_data.len)) {
                    emitRegionViaPread(r, sparse.?, idx.fasta_data.len, &pread_buf, writer);
                    continue;
                }
            }
            emitRegion(r, idx.fasta_data, writer);
        }
    } else if (already_in_offset_order) {
        var release = FastaReleaseCursor{};
        for (resolved) |r| {
            total_bases += r.num_bases;
            emitRegionAndReleaseSequential(r, idx.fasta_data, writer, &release);
        }
        release.flush(idx.fasta_data);
    } else if (!use_sort_buffers) {
        for (resolved) |r| {
            total_bases += r.num_bases;
            emitRegionAndRelease(r, idx.fasta_data, writer);
        }
    } else {
        var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
        defer arena.deinit();
        const sort_allocator = arena.allocator();

        const sorted = sort_allocator.dupe(ResolvedRegion, resolved) catch {
            printErrorAndExit("error: out of memory\n", .{});
        };
        std.mem.sort(ResolvedRegion, sorted, {}, struct {
            fn lessThan(_: void, a: ResolvedRegion, b: ResolvedRegion) bool {
                return a.start_byte < b.start_byte;
            }
        }.lessThan);

        var total_out: u64 = 0;
        for (resolved) |r| total_out += estimateRegionOutputBytes(r);

        var storage = std.Io.Writer.Allocating.init(sort_allocator);
        storage.ensureTotalCapacity(total_out) catch {
            printErrorAndExit("error: out of memory\n", .{});
        };
        const starts = sort_allocator.alloc(usize, requests.len) catch {
            printErrorAndExit("error: out of memory\n", .{});
        };
        const lens = sort_allocator.alloc(usize, requests.len) catch {
            printErrorAndExit("error: out of memory\n", .{});
        };

        var release = FastaReleaseCursor{};
        for (sorted) |r| {
            total_bases += r.num_bases;
            starts[r.original_index] = storage.writer.end;
            emitRegionAndReleaseSequential(r, idx.fasta_data, &storage.writer, &release);
            lens[r.original_index] = storage.writer.end - starts[r.original_index];
        }
        release.flush(idx.fasta_data);

        const bytes = storage.written();
        for (starts, lens) |start, len| {
            writer.writeAll(bytes[start .. start + len]) catch {
                printErrorAndExit("error: write failed\n", .{});
            };
        }
    }

    return .{
        .region_count = requests.len,
        .total_bases = total_bases,
    };
}

fn processBedReaderChunked(
    idx: *LoadedIndex,
    reader: *std.Io.Reader,
    honor_strand: bool,
    global_orientation: Orientation,
    chunk_size: usize,
    annotate_transform: bool,
    writer: anytype,
    sparse: ?SparseFastaSource,
) BatchStats {
    var total = BatchStats{};
    var line_number: usize = 0;
    var reached_end = false;
    var name_arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer name_arena.deinit();
    var last_duplicated_name: ?[]const u8 = null;

    while (true) {
        var chunk_arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
        const chunk_allocator = chunk_arena.allocator();

        var requests = std.ArrayList(ParsedRequest).empty;

        while (requests.items.len < chunk_size) {
            const maybe_line = reader.takeDelimiter('\n') catch |err| switch (err) {
                error.ReadFailed => printErrorAndExit("error: failed to read BED input\n", .{}),
                error.StreamTooLong => printErrorAndExit(
                    "error: BED line {d} exceeds {d}-byte reader buffer (no newline within limit)\n",
                    .{ line_number + 1, bed_line_reader_buffer_bytes },
                ),
            };

            const line = maybe_line orelse {
                reached_end = true;
                break;
            };

            line_number += 1;
            appendBedLineRequest(&requests, line, line_number, honor_strand, global_orientation, chunk_allocator, name_arena.allocator(), true, &last_duplicated_name);
        }

        if (requests.items.len == 0) {
            chunk_arena.deinit();
            break;
        }

        const batch = processParsedRequests(idx, chunk_allocator, requests.items, annotate_transform, writer, false, sparse);
        total.region_count += batch.region_count;
        total.total_bases += batch.total_bases;

        chunk_arena.deinit();

        if (reached_end) break;
    }

    // Full-map DONTNEED only when the mmap was the read path. Positional sparse
    // BED never faulted those pages; a 3 GiB madvise is pure wall cost.
    if (idx.fasta_size > sparse_large_fasta_bytes and sparse == null) {
        index_format.dropFastaSpan(idx.fasta_data, 0, idx.fasta_data.len);
    }

    return total;
}

fn processBedPathChunked(
    io: std.Io,
    idx: *LoadedIndex,
    path: []const u8,
    honor_strand: bool,
    global_orientation: Orientation,
    chunk_size: usize,
    annotate_transform: bool,
    writer: anytype,
    sparse: ?SparseFastaSource,
) BatchStats {
    if (std.mem.eql(u8, path, "-")) {
        var stdin_buf: [bed_line_reader_buffer_bytes]u8 = undefined;
        var stdin_reader = std.Io.File.stdin().reader(io, &stdin_buf);
        return processBedReaderChunked(idx, &stdin_reader.interface, honor_strand, global_orientation, chunk_size, annotate_transform, writer, sparse);
    }

    const file = std.Io.Dir.cwd().openFile(io, path, .{}) catch |err| switch (err) {
        error.FileNotFound => printErrorAndExit("error: file not found: {s}\n", .{path}),
        error.AccessDenied => printErrorAndExit("error: access denied: {s}\n", .{path}),
        else => printErrorAndExit("error: failed to open file: {s}\n", .{path}),
    };
    defer file.close(io);

    var file_buf: [bed_line_reader_buffer_bytes]u8 = undefined;
    var file_reader = file.reader(io, &file_buf);
    return processBedReaderChunked(idx, &file_reader.interface, honor_strand, global_orientation, chunk_size, annotate_transform, writer, sparse);
}

fn processBedPathAllInMemory(
    io: std.Io,
    idx: *LoadedIndex,
    allocator: std.mem.Allocator,
    path: []const u8,
    honor_strand: bool,
    global_orientation: Orientation,
    annotate_transform: bool,
    writer: anytype,
) BatchStats {
    const bed_data = readAllInput(allocator, io, path, .bed_all_in_memory);
    return processBedData(idx, allocator, bed_data, honor_strand, global_orientation, annotate_transform, writer);
}

// ============================================================================
// Public entry point
// ============================================================================

pub fn runGetWithOptions(io: std.Io, fasta_path: []const u8, options: GetOptions) void {
    if (options.chunk_size == 0) {
        printErrorAndExit("error: --chunk-size must be >= 1\n", .{});
    }

    var input_arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer input_arena.deinit();
    const allocator = input_arena.allocator();

    var requests = std.ArrayList(ParsedRequest).empty;
    defer requests.deinit(allocator);

    if (options.names_path) |names_path| {
        const names_data = readAllInput(allocator, io, names_path, .names_file);
        appendNamesRequests(&requests, names_data, options.orientation, allocator);
    }

    appendCliRequests(&requests, options.region_strs, options.orientation, allocator);

    const start_ns = if (options.summary) monotonicNs(io) else 0;

    // Single positional requests skip the name map.
    // Multi-request, names, and BED paths load it once.
    const load_mode: index_format.LoadMode = if (options.bed_path != null or options.names_path != null or requests.items.len > 1) .lookup_full_map else .records_only;
    var idx = index_format.loadIndexWithMode(io, fasta_path, load_mode);
    defer idx.deinit(io);

    const sparse: ?SparseFastaSource = if (idx.fasta_size >= sort_by_offset_min_fasta_bytes)
        .{ .io = io, .file = idx.fasta_map.file }
    else
        null;

    var out_buf: [65536]u8 = undefined;
    var stdout_fw = std.Io.File.Writer.initStreaming(.stdout(), io, &out_buf);
    const writer = &stdout_fw.interface;

    var totals = BatchStats{};

    if (options.bed_path) |bed_path| {
        const batch = if (options.chunk_size == chunk_size_all)
            processBedPathAllInMemory(io, &idx, allocator, bed_path, options.honor_strand, options.orientation, options.annotate_transform, writer)
        else
            processBedPathChunked(io, &idx, bed_path, options.honor_strand, options.orientation, options.chunk_size, options.annotate_transform, writer, sparse);
        totals.region_count += batch.region_count;
        totals.total_bases += batch.total_bases;
    }

    if (requests.items.len > 0) {
        const batch = processParsedRequests(&idx, allocator, requests.items, options.annotate_transform, writer, true, sparse);
        totals.region_count += batch.region_count;
        totals.total_bases += batch.total_bases;
    }

    if (totals.region_count == 0) {
        printErrorAndExit("error: no regions provided\n", .{});
    }

    stdout_fw.flush() catch {
        printErrorAndExit("error: write failed\n", .{});
    };

    if (options.summary) {
        var err_buf: [512]u8 = undefined;
        var stderr_fw = std.Io.File.Writer.initStreaming(.stderr(), io, &err_buf);
        writeSummary(&stderr_fw.interface, totals.region_count, totals.total_bases, monotonicNs(io) - start_ns) catch {
            printErrorAndExit("error: write failed\n", .{});
        };
        stderr_fw.flush() catch {
            printErrorAndExit("error: write failed\n", .{});
        };
    }
}

test "summary writing fails when the destination is full" {
    var buffer: [1]u8 = undefined;
    var writer = std.Io.Writer.fixed(&buffer);

    try std.testing.expectError(error.WriteFailed, writeSummary(&writer, 1, 1, 1));
}

test "processBedReaderChunked matches non-chunked extraction" {
    const test_io = std.testing.io;
    var idx = index_format.loadIndex(test_io, "tests/data/simple.fasta");
    defer idx.deinit(test_io);

    const bed_data =
        "# comment\n" ++
        "seq2\t0\t4\n" ++
        "seq1\t0\t4\n" ++
        "seq1\t4\t8\n";

    var chunk_reader = std.Io.Reader.fixed(bed_data);
    var chunk_writer = std.Io.Writer.Allocating.init(std.testing.allocator);
    defer chunk_writer.deinit();

    const chunked = processBedReaderChunked(&idx, &chunk_reader, false, .{}, 2, false, &chunk_writer.writer, null);

    var batch_writer = std.Io.Writer.Allocating.init(std.testing.allocator);
    defer batch_writer.deinit();

    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const batch = processBedData(&idx, arena.allocator(), bed_data, false, .{}, false, &batch_writer.writer);

    try std.testing.expectEqual(batch.region_count, chunked.region_count);
    try std.testing.expectEqual(batch.total_bases, chunked.total_bases);
    try std.testing.expectEqualStrings(batch_writer.written(), chunk_writer.written());
}

test "processBedReaderChunked preserves strand handling across chunk boundaries" {
    const test_io = std.testing.io;
    var idx = index_format.loadIndex(test_io, "tests/data/simple.fasta");
    defer idx.deinit(test_io);

    const bed_data =
        "seq1\t0\t5\tname\t0\t-\n" ++
        "seq2\t0\t4\tname\t0\t+\n";

    var chunk_reader = std.Io.Reader.fixed(bed_data);
    var chunk_writer = std.Io.Writer.Allocating.init(std.testing.allocator);
    defer chunk_writer.deinit();

    const chunked = processBedReaderChunked(&idx, &chunk_reader, true, .{}, 1, false, &chunk_writer.writer, null);

    try std.testing.expectEqual(@as(usize, 2), chunked.region_count);
    try std.testing.expectEqual(@as(u64, 9), chunked.total_bases);
    try std.testing.expectEqualStrings(
        ">seq1:1-5\nTACGT\n>seq2:1-4\nGGGG\n",
        chunk_writer.written(),
    );
}

test "processBedReaderChunked preserves duplicate chrom cache across chunk boundaries" {
    const test_io = std.testing.io;
    var idx = index_format.loadIndex(test_io, "tests/data/simple.fasta");
    defer idx.deinit(test_io);

    const bed_data =
        "seq1\t0\t4\n" ++
        "seq1\t4\t8\n" ++
        "seq2\t0\t4\n";

    var chunk_reader = std.Io.Reader.fixed(bed_data);
    var chunk_writer = std.Io.Writer.Allocating.init(std.testing.allocator);
    defer chunk_writer.deinit();

    const chunked = processBedReaderChunked(&idx, &chunk_reader, false, .{}, 1, false, &chunk_writer.writer, null);

    var batch_writer = std.Io.Writer.Allocating.init(std.testing.allocator);
    defer batch_writer.deinit();

    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const batch = processBedData(&idx, arena.allocator(), bed_data, false, .{}, false, &batch_writer.writer);

    try std.testing.expectEqual(batch.region_count, chunked.region_count);
    try std.testing.expectEqual(batch.total_bases, chunked.total_bases);
    try std.testing.expectEqualStrings(batch_writer.written(), chunk_writer.written());
}

test "parseRegion matches exhaustive backward-scan reference" {
    const Reference = struct {
        fn findLastColon(input: []const u8, end: usize) ?usize {
            var i = end;
            while (i > 0) {
                i -= 1;
                if (input[i] == ':') return i;
            }
            return null;
        }

        fn parse(input: []const u8) Region {
            var search_end = input.len;
            while (findLastColon(input, search_end)) |colon_pos| {
                if (parseRangeSuffix(input[colon_pos + 1 ..])) |range| {
                    return .{
                        .name = input[0..colon_pos],
                        .start = range.start,
                        .end = range.end,
                        .is_full = false,
                    };
                }
                search_end = colon_pos;
            }
            return .{ .name = input, .start = 1, .end = null, .is_full = true };
        }
    };

    const alphabet = "a:0-1";
    var input: [7]u8 = undefined;
    var combinations: usize = 1;
    var tested: usize = 0;
    for (1..input.len + 1) |length| {
        combinations *= alphabet.len;
        for (0..combinations) |encoded| {
            var value = encoded;
            for (input[0..length]) |*byte| {
                byte.* = alphabet[value % alphabet.len];
                value /= alphabet.len;
            }

            const expected = Reference.parse(input[0..length]);
            const actual = parseRegion(input[0..length]);
            try std.testing.expectEqualStrings(expected.name, actual.name);
            try std.testing.expectEqual(expected.start, actual.start);
            try std.testing.expectEqual(expected.end, actual.end);
            try std.testing.expectEqual(expected.is_full, actual.is_full);
            tested += 1;
        }
    }
    try std.testing.expectEqual(@as(usize, 97_655), tested);
}

test "parseRegion handles coordinate overflow at the rightmost colon" {
    const max = parseRegion("name:18446744073709551615-18446744073709551615");
    try std.testing.expectEqualStrings("name", max.name);
    try std.testing.expectEqual(@as(u64, std.math.maxInt(u64)), max.start);
    try std.testing.expectEqual(@as(?u64, std.math.maxInt(u64)), max.end);
    try std.testing.expect(!max.is_full);

    const full_names = [_][]const u8{
        "name:18446744073709551616-1",
        "name:1-18446744073709551616",
        "name:1-2:bad",
        "name:1-2:18446744073709551616-1",
    };
    for (full_names) |name| {
        const region = parseRegion(name);
        try std.testing.expectEqualStrings(name, region.name);
        try std.testing.expect(region.is_full);
    }
}

test "orientation compose behaves like transform composition" {
    const minus_strand = Orientation.reverseComplement();
    const global_rc = Orientation.reverseComplement();
    const global_reverse = Orientation.reverseOnly();
    const global_complement = Orientation.complementOnly();

    try std.testing.expect(minus_strand.compose(global_rc).isIdentity());
    try std.testing.expectEqual(true, minus_strand.compose(global_reverse).complement);
    try std.testing.expectEqual(false, minus_strand.compose(global_reverse).reverse);
    try std.testing.expectEqual(true, minus_strand.compose(global_complement).reverse);
    try std.testing.expectEqual(false, minus_strand.compose(global_complement).complement);
}

test "processParsedRequests applies complement-only and reverse-only transforms" {
    const test_io = std.testing.io;
    var idx = index_format.loadIndex(test_io, "tests/data/simple.fasta");
    defer idx.deinit(test_io);

    const requests = [_]ParsedRequest{
        .{ .region = parseRegion("seq1:1-5"), .orientation = Orientation.complementOnly() },
        .{ .region = parseRegion("seq1:1-5"), .orientation = Orientation.reverseOnly() },
    };

    var writer = std.Io.Writer.Allocating.init(std.testing.allocator);
    defer writer.deinit();

    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();

    const batch = processParsedRequests(&idx, arena.allocator(), &requests, true, &writer.writer, true, null);

    try std.testing.expectEqual(@as(usize, 2), batch.region_count);
    try std.testing.expectEqualStrings(
        ">seq1:1-5 (complement)\nTGCAT\n>seq1:1-5 (reverse)\nATGCA\n",
        writer.written(),
    );
}

test "detectRecordType classifies nucleotide and protein records" {
    const test_io = std.testing.io;

    var nucleotide_idx = index_format.loadIndex(test_io, "tests/data/simple.fasta");
    defer nucleotide_idx.deinit(test_io);
    try std.testing.expectEqual(
        stats.SequenceType.nucleotide,
        detectRecordType(nucleotide_idx.records[0], nucleotide_idx.fasta_data),
    );

    var protein_idx = index_format.loadIndex(test_io, "tests/data/proteome.fasta");
    defer protein_idx.deinit(test_io);
    try std.testing.expectEqual(
        stats.SequenceType.protein,
        detectRecordType(protein_idx.records[0], protein_idx.fasta_data),
    );
}

test "shouldSortByFileOffset rejects small catalogs and wide gaps" {
    const test_io = std.testing.io;
    var idx = index_format.loadIndex(test_io, "tests/data/simple.fasta");
    defer idx.deinit(test_io);

    var regions: [20]ResolvedRegion = undefined;
    for (&regions, 0..) |*r, i| {
        r.* = .{
            .name = "seq1",
            .start = @intCast(i + 1),
            .display_end = @intCast(i + 1),
            .is_full = false,
            .start_byte = if (i == 0) 0 else 10_000_000 + i,
            .seq_offset = 6,
            .num_bases = 1,
            .line_bases = 24,
            .line_bytes = 25,
            .side_table = &.{},
            .orientation = .{},
            .annotate_transform = false,
            .original_index = i,
        };
    }
    try std.testing.expect(!try shouldSortByFileOffset(&idx, std.testing.allocator, &regions));

    for (&regions, 0..) |*r, i| r.start_byte = 100 + i;
    try std.testing.expect(!try shouldSortByFileOffset(&idx, std.testing.allocator, &regions));
}

test "shouldSortByFileOffsetForBatch requires large FASTA even for large catalogs" {
    var regions: [20]ResolvedRegion = undefined;
    for (&regions, 0..) |*r, i| {
        r.* = .{
            .name = "seq1",
            .start = @intCast(i + 1),
            .display_end = @intCast(i + 1),
            .is_full = false,
            .start_byte = 100 + i,
            .seq_offset = 6,
            .num_bases = 1,
            .line_bases = 24,
            .line_bytes = 25,
            .side_table = &.{},
            .orientation = .{},
            .annotate_transform = false,
            .original_index = i,
        };
    }
    // Proteome-scale catalog on a small FASTA (direct seek wins).
    try std.testing.expect(!try shouldSortByFileOffsetForBatch(
        20_659,
        13 * 1024 * 1024,
        std.testing.allocator,
        &regions,
    ));
    // Transcriptome-scale catalog on a large FASTA with dense offsets.
    try std.testing.expect(try shouldSortByFileOffsetForBatch(
        254_070,
        459 * 1024 * 1024,
        std.testing.allocator,
        &regions,
    ));
    // Large FASTA but sparse median gaps (genome-style within batch).
    for (&regions, 0..) |*r, i| {
        r.start_byte = if (i == 0) 0 else @as(u64, i) * 1024 * 1024;
    }
    try std.testing.expect(!try shouldSortByFileOffsetForBatch(
        254_070,
        459 * 1024 * 1024,
        std.testing.allocator,
        &regions,
    ));
}

test "shouldUseSortBuffers accepts typical 20-region batch" {
    var regions: [20]ResolvedRegion = undefined;
    for (&regions, 0..) |*r, i| {
        r.* = .{
            .name = "seq1",
            .start = @intCast(i + 1),
            .display_end = @intCast(i + 1),
            .is_full = false,
            .start_byte = 100 + i,
            .seq_offset = 6,
            .num_bases = 1,
            .line_bases = 24,
            .line_bytes = 25,
            .side_table = &.{},
            .orientation = .{},
            .annotate_transform = false,
            .original_index = i,
        };
    }
    try std.testing.expect(shouldUseSortBuffers(&regions));
}

test "shouldUseSortBuffers rejects oversized single region" {
    const region = ResolvedRegion{
        .name = "chr1",
        .start = 1,
        .display_end = max_sort_path_region_output_bytes,
        .is_full = false,
        .start_byte = 0,
        .seq_offset = 0,
        .num_bases = max_sort_path_region_output_bytes,
        .line_bases = 60,
        .line_bytes = 61,
        .side_table = &.{},
        .orientation = .{},
        .annotate_transform = false,
        .original_index = 0,
    };
    try std.testing.expect(!shouldUseSortBuffers(&.{region}));
}

test "shouldAdviseSequentialMmap only when sequential scan is active" {
    try std.testing.expect(shouldAdviseSequentialMmap(100, true, false, true));
    try std.testing.expect(!shouldAdviseSequentialMmap(100, false, false, true));
    try std.testing.expect(!shouldAdviseSequentialMmap(100, true, true, true));
    try std.testing.expect(!shouldAdviseSequentialMmap(100, true, false, false));
    try std.testing.expect(!shouldAdviseSequentialMmap(10, true, false, true));
}

test "emitRegion handles non-uniform side-table forward extraction" {
    const fasta = ">seq\nAAA\nCCCC \nGG\n";
    const side_table = [_]index_format.SideTableLine{
        .{ .base_start = 0, .byte_offset = 5, .line_bytes = 4, .line_bases = 3 },
        .{ .base_start = 3, .byte_offset = 9, .line_bytes = 6, .line_bases = 4 },
        .{ .base_start = 7, .byte_offset = 15, .line_bytes = 3, .line_bases = 2 },
    };
    const rec = IndexRecord{
        .name_offset = 1,
        .name_len = 3,
        .seq_offset = 5,
        .seq_len = 9,
        .line_bases = 3,
        .line_bytes = 4,
    };

    var writer = std.Io.Writer.Allocating.init(std.testing.allocator);
    defer writer.deinit();

    emitRegion(.{
        .name = "seq",
        .start = 2,
        .display_end = 8,
        .is_full = false,
        .start_byte = byteOffsetForBase(fasta, rec, &side_table, 1),
        .seq_offset = 5,
        .num_bases = 7,
        .line_bases = 3,
        .line_bytes = 4,
        .side_table = &side_table,
        .orientation = .{},
        .annotate_transform = false,
        .original_index = 0,
    }, fasta, &writer.writer);

    try std.testing.expectEqualStrings(">seq:2-8\nAACCCCG\n", writer.written());
}

test "emitRegion handles non-uniform side-table reverse complement extraction" {
    const fasta = ">seq\nAAA\nCCCC \nGG\n";
    const side_table = [_]index_format.SideTableLine{
        .{ .base_start = 0, .byte_offset = 5, .line_bytes = 4, .line_bases = 3 },
        .{ .base_start = 3, .byte_offset = 9, .line_bytes = 6, .line_bases = 4 },
        .{ .base_start = 7, .byte_offset = 15, .line_bytes = 3, .line_bases = 2 },
    };
    const rec = IndexRecord{
        .name_offset = 1,
        .name_len = 3,
        .seq_offset = 5,
        .seq_len = 9,
        .line_bases = 3,
        .line_bytes = 4,
    };

    var writer = std.Io.Writer.Allocating.init(std.testing.allocator);
    defer writer.deinit();

    emitRegion(.{
        .name = "seq",
        .start = 2,
        .display_end = 8,
        .is_full = false,
        .start_byte = byteOffsetForBase(fasta, rec, &side_table, 1),
        .seq_offset = 5,
        .num_bases = 7,
        .line_bases = 3,
        .line_bytes = 4,
        .side_table = &side_table,
        .orientation = Orientation.reverseComplement(),
        .annotate_transform = false,
        .original_index = 0,
    }, fasta, &writer.writer);

    try std.testing.expectEqualStrings(">seq:2-8\nCGGGGTT\n", writer.written());
}

test "uniform emitters flush complete chunks at the output-buffer boundary" {
    const base_count: usize = 70_000;
    const allocator = std.testing.allocator;
    const fasta = try allocator.alloc(u8, base_count + 1);
    defer allocator.free(fasta);
    @memset(fasta[0..base_count], 'A');
    fasta[base_count] = '\n';

    const orientations = [_]Orientation{ .{}, Orientation.reverseComplement() };
    for (orientations) |orientation| {
        var writer = std.Io.Writer.Allocating.init(allocator);
        defer writer.deinit();

        emitRegion(.{
            .name = "seq",
            .start = 1,
            .display_end = base_count,
            .is_full = true,
            .start_byte = 0,
            .seq_offset = 0,
            .num_bases = base_count,
            .line_bases = base_count,
            .line_bytes = base_count + 1,
            .side_table = &.{},
            .orientation = orientation,
            .annotate_transform = false,
            .original_index = 0,
        }, fasta, &writer.writer);

        var lines = std.mem.splitScalar(u8, writer.written(), '\n');
        try std.testing.expectEqualStrings(">seq", lines.next().?);
        var emitted: usize = 0;
        while (lines.next()) |line| {
            if (line.len == 0) break;
            try std.testing.expectEqual(@min(@as(usize, 60), base_count - emitted), line.len);
            for (line) |base| {
                try std.testing.expectEqual(if (orientation.complement) @as(u8, 'T') else 'A', base);
            }
            emitted += line.len;
        }
        try std.testing.expectEqual(base_count, emitted);
    }
}

test "positional emitter flushes at the output-buffer boundary" {
    const base_count: usize = 65_000;
    const allocator = std.testing.allocator;
    const fasta = try allocator.alloc(u8, base_count + 1);
    defer allocator.free(fasta);
    @memset(fasta[0..base_count], 'A');
    fasta[base_count] = '\n';

    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    const file = try tmp.dir.createFile(std.testing.io, "source.fa", .{ .read = true });
    defer file.close(std.testing.io);
    try std.Io.File.writeStreamingAll(file, std.testing.io, fasta);

    const resolved = ResolvedRegion{
        .name = "seq",
        .start = 1,
        .display_end = base_count,
        .is_full = true,
        .start_byte = 0,
        .seq_offset = 0,
        .num_bases = base_count,
        .line_bases = base_count,
        .line_bytes = base_count + 1,
        .side_table = &.{},
        .orientation = Orientation.complementOnly(),
        .annotate_transform = false,
        .original_index = 0,
    };
    try std.testing.expect(canEmitRegionViaPread(resolved, fasta.len));

    var file_buf: [max_pread_span_bytes]u8 = undefined;
    var writer = std.Io.Writer.Allocating.init(allocator);
    defer writer.deinit();
    emitRegionViaPread(
        resolved,
        .{ .io = std.testing.io, .file = file },
        fasta.len,
        &file_buf,
        &writer.writer,
    );

    const sequence_lines = std.math.divCeil(usize, base_count, 60) catch unreachable;
    try std.testing.expectEqual(">seq\n".len + base_count + sequence_lines, writer.written().len);
    var lines = std.mem.splitScalar(u8, writer.written(), '\n');
    try std.testing.expectEqualStrings(">seq", lines.next().?);
    var emitted: usize = 0;
    while (lines.next()) |line| {
        if (line.len == 0) break;
        try std.testing.expectEqual(@min(@as(usize, 60), base_count - emitted), line.len);
        for (line) |base| try std.testing.expectEqual(@as(u8, 'T'), base);
        emitted += line.len;
    }
    try std.testing.expectEqual(base_count, emitted);
}
