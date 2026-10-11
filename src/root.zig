//! The `z-fasta` package. Only these declarations are supported. Calls fill a caller-owned
//! `diagnostics.Failure` and never print or exit.

const indexed_reader = @import("indexed_reader.zig");

pub const diagnostics = @import("diagnostics");
pub const IndexedReader = indexed_reader.IndexedReader;
pub const IndexSelection = indexed_reader.IndexSelection;
pub const Interval = indexed_reader.Interval;
pub const RecordInfo = indexed_reader.RecordInfo;
