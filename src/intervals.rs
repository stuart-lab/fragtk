use rust_lapper::{IterFind, Lapper};

/// A cursor for [`seek_position`], remembering the last position queried with it.
///
/// Queries made through one cursor should move forward; use a separate cursor for each
/// series of positions (fragment starts and fragment ends are two such series).
pub struct SeekCursor {
    idx: usize,
    last: u32,
}

impl SeekCursor {
    pub fn new() -> Self {
        SeekCursor { idx: 0, last: 0 }
    }

    /// Start over, e.g. on moving to a new chromosome.
    pub fn reset(&mut self) {
        self.idx = 0;
        self.last = 0;
    }
}

/// Find intervals overlapping a single position, advancing `cursor`.
///
/// `Lapper::seek` re-seeds its index only when that index is a valid index into
/// `intervals` AND the interval there starts after the query. Neither guard fires in
/// two cases, and in both the result is silently missed overlaps:
///
/// 1. Once a query runs past the final interval the index is left at
///    `intervals.len()`, so `idx < intervals.len()` is false forever after and every
///    later query on that chromosome returns nothing.
/// 2. If the index sits on an interval that starts before the query but after another
///    interval that also overlaps it, the scan begins too late and returns only some
///    of the matches. This needs overlapping intervals to be reachable at all, so it
///    shows up wherever annotation windows overlap each other (TSS windows of nearby
///    genes) once a previous query pushed the index forward.
///
/// Both are backwards moves relative to where the index was left, so resetting on any
/// backwards query makes `seek` re-run its binary search and return every overlap.
#[inline]
pub fn seek_position<'a>(
    lapper: &'a Lapper<u32, usize>,
    pos: u32,
    cursor: &mut SeekCursor,
) -> IterFind<'a, u32, usize> {
    if pos < cursor.last || cursor.idx >= lapper.intervals.len() {
        cursor.idx = 0;
    }
    cursor.last = pos;
    lapper.seek(pos, pos + 1, &mut cursor.idx)
}
