# 1.7.0

### `fragtk matrix`

- Memory usage improvements.
- Added an `--h5` option to directly output the results as a 10x Genomics-compatible HDF5 file natively.
  Note: This requires compiling with `cargo build --release --features hdf5`. With `--h5`, `--outdir`
  may name the output file itself: a path ending in `.h5` or `.hdf5` is written to directly rather than
  being created as a directory to hold `matrix.h5`.
- Re-architected Matrix Market coordinate outputs so elements are explicitly sorted sequentially by
  genomic position and cell index.
- Changed output format to use `chr:start-end` format for coordinates.
- Now errors on duplicate cell barcodes in `--cells`. Duplicates previously produced a matrix header
  declaring fewer columns than the largest column index written, giving an mtx file that
  `Matrix::readMM` and `scipy.io.mmread` reject.
- Now errors when the fragment file is not sorted by position within a chromosome, which previously
  could silently produce zero counts.

### `fragtk qc`

- Memory usage improvements.
- Added optional `--cells` argument to subset TSS evaluation for only specific cell barcodes.
- Changed TSS enrichment to follow the ENCODE definition: the score is the fold change of the mean
  depth over the central 100bp against the mean depth over the 100bp at each end of the +/-2000bp
  distribution.
- Added a floor on the flanking-region depth used to normalize the TSS enrichment score.
  Cells with very few insertions in the flanking regions were previously divided by
  a near-zero background resulting in a very high score.
- Renamed the `FRiP` column to `FIP` (fraction of insertions in promoters), and changed it to be
  computed over a dedicated set of promoter regions (TSS +/-1000bp), flattened so overlapping
  annotation records cannot count an insertion more than once, and divided by total insertions rather
  than total fragments.

Together these change the value of `TSS_enrichment` and `FIP`: scores are not comparable to earlier
versions and QC thresholds will need recalibrating.

### Other changes

- Fixed missed overlaps in the interval search. The cursor was shared between fragment start and end
  positions and was never rewound, so once a fragment end advanced it past an overlapping interval,
  later queries returned only some of the matching intervals (or none at all, once it ran past the
  last one). This undercounted `TSS_center` wherever TSS windows of nearby genes overlap.
- Added robust magic byte detection for `.gz` file recognition rather than relying on strict `.gz` file
  extensions.
- Memory usage improvements across the `filter` and `count` commands.
- Fixed an issue with compilation on Windows.

# 1.6.0

- Add `fragtk qc` command

# 1.5.0

- Speed and memory improvements

# 1.4.0

- Enable `--cells` file to be gzipped
- Update documentation
- Add option to build man files

# 1.3.0

- Add `--pic` option for paired insertion counting to `fragtk matrix`

# 1.2.0

- Update command line interface
- Allow input BED file for `fragtk matrix` to be gzipped
- Add `--ncell` option to `fragtk count`

# 1.1.0

Performance improvements

# 1.0.0

First release
