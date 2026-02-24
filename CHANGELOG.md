# 1.7.0

- Memory usage improvements across the `matrix`, `filter`, `qc`, `count` commands.
- Added robust magic byte detection for `.gz` file recognition rather than relying on strict `.gz` file extensions.
- Re-architected Matrix Market coordinate outputs in `fragtk matrix` so elements are explicitly sorted sequentially by genomic position and cell index.
- Added optional `--cells` argument to `fragtk qc` to subset TSS evaluation for only specific cell barcodes.
- Changed output format for `fragtk matrix` to use `chr:start-end` format for coordinates.

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
