[![Rust](https://github.com/stuart-lab/fragtk/actions/workflows/rust.yml/badge.svg)](https://github.com/stuart-lab/fragtk/actions/workflows/rust.yml)

# fragtk: fragment file toolkit

`fragtk` is a fast and efficient toolkit implemented in Rust for working with fragment files.

To install `fragtk`, please see the [installation](#installation) section below.

## Usage

### Create region x cell matrix

A region x cell matrix can be created from a fragment file and a peak file:

```
fragtk matrix -f <fragments.tsv.gz> -b <peaks.bed> -c <cells.txt> -o <output>
```

### Count fragments per cell barcode

Select cell barcodes from the fragment file according to their total count:

```
fragtk count -f <fragments.tsv.gz> -o <barcode_counts.tsv> -t <threshold> > barcodes.txt
```

### Filter fragments

Filter fragments according to the cell barcodes:

```
fragtk filter -f <fragments.tsv.gz> -c <barcodes.txt> | bgzip -c > filtered.tsv.gz
```

## Installation

### Using cargo

```
cargo install fragtk
```

### From GitHub

```
git clone git@github.com:stuart-lab/fragtk.git
cd fragtk; cargo install --path .
```

Pre-compiled binaries are also available in the release.
