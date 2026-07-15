# Import association results from STOAT

Read the TSV (potentially gzipped) file produced by STOAT and load it.

## Usage

``` r
import_assoc(assoc_file)
```

## Arguments

- assoc_file:

  Path to STOAT's output (\*assoc.pvalues.tsv.gz)

## Value

A polished data.frame with a subset of the TSV file.
