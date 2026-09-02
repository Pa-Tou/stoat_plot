# Import association results from STOAT

Read the TSV (potentially gzipped) file produced by STOAT and load it.

## Usage

``` r
import_assoc(assoc_file, all_columns = FALSE)
```

## Arguments

- assoc_file:

  Path to STOAT's output (\*assoc.pvalues.tsv.gz)

- all_columns:

  Should all the columns be returned? Default: FALSE

## Value

A polished data.frame with a subset of the TSV file.
