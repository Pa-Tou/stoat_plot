# Path Length Distribution Plots for STOAT

Generate multiple visualizations of path lengths from a TSV file with a
\`TYPE\` column.

## Usage

``` r
path_length_distribution(snarl_file, min = 0, max = Inf, output = NULL)
```

## Arguments

- snarl_file:

  Path to the snarl_file TSV file.

- min:

  Minimum path length to include (default: 0).

- max:

  Maximum path length to include (default: Inf).

- output:

  Path to the output plot (default: "path_length_dot_ecdf.png").

## Value

Saves four plots: dot plot, histogram, ECDF, and boxplot.
