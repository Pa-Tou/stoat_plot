# P-value Histogram for STOAT GWAS Results

Generate histogram of P-values from STOAT GWAS results TSV.

## Usage

``` r
plot_pvalue_hist(
  gwas_file,
  p_column = "P",
  min = 0,
  max = 1,
  bins = 100,
  output = NULL
)
```

## Arguments

- gwas_file:

  Path to the output stoat GWAS TSV file.

- p_column:

  Column name to use for p-values (default: ""). If empty, will use "P"
  or "P_CHI2" if available.

- min:

  Minimun P-value threshold to include in the plot (default: 0).

- max:

  Maximum P-value threshold to include in the plot (default: 1.0).

- bins:

  Number of bins in the histogram (default: 200).

- output:

  Filename to save the output plot (default:
  "pvalue_distribution_plot.png").

## Value

Saves a histogram plot as an image file.
