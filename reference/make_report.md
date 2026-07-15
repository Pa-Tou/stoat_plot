# General a summary report.

An HTML file will be created in the working directory.

## Usage

``` r
make_report(
  assoc_file,
  output_html_file = "stoat.summary.report.html",
  force_overwrite = FALSE
)
```

## Arguments

- assoc_file:

  path to the association file

- output_html_file:

  path for the output HTML file

- force_overwrite:

  should existing Rmd/HTML output file be overwritten?

## Value

a character with the report filename
