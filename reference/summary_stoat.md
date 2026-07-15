# Summary of STOAT's results

Print a few summary metrics from the association results

## Usage

``` r
summary_stoat(assoc)
```

## Arguments

- assoc:

  Either the data.frame imported by \*import_assoc\*, or the path to
  STOAT's output (\*assoc.pvalues.tsv.gz)

## Value

In addition to printing the summary metrics, returns a list with the
values

## Examples

``` r
# prepare the filename
assoc_file = system.file('extdata/stoat.quantitative.assoc.pvalues.tsv.gz', package='StoatPlot')

# print a summary
sum.l = summary_stoat(assoc_file)
#> 
#> 
#> |Metric                   | Value|
#> |:------------------------|-----:|
#> |Total variants           |   680|
#> |Variant PV<0.01          |    31|
#> |Variant adjusted PV<0.01 |     2|
#> |Genomic inflation factor |  1.15|
#> 
#> 
#> |Type | Total| Pv_BH_below_0.01|
#> |:----|-----:|----------------:|
#> |SNP  |   519|                2|
#> |SV   |   161|                0|
#> 
#> 
sum.l 
#> $all
#> # A tibble: 4 × 2
#>   Metric                    Value
#>   <chr>                     <dbl>
#> 1 Total variants           680   
#> 2 Variant PV<0.01           31   
#> 3 Variant adjusted PV<0.01   2   
#> 4 Genomic inflation factor   1.15
#> 
#> $per.type
#> # A tibble: 2 × 3
#>   Type  Total Pv_BH_below_0.01
#>   <fct> <int>            <int>
#> 1 SNP     519                2
#> 2 SV      161                0
#> 
```
