# Using a custom statistical test

``` r

library(StoatPlot)
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
library(knitr) ## for table outputs
library(ggplot2)
library(car) ## to test for specific coeff in regression
#> Loading required package: carData
#> 
#> Attaching package: 'car'
#> The following object is masked from 'package:dplyr':
#> 
#>     recode
```

We use the small test data from the package.

``` r

test_data_dir <- system.file('extdata', package='StoatPlot')
gt_fn = file.path(test_data_dir, 'snarl_genotypes.sorted.tsv.gz')
pheno_fn = file.path(test_data_dir, 'phenotype.quantitative.tsv')
```

## Read the phenotype file

``` r

pheno = readr::read_tsv(pheno_fn)
#> Rows: 40 Columns: 2
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: "\t"
#> chr (1): SAMPLE
#> dbl (1): PHENO
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

pheno |> head() |> kable()
```

| SAMPLE    |      PHENO |
|:----------|-----------:|
| samp_g0_0 | -0.2448110 |
| samp_g0_1 | -0.0251954 |
| samp_g0_2 | -0.1081481 |
| samp_g0_3 |  0.2217557 |
| samp_g0_4 | -1.5727460 |
| samp_g0_5 |  0.1494528 |

## Make two random covariables

``` r

covars = tibble(SAMPLE=pheno$SAMPLE,
                sex=sample(0:1, nrow(pheno), replace=TRUE),
                PC1=rnorm(nrow(pheno)))
```

## Extract snarls in a region

The genotype matrix can get very big and we don’t necessarily want to
load it entirely in memory. Instead, we can load the snarl information
in chunks.

Let’s extract the snarls overlapping `ref:20000-50000`.

``` r

chunk.l = import_genotype_chunk(gt_fn, 'ref', 20000, 50000)
length(chunk.l)
#> [1] 183
```

## Building a test table for one snarl

Let’s take the first snarl info and make a test table by combining with
the phenotype.

``` r

df = combine_chunk_phenotype_covars(chunk.l[[1]]$GT, pheno, covars)

df |> head() |> kable()
```

| sample     | al0 | al1 |  phenotype | covarsex |  covarPC1 |
|:-----------|----:|----:|-----------:|---------:|----------:|
| samp_g0_0  |   1 |   1 | -0.2448110 |        0 | 0.4681544 |
| samp_g0_1  |   0 |   2 | -0.0251954 |        0 | 0.3629513 |
| samp_g0_10 |   1 |   1 | -1.5837474 |        1 | 0.9353632 |
| samp_g0_11 |   0 |   2 | -1.0077583 |        0 | 0.1764886 |
| samp_g0_12 |   1 |   1 |  0.9459124 |        0 | 0.2436855 |
| samp_g0_13 |   0 |   2 |  2.6205592 |        0 | 1.6235489 |

This table could be used for a test, for example fitting a regression
model.

In STOAT, we also preprocess it a bit to remove constant variables,
obvious colinearity. The *sample* column is also not used in the
regression. We’ve created a function to mimick the table preparation
used by STOAT.

``` r

df = preprocess_test_table(df)

df |> head() |> kable()
```

| al1 |  phenotype | covarsex |  covarPC1 |
|----:|-----------:|---------:|----------:|
|   1 | -0.2448110 |        0 | 0.4681544 |
|   2 | -0.0251954 |        0 | 0.3629513 |
|   1 | -1.5837474 |        1 | 0.9353632 |
|   2 | -1.0077583 |        0 | 0.1764886 |
|   1 |  0.9459124 |        0 | 0.2436855 |
|   2 |  2.6205592 |        0 | 1.6235489 |

Now the table only contains:

- a `phenotype` variable
- one or several *allele* variables prefixed by `al`
- one or several *covariables* prefixed by `covar`

In STOAT, we fit a model and test if any of the coefficients for the
*allele* variables are different from 0.

## Compute p-values using linear regression

For example, let’s compute a pvalue by fitting a linear regression model
for each snarl.

``` r

pv.df = lapply(chunk.l, function(snarl_info) {
  ## combine allele count with the phenotype
  df = combine_chunk_phenotype_covars(snarl_info$GT, pheno, covars)
  df = preprocess_test_table(df)

  ## run the regression
  model = lm(phenotype ~ ., data=df)
  ## null hypothesis is all allele coeffs are 0
  null.h = paste0(grep('^al', colnames(df), value=TRUE), '=0')
  lh.o = linearHypothesis(model, null.h)
  pv = lh.o[['Pr(>F)']][2]

  ## return a data.frame with the snarl information and the test result
  return(tibble(START_NODE=snarl_info$START_NODE, END_NODE=snarl_info$END_NODE,
                P=pv))
}) |> bind_rows()

pv.df |> head() |> kable()
```

| START_NODE | END_NODE |         P |
|:-----------|:---------|----------:|
| \>378      | \<381    | 0.6484383 |
| \>381      | \<384    | 0.9463934 |
| \>384      | \<387    | 0.2391300 |
| \>387      | \<389    | 0.4989925 |
| \>389      | \<392    | 0.6818881 |
| \>392      | \<395    | 0.1115229 |
