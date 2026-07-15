# Combine the allele counts for snarl with phenotype and covariables

Combine a phenotype and potentially some covariables with the allele
count information from import_genotype_chunk. If paths are given to
phenotypr or covariables, will try to read the files.

## Usage

``` r
combine_chunk_phenotype_covars(
  snarl_ac,
  phenotype,
  covariables = NULL,
  gene_name = NULL
)
```

## Arguments

- snarl_ac:

  data.frame with the allele counts for the snarl

- phenotype:

  Either a data.frame with the phenotype for each sample, or the path to
  the phenotype file

- covariables:

  Either a data.frame with the covariables for each sample, or the path
  to the covariables file

- gene_name:

  gene of interest, if the phenotype is a matrix of gene expression.

## Value

a data.frame with the new phenotype and covariables columns
