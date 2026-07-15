# Import a chunk of the snarl genotypes for a queried region

Import a slice of the snarl genotype file (sorted, bgzipped and indexed
with Tabix). For each snarl, a table with sample and allele counts is
return in a list, along with other information about the snarl.

## Usage

``` r
import_genotype_chunk(
  genotype_file,
  chrom,
  start_offset,
  end_offset,
  keep_nonvariant = FALSE,
  long_format = FALSE
)
```

## Arguments

- genotype_file:

  path to the snarl genotype file (sorted and indexed with tabix)

- chrom:

  chromosome name

- start_offset:

  start position

- end_offset:

  end position

- keep_nonvariant:

  should we keep variants with only one allele present? Default FALSE.

- long_format:

  should the table we returned in long form. Default FALSE.

## Value

a list of snarls with some information (boundary nodes) and the genotype
table

## Details

The list returned as one element per snarl. Each snarl element has a
\*GT\* data.frame with a \*sample\* column and a column for each alelle
(in the form al0, al1, etc). Each row in this \*GT\* data.frame
correspond to a sample and inform each allele's count. That table could
be combined with the phenotype and used for a regression test, for
example.
