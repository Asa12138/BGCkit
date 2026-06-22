# Read files from bigscape directory

Read files from bigscape directory

## Usage

``` r
read_big_scape_dir(
  big_scape_dir = "2024-05-14_16-00-37_hybrids_glocal",
  cutoff = NULL,
  reassign_GCF = TRUE
)
```

## Arguments

- big_scape_dir:

  Directory of bigscape output such as
  2024-05-14_16-00-37_hybrids_glocal

- cutoff:

  default NULL, set this when you have multiple cutoffs

- reassign_GCF:

  default TRUE. When some BGCs were assigned into different GCFs
  ,reassign them into the GCF with the largest number of BGCs.

## Value

big_scape object
