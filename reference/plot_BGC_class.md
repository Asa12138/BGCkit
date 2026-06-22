# Plot BGC class

Plot BGC class

## Usage

``` r
plot_BGC_class(big_scape_res, mode = 1, rm_mibig = FALSE, ...)
```

## Arguments

- big_scape_res:

  big_scape object from [`read_big_scape_dir()`](read_big_scape_dir.md)

- mode:

  1~2, 1 for doughnut plot, 2 for sankey plot

- rm_mibig:

  logical, remove MIBiG BGCs from the plot

- ...:

  additional parameters for `gghuan` or `my_sankey`, such as `topN=10`

## Value

ggplot
