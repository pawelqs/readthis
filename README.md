
# readthis

<!-- badges: start -->
<!-- badges: end -->

The goal of readthis is to provide comvenient functions for reading results of
algorithms such as CliP.


## Installation

You can install the development version of readthis from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("pawelqs/readthis")
```


## Currently supported algorithms and tools:

- [CliP](https://github.com/wwylab/CliP) is an algorithm for clonal structure identification through penalizing pairwise differences


## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(readthis)

dir <- system.file("CliP", package = "readthis")
read_clip_best_lambda(dir)
```

