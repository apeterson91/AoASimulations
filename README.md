# STAPSimulations 

<!-- badges: start -->
<!-- badges: end -->

The goal of STAPSimulations is to reproduce the tables in the Peterson et al. paper concerning Spatial Temporal Aggregated Predictors.

## Installation

You can install the released version of STAPSimulations from [github](https://github.com) with:

``` r
devtools::install_github("https://github.com/apeterson91/STAPSimulations")
```

## Example

Both table one and table two can be created easily by calling their appropriate functions.
Default values of parameters will be used unless specified otherwise and, unless a filepath is 
specified, tex versions of the tables will be saved to the user's desktop. Be aware
that table 1 in the paper is a combination of both table 1 a and table 1b in this package.

``` r
library(STAPSimulations)
create_table_one()
create_table_two()
```

