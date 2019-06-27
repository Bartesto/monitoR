
<!-- README.md is generated from README.Rmd. Please edit that file -->

# monitoR

<!-- badges: start -->

<!-- badges: end -->

The goal of monitoR is to assist with the visualisation of monitoring
data from the Swan and Canning Rivers.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Bartesto/monitoR")
```

## Example

A basic example:

``` r
library(monitoR)
# To create surfer plots for the SWan River - please see help for parameter descriptions
swan_surfR(path = "C:/some path to sonde spreadsheets", ovit = "green", ocav = "red")
```
