# Matrix of mu2 coefficients for a data frame

Computes pairwise Guttman mu2 coefficients for all numeric columns of a
data frame, returning either a full symmetric correlation matrix or a
lower-triangular distance matrix (`1 - mu2`).

## Usage

``` r
mu2df(df, as_dist = FALSE)
```

## Arguments

- df:

  A data frame (or coercible object) with numeric columns.

- as_dist:

  If `TRUE`, return a lower-triangular
  [`dist`](https://rdrr.io/r/stats/dist.html) object of `1 - mu2`
  values. If `FALSE` (default), return a full symmetric correlation
  matrix.

## Value

A square numeric matrix with `dimnames` set to the column names of `df`
(when `as_dist = FALSE`), or a
[`dist`](https://rdrr.io/r/stats/dist.html) object with `Labels` set to
the column names (when `as_dist = TRUE`).

## Examples

``` r
df <- data.frame(a = 1:5, b = c(2, 3, 1, 5, 4), c = 5:1)
mu2df(df)
#>       a     b     c
#> a  1.00  0.75 -1.00
#> b  0.75  1.00 -0.75
#> c -1.00 -0.75  1.00
mu2df(df, as_dist = TRUE)
#>          a        b
#> b 0.500000         
#> c 1.414214 1.322876
```
