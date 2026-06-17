# Guttman's weak monotonicity coefficient mu2

Computes Guttman's mu2 coefficient measuring the degree of weak monotone
relationship between two numeric vectors. Defined in the HUDAP Manual
(2001, p. 219). Pairs where either value is `NA` are excluded (pairwise
deletion).

## Usage

``` r
mu2(x, y)
```

## Arguments

- x:

  Numeric vector.

- y:

  Numeric vector of the same length as `x`.

## Value

A single numeric value in \[-1, 1\]. Returns `NaN` when all pairs are
tied (zero absolute differences).

## Examples

``` r
mu2(c(1, 2, 3, 4), c(2, 3, 1, 4))
#> [1] 0.5
```
