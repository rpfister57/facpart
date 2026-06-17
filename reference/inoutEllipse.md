# Classify points as inside or outside a given ellipse

Classify points as inside or outside a given ellipse

## Usage

``` r
inoutEllipse(crd, cx, cy, a, b, angle)
```

## Arguments

- crd:

  Numeric matrix or data frame with exactly 2 columns (x, y).

- cx, cy:

  Center of the ellipse.

- a, b:

  Semi-major and semi-minor axes.

- angle:

  Rotation of the major axis in radians (as returned by
  [`ellipseInConfig()`](https://rpfister57.github.io/facpart/reference/ellipseInConfig.md)).

## Value

A character vector of length `nrow(crd)` with values `"inside"` or
`"outside"`.

## Examples

``` r
if (FALSE) { # \dontrun{
crd <- cbind(rnorm(30), rnorm(30))
inoutEllipse(crd, cx = 0, cy = 0, a = 1, b = 1, angle = 0)
} # }
```
