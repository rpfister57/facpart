# Binary radial-elliptic partition (one separating ellipse)

Finds the ellipse minimising misclassification between two groups of 2D
points. All 5 parameters (cx, cy, a, b, angle) are found by Nelder-Mead,
starting from the covariance ellipse of each group in turn; the ordering
with the lower misclassification is retained.

## Usage

``` r
radialEllipse(
  crd,
  group,
  fill = FALSE,
  output = TRUE,
  col = "purple",
  cols = c("steelblue", "tomato"),
  lwd = 2,
  lty = 1,
  add = TRUE
)
```

## Arguments

- crd:

  Numeric matrix or data frame with exactly 2 columns.

- group:

  Factor with exactly 2 levels.

- fill:

  If `TRUE`, shade inner ellipse and outer region.

- output:

  If `TRUE` (default), return results list.

- col:

  Ellipse border colour (default `"purple"`).

- cols:

  Length-2 fill colours (default `c("steelblue", "tomato")`).

- lwd:

  Line width (default `2`).

- lty:

  Line type (default `1`).

- add:

  If `TRUE` (default), add to existing plot; if `FALSE`, call
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) first.

## Value

If `output = TRUE`, a list with `cx`, `cy`, `a`, `b`, `angle` (radians),
`misclass` (list with `n` and `indices`), `misclass_points` (data frame
with columns `x`, `y`, `label` for each misclassified point), `sector`,
and `majority`.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
inner <- cbind(rnorm(20, 0, 0.3), rnorm(20, 0, 0.2))
th    <- runif(20, 0, 2 * pi)
outer <- cbind(2.0 * cos(th), 1.2 * sin(th)) + matrix(rnorm(40, 0, 0.1), 20)
crd <- rbind(inner, outer)
grp <- factor(c(rep("in", 20), rep("out", 20)))
radialEllipse(crd, grp, fill = TRUE, add = FALSE)
} # }
```
