# Binary radial-elliptic partition (one separating ellipse)

Finds the ellipse minimising misclassification between two groups of 2D
points. All 5 parameters (cx, cy, a, b, angle) are found by Nelder-Mead,
starting from the covariance ellipse of each group in turn; the ordering
with the lower misclassification is retained. Because the inner search
is piecewise-constant, Nelder-Mead can stall on a flat plateau around
either start (confirmed on real MDS data: a group's own covariance-based
init can itself be near enough to the configuration centroid to land in
one); if a group's start doesn't reach zero misclassification, a coarse
`n_grid` x `n_grid` scan relocates just the center (holding the fitted
shape and orientation fixed) over the bounding box, padded by half the
data range on each side, and one more full 5-parameter Nelder-Mead run
is seeded from the best cell found.

## Usage

``` r
radialEllipse(
  crd,
  group,
  n_grid = 7L,
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

- n_grid:

  Side length of the coarse fallback grid (`n_grid^2` extra evaluations
  of the cheap inner search per group tried as inner); only used when a
  group's covariance-based start doesn't already reach zero
  misclassification. Default `7L`.

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

If `output = TRUE`, a list with `partition` (`"ellipse"`), `cx`, `cy`,
`a`, `b`, `angle` (radians), `sector`, `majority`, `pcoords` (the input
`crd`), `pgroup` (`group` coerced to factor), `misclass` (list with `n`
and `indices`), and `misclass_points` (data frame with columns `x`, `y`,
`label` for each misclassified point).

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
