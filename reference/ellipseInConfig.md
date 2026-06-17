# Fit and draw an ellipse around a 2D configuration

Three fitting strategies are available:

## Usage

``` r
ellipseInConfig(
  crd,
  labs = NULL,
  mid = "centroid",
  typ = "maxdist",
  output = TRUE,
  col = "black",
  fill = NA,
  lwd = 2,
  lty = 2,
  add = TRUE
)
```

## Arguments

- crd:

  Numeric matrix or data frame with exactly 2 columns (x, y).

- labs:

  Currently unused (kept for backward compatibility).

- mid:

  `"centroid"` (default) places the ellipse center at the data centroid;
  `"nullnull"` forces the center to the origin (0, 0). Only affects
  `"maxdist"`.

- typ:

  Fitting strategy: `"maxdist"` (default), `"bestfit"`, or `"minbound"`.

- output:

  If `TRUE` (default), return a list describing the ellipse.

- col:

  Border colour (default `"black"`).

- fill:

  Fill colour, or `NA` (default) for no fill.

- lwd:

  Line width (default `2`).

- lty:

  Line type (default `2`).

- add:

  If `TRUE` (default), add to existing plot; if `FALSE`, call
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) first.

## Value

If `output = TRUE`, a list with elements `cx`, `cy`, `a`, `b`, `angle`
(in radians).

## Details

- `"maxdist"` — covariance-based ellipse, scaled so all points lie
  inside with a 5% buffer. Semi-axes are `sqrt(eigenvalue) * k`, where
  `k` is the maximum Mahalanobis distance from the center.

- `"bestfit"` — minimises the sum of squared normalised elliptic
  distances `sqrt((u/a)^2 + (v/b)^2) - 1` via BFGS, starting from the
  covariance ellipse.

- `"minbound"` — minimum-area enclosing ellipse via
  [`cluster::ellipsoidhull()`](https://rdrr.io/pkg/cluster/man/ellipsoidhull.html).

The ellipse is drawn on the currently active plot window.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
crd <- cbind(rnorm(30), 0.5 * rnorm(30))
ellipseInConfig(crd, typ = "maxdist", add = FALSE)
ellipseInConfig(crd, typ = "bestfit", col = "blue")
ellipseInConfig(crd, typ = "minbound", col = "red")
} # }
```
