# Replot a stored partition result

Redraws the configuration and partitions from the output list returned
by
[`axialLine()`](https://rpfister57.github.io/facpart/reference/axialLine.md),
[`axialLines()`](https://rpfister57.github.io/facpart/reference/axialLines.md),
[`radialCircle()`](https://rpfister57.github.io/facpart/reference/radialCircle.md),
[`radialCircles()`](https://rpfister57.github.io/facpart/reference/radialCircles.md),
[`radialEllipse()`](https://rpfister57.github.io/facpart/reference/radialEllipse.md),
[`radialEllipses()`](https://rpfister57.github.io/facpart/reference/radialEllipses.md)
or
[`angularPartition()`](https://rpfister57.github.io/facpart/reference/angularPartition.md).
It is a simple replot of the stored partition regions with no
recomputation.

## Usage

``` r
plotPartition(
  Pout,
  fill = FALSE,
  cols = NULL,
  col = NULL,
  lwd = 2,
  lty = 1,
  highlight = FALSE,
  add = FALSE
)
```

## Arguments

- Pout:

  A result list from one of the partition functions above (must have a
  `partition` field naming a known family ("axial", "angular", "circle",
  "ellipse"), plus `pcoords`/`pgroup`).

- fill:

  If `TRUE`, shade the regions between boundaries (default `FALSE`).

- cols:

  Fill colours; auto-generated if `NULL` – `c("steelblue", "tomato")`
  for a 2-region partition, `hcl.colors(k, "Pastel 1")` otherwise,
  matching the palette the originating partition function would have
  used by default. A char vector with k valid colors can be supplied.

- col:

  Boundary line/ray colour; if `NULL` (default), `"darkorange"` for
  `partition == "angular"` and `"purple"` otherwise, matching the
  originating functions' own defaults.

- lwd:

  Line width (default `2`).

- lty:

  Line type (default `1`).

- highlight:

  If `TRUE`, overlay a marker on every misclassified point via
  [`highlightMisclass()`](https://rpfister57.github.io/facpart/reference/highlightMisclass.md)
  (default `FALSE`).

- add:

  If `FALSE` (default), start a fresh plot of `Pout$pcoords` with group
  labels before drawing the boundaries; if `TRUE`, draw only the
  boundaries (and optional fill/highlight) on the already-active plot.

## Value

`invisible(NULL)`, called for its side effect of drawing on the active
plot.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
crd <- rbind(cbind(rnorm(15, -1), rnorm(15)), cbind(rnorm(15, 1), rnorm(15)))
grp <- factor(c(rep("a", 15), rep("b", 15)))
res <- radialCircles(crd, grp, output = TRUE, add = FALSE)
plotPartition(res, fill = TRUE, highlight = TRUE)
} # }
```
