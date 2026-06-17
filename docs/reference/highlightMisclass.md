# Highlight misclassified points on the active plot

Overlays a marker on each misclassified point, using the
`misclass_points` data frame returned by
[`axialLine()`](https://rpfister57.github.io/facpart/reference/axialLine.md),
[`axialLines()`](https://rpfister57.github.io/facpart/reference/axialLines.md),
[`radialCircle()`](https://rpfister57.github.io/facpart/reference/radialCircle.md),
[`radialCircles()`](https://rpfister57.github.io/facpart/reference/radialCircles.md),
[`radialEllipse()`](https://rpfister57.github.io/facpart/reference/radialEllipse.md),
[`radialEllipses()`](https://rpfister57.github.io/facpart/reference/radialEllipses.md),
or
[`angularPartition()`](https://rpfister57.github.io/facpart/reference/angularPartition.md).
The full result list may also be passed directly; `misclass_points` is
extracted automatically.

## Usage

``` r
highlightMisclass(x, col = "red", pch = 4, cex = 1.5, lwd = 2)
```

## Arguments

- x:

  A data frame with columns `x`, `y`, `label` (the `misclass_points`
  element of a partition result), or the full partition result list.

- col:

  Marker colour (default `"red"`).

- pch:

  Point character for the marker (default `4`, an X).

- cex:

  Point size expansion factor (default `1.5`).

- lwd:

  Line width of the marker symbol (default `2`).

## Value

`invisible(NULL)`, called for its side effect of drawing on the active
plot.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
crd <- rbind(cbind(rnorm(15, -1), rnorm(15)),
             cbind(rnorm(15,  1), rnorm(15)))
grp <- factor(c(rep("a", 15), rep("b", 15)))
res <- axialLine(crd, grp, add = FALSE)
highlightMisclass(res)
} # }
```
