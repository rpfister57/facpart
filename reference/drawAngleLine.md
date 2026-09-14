# Draw a line of given angle and length from a starting point

On the currently active plot, draws a straight line starting at `(x, y)`
running in the direction `angle` (in **radians**, counter-clockwise from
the positive x-axis) with length `len`. Optionally marks the endpoint
with a point. Useful for sketching cut rays of an angular partition (see
[`angularPartition()`](https://rpfister57.github.io/facpart/reference/angularPartition.md))
or any other direction indicator.

## Usage

``` r
drawAngleLine(x = 0, y = 0, angle, len = 1, ppoint = FALSE, ...)
```

## Arguments

- x, y:

  Coordinates of the starting point (default `0`, `0`).

- angle:

  Direction of the line in radians.

- len:

  Length of the line (default `1`).

- ppoint:

  Logical: Should the endpoint be marked with a filled point (default
  `FALSE`).

- ...:

  Further graphical parameters passed to
  [`graphics::segments()`](https://rdrr.io/r/graphics/segments.html),
  e.g. `col`, `lwd`, `lty`.

## Value

`invisible(NULL)`, called for its side effect of drawing on the active
plot.

## Examples

``` r
if (FALSE) { # \dontrun{
plot(-2:2, -2:2, type = "n", asp = 1)
drawAngleLine(angle = pi / 4, len = 2, ppoint = TRUE, col = "blue")
drawAngleLine(x = -1, y = -1, angle = 0, len = 1.5, lty = 2)
} # }
```
