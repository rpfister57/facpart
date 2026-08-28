# Project a point (x, y) onto a line (a, b)

On an existing plot including a line with a = intercept and b = slope, a
point with coordinates (x, y) is projected. Optionally, the line as well
as the projection arrow and the distance from intercept to the projected
point are drawn.

## Usage

``` r
projectP2L(Pxy, Lab, add = FALSE, d0 = FALSE, color = "blue", addLine = FALSE)
```

## Arguments

- Pxy:

  A vector c(x,y) with coordinates of a point.

- Lab:

  A vector c(a,b) with intercept (a) and slope (b) of a line.

- add:

  Logical: Should the projection be drawn (default = FALSE).

- d0:

  Logical: Should the distance line from a to projection be drawn
  (default = FALSE).

- color:

  Character - a color name

- addLine:

  Logical: Should the line be drawn.

## Value

A list with the coordinates of the projected point, and the distance
vector of the projected point.

## Examples

``` r
if (FALSE) { # \dontrun{
plot(-2:2, -2:2, type = "n", asp = 1)
abline(a = 0, b = 1)
aPoint <- matrix(c(-1.5, 1), nrow = 1, byrow = TRUE)
points(aPoint, pch = 19)
projectP2L(Pxy = aPoint, Lab = c(0, 1), add = TRUE)
} # }
```
