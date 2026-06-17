# Angular k-way (wedge) partition of a 2D configuration

Finds `k` rays (straight lines) emanating from a center point that
partition a 2D configuration of `n` points into `k` angular sectors,
minimising total misclassification (points whose group differs from the
majority group in their sector), and draws the lines. (see: Shye, S.
(2014). Faceted Smallest Space Analysis (FSSA). In A. Michalos (Ed.),
Encyclopedia of quality of life research (pp. 2129-2133). New York:
Springer.)

## Usage

``` r
angularPartition(
  crd,
  group,
  cx = NULL,
  cy = NULL,
  output = TRUE,
  col = "darkorange",
  lwd = 2,
  lty = 1,
  add = TRUE
)
```

## Arguments

- crd:

  Numeric matrix or numeric data frame with exactly 2 columns; no NAs.

- group:

  Factor with `n >= k >= 2` levels, same length as nrow of crd.

- cx, cy:

  Center; optimised when either is `NULL`.

- output:

  If `TRUE` (default), return list of results.

- col:

  Ray colour (default `"darkorange"`).

- lwd:

  Line width (default `2`).

- lty:

  Line type (default `1`).

- add:

  If `TRUE` (default), add to existing plot, else plot configuration.

## Value

If `output = TRUE`, a list with:

- `cuts` — `numeric[k]`, cut angles of rays (in radians).

- `margin` — minimum distance from any cut ray to nearest point.

- `misclass` — integer, number of misclassified points.

- `misclass_points` — data frame (`x`, `y`, `label`) of misclassified
  points

- `sector` — `integer[n]`, assigned sector `1..k` per point

- `majority` — `character[k]`, majority group per sector

- `center` — `c(cx, cy)` center coordinates of rays

- `pt_angles` — angles (in radians) of points from center

## Details

**Search at fixed center.** Points are sorted by their angle from the
center. A partition is a choice of `k` cut-gaps among the `n` gaps
between consecutive points (`combn(n, k)` candidates); per-arc majority
counts are read from a cumulative count table and total
misclassification is minimised. Cut angles are placed at the arc
midpoint between adjacent points. Tie-breaker: among k-tuples with the
same misclass count, the one with the largest minimum 2D perpendicular
distance from a cut ray to the nearest point wins.

**Search optimal center.** When `cx` and `cy` are `NULL` (default), the
center is optimised by multi-start Nelder-Mead — the brute-force above
runs as the inner objective at each candidate center. Starts are the
data centroid plus the centroid of each non-empty group. `parscale` is
set to the data range. When `cx` and `cy` are supplied, they are used
directly (no optimisation).

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(123)
theta <- rep(c(0, 2 * pi / 3, 4 * pi / 3), each = 12) + rnorm(36, 0, 0.25)
r     <- runif(36, 0.5, 1.5)
crd   <- cbind(r * cos(theta), r * sin(theta))
grp   <- factor(rep(c("a", "b", "c"), each = 12))
plot(crd, asp = 1)
angularPartition(crd, grp)
} # }
```
