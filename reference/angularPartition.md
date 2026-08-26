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
  n_grid = 7L,
  fill = FALSE,
  output = TRUE,
  col = "darkorange",
  cols = NULL,
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

- n_grid:

  Side length of the coarse fallback grid (`n_grid^2` extra evaluations
  of the cheap inner search); only used when the heuristic starts don't
  already reach zero misclassification. Default `7L`.

- fill:

  If `TRUE`, shade the `k` wedge sectors (default `FALSE`).

- output:

  If `TRUE` (default), return list of results.

- col:

  Ray colour (default `"darkorange"`).

- cols:

  Length-`k` fill colours; auto-generated if `NULL`.

- lwd:

  Line width (default `2`).

- lty:

  Line type (default `1`).

- add:

  If `TRUE` (default), add to existing plot, else plot configuration.

## Value

If `output = TRUE`, a list with:

- `partition` — `"angular"`

- `cuts` — `numeric[k]`, cut angles of rays (in radians).

- `margin` — minimum distance from any cut ray to nearest point.

- `sector` — `integer[n]`, assigned sector `1..k` per point

- `majority` — `character[k]`, majority group per sector

- `center` — `c(cx, cy)` center coordinates of rays

- `pt_angles` — angles (in radians) of points from center

- `pcoords` — the input `crd`

- `pgroup` — `group` coerced to factor

- `misclass` — integer, number of misclassified points.

- `misclass_points` — data frame (`x`, `y`, `label`) of misclassified
  points

## Details

**Search at fixed center.** Points are sorted by their angle from the
center. A partition is a choice of `k` cut-gaps among the `n` gaps
between consecutive points (`combn(n, k)` candidates, `M` in total). The
score of a candidate is its exact bijection-constrained
misclassification (best achievable total correct over all `k!`
arc-to-group assignments, each group used exactly once) – the same
region-to-group matching criterion used to derive `sector`/`majority`
below, so the search targets exactly the quantity reported as
`misclass`. Every wedge is non-empty, so a partition leaving a group's
sector empty is not considered (as in
[`axialLines()`](https://rpfister57.github.io/facpart/reference/axialLines.md);
the radial functions do allow it).

All per-candidate quantities are computed as length-`M` vectors indexed
out of a cumulative group-count table, so the loops run over `k` / `k!`
rather than over candidates – analogous to the vectorised cut search in
[`axialLines()`](https://rpfister57.github.io/facpart/reference/axialLines.md).
Scoring proceeds in two stages: a cheap upper bound (each arc
independently taking its own local-majority group, which no bijection
can beat) is evaluated for all `M` candidates, and only those whose
bound is high enough to still win are then scored exactly over all `k!`
assignments. The bound is tight in practice, so typically only a handful
of candidates need the exact scan.

Cut angles are placed at the arc midpoint between adjacent points.
Tie-breaker: among k-tuples with the same misclass count, the one with
the largest minimum 2D perpendicular distance from a cut ray to the
nearest point wins (computed only for the co-minimal candidates).

**Search optimal center.** When `cx` and `cy` are `NULL` (default), the
center is optimised by multi-start Nelder-Mead — the brute-force above
runs as the inner objective at each candidate center. Starts are the
data centroid plus the centroid of each non-empty group and the four
bounding-box corners. `parscale` is set to the data range. Because the
inner objective is piecewise-constant, Nelder-Mead can get stuck on a
flat plateau around any of those starts without ever seeing a better
region; if none of them reaches zero misclassification, a coarse
`n_grid` x `n_grid` scan of the bounding box, padded by half the data
range on each side (the best center is not always inside the convex hull
of the points), locates a cell in a better region and one more
Nelder-Mead run is seeded there. When `cx` and `cy` are supplied, they
are used directly (no optimisation).

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(123)
theta <- rep(c(0, 2 * pi / 3, 4 * pi / 3), each = 12) + rnorm(36, 0, 0.25)
r     <- runif(36, 0.5, 1.5)
crd   <- cbind(r * cos(theta), r * sin(theta))
grp   <- factor(rep(c("a", "b", "c"), each = 12))
plot(crd, asp = 1)
angularPartition(crd, grp, fill = TRUE)
} # }
```
