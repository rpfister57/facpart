# K-group radial partition (nested circles)

Generalises
[`radialCircle()`](https://rpfister57.github.io/facpart/reference/radialCircle.md)
to `k >= 2` groups using `k-1` nested (inclusive) circles. The circles
must be nested (circle `s` contains circle `s-1`). Circles are fitted
sequentially, each minimising misclassification of the innermost `s`
groups (inside) vs the rest (outside) — but **which** groups those are
is searched, not read off the factor levels, so the inside-to-outside
ordering is found from the data and the result does not depend on how
the groups are named.

## Usage

``` r
radialCircles(
  crd,
  group,
  cx = NULL,
  cy = NULL,
  fill = FALSE,
  output = TRUE,
  col = "purple",
  cols = NULL,
  lwd = 2,
  lty = 1,
  .method = "Nelder-Mead",
  add = TRUE
)
```

## Arguments

- crd:

  Numeric matrix or data frame with exactly 2 columns.

- group:

  Factor with `k >= 2` levels. **Factor level order does not matter**:
  the inside-to-outside nesting order is found from the data, so
  relabelling or reordering the levels cannot change the result. Read
  the order off `majority`, which names the group owning each region
  from the innermost outwards.

- cx, cy:

  Optional shared center for all `k-1` circles. When both are supplied
  the circles are concentric at `(cx, cy)`; otherwise each center is
  optimised.

- fill:

  If `TRUE`, shade each ring sector.

- output:

  If `TRUE` (default), return results list.

- col:

  Circle border colour (default `"purple"`).

- cols:

  Length-`k` colour vector; auto-generated if `NULL`.

- lwd:

  Line width (default `2`).

- lty:

  Line type (default `1`).

- .method:

  `"Nelder-Mead"` (default) or `"SANN"`; ignored when `cx` and `cy` are
  both supplied.

- add:

  If `TRUE` (default), add to existing plot; if `FALSE`, call
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) first.

## Value

If `output = TRUE`, a list with `cx`, `cy` (`numeric[k-1]`), `radii`
(innermost to outermost), `misclass`, `misclass_points` (data frame with
columns `x`, `y`, `label` for each misclassified point), `sector`
(`integer[n]`), and `majority` (`character[k]`).

## Details

When `cx` and `cy` are both supplied, all `k-1` circles share that
center (concentric); only the radii are searched, subject to
`r_s >= r_{s-1}`. When either is `NULL` (default), each circle's center
is optimised independently by multi-start Nelder-Mead — the circles are
nested but not generally concentric.

**What is minimised.** Every candidate partition is scored by its exact
bijection-constrained misclassification: the best total correct over all
`k!` ways of matching the `k` regions to the `k` groups, each group used
once. This is the same criterion `sector`/`majority` are derived from,
so the search targets exactly the quantity reported as `misclass`. Radii
that would separate coincident distances are skipped as unrealisable,
and radii leaving a region empty are allowed — on data with no radial
structure such a partition can be the true minimum — but a partition
whose every circle actually splits the points wins any tie.

**Optimality.** Concentric (`cx` and `cy` supplied): all `k-1` radii are
searched *jointly* over every realisable combination, so the result is
globally optimal for that center. If the number of combinations is too
large to enumerate, the function warns and falls back to the sequential
path below. Optimised centers: circles are fitted sequentially, then all
radii are refined by coordinate descent against the criterion above, so
the result is optimal in each individual radius but not jointly in the
centers.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
g1 <- cbind(rnorm(15, 0, 0.2), rnorm(15, 0, 0.2))
th2 <- runif(15, 0, 2 * pi)
g2 <- cbind(1.2 * cos(th2), 1.2 * sin(th2)) + matrix(rnorm(30, 0, 0.1), 15)
th3 <- runif(15, 0, 2 * pi)
g3 <- cbind(2.5 * cos(th3), 2.5 * sin(th3)) + matrix(rnorm(30, 0, 0.1), 15)
crd <- rbind(g1, g2, g3)
grp <- factor(c(rep("a", 15), rep("b", 15), rep("c", 15)))
radialCircles(crd, grp, fill = TRUE, add = FALSE)
} # }
```
