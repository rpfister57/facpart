# Binary radial partition (one separating circle)

Finds the circle minimising misclassification between two groups of 2D
points. The center is found by multi-start Nelder-Mead (starts: data
centroid plus each group's centroid) unless `cx` and `cy` are supplied.
Because the inner search is piecewise-constant, Nelder-Mead can stall on
a flat plateau around any of those starts; if none reaches zero
misclassification, a coarse `n_grid` x `n_grid` scan of the bounding
box, padded by half the data range on each side (the best center is not
always inside the convex hull of the points), locates a cell in a better
region and one more Nelder-Mead run is seeded there.

## Usage

``` r
radialCircle(
  crd,
  group,
  cx = NULL,
  cy = NULL,
  fill = FALSE,
  output = TRUE,
  col = "purple",
  cols = c("steelblue", "tomato"),
  lwd = 2,
  lty = 1,
  .method = "Nelder-Mead",
  n_grid = 7L,
  add = TRUE
)
```

## Arguments

- crd:

  Numeric matrix or data frame with exactly 2 columns.

- group:

  Factor with exactly 2 levels.

- cx, cy:

  Center of the separating circle; optimised when `NULL` (default).

- fill:

  If `TRUE`, shade the inner disc and outer region (default `FALSE`).

- output:

  If `TRUE` (default), return results list.

- col:

  Circle border colour (default `"purple"`).

- cols:

  Length-2 fill colours (default `c("steelblue", "tomato")`).

- lwd:

  Line width (default `2`).

- lty:

  Line type (default `1`).

- .method:

  `"Nelder-Mead"` (default) or `"SANN"`.

- n_grid:

  Side length of the coarse fallback grid (`n_grid^2` extra evaluations
  of the cheap inner search); only used when the heuristic starts don't
  already reach zero misclassification. Default `7L`.

- add:

  If `TRUE` (default), add to existing plot; if `FALSE`, call
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) first.

## Value

If `output = TRUE`, a list with `partition` (`"circle"`), `center`,
`radius`, `sector` (`1` inside / `2` outside per point), `majority`
(`character[2]`), `pcoords` (the input `crd`), `pgroup` (`group` coerced
to factor), `misclass` (integer), and `misclass_points` (data frame with
columns `x`, `y`, `label` for each misclassified point).

## Details

**Search at fixed center.** Points are sorted by distance from the
center and every realisable radius is scanned: the inside/outside split
is scored by its exact bijection-constrained misclassification (the
better of the two ways to match the two regions to the two groups) — the
same criterion used to derive `sector`/`majority` below, so the search
targets exactly the quantity reported as `misclass`, and neither group
is assumed to be the inner one. Radii that separate coincident distances
are skipped as unrealisable. The scan also considers the degenerate
radii that leave one region empty; when the configuration has no radial
structure these can be the true minimum, so they are allowed, but on a
tie a circle that really does split the points is preferred.
[`radialCircles()`](https://rpfister57.github.io/facpart/reference/radialCircles.md)
with `k = 2` and the same center returns the same circle.

## Examples

``` r
set.seed(1)
inner <- cbind(rnorm(20, 0, 0.3), rnorm(20, 0, 0.3))
th <- runif(20, 0, 2 * pi)
outer <- cbind(2 * cos(th), 2 * sin(th)) + matrix(rnorm(40, 0, 0.1), 20)
crd <- rbind(inner, outer)
grp <- factor(c(rep("in", 20), rep("out", 20)))
radialCircle(crd, grp, fill = TRUE, add = FALSE)

#> $partition
#> [1] "circle"
#> 
#> $center
#> [1] -0.03723659 -0.02506890
#> 
#> $radius
#> [1] 1.270447
#> 
#> $sector
#>  [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
#> [39] 2 2
#> 
#> $majority
#> [1] "in"  "out"
#> 
#> $pcoords
#>               [,1]        [,2]
#>  [1,] -0.187936143  0.27569321
#>  [2,]  0.055092997  0.23464089
#>  [3,] -0.250688584  0.02236950
#>  [4,]  0.478584241 -0.59680551
#>  [5,]  0.098852332  0.18594772
#>  [6,] -0.246140515 -0.01683862
#>  [7,]  0.146228716 -0.04673865
#>  [8,]  0.221497412 -0.44122572
#>  [9,]  0.172734405 -0.14344502
#> [10,] -0.091616516  0.12538247
#> [11,]  0.453534351  0.40760387
#> [12,]  0.116952971 -0.03083632
#> [13,] -0.186372174  0.11630148
#> [14,] -0.664409966 -0.01614151
#> [15,]  0.337479275 -0.41311787
#> [16,] -0.013480083 -0.12449837
#> [17,] -0.004857079 -0.11828699
#> [18,]  0.283150863 -0.01779402
#> [19,]  0.246366359  0.33000761
#> [20,]  0.178170396  0.22895272
#> [21,] -1.793995052  0.84577202
#> [22,] -0.527914003 -2.01577754
#> [23,] -1.583880425  1.23670039
#> [24,] -1.024858029  1.68658989
#> [25,]  0.232332669 -2.12338076
#> [26,]  0.783811064  1.94143964
#> [27,] -0.520442649 -1.98495138
#> [28,]  1.338887320  1.38462592
#> [29,]  0.113657376  2.00663066
#> [30,]  1.229075223  1.50820879
#> [31,]  0.370390188  1.93888876
#> [32,]  1.860516813  0.71026419
#> [33,] -1.183586855 -1.44138707
#> [34,]  1.428246617 -1.55524747
#> [35,]  0.287029709 -1.90768980
#> [36,]  0.604663293 -1.87899601
#> [37,] -2.102042587  0.66097942
#> [38,] -1.542666287  1.04034293
#> [39,]  0.761731226 -1.81849758
#> [40,] -1.363562518 -1.19844173
#> 
#> $pgroup
#>  [1] in  in  in  in  in  in  in  in  in  in  in  in  in  in  in  in  in  in  in 
#> [20] in  out out out out out out out out out out out out out out out out out out
#> [39] out out
#> Levels: in out
#> 
#> $misclass
#> [1] 0
#> 
#> $misclass_points
#> [1] x     y     label
#> <0 rows> (or 0-length row.names)
#> 
```
