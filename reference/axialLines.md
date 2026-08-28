# Axial partitions for k \>= 2 groups

Partitions a 2D configuration into `k >= 2` groups using `k-1` parallel
separating lines. Both the slopes (angles) of the lines and their
positions (cuts) are searched to minimise empirical misclassification.

## Usage

``` r
axialLines(
  crd,
  group,
  fill = FALSE,
  output = TRUE,
  col = "purple",
  cols = NULL,
  lwd = 2,
  lty = 1,
  n_angles = 180L,
  add = TRUE
)
```

## Arguments

- crd:

  Numeric matrix or data frame with exactly 2 columns.

- group:

  Factor, character, or integer vector with `k >= 2` levels.

- fill:

  If `TRUE`, shade the `k` strips between consecutive lines.

- output:

  If `TRUE` (default), return a list of results.

- col:

  Separator line colour (default `"purple"`).

- cols:

  Length-`k` fill colours; auto-generated if `NULL`.

- lwd:

  Line width (default `2`).

- lty:

  Line type (default `1`).

- n_angles:

  Number of angles in `[0, pi)` to scan (default `180L`, approximately
  1-degree resolution); LDA's direction is added as an extra candidate.

- add:

  If `TRUE` (default), add to existing plot; if `FALSE`, call
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) first.

## Value

If `output = TRUE`, a list with:

- `partition` — `"axial"`

- `slope` — shared slope of the lines (`Inf` if vertical)

- `intercepts` — `numeric[k-1]`, y-intercepts (or x-positions if
  vertical)

- `angle` — winning direction in radians

- `margin` — minimum perpendicular distance from the cut line(s) to the
  nearest point

- `sector` — `integer[n]`, sector `1..k` per point in input order

- `majority` — `character[k]`, majority group per sector

- `pcoords` — the input `crd`

- `pgroup` — `group` coerced to factor

- `misclass` — total misclassified points

- `misclass_points` — data frame (`x`, `y`, `label`) of misclassified
  points

## Details

A brute-force search over a grid of `n_angles` and all `C(n-1, k-1)`
cuts of n points into k groups finds the `(angle, cuts)` combination
with lowest total misclassification. Tie-breaker: among configurations
with the same misclass count, the one with the largest minimum margin
(perpendicular distance from a cut line to the nearest point) is
preferred.

For `k = 2`, `axialLines()` can differ from
[`axialLine()`](https://rpfister57.github.io/facpart/reference/axialLine.md).

Vertical-line guard: when the separators are vertical, `slope` is
returned as `Inf` and `intercepts` carry the x-positions of the lines.

## See also

[`axialLine()`](https://rpfister57.github.io/facpart/reference/axialLine.md)
for the classical LDA boundary if k=2.

## Examples

``` r
set.seed(1)
crd <- rbind(cbind(rnorm(10, -1), rnorm(10)),
             cbind(rnorm(10,  0), rnorm(10)),
             cbind(rnorm(10,  1), rnorm(10)))
grp <- factor(c(rep("a", 10), rep("b", 10), rep("c", 10)))
axialLines(crd, grp, fill = TRUE, add = FALSE)

#> $partition
#> [1] "axial"
#> 
#> $slope
#> [1] -4.70463
#> 
#> $intercepts
#> [1] -1.162032  2.813375
#> 
#> $angle
#> [1] 0.2094395
#> 
#> $margin
#> [1] 0.002606886
#> 
#> $sector
#>  [1] 1 1 1 2 1 1 1 2 1 1 3 3 2 1 2 2 2 1 2 2 3 3 3 3 3 3 3 3 3 3
#> 
#> $majority
#> [1] "a" "b" "c"
#> 
#> $pcoords
#>              [,1]        [,2]
#>  [1,] -1.62645381  1.51178117
#>  [2,] -0.81635668  0.38984324
#>  [3,] -1.83562861 -0.62124058
#>  [4,]  0.59528080 -2.21469989
#>  [5,] -0.67049223  1.12493092
#>  [6,] -1.82046838 -0.04493361
#>  [7,] -0.51257095 -0.01619026
#>  [8,] -0.26167529  0.94383621
#>  [9,] -0.42421865  0.82122120
#> [10,] -1.30538839  0.59390132
#> [11,]  0.91897737  1.35867955
#> [12,]  0.78213630 -0.10278773
#> [13,]  0.07456498  0.38767161
#> [14,] -1.98935170 -0.05380504
#> [15,]  0.61982575 -1.37705956
#> [16,] -0.05612874 -0.41499456
#> [17,] -0.15579551 -0.39428995
#> [18,] -1.47075238 -0.05931340
#> [19,] -0.47815006  1.10002537
#> [20,]  0.41794156  0.76317575
#> [21,]  0.83547640  0.39810588
#> [22,]  0.74663832 -0.61202639
#> [23,]  1.69696338  0.34111969
#> [24,]  1.55666320 -1.12936310
#> [25,]  0.31124431  1.43302370
#> [26,]  0.29250484  1.98039990
#> [27,]  1.36458196 -0.36722148
#> [28,]  1.76853292 -1.04413463
#> [29,]  0.88765379  0.56971963
#> [30,]  1.88110773 -0.13505460
#> 
#> $pgroup
#>  [1] a a a a a a a a a a b b b b b b b b b b c c c c c c c c c c
#> Levels: a b c
#> 
#> $misclass
#> [1] 6
#> 
#> $misclass_points
#>            x           y label
#> 1  0.5952808 -2.21469989     a
#> 2 -0.2616753  0.94383621     a
#> 3  0.9189774  1.35867955     b
#> 4  0.7821363 -0.10278773     b
#> 5 -1.9893517 -0.05380504     b
#> 6 -1.4707524 -0.05931340     b
#> 
```
