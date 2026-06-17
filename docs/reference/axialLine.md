# Binary LDA separating line for a 2D point configuration

Fits a linear discriminant analysis to two-group 2D data and draws the
classical LDA boundary — perpendicular to LD1 through the midpoint of
the class means — on the current plot.

## Usage

``` r
axialLine(
  crd,
  group,
  fill = FALSE,
  output = TRUE,
  col = "purple",
  cols = c("steelblue", "tomato"),
  lwd = 2,
  lty = 1,
  add = TRUE
)
```

## Arguments

- crd:

  Numeric matrix or data frame with exactly 2 columns (x, y).

- group:

  Factor, character, or integer vector with exactly 2 levels.

- fill:

  If `TRUE`, shade the two half-planes (default `FALSE`).

- output:

  If `TRUE` (default), return a list of line parameters.

- col:

  Separator line colour (default `"purple"`).

- cols:

  Length-2 fill colours; used when `fill = TRUE` (default
  `c("steelblue", "tomato")`).

- lwd:

  Line width (default `2`).

- lty:

  Line type (default `1`).

- add:

  If `TRUE` (default), add to existing plot; if `FALSE`, call
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) first.

## Value

If `output = TRUE`, a list with `slope`, `intercept`, `misclass`
(integer), `misclass_points` (data frame with columns `x`, `y`, `label`
for each misclassified point), and `predicted` (LDA-predicted class
factor).

## Details

When LD1 is along the x-axis (`abs(w[2]) < 1e-10`) the separator is
vertical: `slope` is returned as `Inf` and `intercept` carries the
x-position of the line.

**Not a special case of
[`axialLines()`](https://rpfister57.github.io/facpart/reference/axialLines.md).**
`axialLine()` returns the classical LDA Bayes-rule boundary
(closed-form: LD1 direction, midpoint of class means as cut) — optimal
under multivariate-normal classes with equal covariances, and the
standard reference object in psychometric / Facet Theory pipelines. The
return value includes `predicted`, the LDA class assignment per point.

[`axialLines()`](https://rpfister57.github.io/facpart/reference/axialLines.md)
with `k = 2`, in contrast, searches the full `(angle, cut)` space for
the empirical minimum-misclass parallel line. Its direction need not be
LD1 and its cut need not match the midpoint of means; it can overfit the
training sample when classes are non-Gaussian or have unequal
covariances. Use it when you want the empirically best linear separator;
use `axialLine()` when you want the classical LDA classifier.

## See also

[`axialLines()`](https://rpfister57.github.io/facpart/reference/axialLines.md)
for the empirically-optimal parallel-line partition over both angle and
cut position (any `k >= 2`).

## Examples

``` r
set.seed(1)
crd <- rbind(cbind(rnorm(15, -1), rnorm(15)),
             cbind(rnorm(15,  1), rnorm(15)))
grp <- factor(c(rep("a", 15), rep("b", 15)))
axialLine(crd, grp, fill = TRUE, add = FALSE)

#> $slope
#> [1] -7.867212
#> 
#> $intercept
#>         2 
#> 0.8714933 
#> 
#> $misclass
#> [1] 4
#> 
#> $misclass_points
#>            x           y label
#> 1  0.5952808  0.82122120     a
#> 2  0.5117812 -0.05612874     a
#> 3  0.1249309  0.41794156     a
#> 4 -0.3770596  0.88110773     b
#> 
#> $predicted
#>  [1] a a a b a a a a a a b a a a b b b b b a b b b b b b b b b b
#> Levels: a b
#> 
```
