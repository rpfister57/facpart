# K-group radial-elliptic partition (nested ellipses)

Generalises
[`radialEllipse()`](https://rpfister57.github.io/facpart/reference/radialEllipse.md)
to `k >= 2` groups using `k-1` nested ellipses. The nesting constraint
(ellipse `s-1` lies inside ellipse `s`) is enforced during optimisation.
Groups are ordered by factor level (level 1 = innermost).

## Usage

``` r
radialEllipses(
  crd,
  group,
  ellipse = NULL,
  fill = FALSE,
  output = TRUE,
  col = "purple",
  cols = NULL,
  lwd = 2,
  lty = 1,
  add = TRUE
)
```

## Arguments

- crd:

  Numeric matrix or data frame with exactly 2 columns.

- group:

  Factor with `k >= 2` levels.

- ellipse:

  Optional length-5 numeric vector `(cx, cy, a, b, angle)` specifying
  the innermost ellipse exactly. When supplied, all outer ellipses are
  uniform scalings of this one; only the scale factors are searched.

- fill:

  If `TRUE`, shade each elliptic sector.

- output:

  If `TRUE` (default), return results list.

- col:

  Ellipse border colour (default `"purple"`).

- cols:

  Length-`k` colours; auto-generated if `NULL`.

- lwd:

  Line width (default `2`).

- lty:

  Line type (default `1`).

- add:

  If `TRUE` (default), add to existing plot; if `FALSE`, call
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) first.

## Value

If `output = TRUE`, a list with vectors `cx`, `cy`, `a`, `b`, `angle`
(radians), `misclass` (list with `n` and `indices`), `misclass_points`
(data frame with columns `x`, `y`, `label` for each misclassified
point), `sector`, and `majority`.

## Details

When `ellipse` is supplied as a length-5 vector `(cx, cy, a, b, angle)`,
it defines the innermost ellipse exactly. All `k-1` ellipses then share
that center, orientation, and a:b ratio — outer ellipses are uniform
scalings of the inner one. The scale factors `t_2, ..., t_{k-1}` (with
`t_s > t_{s-1} >= 1`) are found by exact 1D scan over each point's
critical scale `sqrt((u_i/a)^2 + (v_i/b)^2)` in the ellipse's rotated
frame.

When `ellipse` is `NULL` (default), each ellipse has its own
independently optimised center, semi-axes, and rotation, fitted by
Nelder-Mead. Multi-start uses the covariance ellipse of the inner groups
plus an inflated copy of the previous ellipse (a guaranteed- feasible
starting point that avoids stalling on the infeasibility penalty
plateau).

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
g1 <- cbind(rnorm(15, 0, 0.2), rnorm(15, 0, 0.15))
th2 <- runif(15, 0, 2 * pi)
g2 <- cbind(1.4 * cos(th2), 1.0 * sin(th2)) + matrix(rnorm(30, 0, 0.1), 15)
th3 <- runif(15, 0, 2 * pi)
g3 <- cbind(2.8 * cos(th3), 2.0 * sin(th3)) + matrix(rnorm(30, 0, 0.1), 15)
crd <- rbind(g1, g2, g3)
grp <- factor(c(rep("a", 15), rep("b", 15), rep("c", 15)))

# Independent ellipses
radialEllipses(crd, grp, fill = TRUE, add = FALSE)

# Fixed-shape mode: supply the innermost ellipse as a 5-vector
radialEllipses(crd, grp, ellipse = c(0, 0, 0.3, 0.2, 0), fill = TRUE, add = FALSE)
} # }
```
