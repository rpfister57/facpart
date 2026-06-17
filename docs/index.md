# facpart

## Overview

facpart is a package with functions build on Facet Theory. It helps to
partition 2-dimensional point configurations (such as
multidimensional-scaling output) according to typical facet theoretical
partition patterns.

The three types of partion schemes:

- **Axial** — parallel separating lines:
  [`axialLine()`](https://rpfister57.github.io/facpart/reference/axialLine.md),
  [`axialLines()`](https://rpfister57.github.io/facpart/reference/axialLines.md)
- **Radial** — nested circles or ellipses:
  [`radialCircle()`](https://rpfister57.github.io/facpart/reference/radialCircle.md),
  [`radialCircles()`](https://rpfister57.github.io/facpart/reference/radialCircles.md),
  [`radialEllipse()`](https://rpfister57.github.io/facpart/reference/radialEllipse.md),
  [`radialEllipses()`](https://rpfister57.github.io/facpart/reference/radialEllipses.md)
- **Angular** — wedge-shaped sectors:
  [`angularPartition()`](https://rpfister57.github.io/facpart/reference/angularPartition.md)

Plus some utilities:
[`ellipseInConfig()`](https://rpfister57.github.io/facpart/reference/ellipseInConfig.md),
[`inoutEllipse()`](https://rpfister57.github.io/facpart/reference/inoutEllipse.md),
`mu2`

## References

Guttman, R., & Greenbaum, C. W. (1998). Facet theory: Its development
and current status. European Psychologist, 3(1), 13-36.

Shye, S. (2014). Faceted Smallest Space Analysis (FSSA). In A. Michalos
(Ed.), Encyclopedia of quality of life research (pp. 2129-2133). New
York: Springer.

Shye, S. (2015). New directions in facet theory. In S. Shye, E. Solomon,
& I. Borg (Eds.), 15th International Facet Theory Conference
(pp. 147-158). New York City: Fordham University.

## Installation

The current version can be installed from github. Note: `facpart` is in
a very early stage!

``` r

# install.packages("remotes")
remotes::install_github(repo = "https://github.com/rpfister57/facpart.git")
```

## Quick example

``` r

library(facpart)

set.seed(1)
# Three clusters arranged around a common center
theta <- rep(c(0, 2 * pi / 3, 4 * pi / 3), each = 12) + rnorm(36, 0, 0.25)
r     <- runif(36, 0.5, 1.5)
crd   <- cbind(r * cos(theta), r * sin(theta))
grp   <- factor(rep(c("a", "b", "c"), each = 12))

plot(crd, asp = 1)
res <- angularPartition(crd, grp)
res$misclass   # number of misclassified points
res$center     # optimised wedge apex
```

## Dependencies

- `cluster` — minimum-area enclosing ellipse via `ellipsoidhull()`
- `MASS` — `lda()` for axial partitions
- `plotrix` —
  [`draw.circle()`](https://plotrix.github.io/plotrix/reference/draw.circle.html),
  [`draw.ellipse()`](https://plotrix.github.io/plotrix/reference/draw.ellipse.html)

All other functionality uses base R. The package works well with the
`smacof` package.

## License

MIT
