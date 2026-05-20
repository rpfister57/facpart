# facpart

Facet Theory partitioning of 2-dimensional point configurations.

Tools for partitioning 2D configurations (such as multidimensional-scaling
output) into facets, with three types of partion schemes:

- **Axial** — parallel separating lines: `axialLine()`, `axialLines()`
- **Radial** — nested circles or ellipses: `radialCircle()`, `radialCircles()`,
  `radialEllipse()`, `radialEllipses()`
- **Angular** — wedge-shaped sectors: `angularPartition()`, `angularPartition3()`

Plus ellipse-fitting utilities: `ellipseInConfig()`, `inoutEllipse()`.

## Installation

```r
# install.packages("devtools")
devtools::install_local("path/to/facpart")
```

or, from inside the package directory:

```r
devtools::install()
```

To regenerate `man/` documentation and `NAMESPACE` from the inline roxygen
comments:

```r
devtools::document()
```

## Quick example

```r
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
- `plotrix` — `draw.circle()`, `draw.ellipse()`

All other functionality uses base R.

## License

MIT
