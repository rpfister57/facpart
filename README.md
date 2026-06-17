# facpart

## Overview

facpart is a package with functions build on Facet Theory. It helps to partition 2-dimensional point configurations (such as multidimensional-scaling
output) according to typical facet theoretical partition patterns.

The three types of partion schemes:

- **Axial** — parallel separating lines: `axialLine()`, `axialLines()`
- **Radial** — nested circles or ellipses: `radialCircle()`, `radialCircles()`, `radialEllipse()`, `radialEllipses()`
- **Angular** — wedge-shaped sectors: `angularPartition()`

Plus some utilities: `ellipseInConfig()`, `inoutEllipse()`, `mu2`


## References

Guttman, R., & Greenbaum, C. W. (1998). Facet theory: Its development 
   and current status. European Psychologist, 3(1), 13-36.

Shye, S. (2014). Faceted Smallest Space Analysis (FSSA). In A. Michalos (Ed.), 
   Encyclopedia of quality of life research (pp. 2129-2133). 
   New York: Springer.

Shye, S. (2015). New directions in facet theory. In S. Shye, 
   E. Solomon, & I. Borg (Eds.), 15th International Facet Theory 
   Conference (pp. 147-158). New York City: Fordham University.


## Installation
The current version can be installed from github. Note: `facpart` is in a very early stage!
```r
# install.packages("remotes")
remotes::install_github(repo = "https://github.com/rpfister57/facpart.git")
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

All other functionality uses base R. The package works well with the `smacof` package.

## License

MIT
