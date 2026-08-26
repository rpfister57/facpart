# facpart

## Overview

facpart is a package with functions build on Facet Theory. It helps to partition 2-dimensional point configurations (such as multidimensional-scaling output) according to typical facet theoretical partition patterns.

Three types of common partion schemes are supported:

- **Axial** — parallel separating lines: `axialLine()`, `axialLines()`
- **Radial** — nested circles or ellipses: `radialCircle()`, `radialCircles()`, `radialEllipse()`, `radialEllipses()`
- **Angular** — wedge-shaped sectors: `angularPartition()`

Plus some utilities: `ellipseInConfig()`, `inoutEllipse()`, `mu2()`

This package complements similar software for facet partitioning: The R smacof package (Mair, Groenen, & De Leeuw, 2022) has partition functions based on svm; the stand-alone programs HUDAP for Windows (Amar, 2001; Amar & Toledano, 2001), and FSSA for Windows (https://raz-zeevy.github.io/fssa/#) (https://doritalt80.wixsite.com/ftas/copy-of-join-us) provide the same functionality, though outside the R ecosystem.


## References

Amar, R. (2001). Mathematical formulation of regionality in SSA and POSAC/MPOSAC. In: Elizur, D. (ed.). Facet theory: Integrating theory construction with data analysis. Prague: MATFYZPRESS (pp. 63-74).

Amar, R., & Toledano, S. (2001) (2nd ed.). HUDAP Manual. Jerusalem: The Hebrew University of Jerusalem, Computation Authority.

Guttman, R., & Greenbaum, C. W. (1998). Facet theory: Its development 
   and current status. European Psychologist, 3(1), 13-36.

Mair, P., Groenen, P., & De Leeuw, J. (2022). More on multidimensional scaling and unfolding in R: smacof version 2. *Journal of Statistical Software, 102*(10), 1-47. doi:10.18637/jss.v102.i10

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
res2 <- axialLines(crd, grp, fill = TRUE)
res2$misclass
```

## Dependencies

- `cluster` — minimum-area enclosing ellipse via `ellipsoidhull()`
- `MASS` — `lda()` for axial partitions
- `plotrix` — `draw.circle()`, `draw.ellipse()`

All other functionality uses base R. The package works well with the `smacof` package for multidimensional scaling.

## License

MIT

see also the article at: rpfister57.github.io/facpart
