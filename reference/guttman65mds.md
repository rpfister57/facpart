# Guttman 1965 Intelligence Data

Data contain coordinates for 21 items from a MDS analysis of the
original correlation matrix. For each item, the facet assignment is
provided. Data are based on the smacof package.

## Usage

``` r
guttman65mds
```

## Format

A data frame with 21 observations and 3 variables:

- gfacets:

  A factor with 4 levels: analytical, complex, achievement1,
  achievement2

- D1, D2:

  Coordinates from a 2-dimensional MDS

## Source

smacof package, data(Guttman1965)

## Examples

``` r
plot(guttman65mds[ , 2:3], type = "n", asp = 1)
text(guttman65mds[ , 2:3], labels = guttman65mds$gfacets, cex = 0.6)

# check for radial partitions
ellipses_out <- radialEllipses(crd = guttman65mds[ , 2:3], 
   group = guttman65mds$gfacets, fill = TRUE, add = TRUE)
```
