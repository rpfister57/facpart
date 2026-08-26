# Guttman 1991 Intellligence Data

Data are based on the smacof package and on the original article by
Guttman and Levy (1991). The 12 test items are actually 12 subtests from
the WISC-R, and the sample are 2200 U.S. children aged 6.5 to 16.5 years
(see Guttman & Levy, 1991).

## Usage

``` r
gutt91
```

## Format

A list with three components: gutt91_cor, gutt91_dis, gutt91_df:

- gutt91_cor:

  A correlation matrix of 12 intelligence items

- gutt91_dis:

  A distance matrix derived from the correlation matrix

- gutt91_df:

  A data frame with 12 rows and 9 variables:

- items:

  The labels for the tests

- Material:

  The assignment of the Material facet (oral, manual, paper)

- Modality:

  The assignment of the Modality facet (verbal, numerical, figural)

- Rule:

  The assignment of the Rule Task facet (infer, applFact, applSTM)

- X1, X2:

  Coordinates from a 2-dimensional ordinal MDS

- D1, D2, D3:

  Coordinates from a 3-dimensional ordinal MDS

## Source

smacof package, data(Guttman1991)

## References

Guttman, L. & Levy, S. (1991). Two structural laws for intelligence
tests. Intelligence, 15, 79-103.

## Examples

``` r
data(gutt91)
Kor_D <- gutt91$gutt91_dis
Facets <- gutt91$gutt91_df
# MDS
gutt91_mds <- smacof::mds(Kor_D, type = "ordinal")
plot(gutt91_mds, main = "Guttman Levy 1991 Intelligence")

# Angular partition of Modality
angularPartition(crd = gutt91_mds$conf, group = Facets$Modality, add = FALSE)

#> $partition
#> [1] "angular"
#> 
#> $cuts
#> Comprehension        Coding    Arithmetic 
#>     -1.645176      2.189759     -3.076990 
#> 
#> $margin
#> [1] 0.1103395
#> 
#> $sector
#>  [1] 2 2 1 2 2 1 3 3 3 3 3 3
#> 
#> $majority
#> [1] "numerical" "verbal"    "figural"  
#> 
#> $center
#> [1] 0 0
#> 
#> $pt_angles
#>        Information       Similarities         Arithmetic         Vocabulary 
#>         -2.6020094         -2.0511100          2.7780239         -2.6488194 
#>      Comprehension          DigitSpan  PictureCompletion PictureArrangement 
#>         -2.0000848          2.6312076         -1.2902663         -0.9209234 
#>        BlockDesign     ObjectAssembly             Coding              Mazes 
#>          0.2315753         -0.3157876          1.7483111          0.3506980 
#> 
#> $pcoords
#>                            D1          D2
#> Information        -0.2536557 -0.15190509
#> Similarities       -0.1137925 -0.21840762
#> Arithmetic         -0.4410249  0.16780236
#> Vocabulary         -0.2341282 -0.12571647
#> Comprehension      -0.1887449 -0.41232278
#> DigitSpan          -0.8928474  0.49987343
#> PictureCompletion   0.1978224 -0.68657756
#> PictureArrangement  0.3642192 -0.47923328
#> BlockDesign         0.2707838  0.06385235
#> ObjectAssembly      0.5647150 -0.18450416
#> Coding             -0.2123476  1.18363379
#> Mazes               0.9390008  0.34350503
#> 
#> $pgroup
#>  [1] verbal    verbal    numerical verbal    verbal    numerical figural  
#>  [8] figural   figural   figural   figural   figural  
#> Levels: verbal numerical figural
#> 
#> $misclass
#> [1] 0
#> 
#> $misclass_points
#> [1] x     y     label
#> <0 rows> (or 0-length row.names)
#> 
```
