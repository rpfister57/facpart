#' Guttman 1965 Intelligence Data
#' 
#' Data contain coordinates for 21 items from a MDS analysis of
#' the original correlation matrix. For each item, the facet
#' assignment is provided. Data are based on the smacof package.
#' 
#' @format
#' A data frame with 21 observations and 3 variables:
#' \describe{
#'   \item{gfacets}{A factor with 4 levels: analytical, complex, achievement1, achievement2}
#'   \item{D1, D2}{Coordinates from a 2-dimensional MDS}
#' }
#' 
#' @source smacof package, data(Guttman1965)
#' 
#' @examples
#' plot(guttman65mds[ , 2:3], type = "n", asp = 1)
#' text(guttman65mds[ , 2:3], labels = guttman65mds$gfacets, cex = 0.6)
#' 
#' # check for radial partitions
#' ellipses_out <- radialEllipses(crd = guttman65mds[ , 2:3], 
#'    group = guttman65mds$gfacets, fill = TRUE, add = TRUE)
"guttman65mds"


#' Guttman 1991 Intellligence Data
#' 
#' Data are based on the smacof package. It contains the correlations
#' among 12 intelligence items, and the assignment to two facets:
#' Material and Modality.
#' 
#' @format
#' A list with two elements: gutt91_cor, gutt91_var
#' \describe{
#'   \item{gutt91_cor}{A correlation matrix of 12 intelligence items}
#'   \item{gutt91_var}{A data frame with 12 observations and 3 variables:
#'   items = short name of the item, Material = Facet 1 (Oral, Manual, Paper),
#'   Modality = Facet 2 (verbal, numerical, figural)}
#' }
#' 
#' @source smacof package, data(Guttman1991)
#' 
#' @references
#' Guttman, L. & Levy, S. (1991). Two structural laws for intelligence tests. Intelligence, 15, 79-103.
#' 
#' @examples
#' data(gutt91)
#' Kor <- gutt91$gutt91_cor
#' Facets <- gutt91$gutt91_var
#' # convert correlations to distances
#' Kor_D <- smacof::sim2diss(Kor, method = "corr", to.dist = TRUE)
#' # MDS
#' gutt91_mds <- smacof::mds(Kor_D, type = "ordinal")
#' plot(gutt91_mds, main = "Guttman 1991 Intelligence")
#' # Angular partition of Modality
#' angularPartition(crd = gutt91_mds$conf, group = Facets$Modality)
#' 
"gutt91"

