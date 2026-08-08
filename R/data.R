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
#' @source smacof package, data(Guttman1965)
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
#' Data are based on the smacof package and on the original article by 
#' Guttman and Levy (1991). The 12 test items are actually 12 subtests 
#' from the WISC-R, and the sample are 2200 U.S. children 
#' aged 6.5 to 16.5 years (see Guttman & Levy, 1991).
#' 
#' @format
#' A list with three components: gutt91_cor, gutt91_dis, gutt91_df:
#' \describe{
#'   \item{gutt91_cor}{A correlation matrix of 12 intelligence items}
#'   \item{gutt91_dis}{A distance matrix derived from the correlation matrix}
#'   \item{gutt91_df}{A data frame with 12 rows and 9 variables:}
#'   \item{items}{The labels for the tests}
#'   \item{Material}{The assignment of the Material facet (oral, manual, paper)}
#'   \item{Modality}{The assignment of the Modality facet (verbal, numerical, figural)}
#'   \item{Rule}{The assignment of the Rule Task facet (infer, applFact, applSTM)}
#'   \item{X1, X2}{Coordinates from a 2-dimensional ordinal MDS}
#'   \item{D1, D2, D3}{Coordinates from a 3-dimensional ordinal MDS}
#' }
#' @source smacof package, data(Guttman1991)
#' @references
#' Guttman, L. & Levy, S. (1991). Two structural laws for intelligence tests. Intelligence, 15, 79-103.
#' @examples
#' data(gutt91)
#' Kor_D <- gutt91$gutt91_dis
#' Facets <- gutt91$gutt91_df
#' # MDS
#' gutt91_mds <- smacof::mds(Kor_D, type = "ordinal")
#' plot(gutt91_mds, main = "Guttman Levy 1991 Intelligence")
#' # Angular partition of Modality
#' angularPartition(crd = gutt91_mds$conf, group = Facets$Modality, add = FALSE)
#' 
"gutt91"

