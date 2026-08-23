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


#' BIS Intelligence Data
#' 
#' Data are based on analyses in: 
#' Pfister, H.-R., & Beauducel, A. (1993). Data contain the correlation
#' matrix among 12 intelligence tests, and the corresponding 2-facets
#' assignment.
#' 
#' @format
#' A list with two components: Bis1_cor is a correlation matrix, and
#' Bis1_facets are the facet assignments. The correlations are among
#' 12 'cells', each cell is an aggregate variable of four items according
#' to the bimodal Berlin Model of Intelligence Structure (BIS)
#' (Pfister & Jäger, 1992; Süß, 2015).
#' \describe{
#'   \item{Bis1_cor}{A correlation matrix of 12 intelligence tests}
#'   \item{Bis1_facets}{A data frame with 12 rows and 4 variables:}
#'   \item{Operation}{The elements of the 'Operation' facet (B: speed,
#'   M: memory, E: creativity, K: complex cognition)}
#'   \item{Inhalt}{The elements of the 'Content' facet (F: figural,
#'   V: verbal, N: numeric)}
#'   \item{BiMod}{The bimodal combination of Operation and Content facet}
#'   \item{Speed}{A reassignment to speed (B) and power (M, E, K)}
#' }
#' @references
#' Pfister, H.-R., & Beauducel, A. (1993). Stability of operation
#' and content facets: A facet analysis of the Berlin model
#' of intelligence structure BIS. Paper presented at the Fourth
#' International Facet Theory Conference, Prague, 1993.
#' 
#' Pfister, H.-R., & Jäger, A. O. (1992). Topografische Analysen zum
#' Berliner Intelligenzstrukturmodell BIS. Diagnostica, 38(2), 91-115.
#' 
#' Süß, H.-M. (2015). The construct validity of the Berlin Intelligence
#' Structure Model (BIS). In: Roazzi, A., de Souza, B.C., Bilsky, W. (ed.):
#' Facet theory. Recife: Editora UFPE (p. 123-138).
#' @examples
#' data(BIS1)
#' Bis1Cells <- BIS1$Bis1_cor
#' Bis1Facets <- BIS1$Bis1_facets
#' Bis1Cells_D <- smacof::sim2diss(s = Bis1Cells, method = "corr")
#' Bis1Cells_mds <- smacof::mds(delta = Bis1Cells_D, type = "ordinal")
#' Bis1Cells_mds
#' plot(Bis1Cells_mds, main = "BIS 1, Cells", las = 1, asp = 1)
#' Bis1Cells_mds_ax <- axialLines(crd = Bis1Cells_mds$conf,
#'    group = Bis1Facets$Operation, 
#'    col = "steelblue", fill = TRUE)
#' Bis1Cells_mds_ax$misclass
#' 
"BIS1"

#' Life Satisfaction Data
#' 
#' Data are from Levy (1976; see Borg & Groenen, 2005). Life is a list
#' containing the correlations among fifteen items asking about life
#' satisfaction, and the corresponding assignments to two facets.
#' 
#' @format
#' A list with two components: life is a correlation matrix of
#' fifteen items from a survy asking about different areas of life
#' satisfaction; life_facets is a data frame with assignments to
#' two facets: Status (2 elements), and Area (8 elements). For details
#' see Levy (1976) and Borg and Groenen (2005).
#' \describe{
#'   \item{life}{A correlation matrix of 15 items about aspects of life satisfaction}
#'   \item{life_facets}{A data frame with 15 rows and 3 variables:}
#'   \item{Status}{The Status facet with 2 elements}
#'   \item{Area}{The Area facet with 8 elements}
#' }
#' @references
#' Borg, I., & Groenen, P. J. F. (2005) (2nd ed.). Modern multidimensional 
#' scaling. Theory and applications. New York: Springer.
#' 
#' Levy, S. (1976). Use of the mapping sentence for coordinating theory
#' and research: A cross-cultural example. Quality and Quantity, 10, 117-125.
#' @examples
#' data(Life)
#' life <- Life$life
#' life_facets <- Life$life_facets
#' lifeD <- smacof::sim2diss(life, method = "corr")
#' life_mds <- smacof::mds(delta = lifeD, type = "ordinal")
#' life_mds
#' plot(x = life_mds$conf, pch = 19,
#'   xlim = c(-1, 1.5), ylim = c(-1, 1), asp = 1)
#' text(x = life_mds$conf, labels = life_facets$Area, 
#'   cex = 0.7, pos = 4)
#'   angularPartition(crd = life_mds$conf,
#'                    group = life_facets$Area)
#'
"Life"
