#' facpart: Facet Theory Partitioning of 2D Configurations
#'
#' Tools for partitioning 2-dimensional configurations (such as
#' multidimensional-scaling output) into facets, supporting axial (parallel
#' lines), radial (nested circles or ellipses), and angular (wedge-shaped)
#' partitions. Includes utilities for fitting bounding, least-squares, and
#' minimum-area ellipses.
#'
#' @section Partition functions:
#' - [axialLine()], [axialLines()] — parallel-line partitions
#' - [radialCircle()], [radialCircles()] — nested-circle partitions
#' - [radialEllipse()], [radialEllipses()] — nested-ellipse partitions
#' - [angularPartition()], [angularPartition3()] — wedge partitions
#'
#' @section Ellipse utilities:
#' - [ellipseInConfig()] — fit and draw a bounding/best-fit/minimum ellipse
#' - [inoutEllipse()] — classify points as inside/outside a given ellipse
#'
#' @keywords internal
#' @importFrom MASS lda
#' @importFrom cluster ellipsoidhull
#' @importFrom grDevices adjustcolor hcl.colors
#' @importFrom graphics abline par polygon polypath segments
#' @importFrom plotrix draw.circle draw.ellipse
#' @importFrom stats cov optim
#' @importFrom utils combn
"_PACKAGE"
