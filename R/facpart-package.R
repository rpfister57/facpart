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
#' - [angularPartition()] — wedge partitions
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
#' @importFrom stats as.dist cov optim predict
#' @importFrom utils combn
"_PACKAGE"


# ---- Shared internal helpers ----

#' @noRd
# Generate all n! permutations of 1..n as rows of an integer matrix.
.permutations <- function(n) {
    if (n == 1L) return(matrix(1L, nrow = 1L, ncol = 1L))
    prev <- .permutations(n - 1L)
    m    <- nrow(prev)
    out  <- matrix(0L, nrow = n * m, ncol = n)
    for (ins in seq_len(n)) {
        rows <- seq_len(m) + (ins - 1L) * m
        if (ins > 1L) out[rows, seq_len(ins - 1L)]   <- prev[, seq_len(ins - 1L), drop = FALSE]
        out[rows, ins] <- n
        if (ins < n)  out[rows, seq(ins + 1L, n)]    <- prev[, seq(ins, n - 1L),  drop = FALSE]
    }
    out
}


#' @noRd
# Bijective assignment of k groups to k regions that maximises correctly
# classified points. count_mat[g, r] = number of group-g points in region r.
# Returns integer vector assignment where assignment[r] = group index for
# region r. Brute-force over all k! permutations -- fast for small k.
.assign_groups <- function(count_mat) {
    k <- nrow(count_mat)
    if (k == 1L) return(1L)
    perms     <- .permutations(k)
    best_val  <- -1L
    best_perm <- seq_len(k)
    idx_r     <- seq_len(k)
    for (i in seq_len(nrow(perms))) {
        perm <- perms[i, ]
        val  <- sum(count_mat[cbind(perm, idx_r)])
        if (val > best_val) {
            best_val  <- val
            best_perm <- perm
        }
    }
    best_perm
}
