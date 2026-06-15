# ===================================================================
# ==== Shared internal helpers and plotting utilities ====
# ===================================================================


# ---- Internal helpers ----

#' @noRd
# convert degrees to radians, and vice versa
.d2r <- function(degree) {
    return(rad <- degree * (pi/180))
}

#' @noRd
.r2d <- function(rad) {
    return(degree <- rad * (180/pi))
}


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


#' Highlight misclassified points on the active plot
#'
#' Overlays a marker on each misclassified point, using the `misclass_points`
#' data frame returned by [axialLine()], [axialLines()], [radialCircle()],
#' [radialCircles()], [radialEllipse()], [radialEllipses()], or
#' [angularPartition()]. The full result list may also be passed directly;
#' `misclass_points` is extracted automatically.
#'
#' @param x A data frame with columns `x`, `y`, `label` (the `misclass_points`
#'   element of a partition result), or the full partition result list.
#' @param col Marker colour (default `"red"`).
#' @param pch Point character for the marker (default `4`, an X).
#' @param cex Point size expansion factor (default `1.5`).
#' @param lwd Line width of the marker symbol (default `2`).
#'
#' @return `invisible(NULL)`, called for its side effect of drawing on the
#'   active plot.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' crd <- rbind(cbind(rnorm(15, -1), rnorm(15)),
#'              cbind(rnorm(15,  1), rnorm(15)))
#' grp <- factor(c(rep("a", 15), rep("b", 15)))
#' res <- axialLine(crd, grp, add = FALSE)
#' highlightMisclass(res)
#' }
#'
#' @export
highlightMisclass <- function(x,
                               col = "red",
                               pch = 4,
                               cex = 1.5,
                               lwd = 2) {
    if (is.list(x) && !is.data.frame(x)) x <- x$misclass_points
    if (is.null(x) || nrow(x) == 0L) return(invisible(NULL))
    graphics::points(x$x, x$y, col = col, pch = pch, cex = cex, lwd = lwd)
    invisible(NULL)
}
