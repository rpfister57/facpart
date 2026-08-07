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


# Cache for .perms_k(). A package-level environment, mutated in place -- the
# binding itself is never reassigned, so the sealed namespace is not an issue.
.perm_cache <- new.env(parent = emptyenv())


#' @noRd
# All k! permutations of 1..k as rows, memoised by k.
#
# gtools::permutations() is recursive and rebuilds the whole matrix on every
# call. .assign_groups() is called once per candidate radius inside the
# radius-refinement and ordering-search loops, and .bij_best() once per
# candidate angle in axialLines(), so that rebuild dominated the runtime
# (profiled at ~79% of radialCircles() at k = 5). The matrix depends only on
# k, so cache it.
.perms_k <- function(k) {
    key <- as.character(k)
    out <- .perm_cache[[key]]
    if (is.null(out)) {
        out <- gtools::permutations(k, k)
        .perm_cache[[key]] <- out
    }
    return(out)
}


#' @noRd
# Bijective assignment of k groups to r regions (k = r) that maximises
# correctly classified points.
# count_mat[g, r] = number of group-g points in region r.
# Returns integer vector assignment where assignment[r] = group index for
# region r. Brute-force over all k! permutations -- fast for small k.
#
# `lev` (the group level names) only affects *ties*. Several assignments can
# be equally good, and picking the first one in permutation order picks by
# group index -- which depends on the caller's factor level order, so
# `majority` would change when the groups are merely renamed or reordered.
# Given `lev`, the tie goes to the assignment whose region -> name tuple is
# lexicographically smallest, which is intrinsic to the data. Pass it at the
# reporting call sites; omit it in inner scoring loops such as
# .partition_err(), where only the (tie-independent) best value is used and
# the extra key building would be wasted work.
.assign_groups <- function(count_mat, lev = NULL) {
    k <- nrow(count_mat)
    if (k == 1L) return(1L)
    perms <- .perms_k(k)
    n_p   <- nrow(perms)

    # Correct count of every permutation at once. as.vector(perms) runs
    # column-major, so pairing it with each region index repeated n_p times
    # gives a matrix whose [i, r] entry is count_mat[perms[i, r], r]; the row
    # sums are the per-permutation totals. This replaces an R-level loop over
    # all k! permutations and matters because .partition_err() calls this
    # inside the radius-refinement and ordering-search loops.
    hits <- matrix(count_mat[cbind(as.vector(perms),
                                   rep(seq_len(k), each = n_p))],
                   nrow = n_p, ncol = k)
    vals <- .rowSums(hits, n_p, k)

    best_rows <- which(vals == max(vals))
    if (length(best_rows) == 1L || is.null(lev))
        return(perms[best_rows[1L], ])

    keys <- vapply(best_rows,
                   function(i) paste(lev[perms[i, ]], collapse = "\r"),
                   character(1))

    return(perms[best_rows[order(keys)[1L]], ])
}


#' @noRd
# Exact bijection-constrained best-correct classification over M candidate
# segmentations of a sorted point sequence.
#
# cum[i + 1, g] = number of group-g points among sorted points 1..i, so the
# count of group g in a segment spanning sorted positions a+1..b is
# cum[b + 1, g] - cum[a + 1, g]. lo_idx[[s]] / hi_idx[[s]] supply those two
# row indices for segment s, either as a length-M vector or as a scalar when
# the bound is the same for every candidate (which keeps the big vectors down
# to one per segment side).
#
# The score of a candidate is the best total correct over all k! segment-to-
# group assignments with each group used exactly once -- the criterion
# .assign_groups() applies downstream, so a search built on this minimises
# exactly the misclassification its caller goes on to report. Scoring each
# segment by its own local majority instead lets two segments claim the same
# group, and so minimises a total that no bijection can achieve; that was a
# real bug in axialLines(), radialCircle() and angularPartition(). Do not
# reintroduce it.
#
# Scoring is a two-stage bound-and-prune, because the k! scan is the
# expensive part:
#   Stage 1 (all M, cheap) -- `ub`, the *unconstrained* best-correct, where
#     every segment independently keeps its own local-majority group. Any
#     bijection is one particular per-segment choice, so exact <= ub for
#     every candidate: a valid upper bound, at k * (2k - 1) vector ops
#     against k! * (k + 1) for the exact scan.
#   Stage 2 (survivors only) -- walk a threshold down from max(ub), scanning
#     all k! assignments for the candidates that could still win. Once some
#     candidate's exact score reaches the threshold, every unscanned
#     candidate has ub < best and is provably excluded. The bound is tight in
#     practice, which is what makes this much cheaper than scoring all M
#     exactly.
#
# Returns the best achievable correct count and every candidate attaining it
# (`at`), so the caller can apply its own margin tie-break to that pool alone
# rather than to all M.
.bij_best <- function(cum, lo_idx, hi_idx, k, M) {

    # Count of group `g` in segment `s`, over the candidates in `sel`
    # (an index vector into 1..M, or NULL for all M).
    seg_cnt <- function(s, g, sel) {
        ia <- lo_idx[[s]]; ib <- hi_idx[[s]]
        if (!is.null(sel)) {
            if (length(ia) > 1L) ia <- ia[sel]
            if (length(ib) > 1L) ib <- ib[sel]
        }
        cum[ib, g] - cum[ia, g]
    }

    # ---- Stage 1: cheap unconstrained upper bound for all M candidates ----
    ub <- 0L
    for (s in seq_len(k)) {
        r_max <- seg_cnt(s, 1L, NULL)
        for (g in 2L:k) r_max <- pmax(r_max, seg_cnt(s, g, NULL))
        ub <- ub + r_max
    }

    # ---- Stage 2: exact bijection-constrained score, on survivors only ----
    perms <- .perms_k(k)

    exact_correct <- function(sel) {
        cn <- vector("list", k)
        for (s in seq_len(k)) {
            cs <- vector("list", k)
            for (g in seq_len(k)) cs[[g]] <- seg_cnt(s, g, sel)
            cn[[s]] <- cs
        }
        best <- rep(0L, if (is.null(sel)) M else length(sel))
        for (p_idx in seq_len(nrow(perms))) {
            perm    <- perms[p_idx, ]
            correct <- cn[[1L]][[perm[1L]]]
            for (s in 2L:k) correct <- correct + cn[[s]][[perm[s]]]
            best <- pmax(best, correct)
        }
        best
    }

    thresh      <- max(ub)
    max_correct <- -1L
    scanned     <- logical(M)
    repeat {
        sel <- which(ub >= thresh & !scanned)
        if (length(sel) > 0L) {
            max_correct  <- max(max_correct, max(exact_correct(sel)))
            scanned[sel] <- TRUE
        }
        if (max_correct >= thresh) break
        thresh <- thresh - 1L
    }

    # Any candidate attaining max_correct has ub >= max_correct, so this pool
    # is a complete superset of the co-minimal set.
    pool <- which(ub >= max_correct)

    return(list(max_correct = max_correct,
                at          = pool[exact_correct(pool) == max_correct]))
}


#' @noRd
# Candidate nesting orders, each a permutation of the k group indices with the
# innermost group first.
#
# radialCircles() and radialEllipses() fit their boundaries sequentially, and a
# sequential fit has to know which groups belong inside boundary s. That
# ordering must NOT be taken from the factor level order: R's default levels
# are alphabetical, so the result would depend on how the caller happened to
# name the groups, and .assign_groups() cannot repair a nesting the geometry
# search never looked for (the symptom is a `majority` entry owning no points).
# Instead the callers *search* these candidates and keep whichever one
# minimises the bijection-constrained misclassification they go on to report.
#
# Everything here is intrinsic to the configuration, so the returned list --
# and its iteration order -- is invariant under relabelling or reordering of
# the factor levels:
#   * candidate 1 is the data-derived order: groups by mean distance from the
#     configuration centroid, ties broken by level *name*. Callers take strict
#     improvements only, so this wins any tie on misclassification.
#   * the remaining permutations follow in canonical order, sorted by the
#     tuple of level names they induce rather than by level index.
# When k! exceeds `max_perm` only the data-derived order is returned and
# attr(, "exhaustive") is FALSE, so the caller can say so.
.nesting_orders <- function(pcoords, group, max_perm = 120L) {
    lev <- levels(group)
    k   <- length(lev)

    ctr   <- c(mean(pcoords[, 1]), mean(pcoords[, 2]))
    dist_ <- sqrt((pcoords[, 1] - ctr[1])^2 + (pcoords[, 2] - ctr[2])^2)
    mdist <- tapply(dist_, group, mean)
    # an empty level tells us nothing about its radius; park it outermost
    mdist[is.na(mdist)] <- Inf
    # order on (mean distance, level name) so ties never fall back to the
    # level index, which is exactly what must not matter here
    derived <- as.integer(order(mdist, lev))

    if (k < 2L || factorial(k) > max_perm) {
        out <- list(derived)
        attr(out, "exhaustive") <- (k < 2L)
        return(out)
    }

    perms <- .perms_k(k)
    keys  <- apply(perms, 1L, function(p) paste(lev[p], collapse = "\r"))
    perms <- perms[order(keys), , drop = FALSE]
    perms <- perms[!apply(perms, 1L,
                          function(p) identical(as.integer(p), derived)), ,
                   drop = FALSE]

    out <- c(list(derived),
             lapply(seq_len(nrow(perms)), function(i) as.integer(perms[i, ])))
    attr(out, "exhaustive") <- TRUE

    return(out)
}


#' @noRd
# Draw a line from x,y with a given angle (in radians).
# Default is from center (0,0)
.draw_angle_line <- function(x = 0, y = 0, angle, 
                            length = 1, ppoint = FALSE, ...) {
    # x, y: start coordinates
    # angle: in radians
    # Compute endpoint
    x_end <- x + length * cos(angle)
    y_end <- y + length * sin(angle)
    
    # Draw the line
    segments(x, y, x_end, y_end, ...)
    
    if (ppoint) points(x_end, y_end, pch = 19)
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
