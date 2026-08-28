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
# optim()'s Nelder-Mead sizes the initial simplex step for coordinate i as
# fabs(par[i]) when nonzero, and only falls back to a sane parscale-sized
# step when par[i] == 0 exactly (an exact-equality check, not a "close to
# zero" one). Centroids of mean-centered data (e.g. any MDS output) are
# ~1e-16, not exactly 0, so they silently hit the degenerate branch: optim()
# reports convergence after 3 evaluations without moving from the start.
# Snap near-zero start coordinates to exact 0 so the non-degenerate fallback
# step-size kicks in instead. Shared by angularPartition(), radialCircle(s)()
# and radialEllipse(s)(), all of which multi-start Nelder-Mead from a
# centroid-derived point.
.snap_zero <- function(s0, parscale, tol = sqrt(.Machine$double.eps)) {
    ifelse(abs(s0) < parscale * tol, 0, s0)
}


#' @noRd
# Coarse fallback seeds for a 2D center search: fnToOpt() is assumed
# piecewise-constant (an integer misclassification count), so a local,
# gradient-free optim() run can stall on a flat plateau around any start
# without ever finding a genuinely better region. Confirmed on real MDS data
# two ways: a plateau extending well beyond optim()'s own step size but
# still *inside* the data's bounding box, and separately a configuration
# where the true optimum sits *outside* the bounding box entirely. One grid
# widened to reach outside points would be too coarse to land in a narrow
# region still inside the box, so this scans the plain bounding box and,
# separately, one padded by `pad_frac` of the data range on each side, and
# returns the best cell of EACH grid (as a list of two points) for the
# caller to seed one more optim() run from apiece.
#
# Shared by angularPartition(), radialCircle()/radialCircles() and
# radialEllipse()/radialEllipses(); callers trigger this only when their
# heuristic starts didn't already reach the provable optimum (misclass = 0),
# since it costs 2 * n_grid^2 extra evaluations of the cheap inner search.
.grid_seeds <- function(fnToOpt, rng_x, rng_y, n_grid, pad_frac = 0.5) {
    pad_x <- diff(rng_x); if (pad_x == 0) pad_x <- 1
    pad_y <- diff(rng_y); if (pad_y == 0) pad_y <- 1
    pad_x <- pad_x * pad_frac
    pad_y <- pad_y * pad_frac

    grids <- list(
        list(gx = seq(rng_x[1],         rng_x[2],         length.out = n_grid),
             gy = seq(rng_y[1],         rng_y[2],         length.out = n_grid)),
        list(gx = seq(rng_x[1] - pad_x, rng_x[2] + pad_x, length.out = n_grid),
             gy = seq(rng_y[1] - pad_y, rng_y[2] + pad_y, length.out = n_grid))
    )

    lapply(grids, function(grd) {
        best_val <- Inf
        best_p   <- NULL
        for (gx0 in grd$gx) {
            for (gy0 in grd$gy) {
                v <- fnToOpt(c(gx0, gy0))
                if (v < best_val) {
                    best_val <- v
                    best_p   <- c(gx0, gy0)
                }
            }
        }
        best_p
    })
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


#' Project a point (x, y) onto a line (a, b)
#' 
#' On an existing plot including a line with a = intercept and b = slope,
#' a point with coordinates (x, y) is projected. Optionally, the line
#' as well as the projection arrow and the distance from intercept to
#' the projected point are drawn.
#' 
#' @param Pxy A vector c(x,y) with coordinates of a point.
#' @param Lab A vector c(a,b) with intercept (a) and slope (b) of a line.
#' @param add Logical: Should the projection be drawn (default = FALSE).
#' @param d0 Logical: Should the distance line from a to 
#'     projection be drawn (default = FALSE).
#' @param color Character - a color name
#' @param addLine Logical: Should the line be drawn.
#' 
#' @return A list with the coordinates of the projected point, and 
#'    the distance vector of the projected point.
#'    
#' @examples
#' \dontrun{
#' plot(-2:2, -2:2, type = "n", asp = 1)
#' abline(a = 0, b = 1)
#' aPoint <- matrix(c(-1.5, 1), nrow = 1, byrow = TRUE)
#' points(aPoint, pch = 19)
#' projectP2L(Pxy = aPoint, Lab = c(0, 1), add = TRUE)
#' }
#' 
#' @export
projectP2L <- function(Pxy, Lab, 
                       add = FALSE, 
                       d0 = FALSE,
                       color = "blue",
                       addLine = FALSE){
    x0 <- Pxy[1]
    y0 <- Pxy[2]
    a <- Lab[1]
    b <- Lab[2]
    
    xL <- (x0 + b * (y0 - a)) / (1 + b^2)
    yL <- a + b*xL
    Pproj <- c(xL, yL)
    names(Pproj) <- c("x", "y")
    
    if (addLine) graphics::abline(a, b)
    
    if (add) {
        graphics::points(x0, y0, col = color)
        graphics::points(xL, yL, col = color, pch = 19)
        graphics::arrows(x0, y0, xL, yL, col = color,
               length = 0.1, angle = 20)
    }
    
    # convert slope b to theta angle
    theta <- atan(b)
    w_t <- c(cos(theta), sin(theta))
    td <- Pxy %*% w_t
    dist0 <- sign(td) * sqrt(sum(Pproj^2))
    if (d0) graphics::arrows(0, a, xL, yL, 
                   length = 0.1, angle = 20, lwd = 2, col = color)
    
    return(list(projection = Pproj,
                distance = dist0))
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


#' @noRd
# Default fill palette when the caller doesn't supply `cols`: the binary
# c("steelblue", "tomato") pair used by axialLine()/radialCircle()/
# radialEllipse(), or hcl.colors(k, "Pastel 1") used by axialLines()/
# radialCircles()/radialEllipses()/angularPartition() for k > 2. Shared so
# plotPartition() reproduces whichever palette the originating call would
# have used by default.
.default_fill_cols <- function(k) {
    if (k == 2L) c("steelblue", "tomato") else hcl.colors(k, palette = "Pastel 1")
}


#' @noRd
# Redraw the k-1 nested circles of a radialCircle()/radialCircles() result.
# radialCircle() (binary) stores a single `center`/`radius` pair rather than
# the `cx`/`cy`/`radii` vectors radialCircles() stores; both are normalised
# here to vectors of length k - 1 before drawing. Fill logic mirrors
# radialCircles()'s own (radial-circle.R): even-odd polypath rings from the
# outside in, with the innermost disc drawn directly.
.plot_circle_partition <- function(Pout, fill, cols, col, lwd, lty) {
    if (!is.null(Pout$radii)) {
        cx_vec <- Pout$cx; cy_vec <- Pout$cy; radii_vec <- Pout$radii
    } else {
        cx_vec <- Pout$center[1]; cy_vec <- Pout$center[2]; radii_vec <- Pout$radius
    }
    n <- length(radii_vec)
    k <- n + 1L

    if (fill) {
        if (is.null(cols)) cols <- .default_fill_cols(k)
        usr    <- par("usr")
        n_c    <- 200L
        rect_x <- c(usr[1], usr[2], usr[2], usr[1])
        rect_y <- c(usr[3], usr[3], usr[4], usr[4])

        outer_c <- .circle_pts(cx_vec[n], cy_vec[n], radii_vec[n], n_c)
        polypath(x = c(rect_x, NA, outer_c[, 1]), y = c(rect_y, NA, outer_c[, 2]),
                 rule = "evenodd", col = adjustcolor(cols[k], alpha.f = 0.15), border = NA)

        if (n >= 2L) {
            for (s in n:2L) {
                outer_c <- .circle_pts(cx_vec[s],      cy_vec[s],      radii_vec[s],      n_c)
                inner_c <- .circle_pts(cx_vec[s - 1L], cy_vec[s - 1L], radii_vec[s - 1L], n_c)
                polypath(x = c(outer_c[, 1], NA, inner_c[, 1]),
                         y = c(outer_c[, 2], NA, inner_c[, 2]),
                         rule = "evenodd", col = adjustcolor(cols[s], alpha.f = 0.15), border = NA)
            }
        }

        draw.circle(cx_vec[1L], cy_vec[1L], radii_vec[1L],
                    col = adjustcolor(cols[1L], alpha.f = 0.15), border = NA)
    }

    for (s in seq_len(n))
        draw.circle(cx_vec[s], cy_vec[s], radii_vec[s], border = col, lwd = lwd, lty = lty)
}


#' @noRd
# Redraw the k-1 nested ellipses of a radialEllipse()/radialEllipses()
# result. Both functions already share the same field names (`cx`, `cy`,
# `a`, `b`, `angle`) -- scalars for radialEllipse(), length-(k-1) vectors for
# radialEllipses() -- so no per-caller branching is needed here, unlike the
# circle case above. Fill logic mirrors radialEllipses()'s own
# (radial-ellipse.R).
.plot_ellipse_partition <- function(Pout, fill, cols, col, lwd, lty) {
    cx_vec    <- Pout$cx
    cy_vec    <- Pout$cy
    a_vec     <- Pout$a
    b_vec     <- Pout$b
    angle_vec <- Pout$angle
    n <- length(cx_vec)
    k <- n + 1L

    if (fill) {
        if (is.null(cols)) cols <- .default_fill_cols(k)
        usr    <- par("usr")
        rect_x <- c(usr[1], usr[2], usr[2], usr[1])
        rect_y <- c(usr[3], usr[3], usr[4], usr[4])

        outer_bnd <- .ellipse_pts(cx_vec[n], cy_vec[n], a_vec[n], b_vec[n], angle_vec[n])
        polypath(x = c(rect_x, NA, outer_bnd[, 1]), y = c(rect_y, NA, outer_bnd[, 2]),
                 rule = "evenodd", col = adjustcolor(cols[k], alpha.f = 0.15), border = NA)

        if (n >= 2L) {
            for (s in n:2L) {
                outer_bnd <- .ellipse_pts(cx_vec[s], cy_vec[s], a_vec[s], b_vec[s], angle_vec[s])
                inner_bnd <- .ellipse_pts(cx_vec[s - 1L], cy_vec[s - 1L],
                                          a_vec[s - 1L], b_vec[s - 1L], angle_vec[s - 1L])
                polypath(x = c(outer_bnd[, 1], NA, inner_bnd[, 1]),
                         y = c(outer_bnd[, 2], NA, inner_bnd[, 2]),
                         rule = "evenodd", col = adjustcolor(cols[s], alpha.f = 0.15), border = NA)
            }
        }

        draw.ellipse(cx_vec[1L], cy_vec[1L], a = a_vec[1L], b = b_vec[1L],
                     angle = angle_vec[1L] * 180 / pi,
                     col = adjustcolor(cols[1L], alpha.f = 0.15), border = NA)
    }

    for (s in seq_len(n))
        draw.ellipse(cx_vec[s], cy_vec[s], a = a_vec[s], b = b_vec[s],
                     angle = angle_vec[s] * 180 / pi, border = col, lwd = lwd, lty = lty)
}


#' @noRd
# Redraw the k-1 parallel separator(s) of an axialLine()/axialLines() result.
# axialLine() (binary LDA) stores a single `intercept`; axialLines() stores
# `intercepts` (length k - 1) plus the exact search `angle`. Both need the
# projection direction `w` that `intercepts` was expressed in: axialLines()'s
# `angle` reconstructs it exactly, but axialLine() has no such field, so its
# `w` is rebuilt from `slope` up to an orientation choice -- that choice only
# relabels which end of the split gets cols[1] in the fill below, it never
# changes the drawn line itself (`abline()` uses `slope`/`intercept` directly).
.plot_axial_partition <- function(Pout, fill, cols, col, lwd, lty) {
    slope      <- Pout$slope
    vertical   <- is.infinite(slope)
    intercepts <- if (!is.null(Pout$intercepts)) Pout$intercepts else Pout$intercept
    n <- length(intercepts)
    k <- n + 1L

    if (!is.null(Pout$angle)) {
        w <- c(cos(Pout$angle), sin(Pout$angle))
    } else if (vertical) {
        w <- c(1, 0)
    } else {
        w2 <- 1 / sqrt(1 + slope^2)
        w  <- c(-slope * w2, w2)
    }
    cuts <- if (vertical) intercepts * w[1] else intercepts * w[2]

    if (fill) {
        if (is.null(cols)) cols <- .default_fill_cols(k)
        usr    <- par("usr")
        rect_x <- c(usr[1], usr[2], usr[2], usr[1])
        rect_y <- c(usr[3], usr[3], usr[4], usr[4])

        cuts_sorted <- sort(cuts)
        cuts_ext    <- c(-Inf, cuts_sorted, Inf)

        for (s in seq_len(k)) {
            lo   <- cuts_ext[s]
            hi   <- cuts_ext[s + 1L]
            poly <- list(x = rect_x, y = rect_y)
            if (is.finite(hi)) poly <- .clip_halfplane(poly$x, poly$y, w, hi, keep_le = TRUE)
            if (is.finite(lo)) poly <- .clip_halfplane(poly$x, poly$y, w, lo, keep_le = FALSE)
            if (length(poly$x) >= 3L)
                polygon(poly$x, poly$y, col = adjustcolor(cols[s], alpha.f = 0.15), border = NA)
        }
    }

    for (ic in intercepts) {
        if (vertical) abline(v = ic, col = col, lwd = lwd, lty = lty)
        else          abline(a = ic, b = slope, col = col, lwd = lwd, lty = lty)
    }
}


#' @noRd
# Redraw the k rays of an angularPartition() result. Mirrors
# angularPartition()'s own drawing code (angular.R): a pie-slice polygon per
# wedge for the fill, relying on default plot-region clipping, then the cut
# rays themselves.
.plot_angular_partition <- function(Pout, fill, cols, col, lwd, lty) {
    cx <- Pout$center[1]; cy <- Pout$center[2]
    k  <- length(Pout$cuts)
    sc <- sort(Pout$cuts %% (2 * pi))

    usr     <- par("usr")
    ray_len <- 2 * sqrt((usr[2] - usr[1])^2 + (usr[4] - usr[3])^2)

    if (fill) {
        if (is.null(cols)) cols <- hcl.colors(k, palette = "Pastel 1")

        wedge_starts <- sc
        wedge_ends   <- c(sc[-1L], sc[1L] + 2 * pi)
        n_arc        <- 100L

        for (r in seq_len(k)) {
            arc_ang <- seq(wedge_starts[r], wedge_ends[r], length.out = n_arc)
            polygon(c(cx, cx + ray_len * cos(arc_ang)),
                    c(cy, cy + ray_len * sin(arc_ang)),
                    col = adjustcolor(cols[r], alpha.f = 0.15), border = NA)
        }
    }

    for (ang in Pout$cuts)
        segments(cx, cy, cx + ray_len * cos(ang), cy + ray_len * sin(ang),
                 col = col, lwd = lwd, lty = lty)
}


#' Replot a stored partition result
#'
#' Redraws the configuration and partitions from the output list returned by
#' [axialLine()], [axialLines()], [radialCircle()], [radialCircles()],
#' [radialEllipse()], [radialEllipses()] or [angularPartition()]. It is a 
#' simple replot of the stored partition regions with no recomputation.
#'
#' @param Pout A result list from one of the partition functions above (must
#'   have a `partition` field naming a known family 
#'   ("axial", "angular", "circle", "ellipse"), plus `pcoords`/`pgroup`).
#' @param fill If `TRUE`, shade the regions between boundaries (default
#'   `FALSE`).
#' @param cols Fill colours; auto-generated if `NULL` -- `c("steelblue",
#'   "tomato")` for a 2-region partition, `hcl.colors(k, "Pastel 1")`
#'   otherwise, matching the palette the originating partition function would
#'   have used by default. A char vector with k valid colors can be supplied.
#' @param col Boundary line/ray colour; if `NULL` (default), `"darkorange"`
#'   for `partition == "angular"` and `"purple"` otherwise, matching the
#'   originating functions' own defaults.
#' @param lwd Line width (default `2`).
#' @param lty Line type (default `1`).
#' @param highlight If `TRUE`, overlay a marker on every misclassified point
#'   via [highlightMisclass()] (default `FALSE`).
#' @param add If `FALSE` (default), start a fresh plot of `Pout$pcoords` with
#'   group labels before drawing the boundaries; if `TRUE`, draw only the
#'   boundaries (and optional fill/highlight) on the already-active plot.
#'
#' @return `invisible(NULL)`, called for its side effect of drawing on the
#'   active plot.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' crd <- rbind(cbind(rnorm(15, -1), rnorm(15)), cbind(rnorm(15, 1), rnorm(15)))
#' grp <- factor(c(rep("a", 15), rep("b", 15)))
#' res <- radialCircles(crd, grp, output = TRUE, add = FALSE)
#' plotPartition(res, fill = TRUE, highlight = TRUE)
#' }
#'
#' @export
plotPartition <- function(Pout,
                           fill      = FALSE,
                           cols      = NULL,
                           col       = NULL,
                           lwd       = 2,
                           lty       = 1,
                           highlight = FALSE,
                           add       = FALSE) {

    partitions_set <- c("axial", "angular", "circle", "ellipse")
    if (!is.list(Pout) || is.null(Pout$partition) ||
        !(Pout$partition %in% partitions_set) ||
        is.null(Pout$pcoords) || is.null(Pout$pgroup))
        stop("Not a valid partition object")

    pp <- Pout$partition

    if (is.null(col)) col <- if (pp == "angular") "darkorange" else "purple"

    if (!add) {
        plot(Pout$pcoords, las = 1, asp = 1)
        graphics::text(Pout$pcoords, labels = Pout$pgroup, cex = 0.7, pos = 4)
    }

    switch(pp,
           axial   = .plot_axial_partition(Pout, fill, cols, col, lwd, lty),
           circle  = .plot_circle_partition(Pout, fill, cols, col, lwd, lty),
           ellipse = .plot_ellipse_partition(Pout, fill, cols, col, lwd, lty),
           angular = .plot_angular_partition(Pout, fill, cols, col, lwd, lty))

    if (highlight) highlightMisclass(Pout)

    invisible(NULL)
}
