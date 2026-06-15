# ===================================================================
# ==== Functions for Angular Partitions ====
# ===================================================================


# ---- Internal helpers ----

#' @noRd
.arc_mid <- function(a, b) {
    # given two angles a < b, find the middle one
    # e.g.: two radians a < b, delta is the diff b - a
    # returns the radiant of the line bisecting b - a
    # delta (also if cycled around 2pi):
    delta <- (b - a) %% (2 * pi)
    return(atan2(sin(a + delta / 2), 
                 cos(a + delta / 2)) )
}


#' @noRd
.angular_search <- function(coords, 
                            grp_int,
                            cx, cy,
                            k, n_pts, combos,
                            full = TRUE) {
    # Given center cx, cy: find the best angular k-partition by brute force.
    # A candidate partition is a choice of `k` cut-gaps out of the `n` gaps
    # between angle-sorted points; `combos` supplies these as `combn(n, k)`.
    # Per-arc majority counts are read from a cumulative group-count table 
    # per arc, without re-tabulating. With `full = FALSE` only the minimum
    # misclassification is returned (margin and cut angles skipped) for use as
    # the optimiser objective; `full = TRUE` additionally returns the 2D margin
    # and the cut angles of the best partition.

    # order points according to their angles (in radians)
    # Note: atan2(y, x) has y-coordinate first!
    pt_angles <- atan2(coords[, 2] - cy, coords[, 1] - cx)
    ord       <- order(pt_angles)
    s_ang     <- pt_angles[ord]
    s_grp     <- grp_int[ord]

    # Cumulative per-group counts over the angle-sorted sequence:
    # cum[i + 1, g] = number of group-g points among sorted points 1..i.
    # The count of group g on an arc a..b is cum[b + 1, g] - cum[a, g].
    cum <- matrix(0L, nrow = n_pts + 1, ncol = k)
    # loop over all k columns: cumulate group number
    for (g in seq_len(k)) cum[-1L, g] <- cumsum(s_grp == g)
    total <- cum[n_pts + 1L, ]

    # ordered radius lengths
    if (full) {
        s_rad <- sqrt((coords[, 1] - cx)^2 + (coords[, 2] - cy)^2)[ord]
    }

    best_err    <- n_pts + 1L
    best_margin <- -Inf
    best_gaps   <- NULL

    # loop over all ci = 1..combn(n, k) cut-gaps from combos:
    for (ci in seq_len(ncol(combos))) {
        # g is one specific partitioning,
        # increasing gap indices g[1] < .. < g[k]
        g <- combos[, ci]   

        # interior arcs (g[s]+1 .. g[s+1]); 
        #   branch-and-bound on partial error.
        # Strict (> best_err) when full, 
        #   so co-minimal partitions survive for the margin tie-break; 
        # loose (>= best_err) otherwise, since only the
        #   minimum count is needed.
        err  <- 0L
        drop <- FALSE
        
        # loop over groups s = 1..k-1 in a given g cut (from combos):
        for (s in seq_len(k - 1L)) {
            # cnt: group frequencies in a given arc g[s]..g[s+1]
            cnt <- cum[g[s + 1L] + 1L, ] - cum[g[s] + 1L, ]
            # err: misclassification = n points in cut - maximum group
            #      cumulates over all arcs
            err <- err + (g[s + 1L] - g[s]) - max(cnt)
            if (if (full) err > best_err else err >= best_err) {
                drop <- TRUE
                break
            }
        }
        if (drop) next  # something went wrong ...

        # wrap-around arc (g[k]+1 .. n, 1 .. g[1]): from last to first cut
        cnt_wrap <- (total - cum[g[k] + 1L, ]) + cum[g[1L] + 1L, ]
        # add wrap-around to err
        err      <- err + (n_pts - g[k] + g[1L]) - max(cnt_wrap)

        if (!full) {
            if (err < best_err) {
                best_err  <- err
                best_gaps <- g
                if (best_err == 0L) break        # cannot improve on zero
            }
            next
        }

        # full path: keep min error, break ties by largest 2D margin.
        # 2D margin at a gap: min(r at its two points) * sin(angular gap / 2).
        if (err > best_err) next
        
        # next cut after current g:
        nx      <- g + 1L
        nx[g == n_pts] <- 1L
        
        ang_gap <- (s_ang[nx] - s_ang[g]) %% (2 * pi)
        margin  <- min(pmin(s_rad[g], s_rad[nx]) * sin(ang_gap / 2))
        if (err < best_err || margin > best_margin) {
            best_err    <- err
            best_margin <- margin
            best_gaps   <- g
        }
    }  # end of loop over all cuts in combos

    
    if (!full) {
        return(list(misclass = best_err, pt_angles = pt_angles))
    }

    # cut angles placed at the arc midpoint of each chosen gap
    g  <- best_gaps
    nx <- g + 1L; nx[g == n_pts] <- 1L
    list(misclass  = best_err,
         margin    = best_margin,
         cuts      = .arc_mid(s_ang[g], s_ang[nx]),
         pt_angles = pt_angles)
}


#' Angular k-way (wedge) partition of a 2D configuration
#'
#' Finds `k` rays (straight lines) emanating from a center point that partition 
#' a 2D configuration of `n` points into `k` angular sectors, minimising total
#' misclassification (points whose group differs from the majority group
#' in their sector), and draws the lines. (see: Shye, S. (2014). Faceted Smallest Space Analysis (FSSA). 
#' In A. Michalos (Ed.), Encyclopedia of quality of life research (pp. 2129-2133). 
#' New York: Springer.)
#'
#' **Search at fixed center.** Points are sorted by their angle from the
#' center. A partition is a choice of `k` cut-gaps among the `n` gaps
#' between consecutive points (`combn(n, k)` candidates); per-arc majority
#' counts are read from a cumulative count table and total misclassification
#' is minimised. Cut angles are placed at the arc midpoint between adjacent
#' points. Tie-breaker: among k-tuples with the same misclass count, the one
#' with the largest minimum 2D perpendicular distance from a cut ray to the
#' nearest point wins.
#'
#' **Search optimal center.** When `cx` and `cy` are `NULL` (default), the center
#' is optimised by multi-start Nelder-Mead — the brute-force above runs
#' as the inner objective at each candidate center. Starts are the data
#' centroid plus the centroid of each non-empty group. `parscale` is set
#' to the data range. When `cx` and `cy` are supplied, they are used
#' directly (no optimisation).
#'
#' @param crd Numeric matrix or numeric data frame with exactly 2 columns; no NAs.
#' @param group Factor with `n >= k >= 2` levels, same length as nrow of crd.
#' @param cx,cy Center; optimised when either is `NULL`.
#' @param output If `TRUE` (default), return list of results.
#' @param col Ray colour (default `"darkorange"`).
#' @param lwd Line width (default `2`).
#' @param lty Line type (default `1`).
#' @param add If `TRUE` (default), add to existing plot, else plot configuration.
#'
#' @return If `output = TRUE`, a list with:
#'   - `cuts` — `numeric[k]`, cut angles of rays (in radians).
#'   - `margin` — minimum distance from any cut ray to nearest point.
#'   - `misclass` — integer, number of misclassified points.
#'   - `misclass_points` — data frame (`x`, `y`, `label`) of misclassified points
#'   - `sector` — `integer[n]`, assigned sector `1..k` per point
#'   - `majority` — `character[k]`, majority group per sector
#'   - `center` — `c(cx, cy)` center coordinates of rays
#'   - `pt_angles` — angles (in radians) of points from center
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' theta <- rep(c(0, 2 * pi / 3, 4 * pi / 3), each = 12) + rnorm(36, 0, 0.25)
#' r     <- runif(36, 0.5, 1.5)
#' crd   <- cbind(r * cos(theta), r * sin(theta))
#' grp   <- factor(rep(c("a", "b", "c"), each = 12))
#' plot(crd, asp = 1)
#' angularPartition(crd, grp)
#' }
#'
#' @export
angularPartition <- function(crd,
                             group,
                             cx = NULL,
                             cy = NULL,
                             output = TRUE,
                             col = "darkorange",
                             lwd = 2,
                             lty = 1,
                             add = TRUE) {

    # ---- Input validation ----
    
    if (!all(is.numeric(as.matrix(crd))))  stop("Input matrix crd must be numeric!")
    if (!is.factor(group)) stop("group must be a factor!")
    if (any(is.na(crd)))            stop("No NA allowed in input!")
    if (length(dim(crd)) != 2)      stop("Input data must have two dimensions!")
    if (dim(crd)[2] != 2)           stop("Coordinates must be 2-dimensional!")
    if (nrow(crd) != length(group)) stop("nrow(crd) must equal length(group)!")
    
    group <- factor(group, exclude = NA)
    coords  <- as.matrix(crd)
    
    k <- nlevels(group)
    if (k < 2)         stop("group must have at least 2 levels!")
    if (nrow(crd) < k) stop("Number of points must be >= number of groups!")

    
    # ---- Basic parameters ----
    
    n_pts   <- nrow(coords)
    grp_int <- as.integer(group)

    # Generate all possible partitions = cut-gap combinations: 
    # choose k of the n gaps between all points. 
    # Each angular k-partition corresponds to exactly
    # one such choice (see .angular_search).
    combos <- combn(n_pts, k)
    

    # ---- Center: optimise if cx or cy is NULL, use as given otherwise ----
    
    # Optimise center
    if (is.null(cx) || is.null(cy)) {
        x_range <- diff(range(coords[, 1]))
        y_range <- diff(range(coords[, 2]))
        if (x_range == 0) x_range <- 1
        if (y_range == 0) y_range <- 1
        parscale <- c(x_range, y_range)

        # start values for center: overall mean, group means
        starts <- list(c(mean(coords[, 1]), mean(coords[, 2])))
        for (g in seq_len(k)) {
            in_g <- which(grp_int == g)
            if (length(in_g) > 0L) {
                starts <- c(starts,
                            list(c(mean(coords[in_g, 1]),
                                   mean(coords[in_g, 2]))))
            }
        }

        # function to optimize: n of misclassification from .angular_search()
        # parameter to optimize: p = center coordinates
        fnToOpt <- function(p) {
            .angular_search(coords, grp_int, p[1], p[2],
                            k, n_pts, combos, full = FALSE)$misclass
        }

        best <- NULL
        
        # minimization loops over all starting values from starts:
        for (s0 in starts) {
            
            opt <- optim(par     = s0,
                         fn      = fnToOpt,
                         method  = "Nelder-Mead",
                         control = list(reltol = 1e-8, 
                                        maxit = 2000,
                                        parscale = parscale))
            
            if (is.null(best) || opt$value < best$value) best <- opt
            if (best$value == 0) break          # zero misclass is optimal
        }
        cx <- best$par[1]
        cy <- best$par[2]
    }

    # ---- Final search at chosen center ----
    res         <- .angular_search(coords, grp_int,
                                   cx, cy,
                                   k, n_pts, combos, full = TRUE)
    
    best_err    <- res$misclass
    best_cuts   <- res$cuts
    best_margin <- res$margin
    pt_angles   <- res$pt_angles
    
    
    # ---- Plot points if add = FALSE ----
    if (!add) {
        plot(coords, asp = 1)
        graphics::text(coords, labels = group, cex = 0.7, pos = 4)
    }
    
    # ---- Draw rays from center at the optimal cut angles ----
    usr     <- par("usr")
    ray_len <- 2 * sqrt( (usr[2] - usr[1])^2 + (usr[4] - usr[3])^2 )

    for (ang in best_cuts) {
        segments(cx, cy,
                 cx + ray_len * cos(ang),
                 cy + ray_len * sin(ang),
                 col = col, lwd = lwd, lty = lty)
    }

    # ---- Sector assignment for each point (in original input order) ----
    sc       <- sort(best_cuts %% (2 * pi))
    norm_pts <- pt_angles %% (2 * pi)
    fi       <- findInterval(norm_pts, sc)
    sector   <- ifelse(fi == 0L | fi == k, k, fi)

    # ---- Unique group assignment ----
    count_mat <- matrix(0L, nrow = k, ncol = k)
    for (r in seq_len(k)) {
        pts_r <- grp_int[sector == r]
        if (length(pts_r) > 0L)
            count_mat[, r] <- tabulate(pts_r, nbins = k)
    }
    assignment <- .assign_groups(count_mat)
    majority   <- levels(group)[assignment]

    misclass_idx    <- which(as.character(group) != majority[sector])
    misclass_points <- data.frame(
        x     = coords[misclass_idx, 1],
        y     = coords[misclass_idx, 2],
        label = group[misclass_idx]
    )

    if (!output) return(invisible(NULL))
    else return(list(
        cuts            = best_cuts,
        margin          = best_margin,
        misclass        = length(misclass_idx),
        misclass_points = misclass_points,
        sector          = sector,
        majority        = majority,
        center          = c(cx, cy),
        pt_angles       = pt_angles))

}

