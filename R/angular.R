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
    # between angle-sorted points; `combos` supplies these as `combn(n, k)`,
    # M = ncol(combos) candidates in total. The misclassification of a
    # candidate is the bijection-constrained error (best achievable total
    # correct over all k! arc-to-group assignments, each group used exactly
    # once) -- the same criterion `.assign_groups()` applies to the final
    # sector/majority mapping in angularPartition(), so the search targets
    # exactly the quantity the function ultimately reports.
    #
    # All M candidates are scored with vector arithmetic over the candidate
    # axis, analogous to the vectorised cut search in axialLines(): every
    # quantity below is a length-M vector built by indexing the cumulative
    # count table, so the loops run over k / k! (small, fixed, independent
    # of M) and never over candidates.
    #
    # Scoring is a two-stage bound-and-prune, because the exact k!
    # assignment scan is the expensive part:
    #
    #   Stage 1 (all M, cheap) -- `ub`, the *unconstrained* best-correct:
    #     each arc independently keeps its own local-majority group. Since
    #     any bijection is one particular choice of group per arc,
    #     exact(c) <= ub(c) for every candidate, so `ub` is a valid upper
    #     bound. It costs k * (2k - 1) vector ops versus k! * (k + 1) for
    #     the exact scan (~7x cheaper at k = 4) and needs no M x k matrix.
    #
    #   Stage 2 (tiny subset, exact) -- only candidates whose bound is high
    #     enough to still win are scanned over all k! assignments. Any
    #     candidate with ub < the best exact score found cannot beat it and
    #     is skipped. The bound is tight in practice (typically only a
    #     handful of the M candidates survive), which is what makes this
    #     much faster than scoring every candidate exactly.
    #
    # `full = TRUE` additionally needs the 2D margin, but only to break
    # ties among the co-minimal candidates -- so the margin is computed for
    # those alone, not for all M.

    # order points according to their angles (in radians)
    # Note: atan2(y, x) has y-coordinate first!
    pt_angles <- atan2(coords[, 2] - cy, coords[, 1] - cx)
    ord       <- order(pt_angles)
    s_ang     <- pt_angles[ord]
    s_grp     <- grp_int[ord]

    # Cumulative per-group counts over the angle-sorted sequence:
    # cum[i + 1, g] = number of group-g points among sorted points 1..i.
    # The count of group g on an arc a..b is cum[b + 1, g] - cum[a, g].
    cum <- matrix(0L, nrow = n_pts + 1L, ncol = k)
    # loop over all k columns: cumulate group number
    for (g in seq_len(k)) cum[-1L, g] <- cumsum(s_grp == g)
    total <- cum[n_pts + 1L, ]

    M <- ncol(combos)

    # Row indices into `cum` for each cut: end_idx[[s]] = combos[s, ] + 1L.
    # Interior arc s (s < k) spans sorted positions combos[s, ] + 1 ..
    # combos[s + 1, ]; arc k wraps combos[k, ] + 1 .. n, 1 .. combos[1, ].
    end_idx <- vector("list", k)
    for (s in seq_len(k)) end_idx[[s]] <- combos[s, ] + 1L

    # Counts of group `g` on arc `s`, as a vector over the candidates in
    # `sel` (an index vector into 1..M, or NULL for all M).
    arc_cnt <- function(s, g, sel) {
        if (s < k) {
            ia <- end_idx[[s]]; ib <- end_idx[[s + 1L]]
            if (!is.null(sel)) { ia <- ia[sel]; ib <- ib[sel] }
            cum[ib, g] - cum[ia, g]
        } else {
            ia <- end_idx[[k]]; ib <- end_idx[[1L]]
            if (!is.null(sel)) { ia <- ia[sel]; ib <- ib[sel] }
            total[g] - cum[ia, g] + cum[ib, g]
        }
    }

    # ---- Stage 1: cheap unconstrained upper bound for all M candidates ----
    ub <- 0L
    for (s in seq_len(k)) {
        arc_max <- arc_cnt(s, 1L, NULL)
        for (g in 2L:k) arc_max <- pmax(arc_max, arc_cnt(s, g, NULL))
        ub <- ub + arc_max
    }

    # ---- Stage 2: exact bijection-constrained score, on survivors only ----
    perms <- .perms_k(k)

    # Exact best-correct over all k! arc-to-group assignments (each group
    # used exactly once) for the candidates in `sel` -- the same criterion
    # .assign_groups() applies to the winning combo downstream.
    exact_correct <- function(sel) {
        cn <- vector("list", k)
        for (s in seq_len(k)) {
            cs <- vector("list", k)
            for (g in seq_len(k)) cs[[g]] <- arc_cnt(s, g, sel)
            cn[[s]] <- cs
        }
        best <- rep(0L, length(sel))
        for (p_idx in seq_len(nrow(perms))) {
            perm    <- perms[p_idx, ]
            correct <- cn[[1L]][[perm[1L]]]
            for (s in 2L:k) correct <- correct + cn[[s]][[perm[s]]]
            best <- pmax(best, correct)
        }
        best
    }

    # Walk the bound down until a candidate's exact score reaches it: at
    # that point every unscanned candidate has ub < best and cannot win.
    thresh       <- max(ub)
    max_correct  <- -1L
    scanned      <- logical(M)
    repeat {
        sel <- which(ub >= thresh & !scanned)
        if (length(sel) > 0L) {
            max_correct   <- max(max_correct, max(exact_correct(sel)))
            scanned[sel]  <- TRUE
        }
        if (max_correct >= thresh) break
        thresh <- thresh - 1L
    }
    m <- n_pts - max_correct

    if (!full) {
        return(list(misclass = m, pt_angles = pt_angles))
    }

    # ---- Co-minimal candidates, then the margin tie-break among them ----
    # Any candidate attaining max_correct must have ub >= max_correct, so
    # this pool is a complete superset of the co-minimal set.
    pool <- which(ub >= max_correct)
    at   <- pool[exact_correct(pool) == max_correct]

    # ordered radius lengths, needed only for the 2D margin tie-break
    s_rad <- sqrt((coords[, 1] - cx)^2 + (coords[, 2] - cy)^2)[ord]

    # Margin of a candidate: min over its k gaps of
    # min(r at the gap's two adjacent points) * sin(angular gap / 2).
    A      <- length(at)
    g_mat  <- t(combos[, at, drop = FALSE])
    nx_mat <- g_mat + 1L
    nx_mat[g_mat == n_pts] <- 1L

    ang_gap_mat <- (matrix(s_ang[nx_mat], nrow = A) -
                    matrix(s_ang[g_mat],  nrow = A)) %% (2 * pi)
    margin_mat  <- pmin(matrix(s_rad[g_mat],  nrow = A),
                        matrix(s_rad[nx_mat], nrow = A)) * sin(ang_gap_mat / 2)
    margin_at   <- do.call(pmin, as.data.frame(margin_mat))

    # Ties broken by largest margin. which.max returns the first maximum and
    # `at` is increasing, so the earliest candidate in scan order wins ties,
    # matching axialLines()'s convention.
    w <- which.max(margin_at)
    j <- at[w]

    # cut angles placed at the arc midpoint of each chosen gap
    g  <- combos[, j]
    nx <- g + 1L; nx[g == n_pts] <- 1L
    list(misclass  = m,
         margin    = margin_at[w],
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
#' between consecutive points (`combn(n, k)` candidates, `M` in total). The
#' score of a candidate is its exact bijection-constrained misclassification
#' (best achievable total correct over all `k!` arc-to-group assignments,
#' each group used exactly once) -- the same region-to-group matching
#' criterion used to derive `sector`/`majority` below, so the search targets
#' exactly the quantity reported as `misclass`. Every wedge is non-empty, so a
#' partition leaving a group's sector empty is not considered (as in
#' [axialLines()]; the radial functions do allow it).
#'
#' All per-candidate quantities are computed as length-`M` vectors indexed
#' out of a cumulative group-count table, so the loops run over `k` / `k!`
#' rather than over candidates -- analogous to the vectorised cut search in
#' [axialLines()]. Scoring proceeds in two stages: a cheap upper bound (each
#' arc independently taking its own local-majority group, which no bijection
#' can beat) is evaluated for all `M` candidates, and only those whose bound
#' is high enough to still win are then scored exactly over all `k!`
#' assignments. The bound is tight in practice, so typically only a handful
#' of candidates need the exact scan.
#'
#' Cut angles are placed at the arc midpoint between adjacent points.
#' Tie-breaker: among k-tuples with the same misclass count, the one with the
#' largest minimum 2D perpendicular distance from a cut ray to the nearest
#' point wins (computed only for the co-minimal candidates).
#'
#' **Search optimal center.** When `cx` and `cy` are `NULL` (default), the center
#' is optimised by multi-start Nelder-Mead — the brute-force above runs
#' as the inner objective at each candidate center. Starts are the data
#' centroid plus the centroid of each non-empty group and the four
#' bounding-box corners. `parscale` is set to the data range. Because
#' the inner objective is piecewise-constant, Nelder-Mead can get stuck
#' on a flat plateau around any of those starts without ever seeing a
#' better region; if none of them reaches zero misclassification, a
#' coarse `n_grid` x `n_grid` scan of the bounding box, padded by half the
#' data range on each side (the best center is not always inside the
#' convex hull of the points), locates a cell in a better region and one
#' more Nelder-Mead run is seeded there. When `cx` and `cy` are supplied,
#' they are used directly (no optimisation).
#'
#' @param crd Numeric matrix or numeric data frame with exactly 2 columns; no NAs.
#' @param group Factor with `n >= k >= 2` levels, same length as nrow of crd.
#' @param cx,cy Center; optimised when either is `NULL`.
#' @param n_grid Side length of the coarse fallback grid (`n_grid^2` extra
#'   evaluations of the cheap inner search); only used when the heuristic
#'   starts don't already reach zero misclassification. Default `7L`.
#' @param fill If `TRUE`, shade the `k` wedge sectors (default `FALSE`).
#' @param output If `TRUE` (default), return list of results.
#' @param col Ray colour (default `"darkorange"`).
#' @param cols Length-`k` fill colours; auto-generated if `NULL`.
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
#' angularPartition(crd, grp, fill = TRUE)
#' }
#'
#' @export
angularPartition <- function(crd,
                             group,
                             cx = NULL,
                             cy = NULL,
                             n_grid = 7L,
                             fill = FALSE,
                             output = TRUE,
                             col = "darkorange",
                             cols = NULL,
                             lwd = 2,
                             lty = 1,
                             add = TRUE) {

    # ---- Input validation ----
    
    if (!all(is.numeric(as.matrix(crd))))  stop("Coordinate data must be numeric!")
    if (any(is.na(crd)))            stop("No NAs allowed in crd!")
    if (any(is.na(group)))          stop("No NAs allowed in group!")
    if (length(dim(crd)) != 2)      stop("Coordinates must have two dimensions!")
    if (dim(crd)[2] != 2)           stop("Coordinates must have 2 columns!")
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
        # parscale <- c(1, 1)

        # start values for center: overall mean, group means, max, min
        starts <- list(c(mean(coords[, 1]), mean(coords[, 2])))
        
        for (g in seq_len(k)) {
            in_g <- which(grp_int == g)
            if (length(in_g) > 0L) {
                starts <- c(starts,
                            list(c(mean(coords[in_g, 1]),
                                   mean(coords[in_g, 2]))))
            }
        }
        
        starts <- c(starts, list(c(min(coords[ , 1]), 
                                   min(coords[ , 2]))))
        starts <- c(starts, list(c(max(coords[ , 1]), 
                                   max(coords[ , 2]))))
        starts <- c(starts, list(c(min(coords[ , 1]), 
                                   max(coords[ , 2]))))
        starts <- c(starts, list(c(max(coords[ , 1]), 
                                   min(coords[ , 2]))))
        

        # function to optimize: n of misclassification from .angular_search()
        # parameter to optimize: p = center coordinates
        fnToOpt <- function(p) {
            .angular_search(coords, grp_int, p[1], p[2],
                            k, n_pts, combos, full = FALSE)$misclass
        }

        best <- NULL

        # minimization loops over all starting values from starts:
        for (s0 in starts) {
            s0 <- .snap_zero(s0, parscale)

            opt <- optim(par     = s0,
                         fn      = fnToOpt,
                         method  = "Nelder-Mead",
                         control = list(reltol = 1e-8,
                                        maxit = 3000,
                                        parscale = parscale))

            if ((is.null(best)) || (opt$value < best$value)) best <- opt
            if (best$value == 0) break          # zero misclass is optimal
        }

        # ---- Coarse grid fallback ----
        # fnToOpt() is piecewise-constant, so Nelder-Mead's local,
        # gradient-free search can get stuck on a flat plateau around any
        # heuristic start above without ever seeing a better region --
        # even past the exact-zero degeneracy .snap_zero() rules out, its
        # fallback step is still only a fixed fraction of parscale, and a
        # plateau can be wider than that (confirmed on real MDS data: a
        # plateau of misclass = 1 extending +-0.15 around the centroid,
        # with the true misclass = 0 region starting only around 0.3-0.6
        # away -- outside that fallback step, so every heuristic start
        # stalled on it). Only pay for the scan when the heuristics above
        # didn't already reach the provable optimum (misclass = 0); it
        # gives Nelder-Mead a non-degenerate foothold inside the region
        # the scan found, rather than asking it to escape a plateau blind.
        # `.grid_seeds()` (utils.R) scans both the plain bounding box and
        # one padded by half the data range on each side -- see its
        # comment for why one widened grid isn't enough. Shared with
        # radialCircle(s)() and radialEllipse(s)(), which hit the same
        # plateau-stall pathology on their own center searches.
        if (best$value > 0) {
            for (p in .grid_seeds(fnToOpt, range(coords[, 1]), range(coords[, 2]), n_grid)) {
                opt <- optim(par     = .snap_zero(p, parscale),
                             fn      = fnToOpt,
                             method  = "Nelder-Mead",
                             control = list(reltol = 1e-8,
                                            maxit = 3000,
                                            parscale = parscale))

                if (opt$value < best$value) best <- opt
                if (best$value == 0) break          # zero misclass is optimal
            }
        }

        cx_best <- best$par[1]
        cy_best <- best$par[2]
    }

    # ---- Final search at chosen center ----
    if (is.null(cx) || is.null(cy)) {
        cx <- cx_best
        cy <- cy_best
    }
    
    res         <- .angular_search(coords, grp_int,
                                   cx = cx, cy = cy,
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
    
    # ---- Sorted cut angles: shared by the fill and the sector assignment ----
    sc <- sort(best_cuts %% (2 * pi))

    usr     <- par("usr")
    ray_len <- 2 * sqrt( (usr[2] - usr[1])^2 + (usr[4] - usr[3])^2 )

    # ---- Optional wedge shading ----
    # Each sector r spans (sc[r], sc[r + 1]), with sector k wrapping from
    # sc[k] back to sc[1]. Filled as a pie slice out to ray_len; R clips
    # polygon() to the plot region (par("xpd") is FALSE by default), so the
    # slice need not be trimmed to the visible rectangle by hand.
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

    # ---- Draw rays from center at the optimal cut angles ----
    for (ang in best_cuts) {
        segments(cx, cy,
                 cx + ray_len * cos(ang),
                 cy + ray_len * sin(ang),
                 col = col, lwd = lwd, lty = lty)
    }

    # ---- Sector assignment for each point (in original input order) ----
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
    assignment <- .assign_groups(count_mat, levels(group))
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
        sector          = sector,
        majority        = majority,
        center          = c(cx, cy),
        pt_angles       = pt_angles,
        misclass        = length(misclass_idx),
        misclass_points = misclass_points))

}

