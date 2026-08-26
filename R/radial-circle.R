# ===================================================================
# ==== Functions for Circular Radial Partitions ====
# ===================================================================


# ---- Internal helpers ----

#' @noRd
.circle_pts <- function(cx, cy, r, n = 200L) {
    # draw points of a circle
    theta <- seq(0, 2 * pi, length.out = n + 1L)[-1L]
    cbind(cx + r * cos(theta), cy + r * sin(theta))
}


#' @noRd
.best_radius <- function(s_dist, s_flag, r_min = 0) {
    n        <- length(s_dist)
    cs_not   <- cumsum(!s_flag)
    total_in <- sum(s_flag)
    cs_in    <- cumsum(s_flag)
    r_cands  <- (s_dist[-n] + s_dist[-1L]) / 2
    errs     <- cs_not[-n] + (total_in - cs_in[-n])

    # A split at index i (points 1..i inside) is realised by any radius in
    # [max(s_dist[i], r_min), s_dist[i+1]) -- note the inside test is
    # `dist <= r`. Prefer the midpoint, clamped up to r_min when the nesting
    # constraint binds: testing r_cands[i] >= r_min instead would drop splits
    # that are still reachable at r = r_min. A split with
    # s_dist[i] == s_dist[i+1] is not realisable at all (no radius separates
    # coincident distances), so it must not be scored.
    sp_ok    <- s_dist[-n] < s_dist[-1L] & s_dist[-1L] > r_min
    r_try    <- pmax(r_cands, r_min)
    best_err <- 2L * n
    best_r   <- NA_real_

    if (any(sp_ok)) {
        idx      <- which(sp_ok)
        best_i   <- idx[which.min(errs[idx])]
        best_err <- errs[best_i]
        best_r   <- r_try[best_i]
    }

    r_all_in   <- max(s_dist) * 1.05 + 1e-6
    if (r_all_in < r_min) r_all_in <- r_min * 1.05 + 1e-6
    err_all_in <- sum(!s_flag)
    if (err_all_in < best_err || is.na(best_r)) {
        best_err <- err_all_in
        best_r   <- r_all_in
    }

    if (r_min < s_dist[1L]) {
        err_all_out <- sum(s_flag)
        if (err_all_out < best_err) {
            best_err <- err_all_out
            best_r   <- (r_min + s_dist[1L]) / 2
        }
    }

    return(list(
        r = best_r,
        misclass = best_err))
}


#' @noRd
.optimize_circle <- function(pcoords,
                             inner_flag,
                             prev_cx, prev_cy, prev_r,
                             starts,
                             meth = "Nelder-Mead",
                             n_grid = 7L) {
    has_prev <- !is.null(prev_r)

    parscale <- c(diff(range(pcoords[, 1])), diff(range(pcoords[, 2])))
    parscale[parscale == 0] <- 1

    eval_ctr <- function(ctr) {
        cx_  <- ctr[1]; cy_  <- ctr[2]
        d    <- sqrt((pcoords[, 1] - cx_)^2 + (pcoords[, 2] - cy_)^2)
        rmin <- if (has_prev) sqrt((cx_ - prev_cx)^2 + (cy_ - prev_cy)^2) + prev_r else 0
        ord  <- order(d)
        return(
            .best_radius(d[ord], inner_flag[ord], rmin)$misclass)
    }

    # multi-start: keep the optim() result with lowest misclassification.
    # eval_ctr() is piecewise-constant, so Nelder-Mead can stall on a flat
    # plateau around any of `starts` -- including the centroid, which is
    # ~1e-16 for mean-centered data (e.g. any MDS output) and hits optim()'s
    # near-zero initial-simplex degeneracy on top of that (see .snap_zero(),
    # utils.R). Confirmed on real MDS data: every single-group-vs-rest
    # circle fit here got stuck several misclassifications above the true
    # optimum when started from the (near-zero) overall centroid alone.
    best <- NULL
    for (s0 in starts) {
        opt <- optim(par = .snap_zero(s0, parscale),
                     fn = eval_ctr,
                     method = meth,
                     control = list(reltol = 1e-8, maxit = 2000,
                                    parscale = parscale))
        if (is.null(best) || opt$value < best$value) best <- opt
        if (best$value == 0) break          # zero misclass is optimal
    }

    # ---- Coarse grid fallback ----
    # Only pay for it when the starts above didn't already reach the
    # provable optimum; see .grid_seeds() (utils.R) for why two grids are
    # scanned rather than one.
    if (best$value > 0) {
        for (p in .grid_seeds(eval_ctr, range(pcoords[, 1]), range(pcoords[, 2]), n_grid)) {
            opt <- optim(par = .snap_zero(p, parscale),
                         fn = eval_ctr,
                         method = meth,
                         control = list(reltol = 1e-8, maxit = 2000,
                                        parscale = parscale))
            if (opt$value < best$value) best <- opt
            if (best$value == 0) break      # zero misclass is optimal
        }
    }

    cx_opt <- best$par[1]; cy_opt <- best$par[2]
    d_opt  <- sqrt((pcoords[, 1] - cx_opt)^2 + (pcoords[, 2] - cy_opt)^2)
    r_min  <- if (has_prev) sqrt((cx_opt - prev_cx)^2 + (cy_opt - prev_cy)^2) + prev_r else 0
    ord    <- order(d_opt)
    res    <- .best_radius(d_opt[ord], inner_flag[ord], r_min)

    return(list(
        cx = cx_opt,
        cy = cy_opt,
        r = res$r,
        misclass = res$misclass))
}


#' @noRd
.cut_radii <- function(s_dist) {
    # Realisable cut positions along the distance-sorted sequence, with the
    # radius and margin realising each. Position p means "the p nearest
    # points are inside", realised by any radius in [s_dist[p], s_dist[p+1])
    # -- the inside test is `dist <= r`. So:
    #   p = 0 (nothing inside) needs r < s_dist[1], possible only when
    #     s_dist[1] > 0;
    #   0 < p < n needs s_dist[p] < s_dist[p+1] -- no radius separates
    #     coincident distances, so such a p is NOT realisable;
    #   p = n (everything inside) always is.
    # The margin is the half-gap the boundary sits in, and is 0 for p = 0
    # and p = n: those boundaries separate no points, so on a tie the
    # margin tie-break prefers a partition whose every circle really does
    # split the configuration.
    n   <- length(s_dist)
    pos <- integer(0); rad <- numeric(0); mrg <- numeric(0)
    if (s_dist[1L] > 0) {
        pos <- 0L
        rad <- s_dist[1L] / 2
        mrg <- 0
    }
    if (n >= 2L) {
        ok <- which(s_dist[-n] < s_dist[-1L])
        if (length(ok) > 0L) {
            pos <- c(pos, ok)
            rad <- c(rad, (s_dist[ok] + s_dist[ok + 1L]) / 2)
            mrg <- c(mrg, (s_dist[ok + 1L] - s_dist[ok]) / 2)
        }
    }
    pos <- c(pos, n)
    rad <- c(rad, s_dist[n] * 1.05 + 1e-6)
    mrg <- c(mrg, 0)

    return(list(pos = pos, radius = rad, margin = mrg))
}


#' @noRd
.radial_search <- function(s_dist, s_grp, k, full = TRUE) {
    # Given points already sorted by distance from the center, find the best
    # nested k-partition of the distance axis by brute force. A candidate is
    # a choice of k-1 cut positions p_1 <= ... <= p_{k-1} out of the
    # realisable positions from .cut_radii(); region s spans sorted
    # positions p_{s-1}+1 .. p_s (p_0 = 0, p_k = n). Repeats are allowed, so
    # a region may be empty -- when the data have no radial structure an
    # empty region really is the optimum, and forbidding it would make the
    # reported misclass worse than the criterion's true minimum.
    #
    # Candidates are scored by .bij_best() (utils.R): the exact
    # bijection-constrained misclassification, i.e. the best total correct
    # over all k! region-to-group assignments with each group used exactly
    # once. That is the same criterion .assign_groups() applies to derive
    # sector/majority downstream, so the search targets exactly the quantity
    # reported as `misclass`.
    n_pts <- length(s_dist)
    cp    <- .cut_radii(s_dist)
    m     <- length(cp$pos)

    # Non-decreasing (k-1)-tuples of cut positions = combinations with
    # repetition: take the strictly increasing (k-1)-subsets of
    # 1..(m + k - 2) and shift row i down by i - 1.
    pidx <- combn(m + k - 2L, k - 1L) - (0L:(k - 2L))
    M    <- ncol(pidx)
    pos  <- matrix(cp$pos[pidx], nrow = k - 1L)

    # cum[i + 1, g] = number of group-g points among sorted points 1..i, so
    # the count of group g on positions a+1..b is cum[b + 1, g] - cum[a + 1, g].
    cum <- matrix(0L, nrow = n_pts + 1L, ncol = k)
    for (g in seq_len(k)) cum[-1L, g] <- cumsum(s_grp == g)

    # Row indices into `cum` bounding region s; scalars where constant over
    # candidates (region 1 always starts at position 0, region k always
    # ends at n), which keeps the big vectors down to k - 1 per side.
    lo_idx <- vector("list", k)
    hi_idx <- vector("list", k)
    for (s in seq_len(k)) {
        lo_idx[[s]] <- if (s == 1L) 1L         else pos[s - 1L, ] + 1L
        hi_idx[[s]] <- if (s == k)  n_pts + 1L else pos[s, ]      + 1L
    }

    bb <- .bij_best(cum, lo_idx, hi_idx, k, M)

    if (!full) return(list(misclass = n_pts - bb$max_correct))

    # ---- Margin tie-break among the co-minimal candidates ----
    max_correct <- bb$max_correct
    at          <- bb$at
    mrg_mat <- matrix(cp$margin[pidx[, at, drop = FALSE]], nrow = k - 1L)
    mrg_at  <- if (k == 2L) mrg_mat[1L, ] else apply(mrg_mat, 2L, min)

    # which.max returns the first maximum and `at` is increasing, so the
    # earliest candidate in scan order wins ties, matching the convention in
    # axialLines() and .angular_search().
    w <- which.max(mrg_at)
    j <- at[w]

    return(list(misclass = n_pts - max_correct,
                margin   = mrg_at[w],
                pos      = pos[, j],
                radii    = cp$radius[pidx[, j]]))
}


#' @noRd
.nested_sector <- function(pcoords, cx_vec, cy_vec, radii, k) {
    # Region of each point: the innermost circle containing it, else k.
    sector <- rep(k, nrow(pcoords))
    for (s in (k - 1L):1L) {
        d_s <- sqrt((pcoords[, 1] - cx_vec[s])^2 + (pcoords[, 2] - cy_vec[s])^2)
        sector[d_s <= radii[s]] <- s
    }
    return(sector)
}


#' @noRd
.partition_err <- function(sector, grp_int, k) {
    # Bijection-constrained misclassification of a region assignment -- the
    # exact quantity radialCircle()/radialCircles() report as `misclass`.
    count_mat <- matrix(0L, nrow = k, ncol = k)
    for (r in seq_len(k)) {
        pts_r <- grp_int[sector == r]
        if (length(pts_r) > 0L)
            count_mat[, r] <- tabulate(pts_r, nbins = k)
    }
    return(sum(grp_int != .assign_groups(count_mat)[sector]))
}


#' @noRd
.refine_radii <- function(pcoords, grp_int, cx_vec, cy_vec, radii, k) {
    # Coordinate descent on the radii against the *reported* objective.
    # The sequential fit scores each circle in isolation (groups 1..s inside
    # vs s+1..k outside), which is not the quantity radialCircles() reports,
    # so a radius optimal for the sequential proxy need not be optimal for
    # the final k-region partition. Each sweep rescans one radius over all
    # realisable values -- honouring the nesting constraint on both sides,
    # since circle s must contain circle s-1 and fit inside circle s+1 --
    # and keeps a change only when the total strictly improves. The error is
    # a bounded integer that never increases, so this terminates.
    n_pts <- nrow(pcoords)
    d_mat <- vapply(seq_len(k - 1L),
                    function(s) sqrt((pcoords[, 1] - cx_vec[s])^2 +
                                     (pcoords[, 2] - cy_vec[s])^2),
                    numeric(n_pts))
    cand <- lapply(seq_len(k - 1L),
                   function(s) .cut_radii(sort(d_mat[, s]))$radius)

    # The centers are fixed throughout, so score off the precomputed distances
    # rather than calling .nested_sector(), which would recompute all k - 1
    # distance vectors on every one of the many thousands of evaluations below.
    sector_at <- function(rr) {
        sec <- rep(k, n_pts)
        for (s in (k - 1L):1L) sec[d_mat[, s] <= rr[s]] <- s
        sec
    }

    cur <- .partition_err(sector_at(radii), grp_int, k)

    repeat {
        improved <- FALSE
        for (s in seq_len(k - 1L)) {
            lo <- if (s == 1L) 0 else
                sqrt((cx_vec[s] - cx_vec[s - 1L])^2 +
                     (cy_vec[s] - cy_vec[s - 1L])^2) + radii[s - 1L]
            hi <- if (s == k - 1L) Inf else
                radii[s + 1L] - sqrt((cx_vec[s + 1L] - cx_vec[s])^2 +
                                     (cy_vec[s + 1L] - cy_vec[s])^2)
            for (r_try in cand[[s]][cand[[s]] >= lo & cand[[s]] <= hi]) {
                rr    <- radii
                rr[s] <- r_try
                err   <- .partition_err(sector_at(rr), grp_int, k)
                if (err < cur) {
                    cur      <- err
                    radii    <- rr
                    improved <- TRUE
                }
            }
        }
        if (!improved) break
    }

    return(radii)
}


#' Binary radial partition (one separating circle)
#'
#' Finds the circle minimising misclassification between two groups of 2D
#' points. The center is found by multi-start Nelder-Mead (starts: data
#' centroid plus each group's centroid) unless `cx` and `cy` are supplied.
#' Because the inner search is piecewise-constant, Nelder-Mead can stall on
#' a flat plateau around any of those starts; if none reaches zero
#' misclassification, a coarse `n_grid` x `n_grid` scan of the bounding
#' box, padded by half the data range on each side (the best center is not
#' always inside the convex hull of the points), locates a cell in a
#' better region and one more Nelder-Mead run is seeded there.
#'
#' **Search at fixed center.** Points are sorted by distance from the center
#' and every realisable radius is scanned: the inside/outside split is scored
#' by its exact bijection-constrained misclassification (the better of the
#' two ways to match the two regions to the two groups) — the same criterion
#' used to derive `sector`/`majority` below, so the search targets exactly
#' the quantity reported as `misclass`, and neither group is assumed to be
#' the inner one. Radii that separate coincident distances are skipped as
#' unrealisable. The scan also considers the degenerate radii that leave one
#' region empty; when the configuration has no radial structure these can be
#' the true minimum, so they are allowed, but on a tie a circle that really
#' does split the points is preferred. `radialCircles()` with `k = 2` and the
#' same center returns the same circle.
#'
#' @param crd Numeric matrix or data frame with exactly 2 columns.
#' @param group Factor with exactly 2 levels.
#' @param cx,cy Center of the separating circle; optimised when `NULL` (default).
#' @param fill If `TRUE`, shade the inner disc and outer region (default `FALSE`).
#' @param output If `TRUE` (default), return results list.
#' @param col Circle border colour (default `"purple"`).
#' @param cols Length-2 fill colours (default `c("steelblue", "tomato")`).
#' @param lwd Line width (default `2`).
#' @param lty Line type (default `1`).
#' @param .method `"Nelder-Mead"` (default) or `"SANN"`.
#' @param n_grid Side length of the coarse fallback grid (`n_grid^2` extra
#'   evaluations of the cheap inner search); only used when the heuristic
#'   starts don't already reach zero misclassification. Default `7L`.
#' @param add If `TRUE` (default), add to existing plot; if `FALSE`, call
#'   `plot()` first.
#'
#' @return If `output = TRUE`, a list with `partition` (`"circle"`), `center`,
#'   `radius`, `sector` (`1` inside / `2` outside per point), `majority`
#'   (`character[2]`), `pcoords` (the input `crd`), `pgroup` (`group` coerced
#'   to factor), `misclass` (integer), and `misclass_points` (data frame with
#'   columns `x`, `y`, `label` for each misclassified point).
#'
#' @examples
#' set.seed(1)
#' inner <- cbind(rnorm(20, 0, 0.3), rnorm(20, 0, 0.3))
#' th <- runif(20, 0, 2 * pi)
#' outer <- cbind(2 * cos(th), 2 * sin(th)) + matrix(rnorm(40, 0, 0.1), 20)
#' crd <- rbind(inner, outer)
#' grp <- factor(c(rep("in", 20), rep("out", 20)))
#' radialCircle(crd, grp, fill = TRUE, add = FALSE)
#'
#' @export
radialCircle <- function(crd,
                       group,
                       cx = NULL,
                       cy = NULL,
                       fill = FALSE,
                       output = TRUE,
                       col = "purple",
                       cols = c("steelblue", "tomato"),
                       lwd = 2,
                       lty = 1,
                       .method = "Nelder-Mead",
                       n_grid = 7L,
                       add = TRUE) {

    # ---- Input validation ----
    if (!is.numeric(as.matrix(crd))) stop("Coordinate data must be numeric!")
    if (any(is.na(crd)))            stop("No NAs allowed in crd!")
    if (any(is.na(group)))          stop("No NAs allowed in group!")
    if (length(dim(crd)) != 2)      stop("Coordinates must have two dimensions!")
    if (dim(crd)[2] != 2)           stop("Coordinates must have 2 columns!")
    if (nrow(crd) != length(group)) stop("nrow(crd) must equal length(group)!")
    if (!(.method %in% c("Nelder-Mead", "SANN"))) stop("Method not available!")

    group <- as.factor(group)
    if (nlevels(group) != 2) stop("group must have exactly 2 levels!")

    pcoords  <- as.matrix(crd)
    grp_int <- as.integer(group)

    # ---- Optimise center if not given ----
    # Multi-start: try overall centroid + each group's centroid. The
    # objective (integer misclassification) is a step function, so
    # Nelder-Mead is prone to stalling on plateaus from a single start.
    
    if (is.null(cx) || is.null(cy)) {
        # parscale sizes the simplex to the data range:
        parscale <- c(diff(range(pcoords[, 1])), diff(range(pcoords[, 2])))
        parscale[parscale == 0] <- 1

        # start centers (overall, group means of the 2 groups):
        starts <- list(
            c(mean(pcoords[, 1]),              mean(pcoords[, 2])),
            c(mean(pcoords[grp_int == 1L, 1]), mean(pcoords[grp_int == 1L, 2])),
            c(mean(pcoords[grp_int == 2L, 1]), mean(pcoords[grp_int == 2L, 2]))
        )

        # define function to optimize: fewest misclassified points
        fnToOpt <- function(p) {
            d <- sqrt((pcoords[, 1] - p[1])^2 + (pcoords[, 2] - p[2])^2)
            o <- order(d)
            .radial_search(d[o], grp_int[o], 2L, full = FALSE)$misclass}

        # loop over start values, find best center coordinates cx,cy:
        best <- NULL
        for (s0 in starts) {
            opt <- optim(
                par     = .snap_zero(s0, parscale),
                fn      = fnToOpt,
                method  = .method,
                control = list(reltol = 1e-8, maxit = 2000,
                               parscale = parscale)
            )
            if (is.null(best) || opt$value < best$value) best <- opt
            if (best$value == 0) break          # zero misclass is optimal
        }

        # ---- Coarse grid fallback ----
        # Only pay for it when the starts above didn't already reach the
        # provable optimum; see .grid_seeds() (utils.R) for why two grids
        # (plain bounding box, and one padded by half the data range) are
        # scanned rather than one.
        if (best$value > 0) {
            for (p in .grid_seeds(fnToOpt, range(pcoords[, 1]), range(pcoords[, 2]), n_grid)) {
                opt <- optim(
                    par     = .snap_zero(p, parscale),
                    fn      = fnToOpt,
                    method  = .method,
                    control = list(reltol = 1e-8, maxit = 2000,
                                   parscale = parscale)
                )
                if (opt$value < best$value) best <- opt
                if (best$value == 0) break      # zero misclass is optimal
            }
        }

        cx <- best$par[1]
        cy <- best$par[2]
    }

    # ---- Radius and sectors at chosen center ----
    dists  <- sqrt((pcoords[, 1] - cx)^2 + (pcoords[, 2] - cy)^2)
    ord    <- order(dists)
    res    <- .radial_search(dists[ord], grp_int[ord], 2L)
    radius <- res$radii[1L]
    sector <- ifelse(dists <= radius, 1L, 2L)

    # ---- Unique group assignment ----
    count_mat <- matrix(0L, nrow = 2L, ncol = 2L)
    for (r in 1:2) {
        pts_r <- grp_int[sector == r]
        if (length(pts_r) > 0L)
            count_mat[, r] <- tabulate(pts_r, nbins = 2L)
    }
    assignment <- .assign_groups(count_mat, levels(group))
    majority   <- levels(group)[assignment]

    if (!add) {
        plot(pcoords, asp = 1)
        graphics::text(pcoords, labels = group, cex = 0.7, pos = 4)
    }

    # ---- Optional fill ----
    if (fill) {
        usr    <- par("usr")
        n_c    <- 200L
        inner  <- .circle_pts(cx, cy, radius, n_c)

        draw.circle(cx, cy, radius,
                    col = adjustcolor(cols[1], alpha.f = 0.15), border = NA)

        rect_x <- c(usr[1], usr[2], usr[2], usr[1])
        rect_y <- c(usr[3], usr[3], usr[4], usr[4])
        polypath(
            x    = c(rect_x, NA, inner[, 1]),
            y    = c(rect_y, NA, inner[, 2]),
            rule = "evenodd",
            col  = adjustcolor(cols[2], alpha.f = 0.15), border = NA
        )
    }

    # ---- Draw circle ----
    draw.circle(cx, cy, radius, border = col, lwd = lwd, lty = lty)

    misclass_idx    <- which(as.character(group) != majority[sector])
    misclass_points <- data.frame(
        x     = pcoords[misclass_idx, 1],
        y     = pcoords[misclass_idx, 2],
        label = group[misclass_idx]
    )

    if (!output) return(invisible(NULL))

    list(
        partition       = "circle",
        center          = c(cx, cy),
        radius          = radius,
        sector          = sector,
        majority        = majority,
        pcoords         = crd,
        pgroup          = group,
        misclass        = length(misclass_idx),
        misclass_points = misclass_points
    )
}


#' K-group radial partition (nested circles)
#'
#' Generalises [radialCircle()] to `k >= 2` groups using `k-1` nested
#' (inclusive) circles. The circles must be nested (circle `s` contains
#' circle `s-1`). Circles are fitted sequentially, each minimising
#' misclassification of the innermost `s` groups (inside) vs the rest
#' (outside) — but **which** groups those are is searched, not read off the
#' factor levels, so the inside-to-outside ordering is found from the data
#' and the result does not depend on how the groups are named.
#'
#' When `cx` and `cy` are both supplied, all `k-1` circles share that
#' center (concentric); only the radii are searched, subject to
#' `r_s >= r_{s-1}`. When either is `NULL` (default), each circle's
#' center is optimised independently by multi-start Nelder-Mead — the
#' circles are nested but not generally concentric.
#'
#' **What is minimised.** Every candidate partition is scored by its exact
#' bijection-constrained misclassification: the best total correct over all
#' `k!` ways of matching the `k` regions to the `k` groups, each group used
#' once. This is the same criterion `sector`/`majority` are derived from, so
#' the search targets exactly the quantity reported as `misclass`. Radii that
#' would separate coincident distances are skipped as unrealisable, and radii
#' leaving a region empty are allowed — on data with no radial structure such
#' a partition can be the true minimum — but a partition whose every circle
#' actually splits the points wins any tie.
#'
#' **Optimality.** Concentric (`cx` and `cy` supplied): all `k-1` radii are
#' searched *jointly* over every realisable combination, so the result is
#' globally optimal for that center. If the number of combinations is too
#' large to enumerate, the function warns and falls back to the sequential
#' path below. Optimised centers: circles are fitted sequentially, then all
#' radii are refined by coordinate descent against the criterion above, so
#' the result is optimal in each individual radius but not jointly in the
#' centers. Each circle's own center search is multi-start Nelder-Mead
#' (starts: that circle's inner-group centroid, the overall centroid, and
#' the previous circle's center); as in [radialCircle()], a coarse grid
#' fallback seeds one more run when none of those reaches zero
#' misclassification for that circle.
#'
#' @param crd Numeric matrix or data frame with exactly 2 columns.
#' @param group Factor with `k >= 2` levels. **Factor level order does not
#'   matter**: the inside-to-outside nesting order is found from the data, so
#'   relabelling or reordering the levels cannot change the result. Read the
#'   order off `majority`, which names the group owning each region from the
#'   innermost outwards.
#' @param cx,cy Optional shared center for all `k-1` circles. When both
#'   are supplied the circles are concentric at `(cx, cy)`; otherwise
#'   each center is optimised.
#' @param fill If `TRUE`, shade each ring sector.
#' @param output If `TRUE` (default), return results list.
#' @param col Circle border colour (default `"purple"`).
#' @param cols Length-`k` colour vector; auto-generated if `NULL`.
#' @param lwd Line width (default `2`).
#' @param lty Line type (default `1`).
#' @param .method `"Nelder-Mead"` (default) or `"SANN"`; ignored when
#'   `cx` and `cy` are both supplied.
#' @param n_grid Side length of the coarse fallback grid used when
#'   optimising a circle's center (`n_grid^2` extra evaluations of the
#'   cheap inner search per circle); only used when a circle's heuristic
#'   starts don't already reach zero misclassification, and ignored when
#'   `cx` and `cy` are both supplied. Default `7L`.
#' @param add If `TRUE` (default), add to existing plot; if `FALSE`, call
#'   `plot()` first.
#'
#' @return If `output = TRUE`, a list with `partition` (`"circle"`), `cx`,
#'   `cy` (`numeric[k-1]`), `radii` (innermost to outermost), `sector`
#'   (`integer[n]`), `majority` (`character[k]`), `pcoords` (the input
#'   `crd`), `pgroup` (`group` coerced to factor), `misclass`, and
#'   `misclass_points` (data frame with columns `x`, `y`, `label` for each
#'   misclassified point).
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' g1 <- cbind(rnorm(15, 0, 0.2), rnorm(15, 0, 0.2))
#' th2 <- runif(15, 0, 2 * pi)
#' g2 <- cbind(1.2 * cos(th2), 1.2 * sin(th2)) + matrix(rnorm(30, 0, 0.1), 15)
#' th3 <- runif(15, 0, 2 * pi)
#' g3 <- cbind(2.5 * cos(th3), 2.5 * sin(th3)) + matrix(rnorm(30, 0, 0.1), 15)
#' crd <- rbind(g1, g2, g3)
#' grp <- factor(c(rep("a", 15), rep("b", 15), rep("c", 15)))
#' radialCircles(crd, grp, fill = TRUE, add = FALSE)
#' }
#'
#' @export
radialCircles <- function(crd,
                        group,
                        cx = NULL,
                        cy = NULL,
                        fill = FALSE,
                        output = TRUE,
                        col = "purple",
                        cols = NULL,
                        lwd = 2,
                        lty = 1,
                        .method = "Nelder-Mead",
                        n_grid = 7L,
                        add = TRUE) {

    # ---- Input validation ----
    if (!is.numeric(as.matrix(crd))) stop("Coordinate data must be numeric!")
    if (any(is.na(crd)))            stop("No NAs allowed in crd!")
    if (any(is.na(group)))          stop("No NAs allowed in group!")
    if (length(dim(crd)) != 2)      stop("Coordinates must have two dimensions!")
    if (dim(crd)[2] != 2)           stop("Coordinates must have 2 columns!")
    if (nrow(crd) != length(group)) stop("nrow(crd) must equal length(group)!")
    if (!(.method %in% c("Nelder-Mead", "SANN"))) stop("Method not available!")

    group <- as.factor(group)
    k     <- nlevels(group)
    if (k < 2L) stop("group must have at least 2 levels!")

    pcoords  <- as.matrix(crd)
    grp_int <- as.integer(group)
    n_pts   <- nrow(pcoords)

    if (is.null(cols)) cols <- hcl.colors(k, palette = "Pastel 1")

    # ---- Sequential circle fitting ----
    # If cx/cy are both supplied, all circles share that center
    # (concentric); only radii are searched, subject to r_s >= r_{s-1}.
    # Otherwise each circle's center is optimised independently.
    fixed_center <- !(is.null(cx) || is.null(cy))

    circles   <- vector("list", k - 1L)
    exact_fit <- FALSE

    if (fixed_center) {
        # Concentric branch: shared (cx, cy) for all k-1 circles. All k-1
        # radii are searched *jointly* over every realisable combination of
        # cut positions, scored by the bijection-constrained criterion the
        # function reports -- so this branch is exactly optimal, not greedy.
        d        <- sqrt((pcoords[, 1] - cx)^2 + (pcoords[, 2] - cy)^2)
        ord      <- order(d)
        d_sorted <- d[ord]

        # The joint search builds (k-1) x M integer index matrices, so its
        # cost grows as choose(m + k - 2, k - 1); past a few tens of
        # millions of cells that is gigabytes. Beyond the cap, fall back to
        # the same sequential fit plus radius refinement the
        # independent-centers branch uses -- locally, not globally, optimal,
        # so say so rather than silently returning a weaker answer.
        m_pos <- length(.cut_radii(d_sorted)$pos)
        if (choose(m_pos + k - 2L, k - 1L) * (k - 1L) <= 3e7) {
            res <- .radial_search(d_sorted, grp_int[ord], k)
            for (s in seq_len(k - 1L))
                circles[[s]] <- list(cx = cx, cy = cy, r = res$radii[s])
            exact_fit <- TRUE
        } else {
            warning("too many radius combinations for an exact concentric ",
                    "search (k = ", k, ", n = ", n_pts, "); falling back to ",
                    "the sequential fit with radius refinement, which is ",
                    "not guaranteed to be globally optimal", call. = FALSE)
        }
    }

    if (!exact_fit) {
        # ---- Sequential fit, searched over candidate nesting orders ----
        # A sequential fit needs to know which groups belong inside circle s.
        # Taking that from the factor level order would make the result depend
        # on how the groups happen to be named, so instead every candidate
        # order from .nesting_orders() is fitted and the one whose final
        # geometry minimises the *reported* criterion is kept. Both the
        # candidate set and its iteration order are intrinsic to the
        # configuration, and improvements are strict, so the outcome is
        # invariant under any relabelling or reordering of the levels.
        fit_one <- function(nest_ord) {
            rank_of           <- integer(k)
            rank_of[nest_ord] <- seq_len(k)
            rk                <- rank_of[grp_int]

            out     <- vector("list", k - 1L)
            prev_cx <- NULL; prev_cy <- NULL; prev_r <- NULL

            if (fixed_center) {
                for (s in seq_len(k - 1L)) {
                    r_min    <- if (is.null(prev_r)) 0 else prev_r
                    res      <- .best_radius(d_sorted, (rk <= s)[ord], r_min)
                    out[[s]] <- list(cx = cx, cy = cy, r = res$r)
                    prev_r   <- res$r
                }
            } else {
                overall_ctr <- c(mean(pcoords[, 1]), mean(pcoords[, 2]))
                for (s in seq_len(k - 1L)) {
                    inner_flag <- rk <= s
                    inner_ctr  <- c(mean(pcoords[inner_flag, 1]),
                                    mean(pcoords[inner_flag, 2]))

                    starts <- list(inner_ctr, overall_ctr)
                    if (!is.null(prev_cx))
                        starts <- c(starts, list(c(prev_cx, prev_cy)))

                    circ     <- .optimize_circle(pcoords, inner_flag,
                                                 prev_cx, prev_cy, prev_r,
                                                 starts,
                                                 meth = .method,
                                                 n_grid = n_grid)
                    out[[s]] <- circ
                    prev_cx  <- circ$cx; prev_cy <- circ$cy; prev_r <- circ$r
                }
            }

            out
        }

        cands <- .nesting_orders(pcoords, group)
        if (!isTRUE(attr(cands, "exhaustive")))
            warning("too many group orderings to search exhaustively (k = ", k,
                    "); using the nesting order implied by mean distance from ",
                    "the configuration centroid, which may not be optimal",
                    call. = FALSE)

        # Every candidate is refined before being scored. Refinement is the
        # expensive half, so it is tempting to rank the orderings on their
        # unrefined error and refine only the leaders -- but that measurably
        # loses quality (refinement re-ranks orderings, and no slack band
        # recovers it), so all candidates get the full treatment.
        best <- NULL
        for (nest_ord in cands) {
            cc  <- fit_one(nest_ord)
            cxv <- vapply(cc, function(z) z$cx, numeric(1))
            cyv <- vapply(cc, function(z) z$cy, numeric(1))
            # The sequential fit scores each circle in isolation, which is not
            # the quantity reported below, so refine before scoring.
            rv  <- .refine_radii(pcoords, grp_int, cxv, cyv,
                                 vapply(cc, function(z) z$r, numeric(1)), k)
            err <- .partition_err(
                .nested_sector(pcoords, cxv, cyv, rv, k), grp_int, k)
            if (is.null(best) || err < best$err)
                best <- list(cx = cxv, cy = cyv, r = rv, err = err)
        }

        for (s in seq_len(k - 1L))
            circles[[s]] <- list(cx = best$cx[s], cy = best$cy[s], r = best$r[s])
    }

    cx_vec    <- vapply(circles, function(cc) cc$cx, numeric(1))
    cy_vec    <- vapply(circles, function(cc) cc$cy, numeric(1))
    radii_vec <- vapply(circles, function(cc) cc$r,  numeric(1))

    # ---- Sector assignment ----
    sector <- .nested_sector(pcoords, cx_vec, cy_vec, radii_vec, k)

    # ---- Unique group assignment ----
    count_mat <- matrix(0L, nrow = k, ncol = k)
    for (r in seq_len(k)) {
        pts_r <- grp_int[sector == r]
        if (length(pts_r) > 0L)
            count_mat[, r] <- tabulate(pts_r, nbins = k)
    }
    assignment <- .assign_groups(count_mat, levels(group))
    majority   <- levels(group)[assignment]

    if (!add) {
        plot(pcoords, asp = 1)
        graphics::text(pcoords, labels = group, cex = 0.7, pos = 4)
    }

    # ---- Optional fill ----
    if (fill) {
        usr <- par("usr")
        n_c <- 200L

        outer_circ <- .circle_pts(cx_vec[k - 1L], cy_vec[k - 1L], radii_vec[k - 1L], n_c)
        rect_x     <- c(usr[1], usr[2], usr[2], usr[1])
        rect_y     <- c(usr[3], usr[3], usr[4], usr[4])
        polypath(
            x    = c(rect_x, NA, outer_circ[, 1]),
            y    = c(rect_y, NA, outer_circ[, 2]),
            rule = "evenodd",
            col  = adjustcolor(cols[k], alpha.f = 0.15), border = NA
        )

        if (k >= 3L) {
            for (s in (k - 1L):2L) {
                outer_c <- .circle_pts(cx_vec[s],      cy_vec[s],      radii_vec[s],      n_c)
                inner_c <- .circle_pts(cx_vec[s - 1L], cy_vec[s - 1L], radii_vec[s - 1L], n_c)
                polypath(
                    x    = c(outer_c[, 1], NA, inner_c[, 1]),
                    y    = c(outer_c[, 2], NA, inner_c[, 2]),
                    rule = "evenodd",
                    col  = adjustcolor(cols[s], alpha.f = 0.15), border = NA
                )
            }
        }

        draw.circle(cx_vec[1L], cy_vec[1L], radii_vec[1L],
                    col = adjustcolor(cols[1L], alpha.f = 0.15), border = NA)
    }

    # ---- Draw k-1 circles ----
    for (s in seq_len(k - 1L)) {
        draw.circle(cx_vec[s], cy_vec[s], radii_vec[s], border = col, lwd = lwd, lty = lty)
    }

    misclass_idx    <- which(as.character(group) != majority[sector])
    misclass_points <- data.frame(
        x     = pcoords[misclass_idx, 1],
        y     = pcoords[misclass_idx, 2],
        label = group[misclass_idx]
    )

    if (!output) return(invisible(NULL))

    list(
        partition       = "circle",
        cx              = cx_vec,
        cy              = cy_vec,
        radii           = radii_vec,
        sector          = sector,
        majority        = majority,
        pcoords         = crd,
        pgroup          = group,
        misclass        = length(misclass_idx),
        misclass_points = misclass_points
    )
}
