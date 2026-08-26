# ===================================================================
# ==== Functions for Radial Elliptic Partitions ====
# ===================================================================


# ---- Internal helpers ----

#' @noRd
.ellipse_pts <- function(cx, cy, a, b, angle_rad, n = 200L) {
    # compute n coordinates for ellipse cx,cy,a,b,angle_rad
    theta   <- seq(0, 2 * pi, length.out = n + 1L)[-1L]
    cos_a   <- cos(angle_rad)
    sin_a   <- sin(angle_rad)
    x_local <- a * cos(theta)
    y_local <- b * sin(theta)
    cbind(cx + cos_a * x_local - sin_a * y_local,
          cy + sin_a * x_local + cos_a * y_local)
}


#' @noRd
.in_ellipse <- function(coords, cx, cy, a, b, angle_rad) {
    # check if coords are inside ellipse: TRUE/FALSE
    cos_a <- cos(angle_rad)
    sin_a <- sin(angle_rad)
    dx    <- coords[, 1] - cx
    dy    <- coords[, 2] - cy
    u     <-  cos_a * dx + sin_a * dy
    v     <- -sin_a * dx + cos_a * dy
    (u / a)^2 + (v / b)^2 <= 1
}


#' @noRd
.eval_ellipse <- function(params, coords, 
                          inner_flag, prev_bnd = NULL,
                          penalty = 1e6) {
    cx    <- params[1]
    cy    <- params[2]
    a     <- abs(params[3])
    b     <- abs(params[4])
    angle <- params[5]

    if (a < 1e-6 || b < 1e-6) return(penalty)

    cos_a <- cos(angle)
    sin_a <- sin(angle)

    if (!is.null(prev_bnd)) {
        dx <- prev_bnd[, 1] - cx;  dy <- prev_bnd[, 2] - cy
        u  <-  cos_a * dx + sin_a * dy
        v  <- -sin_a * dx + cos_a * dy
        if (any((u / a)^2 + (v / b)^2 > 1)) return(penalty)
    }

    dx     <- coords[, 1] - cx;  dy <- coords[, 2] - cy
    u      <-  cos_a * dx + sin_a * dy
    v      <- -sin_a * dx + cos_a * dy
    inside <- (u / a)^2 + (v / b)^2 <= 1
    sum(inner_flag & !inside) + sum(!inner_flag & inside)
}


#' @noRd
.init_ellipse_params <- function(coords, inner_flag) {
    inner_pts <- coords[inner_flag, , drop = FALSE]
    cx        <- mean(inner_pts[, 1])
    cy        <- mean(inner_pts[, 2])

    if (nrow(inner_pts) < 3L) {
        return(c(cx, cy, 0.1, 0.1, 0))
    }

    S     <- cov(inner_pts)
    S_reg <- S + diag(1e-10, 2)
    eig   <- eigen(S_reg)
    angle <- atan2(eig$vectors[2, 1], eig$vectors[1, 1])

    S_inv    <- solve(S_reg)
    centered <- sweep(inner_pts, 2, c(cx, cy))
    mah_dist <- sqrt(rowSums((centered %*% S_inv) * centered))

    k <- max(mah_dist) * 1.05 + 1e-6
    a <- sqrt(abs(eig$values[1])) * k
    b <- sqrt(abs(eig$values[2])) * k

    c(cx, cy, a, b, angle)
}


#' @noRd
.optimize_ellipse <- function(coords, inner_flag, 
                              prev = NULL, starts,
                              n_grid = 7L) {
    x_range   <- diff(range(coords[, 1]))
    y_range   <- diff(range(coords[, 2]))
    if (x_range == 0) x_range <- 1
    if (y_range == 0) y_range <- 1
    max_range <- max(x_range, y_range)
    parscale  <- c(x_range, y_range, max_range, max_range, pi)

    prev_bnd <- if (!is.null(prev))
        .ellipse_pts(prev$cx, prev$cy, prev$a, prev$b, prev$angle, 
                     n = 200L)
    else
        NULL

    # .eval_ellipse()'s misclass is piecewise-constant, so Nelder-Mead can
    # stall on a flat plateau around any of `starts` -- including a start
    # whose (cx, cy) happens to be near zero (e.g. a group centered close
    # to the configuration centroid of mean-centered MDS data), which also
    # hits optim()'s near-zero initial-simplex degeneracy (see
    # .snap_zero(), utils.R). Confirmed on real MDS data: one group's own
    # covariance-based init landed inside such a plateau and stalled
    # several misclassifications above the true optimum.
    best <- NULL
    for (init in starts) {
        init[1:2] <- .snap_zero(init[1:2], parscale[1:2])
        opt <- optim(
            par        = init,
            fn         = .eval_ellipse,
            coords     = coords,
            inner_flag = inner_flag,
            prev_bnd   = prev_bnd,
            method     = "Nelder-Mead",
            control    = list(reltol = 1e-8, maxit = 5000,
                              parscale = parscale)
        )
        if (is.null(best) || opt$value < best$value) best <- opt
        if (best$value == 0) break          # zero misclass is optimal
    }

    # ---- Coarse grid fallback (center only) ----
    # A full 5D grid is infeasible, so this relocates just (cx, cy):
    # holding (a, b, angle) at whatever the starts above already found,
    # scan for a better center, then re-optimise all 5 parameters together
    # from there. Only pay for it when the starts didn't already reach the
    # provable optimum; see .grid_seeds() (utils.R) for why two grids
    # (plain bounding box, and one padded by half the data range) are
    # scanned rather than one.
    if (best$value > 0) {
        shape    <- best$par[3:5]
        eval_ctr <- function(ctr) {
            .eval_ellipse(c(ctr[1], ctr[2], shape),
                         coords, inner_flag, prev_bnd)
        }
        for (p in .grid_seeds(eval_ctr, range(coords[, 1]), range(coords[, 2]), n_grid)) {
            opt <- optim(
                par        = c(.snap_zero(p, parscale[1:2]), shape),
                fn         = .eval_ellipse,
                coords     = coords,
                inner_flag = inner_flag,
                prev_bnd   = prev_bnd,
                method     = "Nelder-Mead",
                control    = list(reltol = 1e-8, maxit = 5000,
                                  parscale = parscale)
            )
            if (opt$value < best$value) best <- opt
            if (best$value == 0) break      # zero misclass is optimal
        }
    }

    p <- best$par
    list(cx       = p[1],
         cy       = p[2],
         a        = abs(p[3]),
         b        = abs(p[4]),
         angle    = p[5],
         misclass = as.integer(round(best$value)))
}


#' @noRd
.elliptic_cuts_2 <- function(coords, grp_int, n_grid = 7L) {
    best_res <- NULL
    best_mc  <- .Machine$integer.max

    for (which_inner in 1L:2L) {
        inner_flag <- grp_int == which_inner
        init       <- .init_ellipse_params(coords, inner_flag)
        res        <- .optimize_ellipse(coords, inner_flag,
                                        prev = NULL, starts = list(init),
                                        n_grid = n_grid)
        if (res$misclass < best_mc) {
            best_mc  <- res$misclass
            best_res <- res
        }
    }

    inside <- .in_ellipse(coords, best_res$cx, best_res$cy,
                          best_res$a, best_res$b, best_res$angle)
    list(ellipse = best_res, sector = ifelse(inside, 1L, 2L))
}


#' Binary radial-elliptic partition (one separating ellipse)
#'
#' Finds the ellipse minimising misclassification between two groups of 2D
#' points. All 5 parameters (cx, cy, a, b, angle) are found by Nelder-Mead,
#' starting from the covariance ellipse of each group in turn; the ordering
#' with the lower misclassification is retained. Because the inner search is
#' piecewise-constant, Nelder-Mead can stall on a flat plateau around either
#' start (confirmed on real MDS data: a group's own covariance-based init
#' can itself be near enough to the configuration centroid to land in one);
#' if a group's start doesn't reach zero misclassification, a coarse
#' `n_grid` x `n_grid` scan relocates just the center (holding the fitted
#' shape and orientation fixed) over the bounding box, padded by half the
#' data range on each side, and one more full 5-parameter Nelder-Mead run is
#' seeded from the best cell found.
#'
#' @param crd Numeric matrix or data frame with exactly 2 columns.
#' @param group Factor with exactly 2 levels.
#' @param n_grid Side length of the coarse fallback grid (`n_grid^2` extra
#'   evaluations of the cheap inner search per group tried as inner);
#'   only used when a group's covariance-based start doesn't already reach
#'   zero misclassification. Default `7L`.
#' @param fill If `TRUE`, shade inner ellipse and outer region.
#' @param output If `TRUE` (default), return results list.
#' @param col Ellipse border colour (default `"purple"`).
#' @param cols Length-2 fill colours (default `c("steelblue", "tomato")`).
#' @param lwd Line width (default `2`).
#' @param lty Line type (default `1`).
#' @param add If `TRUE` (default), add to existing plot; if `FALSE`, call
#'   `plot()` first.
#'
#' @return If `output = TRUE`, a list with `partition` (`"ellipse"`), `cx`,
#'   `cy`, `a`, `b`, `angle` (radians), `sector`, `majority`, `pcoords` (the
#'   input `crd`), `pgroup` (`group` coerced to factor), `misclass` (list
#'   with `n` and `indices`), and `misclass_points` (data frame with columns
#'   `x`, `y`, `label` for each misclassified point).
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' inner <- cbind(rnorm(20, 0, 0.3), rnorm(20, 0, 0.2))
#' th    <- runif(20, 0, 2 * pi)
#' outer <- cbind(2.0 * cos(th), 1.2 * sin(th)) + matrix(rnorm(40, 0, 0.1), 20)
#' crd <- rbind(inner, outer)
#' grp <- factor(c(rep("in", 20), rep("out", 20)))
#' radialEllipse(crd, grp, fill = TRUE, add = FALSE)
#' }
#'
#' @export
radialEllipse <- function(crd,
                           group,
                           n_grid = 7L,
                           fill = FALSE,
                           output = TRUE,
                           col = "purple",
                           cols = c("steelblue", "tomato"),
                           lwd = 2,
                           lty = 1,
                           add = TRUE) {

    # ---- Input validation ----
    if (!is.numeric(as.matrix(crd))) stop("Coordinate data must be numeric!")
    if (any(is.na(crd)))            stop("No NAs allowed in crd!")
    if (any(is.na(group)))          stop("No NAs allowed in group!")
    if (length(dim(crd)) != 2)      stop("Coordinates must have two dimensions!")
    if (dim(crd)[2] != 2)           stop("Coordinates must have 2 columns!")
    if (nrow(crd) != length(group)) stop("nrow(crd) must equal length(group)!")

    group <- as.factor(group)
    if (nlevels(group) != 2) stop("group must have exactly 2 levels!")

    coords  <- as.matrix(crd)
    grp_int <- as.integer(group)

    # ---- Optimise ellipse ----
    res    <- .elliptic_cuts_2(coords, grp_int, n_grid = n_grid)
    ell    <- res$ellipse
    sector <- res$sector

    # ---- Majority labels and misclassification ----
    count_mat <- matrix(0L, nrow = 2L, ncol = 2L)
    for (r in 1L:2L) {
        pts_r <- grp_int[sector == r]
        if (length(pts_r) > 0L)
            count_mat[, r] <- tabulate(pts_r, nbins = 2L)
    }
    assignment      <- .assign_groups(count_mat, levels(group))
    majority        <- levels(group)[assignment]
    misclass_idx    <- which(levels(group)[grp_int] != majority[sector])
    misclass        <- list(n = length(misclass_idx), indices = misclass_idx)
    misclass_points <- data.frame(
        x     = coords[misclass_idx, 1],
        y     = coords[misclass_idx, 2],
        label = group[misclass_idx]
    )

    if (!add) {
        plot(coords, asp = 1)
        graphics::text(coords, labels = group, cex = 0.7, pos = 4)
    }

    # ---- Optional fill ----
    if (fill) {
        usr <- par("usr")
        bnd <- .ellipse_pts(ell$cx, ell$cy, ell$a, ell$b, ell$angle)

        draw.ellipse(ell$cx, ell$cy, a = ell$a, b = ell$b,
                     angle = ell$angle * 180 / pi,
                     col = adjustcolor(cols[1], alpha.f = 0.15), border = NA)

        rect_x <- c(usr[1], usr[2], usr[2], usr[1])
        rect_y <- c(usr[3], usr[3], usr[4], usr[4])
        polypath(
            x    = c(rect_x, NA, bnd[, 1]),
            y    = c(rect_y, NA, bnd[, 2]),
            rule = "evenodd",
            col  = adjustcolor(cols[2], alpha.f = 0.15), border = NA
        )
    }

    # ---- Draw ellipse ----
    draw.ellipse(ell$cx, ell$cy, a = ell$a, b = ell$b,
                 angle = ell$angle * 180 / pi,
                 border = col, lwd = lwd, lty = lty)

    if (!output) return(invisible(NULL))

    list(
        partition       = "ellipse",
        cx              = ell$cx,
        cy              = ell$cy,
        a               = ell$a,
        b               = ell$b,
        angle           = ell$angle,
        sector          = sector,
        majority        = majority,
        pcoords         = crd,
        pgroup          = group,
        misclass        = misclass,
        misclass_points = misclass_points
    )
}


#' K-group radial-elliptic partition (nested ellipses)
#'
#' Generalises [radialEllipse()] to `k >= 2` groups using `k-1` nested
#' ellipses. The nesting constraint (ellipse `s-1` lies inside ellipse
#' `s`) is enforced during optimisation. The inside-to-outside ordering of
#' the groups is searched, not read off the factor levels, so the result
#' does not depend on how the groups are named.
#'
#' When `ellipse` is supplied as a length-5 vector `(cx, cy, a, b, angle)`,
#' it defines the innermost ellipse exactly. All `k-1` ellipses then share
#' that center, orientation, and a:b ratio — outer ellipses are uniform
#' scalings of the inner one. The scale factors `t_2, ..., t_{k-1}` (with
#' `t_s > t_{s-1} >= 1`) are found by exact 1D scan over each point's
#' critical scale `sqrt((u_i/a)^2 + (v_i/b)^2)` in the ellipse's rotated
#' frame.
#'
#' When `ellipse` is `NULL` (default), each ellipse has its own
#' independently optimised center, semi-axes, and rotation, fitted by
#' Nelder-Mead. Multi-start uses the covariance ellipse of the inner
#' groups plus an inflated copy of the previous ellipse (a guaranteed-
#' feasible starting point that avoids stalling on the infeasibility
#' penalty plateau). As in [radialEllipse()], the inner search is
#' piecewise-constant and can stall on a flat plateau around either start;
#' when that happens for a given ellipse, a coarse `n_grid` x `n_grid`
#' fallback relocates its center before one more full refit (ignored in
#' the fixed-`ellipse` mode, which has no Nelder-Mead step).
#'
#' @param crd Numeric matrix or data frame with exactly 2 columns.
#' @param group Factor with `k >= 2` levels. **Factor level order does not
#'   matter**: the inside-to-outside nesting order is found from the data, so
#'   relabelling or reordering the levels cannot change the result. Read the
#'   order off `majority`, which names the group owning each region from the
#'   innermost outwards.
#' @param ellipse Optional length-5 numeric vector `(cx, cy, a, b, angle)`
#'   specifying the innermost ellipse exactly. When supplied, all outer
#'   ellipses are uniform scalings of this one; only the scale factors
#'   are searched.
#' @param n_grid Side length of the coarse fallback grid used when
#'   optimising an ellipse's center (`n_grid^2` extra evaluations of the
#'   cheap inner search per ellipse); only used when an ellipse's
#'   heuristic starts don't already reach zero misclassification, and
#'   ignored when `ellipse` is supplied. Default `7L`.
#' @param fill If `TRUE`, shade each elliptic sector.
#' @param output If `TRUE` (default), return results list.
#' @param col Ellipse border colour (default `"purple"`).
#' @param cols Length-`k` colours; auto-generated if `NULL`.
#' @param lwd Line width (default `2`).
#' @param lty Line type (default `1`).
#' @param add If `TRUE` (default), add to existing plot; if `FALSE`, call
#'   `plot()` first.
#'
#' @return If `output = TRUE`, a list with `partition` (`"ellipse"`), vectors
#'   `cx`, `cy`, `a`, `b`, `angle` (radians), `sector`, `majority`, `pcoords`
#'   (the input `crd`), `pgroup` (`group` coerced to factor), `misclass`
#'   (list with `n` and `indices`), and `misclass_points` (data frame with
#'   columns `x`, `y`, `label` for each misclassified point).
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' g1 <- cbind(rnorm(15, 0, 0.2), rnorm(15, 0, 0.15))
#' th2 <- runif(15, 0, 2 * pi)
#' g2 <- cbind(1.4 * cos(th2), 1.0 * sin(th2)) + matrix(rnorm(30, 0, 0.1), 15)
#' th3 <- runif(15, 0, 2 * pi)
#' g3 <- cbind(2.8 * cos(th3), 2.0 * sin(th3)) + matrix(rnorm(30, 0, 0.1), 15)
#' crd <- rbind(g1, g2, g3)
#' grp <- factor(c(rep("a", 15), rep("b", 15), rep("c", 15)))
#'
#' # Independent ellipses
#' radialEllipses(crd, grp, fill = TRUE, add = FALSE)
#'
#' # Fixed-shape mode: supply the innermost ellipse as a 5-vector
#' radialEllipses(crd, grp, ellipse = c(0, 0, 0.3, 0.2, 0), fill = TRUE, add = FALSE)
#' }
#'
#' @export
radialEllipses <- function(crd,
                            group,
                            ellipse = NULL,
                            n_grid = 7L,
                            fill = FALSE,
                            output = TRUE,
                            col = "purple",
                            cols = NULL,
                            lwd = 2,
                            lty = 1,
                            add = TRUE) {

    # ---- Input validation ----
    if (!is.numeric(as.matrix(crd))) stop("Coordinate data must be numeric!")
    if (any(is.na(crd)))            stop("No NAs allowed in crd!")
    if (any(is.na(group)))          stop("No NAs allowed in group!")
    if (length(dim(crd)) != 2)      stop("Coordinates must have two dimensions!")
    if (dim(crd)[2] != 2)           stop("Coordinates must have 2 columns!")
    if (nrow(crd) != length(group)) stop("nrow(crd) must equal length(group)!")
    if (!is.null(ellipse)) {
        if (!is.numeric(ellipse) || length(ellipse) != 5L) {
            stop("ellipse must be a length-5 numeric vector (cx, cy, a, b, angle)!")
        }
        if (abs(ellipse[3]) < 1e-6 || abs(ellipse[4]) < 1e-6) {
            stop("ellipse semi-axes (a, b) must be positive!")
        }
    }

    group <- as.factor(group)
    k     <- nlevels(group)
    if (k < 2L) stop("group must have at least 2 levels!")

    coords  <- as.matrix(crd)
    grp_int <- as.integer(group)
    n_pts   <- nrow(coords)

    if (is.null(cols)) cols <- hcl.colors(k, palette = "Pastel 1")

    # ---- Sequential ellipse fitting ----
    # If `ellipse` is supplied, the k-1 ellipses are uniform scalings of
    # this template: innermost = supplied ellipse (scale t_1 = 1), outer
    # scales t_s > t_{s-1} found by exact 1D scan over the critical scales
    # c_i = sqrt((u_i/a)^2 + (v_i/b)^2). Otherwise each ellipse has its
    # own independently optimised center, semi-axes, and rotation.
    fixed_shape <- !is.null(ellipse)

    if (fixed_shape) {
        cx_fix    <- ellipse[1]
        cy_fix    <- ellipse[2]
        a_fix     <- abs(ellipse[3])
        b_fix     <- abs(ellipse[4])
        angle_fix <- ellipse[5]

        # Project points to the (u, v) frame of the fixed ellipse and
        # compute each point's critical scale c_i. Point i lies inside
        # the scaled ellipse t*ellipse iff t >= c_i.
        cos_a    <- cos(angle_fix)
        sin_a    <- sin(angle_fix)
        dx       <- coords[, 1] - cx_fix
        dy       <- coords[, 2] - cy_fix
        u        <- cos_a * dx + sin_a * dy
        v        <- -sin_a * dx + cos_a * dy
        crit_t   <- sqrt((u / a_fix)^2 + (v / b_fix)^2)
        ord      <- order(crit_t)
        t_sorted <- crit_t[ord]
    }

    # Fit all k-1 ellipses for one candidate nesting order (innermost group
    # first). Both paths need to know which groups belong inside ellipse s;
    # taking that from the factor level order would make the result depend on
    # how the groups happen to be named, so the order is searched below
    # instead.
    fit_one <- function(nest_ord) {
        rank_of           <- integer(k)
        rank_of[nest_ord] <- seq_len(k)
        rk                <- rank_of[grp_int]

        out <- vector("list", k - 1L)

        if (fixed_shape) {
            # Innermost ellipse: exactly as supplied (t_1 = 1)
            out[[1]] <- list(cx = cx_fix, cy = cy_fix,
                             a = a_fix, b = b_fix, angle = angle_fix)

            # Outer ellipses: search t_s subject to t_s > t_{s-1}, reusing
            # .best_radius() with the critical scales playing the role of
            # distances and t playing the role of radius.
            if (k >= 3L) {
                prev_t <- 1
                for (s in 2L:(k - 1L)) {
                    res <- .best_radius(t_sorted, (rk <= s)[ord],
                                        r_min = prev_t)
                    t_s <- res$r
                    out[[s]] <- list(cx = cx_fix, cy = cy_fix,
                                     a = a_fix * t_s,
                                     b = b_fix * t_s,
                                     angle = angle_fix)
                    prev_t <- t_s
                }
            }
        } else {
            # Independent-ellipse branch (original behaviour).
            prev <- NULL
            for (s in seq_len(k - 1L)) {
                inner_flag <- rk <= s
                starts     <- list(.init_ellipse_params(coords, inner_flag))
                if (!is.null(prev)) {
                    prev_inflated <- c(prev$cx, prev$cy,
                                       prev$a * 1.2, prev$b * 1.2,
                                       prev$angle)
                    starts <- c(starts, list(prev_inflated))
                }
                ell      <- .optimize_ellipse(coords, inner_flag,
                                              prev = prev, starts = starts,
                                              n_grid = n_grid)
                out[[s]] <- ell
                prev     <- ell
            }
        }

        out
    }

    # Region of each point: the innermost ellipse containing it, else k.
    sector_of <- function(ell) {
        sec <- rep(k, n_pts)
        if (fixed_shape) {
            t_vec <- vapply(ell, function(e) e$a, numeric(1)) / a_fix
            for (s in (k - 1L):1L) sec[crit_t <= t_vec[s]] <- s
        } else {
            for (s in (k - 1L):1L) {
                inside      <- .in_ellipse(coords, ell[[s]]$cx, ell[[s]]$cy,
                                           ell[[s]]$a, ell[[s]]$b,
                                           ell[[s]]$angle)
                sec[inside] <- s
            }
        }
        sec
    }

    # ---- Search over candidate nesting orders ----
    # Both the candidate set and its iteration order are intrinsic to the
    # configuration, and improvements are strict, so the result is invariant
    # under any relabelling or reordering of the factor levels.
    cands <- .nesting_orders(coords, group)
    if (!isTRUE(attr(cands, "exhaustive")))
        warning("too many group orderings to search exhaustively (k = ", k,
                "); using the nesting order implied by mean distance from the ",
                "configuration centroid, which may not be optimal",
                call. = FALSE)

    best <- NULL
    for (nest_ord in cands) {
        ell <- fit_one(nest_ord)
        sec <- sector_of(ell)
        err <- .partition_err(sec, grp_int, k)
        if (is.null(best) || err < best$err)
            best <- list(ell = ell, sector = sec, err = err)
    }

    ellipses <- best$ell
    sector   <- best$sector

    cx_vec    <- vapply(ellipses, function(e) e$cx,    numeric(1))
    cy_vec    <- vapply(ellipses, function(e) e$cy,    numeric(1))
    a_vec     <- vapply(ellipses, function(e) e$a,     numeric(1))
    b_vec     <- vapply(ellipses, function(e) e$b,     numeric(1))
    angle_vec <- vapply(ellipses, function(e) e$angle, numeric(1))

    # ---- Majority labels and misclassification ----
    count_mat <- matrix(0L, nrow = k, ncol = k)
    for (r in seq_len(k)) {
        pts_r <- grp_int[sector == r]
        if (length(pts_r) > 0L)
            count_mat[, r] <- tabulate(pts_r, nbins = k)
    }
    assignment      <- .assign_groups(count_mat, levels(group))
    majority        <- levels(group)[assignment]
    misclass_idx    <- which(levels(group)[grp_int] != majority[sector])
    misclass        <- list(n = length(misclass_idx), 
                            indices = misclass_idx)
    misclass_points <- data.frame(
        x     = coords[misclass_idx, 1],
        y     = coords[misclass_idx, 2],
        label = group[misclass_idx]
    )

    if (!add) {
        plot(coords, las = 1, asp = 1)
        graphics::text(coords, labels = group, cex = 0.7, pos = 4)
    }

    # ---- Optional fill ----
    if (fill) {
        usr <- par("usr")

        outer_bnd <- .ellipse_pts(cx_vec[k - 1L], cy_vec[k - 1L],
                                   a_vec[k - 1L], b_vec[k - 1L],
                                   angle_vec[k - 1L])
        rect_x <- c(usr[1], usr[2], usr[2], usr[1])
        rect_y <- c(usr[3], usr[3], usr[4], usr[4])
        polypath(
            x    = c(rect_x, NA, outer_bnd[, 1]),
            y    = c(rect_y, NA, outer_bnd[, 2]),
            rule = "evenodd",
            col  = adjustcolor(cols[k], alpha.f = 0.15), 
            border = NA
        )

        if (k >= 3L) {
            for (s in (k - 1L):2L) {
                outer_bnd <- .ellipse_pts(cx_vec[s], cy_vec[s],
                                          a_vec[s], b_vec[s], 
                                          angle_vec[s])
                inner_bnd <- .ellipse_pts(cx_vec[s - 1L], cy_vec[s - 1L],
                                           a_vec[s - 1L], b_vec[s - 1L],
                                           angle_vec[s - 1L])
                polypath(
                    x    = c(outer_bnd[, 1], NA, inner_bnd[, 1]),
                    y    = c(outer_bnd[, 2], NA, inner_bnd[, 2]),
                    rule = "evenodd",
                    col  = adjustcolor(cols[s], alpha.f = 0.15), 
                    border = NA
                )
            }
        }

        draw.ellipse(cx_vec[1L], cy_vec[1L], 
                     a = a_vec[1L], b = b_vec[1L],
                     angle = angle_vec[1L] * 180 / pi,
                     col = adjustcolor(cols[1L], alpha.f = 0.15), 
                     border = NA)
    }

    # ---- Draw k-1 ellipses ----
    for (s in seq_len(k - 1L)) {
        draw.ellipse(cx_vec[s], cy_vec[s], 
                     a = a_vec[s], b = b_vec[s],
                     angle = angle_vec[s] * 180 / pi,
                     border = col, lwd = lwd, lty = lty)
    }

    if (!output) return(invisible(NULL))

    list(
        partition       = "ellipse",
        cx              = cx_vec,
        cy              = cy_vec,
        a               = a_vec,
        b               = b_vec,
        angle           = angle_vec,
        sector          = sector,
        majority        = majority,
        pcoords         = crd,
        pgroup          = group,
        misclass        = misclass,
        misclass_points = misclass_points
    )
}
