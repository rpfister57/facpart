# ===================================================================
# ==== Functions for Radial (Elliptic) Partitions ====
# ===================================================================


# ---- Internal helpers ----

#' @noRd
.ellipse_pts <- function(cx, cy, a, b, angle_rad, n = 200L) {
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
    cos_a <- cos(angle_rad)
    sin_a <- sin(angle_rad)
    dx    <- coords[, 1] - cx
    dy    <- coords[, 2] - cy
    u     <-  cos_a * dx + sin_a * dy
    v     <- -sin_a * dx + cos_a * dy
    (u / a)^2 + (v / b)^2 <= 1
}


#' @noRd
.eval_ellipse <- function(params, coords, inner_flag, prev = NULL,
                          penalty = 1e6) {
    cx    <- params[1]
    cy    <- params[2]
    a     <- abs(params[3])
    b     <- abs(params[4])
    angle <- params[5]

    if (a < 1e-6 || b < 1e-6) return(penalty)

    if (!is.null(prev)) {
        bnd <- .ellipse_pts(prev$cx, prev$cy, prev$a, prev$b, prev$angle,
                            n = 200L)
        if (!all(.in_ellipse(bnd, cx, cy, a, b, angle))) return(penalty)
    }

    inside <- .in_ellipse(coords, cx, cy, a, b, angle)
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
.optimize_ellipse <- function(coords, inner_flag, prev = NULL, starts) {
    x_range   <- diff(range(coords[, 1]))
    y_range   <- diff(range(coords[, 2]))
    if (x_range == 0) x_range <- 1
    if (y_range == 0) y_range <- 1
    max_range <- max(x_range, y_range)
    parscale  <- c(x_range, y_range, max_range, max_range, pi)

    best <- NULL
    for (init in starts) {
        opt <- optim(
            par        = init,
            fn         = .eval_ellipse,
            coords     = coords,
            inner_flag = inner_flag,
            prev       = prev,
            method     = "Nelder-Mead",
            control    = list(reltol = 1e-8, maxit = 5000,
                              parscale = parscale)
        )
        if (is.null(best) || opt$value < best$value) best <- opt
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
.elliptic_cuts_2 <- function(coords, grp_int) {
    best_res <- NULL
    best_mc  <- .Machine$integer.max

    for (which_inner in 1L:2L) {
        inner_flag <- grp_int == which_inner
        init       <- .init_ellipse_params(coords, inner_flag)
        res        <- .optimize_ellipse(coords, inner_flag,
                                        prev = NULL, starts = list(init))
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
#' with the lower misclassification is retained.
#'
#' @param crd Numeric matrix or data frame with exactly 2 columns.
#' @param group Factor with exactly 2 levels.
#' @param fill If `TRUE`, shade inner ellipse and outer region.
#' @param output If `TRUE` (default), return results list.
#' @param col Ellipse border colour (default `"purple"`).
#' @param cols Length-2 fill colours (default `c("steelblue", "tomato")`).
#' @param lwd Line width (default `2`).
#' @param lty Line type (default `1`).
#'
#' @return If `output = TRUE`, a list with `cx`, `cy`, `a`, `b`, `angle`
#'   (radians), `misclass` (list with `n` and `indices`), `sector`, and
#'   `majority`.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' inner <- cbind(rnorm(20, 0, 0.3), rnorm(20, 0, 0.2))
#' th    <- runif(20, 0, 2 * pi)
#' outer <- cbind(2.0 * cos(th), 1.2 * sin(th)) + matrix(rnorm(40, 0, 0.1), 20)
#' crd <- rbind(inner, outer)
#' grp <- factor(c(rep("in", 20), rep("out", 20)))
#' plot(crd, asp = 1)
#' radialEllipse(crd, grp, fill = TRUE)
#' }
#'
#' @export
radialEllipse <- function(crd,
                           group,
                           fill = FALSE,
                           output = TRUE,
                           col = "purple",
                           cols = c("steelblue", "tomato"),
                           lwd = 2,
                           lty = 1) {

    # ---- Input validation ----
    if (!is.numeric(crd))           stop("Coordinate data must be numeric!")
    if (length(dim(crd)) != 2)      stop("Coordinates must have two dimensions!")
    if (dim(crd)[2] != 2)           stop("Coordinates must be 2-dimensional!")
    if (nrow(crd) != length(group)) stop("nrow(crd) must equal length(group)!")

    group <- as.factor(group)
    if (nlevels(group) != 2) stop("group must have exactly 2 levels!")

    coords  <- as.matrix(crd)
    grp_int <- as.integer(group)

    # ---- Optimise ellipse ----
    res    <- .elliptic_cuts_2(coords, grp_int)
    ell    <- res$ellipse
    sector <- res$sector

    # ---- Majority labels and misclassification ----
    majority <- character(2L)
    for (s in 1L:2L) {
        pts_s       <- grp_int[sector == s]
        majority[s] <- levels(group)[which.max(tabulate(pts_s, nbins = 2L))]
    }
    misclass_idx <- which(levels(group)[grp_int] != majority[sector])
    misclass     <- list(n = length(misclass_idx), indices = misclass_idx)

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
        cx       = ell$cx,
        cy       = ell$cy,
        a        = ell$a,
        b        = ell$b,
        angle    = ell$angle,
        misclass = misclass,
        sector   = sector,
        majority = majority
    )
}


#' K-group radial-elliptic partition (nested ellipses)
#'
#' Generalises [radialEllipse()] to `k >= 2` groups using `k-1` nested
#' ellipses. The nesting constraint (ellipse `s-1` lies inside ellipse
#' `s`) is enforced during optimisation. Groups are ordered by factor
#' level (level 1 = innermost).
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
#' penalty plateau).
#'
#' @param crd Numeric matrix or data frame with exactly 2 columns.
#' @param group Factor with `k >= 2` levels.
#' @param ellipse Optional length-5 numeric vector `(cx, cy, a, b, angle)`
#'   specifying the innermost ellipse exactly. When supplied, all outer
#'   ellipses are uniform scalings of this one; only the scale factors
#'   are searched.
#' @param fill If `TRUE`, shade each elliptic sector.
#' @param output If `TRUE` (default), return results list.
#' @param col Ellipse border colour (default `"purple"`).
#' @param cols Length-`k` colours; auto-generated if `NULL`.
#' @param lwd Line width (default `2`).
#' @param lty Line type (default `1`).
#'
#' @return If `output = TRUE`, a list with vectors `cx`, `cy`, `a`, `b`,
#'   `angle` (radians), `misclass` (list with `n` and `indices`),
#'   `sector`, and `majority`.
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
#' plot(crd, asp = 1)
#'
#' # Independent ellipses
#' radialEllipses(crd, grp, fill = TRUE)
#'
#' # Fixed-shape mode: supply the innermost ellipse as a 5-vector
#' radialEllipses(crd, grp, ellipse = c(0, 0, 0.3, 0.2, 0), fill = TRUE)
#' }
#'
#' @export
radialEllipses <- function(crd,
                            group,
                            ellipse = NULL,
                            fill = FALSE,
                            output = TRUE,
                            col = "purple",
                            cols = NULL,
                            lwd = 2,
                            lty = 1) {

    # ---- Input validation ----
    if (!is.numeric(crd))           stop("Coordinate data must be numeric!")
    if (length(dim(crd)) != 2)      stop("Coordinates must have two dimensions!")
    if (dim(crd)[2] != 2)           stop("Coordinates must be 2-dimensional!")
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

    ellipses <- vector("list", k - 1L)

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

        # Innermost ellipse: exactly as supplied (t_1 = 1)
        ellipses[[1]] <- list(cx = cx_fix, cy = cy_fix,
                              a = a_fix, b = b_fix, angle = angle_fix)

        # Outer ellipses: search t_s subject to t_s > t_{s-1}, reusing
        # .best_radius() with the critical scales playing the role of
        # distances and t playing the role of radius.
        if (k >= 3L) {
            prev_t <- 1
            for (s in 2L:(k - 1L)) {
                inner_flag <- grp_int <= s
                res <- .best_radius(t_sorted, inner_flag[ord], r_min = prev_t)
                t_s <- res$r
                ellipses[[s]] <- list(cx = cx_fix, cy = cy_fix,
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
            inner_flag <- grp_int <= s
            starts     <- list(.init_ellipse_params(coords, inner_flag))
            if (!is.null(prev)) {
                prev_inflated <- c(prev$cx, prev$cy,
                                   prev$a * 1.2, prev$b * 1.2,
                                   prev$angle)
                starts <- c(starts, list(prev_inflated))
            }
            ell <- .optimize_ellipse(coords, inner_flag,
                                     prev = prev, starts = starts)
            ellipses[[s]] <- ell
            prev          <- ell
        }
    }

    cx_vec    <- vapply(ellipses, function(e) e$cx,    numeric(1))
    cy_vec    <- vapply(ellipses, function(e) e$cy,    numeric(1))
    a_vec     <- vapply(ellipses, function(e) e$a,     numeric(1))
    b_vec     <- vapply(ellipses, function(e) e$b,     numeric(1))
    angle_vec <- vapply(ellipses, function(e) e$angle, numeric(1))

    # ---- Sector assignment ----
    sector <- rep(k, n_pts)
    for (s in (k - 1L):1L) {
        inside         <- .in_ellipse(coords, cx_vec[s], cy_vec[s],
                                      a_vec[s], b_vec[s], angle_vec[s])
        sector[inside] <- s
    }

    # ---- Majority labels and misclassification ----
    majority <- character(k)
    for (s in seq_len(k)) {
        pts_s <- grp_int[sector == s]
        if (length(pts_s) == 0L) {
            majority[s] <- NA_character_
        } else {
            majority[s] <- levels(group)[which.max(tabulate(pts_s, nbins = k))]
        }
    }
    misclass_idx <- which(levels(group)[grp_int] != majority[sector])
    misclass     <- list(n = length(misclass_idx), indices = misclass_idx)

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
            col  = adjustcolor(cols[k], alpha.f = 0.15), border = NA
        )

        if (k >= 3L) {
            for (s in (k - 1L):2L) {
                outer_bnd <- .ellipse_pts(cx_vec[s], cy_vec[s],
                                           a_vec[s], b_vec[s], angle_vec[s])
                inner_bnd <- .ellipse_pts(cx_vec[s - 1L], cy_vec[s - 1L],
                                           a_vec[s - 1L], b_vec[s - 1L],
                                           angle_vec[s - 1L])
                polypath(
                    x    = c(outer_bnd[, 1], NA, inner_bnd[, 1]),
                    y    = c(outer_bnd[, 2], NA, inner_bnd[, 2]),
                    rule = "evenodd",
                    col  = adjustcolor(cols[s], alpha.f = 0.15), border = NA
                )
            }
        }

        draw.ellipse(cx_vec[1L], cy_vec[1L], a = a_vec[1L], b = b_vec[1L],
                     angle = angle_vec[1L] * 180 / pi,
                     col = adjustcolor(cols[1L], alpha.f = 0.15), border = NA)
    }

    # ---- Draw k-1 ellipses ----
    for (s in seq_len(k - 1L)) {
        draw.ellipse(cx_vec[s], cy_vec[s], a = a_vec[s], b = b_vec[s],
                     angle = angle_vec[s] * 180 / pi,
                     border = col, lwd = lwd, lty = lty)
    }

    if (!output) return(invisible(NULL))

    list(
        cx       = cx_vec,
        cy       = cy_vec,
        a        = a_vec,
        b        = b_vec,
        angle    = angle_vec,
        misclass = misclass,
        sector   = sector,
        majority = majority
    )
}
