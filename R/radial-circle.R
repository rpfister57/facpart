# ===================================================================
# ==== Functions for Circular Radial Partitions ====
# ===================================================================


# ---- Internal helpers ----

#' @noRd
.circle_pts <- function(cx, cy, r, n = 200L) {
    theta <- seq(0, 2 * pi, length.out = n + 1L)[-1L]
    cbind(cx + r * cos(theta), cy + r * sin(theta))
}


#' @noRd
.best_radius <- function(s_dist, s_flag, r_min = 0) {
    n        <- length(s_dist)
    best_err <- 2*n
    best_r   <- NA

    for (sp in 1L:(n - 1L)) {
        r_cand <- (s_dist[sp] + s_dist[sp + 1L]) / 2
        if (r_cand < r_min) next
        err <- sum(!s_flag[1L:sp]) + sum(s_flag[(sp + 1L):n])
        if (err < best_err) {
            best_err <- err
            best_r   <- r_cand
        }
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
.optimize_circle <- function(coords,
                             inner_flag,
                             prev_cx, prev_cy, prev_r,
                             starts,
                             meth = "Nelder-Mead") {
    has_prev <- !is.null(prev_r)

    parscale <- c(diff(range(coords[, 1])), diff(range(coords[, 2])))
    parscale[parscale == 0] <- 1

    eval_ctr <- function(ctr) {
        cx_  <- ctr[1]; cy_  <- ctr[2]
        d    <- sqrt((coords[, 1] - cx_)^2 + (coords[, 2] - cy_)^2)
        rmin <- if (has_prev) sqrt((cx_ - prev_cx)^2 + (cy_ - prev_cy)^2) + prev_r else 0
        ord  <- order(d)
        return(
            .best_radius(d[ord], inner_flag[ord], rmin)$misclass)
    }

    # multi-start: keep the optim() result with lowest misclassification
    best <- NULL
    for (s0 in starts) {
        opt <- optim(par = s0,
                     fn = eval_ctr,
                     method = meth,
                     control = list(reltol = 1e-8, maxit = 2000,
                                    parscale = parscale))
        if (is.null(best) || opt$value < best$value) best <- opt
    }

    cx_opt <- best$par[1]; cy_opt <- best$par[2]
    d_opt  <- sqrt((coords[, 1] - cx_opt)^2 + (coords[, 2] - cy_opt)^2)
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
.radial_cuts_2 <- function(coords, grp_int, cx, cy) {
    dists  <- sqrt((coords[, 1] - cx)^2 + (coords[, 2] - cy)^2)
    ord    <- order(dists)
    s_dist <- dists[ord]
    s_grp  <- grp_int[ord]
    n_pts  <- length(dists)

    best_err <- n_pts + 1
    best_r   <- NA

    for (sp in 1L:(n_pts - 1L)) {
        seg1 <- s_grp[1L:sp]
        seg2 <- s_grp[(sp + 1L):n_pts]
        err  <- (sp         - max(tabulate(seg1, nbins = 2L))) +
                (n_pts - sp - max(tabulate(seg2, nbins = 2L)))
        if (err < best_err) {
            best_err <- err
            best_r   <- (s_dist[sp] + s_dist[sp + 1L]) / 2
        }
    }

    sector <- ifelse(dists <= best_r, 1L, 2L)
    return(
        list(radius = best_r, misclass = best_err, sector = sector))
}


#' Binary radial partition (one separating circle)
#'
#' Finds the circle minimising misclassification between two groups of 2D
#' points. The center is found by multi-start Nelder-Mead (starts: data
#' centroid plus each group's centroid) unless `cx` and `cy` are supplied.
#' For a given center, the optimal radius is found by a linear scan over
#' all `n-1` candidate midpoints without assuming which group is inner.
#'
#' @param crd Numeric matrix or data frame with exactly 2 columns.
#' @param group Factor, character, or integer vector with exactly 2 levels.
#' @param cx,cy Center of the separating circle; optimised when `NULL`.
#' @param fill If `TRUE`, shade the inner disc and outer region.
#' @param output If `TRUE` (default), return results list.
#' @param col Circle border colour (default `"purple"`).
#' @param cols Length-2 fill colours (default `c("steelblue", "tomato")`).
#' @param lwd Line width (default `2`).
#' @param lty Line type (default `1`).
#' @param .method `"Nelder-Mead"` (default) or `"SANN"`.
#' @param add If `TRUE` (default), add to existing plot; if `FALSE`, call
#'   `plot()` first.
#'
#' @return If `output = TRUE`, a list with `center`, `radius`, `misclass`
#'   (integer), `misclass_points` (data frame with columns `x`, `y`, `label`
#'   for each misclassified point), `sector` (`1` inside / `2` outside per
#'   point), and `majority` (`character[2]`).
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' inner <- cbind(rnorm(20, 0, 0.3), rnorm(20, 0, 0.3))
#' th <- runif(20, 0, 2 * pi)
#' outer <- cbind(2 * cos(th), 2 * sin(th)) + matrix(rnorm(40, 0, 0.1), 20)
#' crd <- rbind(inner, outer)
#' grp <- factor(c(rep("in", 20), rep("out", 20)))
#' radialCircle(crd, grp, fill = TRUE, add = FALSE)
#' }
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
                       add = TRUE) {

    # ---- Input validation ----
    if (!is.numeric(crd))           stop("Coordinate data must be numeric!")
    if (length(dim(crd)) != 2)      stop("Coordinates must have two dimensions!")
    if (dim(crd)[2] != 2)           stop("Coordinates must have 2 columns!")
    if (nrow(crd) != length(group)) stop("nrow(crd) must equal length(group)!")
    if (!(.method %in% c("Nelder-Mead", "SANN"))) stop("Method not available!")

    group <- as.factor(group)
    if (nlevels(group) != 2) stop("group must have exactly 2 levels!")

    coords  <- as.matrix(crd)
    grp_int <- as.integer(group)

    # ---- Optimise center if not given ----
    # Multi-start: try overall centroid + each group's centroid. The
    # objective (integer misclassification) is a step function, so
    # Nelder-Mead is prone to stalling on plateaus from a single start.
    # parscale sizes the simplex to the data range.
    if (is.null(cx) || is.null(cy)) {
        parscale <- c(diff(range(coords[, 1])), diff(range(coords[, 2])))
        parscale[parscale == 0] <- 1

        starts <- list(
            c(mean(coords[, 1]),              mean(coords[, 2])),
            c(mean(coords[grp_int == 1L, 1]), mean(coords[grp_int == 1L, 2])),
            c(mean(coords[grp_int == 2L, 1]), mean(coords[grp_int == 2L, 2]))
        )

        fnToOpt <- function(p) {
            .radial_cuts_2(coords, grp_int, p[1], p[2])$misclass}

        best <- NULL
        for (s0 in starts) {
            opt <- optim(
                par     = s0,
                fn      = fnToOpt,
                method  = .method,
                control = list(reltol = 1e-8, maxit = 2000,
                               parscale = parscale)
            )
            if (is.null(best) || opt$value < best$value) best <- opt
        }
        cx <- best$par[1]
        cy <- best$par[2]
    }

    # ---- Radius and sectors at chosen center ----
    res    <- .radial_cuts_2(coords, grp_int, cx, cy)
    radius <- res$radius
    sector <- res$sector

    # ---- Unique group assignment ----
    count_mat <- matrix(0L, nrow = 2L, ncol = 2L)
    for (r in 1L:2L) {
        pts_r <- grp_int[sector == r]
        if (length(pts_r) > 0L)
            count_mat[, r] <- tabulate(pts_r, nbins = 2L)
    }
    assignment <- .assign_groups(count_mat)
    majority   <- levels(group)[assignment]

    if (!add) {
        plot(coords, asp = 1)
        graphics::text(coords, labels = group, cex = 0.7, pos = 4)
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
        x     = coords[misclass_idx, 1],
        y     = coords[misclass_idx, 2],
        label = group[misclass_idx]
    )

    if (!output) return(invisible(NULL))

    list(
        center          = c(cx, cy),
        radius          = radius,
        misclass        = length(misclass_idx),
        misclass_points = misclass_points,
        sector          = sector,
        majority        = majority
    )
}


#' K-group radial partition (nested circles)
#'
#' Generalises [radialCircle()] to `k >= 2` groups using `k-1` nested
#' (inclusive) circles. The circles must be nested (circle `s` contains
#' circle `s-1`). Groups are ordered by factor level (level 1 = innermost).
#' Circles are fitted sequentially, each minimising misclassification of
#' groups `1..s` (inside) vs groups `s+1..k` (outside).
#'
#' When `cx` and `cy` are both supplied, all `k-1` circles share that
#' center (concentric); only the radii are searched, subject to
#' `r_s >= r_{s-1}`. When either is `NULL` (default), each circle's
#' center is optimised independently by multi-start Nelder-Mead — the
#' circles are nested but not generally concentric.
#'
#' @param crd Numeric matrix or data frame with exactly 2 columns.
#' @param group Factor with `k >= 2` levels (factor levels define the
#'   inside-to-outside ordering).
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
#' @param add If `TRUE` (default), add to existing plot; if `FALSE`, call
#'   `plot()` first.
#'
#' @return If `output = TRUE`, a list with `cx`, `cy` (`numeric[k-1]`),
#'   `radii` (innermost to outermost), `misclass`, `misclass_points` (data
#'   frame with columns `x`, `y`, `label` for each misclassified point),
#'   `sector` (`integer[n]`), and `majority` (`character[k]`).
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
                        add = TRUE) {

    # ---- Input validation ----
    if (!is.numeric(crd))           stop("Coordinate data must be numeric!")
    if (length(dim(crd)) != 2)      stop("Coordinates must have two dimensions!")
    if (dim(crd)[2] != 2)           stop("Coordinates must have 2 columns")
    if (nrow(crd) != length(group)) stop("nrow(crd) must equal length(group)!")
    if (!(.method %in% c("Nelder-Mead", "SANN"))) stop("Method not available!")

    group <- as.factor(group)
    k     <- nlevels(group)
    if (k < 2L) stop("group must have at least 2 levels!")

    coords  <- as.matrix(crd)
    grp_int <- as.integer(group)
    n_pts   <- nrow(coords)

    if (is.null(cols)) cols <- hcl.colors(k, palette = "Pastel 1")

    # ---- Sequential circle fitting ----
    # If cx/cy are both supplied, all circles share that center
    # (concentric); only radii are searched, subject to r_s >= r_{s-1}.
    # Otherwise each circle's center is optimised independently.
    fixed_center <- !(is.null(cx) || is.null(cy))

    circles <- vector("list", k - 1L)
    prev_cx <- NULL; prev_cy <- NULL; prev_r <- NULL

    if (fixed_center) {
        # Concentric branch: shared (cx, cy) for all k-1 circles.
        d        <- sqrt((coords[, 1] - cx)^2 + (coords[, 2] - cy)^2)
        ord      <- order(d)
        d_sorted <- d[ord]

        for (s in seq_len(k - 1L)) {
            inner_flag <- grp_int <= s
            r_min      <- if (is.null(prev_r)) 0 else prev_r
            res        <- .best_radius(d_sorted, inner_flag[ord], r_min)
            circles[[s]] <- list(cx = cx, cy = cy, r = res$r,
                                 misclass = res$misclass)
            prev_r <- res$r
        }
    } else {
        # Independent-centers branch (original behaviour).
        overall_ctr <- c(mean(coords[, 1]), mean(coords[, 2]))
        for (s in seq_len(k - 1L)) {
            inner_flag <- grp_int <= s
            inner_ctr  <- c(mean(coords[inner_flag, 1]), mean(coords[inner_flag, 2]))

            starts <- list(inner_ctr, overall_ctr)
            if (!is.null(prev_cx)) starts <- c(starts, list(c(prev_cx, prev_cy)))

            circ <- .optimize_circle(coords, inner_flag,
                                     prev_cx, prev_cy, prev_r,
                                     starts,
                                     meth = .method)
            circles[[s]] <- circ
            prev_cx <- circ$cx; prev_cy <- circ$cy; prev_r <- circ$r
        }
    }

    cx_vec    <- vapply(circles, function(c) c$cx, numeric(1))
    cy_vec    <- vapply(circles, function(c) c$cy, numeric(1))
    radii_vec <- vapply(circles, function(c) c$r,  numeric(1))

    # ---- Sector assignment ----
    sector <- rep(k, n_pts)
    for (s in (k - 1L):1L) {
        d <- sqrt((coords[, 1] - cx_vec[s])^2 + (coords[, 2] - cy_vec[s])^2)
        sector[d <= radii_vec[s]] <- s
    }

    # ---- Unique group assignment ----
    count_mat <- matrix(0L, nrow = k, ncol = k)
    for (r in seq_len(k)) {
        pts_r <- grp_int[sector == r]
        if (length(pts_r) > 0L)
            count_mat[, r] <- tabulate(pts_r, nbins = k)
    }
    assignment <- .assign_groups(count_mat)
    majority   <- levels(group)[assignment]

    if (!add) {
        plot(coords, asp = 1)
        graphics::text(coords, labels = group, cex = 0.7, pos = 4)
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
        x     = coords[misclass_idx, 1],
        y     = coords[misclass_idx, 2],
        label = group[misclass_idx]
    )

    if (!output) return(invisible(NULL))

    list(
        cx              = cx_vec,
        cy              = cy_vec,
        radii           = radii_vec,
        misclass        = length(misclass_idx),
        misclass_points = misclass_points,
        sector          = sector,
        majority        = majority
    )
}
