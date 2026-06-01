# ===================================================================
# ==== Functions for Angular Partitions ====
# ===================================================================


# ---- Internal helpers ----

#' @noRd
#' e.g.: two radians b > a, delta is the diff b - a
#' returns the radiant of the line bisecting b - a
.arc_mid <- function(a, b) {
    delta <- (b - a) %% (2 * pi)
    return(atan2(sin(a + delta / 2), cos(a + delta / 2)))
}


#' @noRd
#' given: center cx, cy; find: best partitions by brute force search
.angular_search <- function(coords, grp_int, 
                            cx, cy, 
                            k, n_pts, combos) {
    
    # order points according to their angles
    pt_angles <- atan2(coords[, 2] - cy, coords[, 1] - cx)
    pt_radii  <- sqrt((coords[, 1] - cx)^2 + (coords[, 2] - cy)^2)
    ord       <- order(pt_angles)
    s_ang     <- pt_angles[ord]
    s_grp     <- grp_int[ord]
    s_rad     <- pt_radii[ord]

    best_err    <- Inf
    best_margin <- -Inf
    best_cuts   <- NULL

    # loop over point sequences:
    for (start in 0L:(n_pts - 1L)) {
        ri <- ((seq_len(n_pts) - 1L + start) %% n_pts) + 1L
        rg <- s_grp[ri]
        ra <- s_ang[ri]
        rr <- s_rad[ri]

        # loop over all k-1 combinations of n-1 points:
        for (ci in seq_len(ncol(combos))) {
            sp     <- combos[, ci]
            bounds <- c(0, sp, n_pts)

            err <- 0
            
            # loop over 1..k partitions:
            for (s in seq_len(k)) {
                seg <- rg[seq(bounds[s] + 1L, bounds[s + 1L])]
                err <- err + length(seg) - max(tabulate(seg, nbins = k))
            }

            # 2D margin at each cut: r_min * sin(angular gap / 2).
            gaps_int   <- (ra[sp + 1L] - ra[sp]) %% (2 * pi)
            r_int      <- pmin(rr[sp], rr[sp + 1L])
            margin_int <- r_int * sin(gaps_int / 2)

            gap_wrap    <- (ra[1L] - ra[n_pts]) %% (2 * pi)
            r_wrap      <- min(rr[n_pts], rr[1L])
            margin_wrap <- r_wrap * sin(gap_wrap / 2)

            margin <- min(margin_int, margin_wrap)

            if (err < best_err ||
                (err == best_err && margin > best_margin)) {
                best_err    <- err
                best_margin <- margin
                best_cuts   <- c(
                    vapply(seq_len(k - 1L), function(s)
                        .arc_mid(ra[sp[s]], ra[sp[s] + 1L]), numeric(1)),
                    .arc_mid(ra[n_pts], ra[1L])
                )
            }
        }
    }

    list(misclass  = best_err,
         margin    = best_margin,
         cuts      = best_cuts,
         pt_angles = pt_angles)
}


#' Angular k-way (wedge) partition of a 2D configuration
#'
#' Finds `k` rays emanating from a center point that partition a 2D
#' configuration of `n` points into `k` angular sectors, minimising total
#' misclassification (points whose group differs from the majority group
#' in their sector).
#'
#' **Search at fixed center.** Points are sorted by their angle from
#' the center. For each circular rotation (n choices, implicitly placing
#' the wrap-around cut) and each of `combn(n-1, k-1)` ways to choose the
#' internal split positions, total misclassification is computed.
#' Cut angles are placed at the arc midpoint between adjacent points.
#' Tie-breaker: among k-tuples with the same misclass count, the one with
#' the largest minimum 2D perpendicular distance from a cut ray to the
#' nearest point wins.
#'
#' **Search optimal center.** When `cx` and `cy` are `NULL` (default), the center
#' is optimised by multi-start Nelder-Mead — the brute-force above runs
#' as the inner objective at each candidate center. Starts are the data
#' centroid plus the centroid of each non-empty group. `parscale` is set
#' to the data range. When `cx` and `cy` are supplied, they are used
#' directly (no optimisation).
#'
#' @param crd Numeric matrix or data frame with exactly 2 columns; no NAs.
#' @param group Factor with `k >= 2` levels, same nrow as crd.
#' @param cx,cy Center; optimised when either is `NULL`.
#' @param output If `TRUE` (default), return list of results.
#' @param col Ray colour (default `"darkorange"`).
#' @param lwd Line width (default `2`).
#' @param lty Line type (default `1`).
#' @param add If `TRUE` (default), add to existing plot, else plot configuration.
#'
#' @return If `output = TRUE`, a list with:
#'   - `cuts` — `numeric[k]`, cut angles of rays in radians.
#'   - `margin` — minimum distance from any cut ray to nearest point.
#'   - `misclass` — integer, number of misclassified points.
#'   - `sector` — `integer[n]`, sector `1..k` per point
#'   - `majority` — `character[k]`, majority group per sector
#'   - `center` — `c(cx, cy)`
#'
#' @examples
#' \dontrun{
#' set.seed(1)
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
    if (length(dim(crd)) != 2)      stop("Coordinates must have two dimensions!")
    if (dim(crd)[2] != 2)           stop("Coordinates must be 2-dimensional!")
    if (nrow(crd) != length(group)) stop("nrow(crd) must equal length(group)!")
    if (!is.numeric(as.matrix(crd)))  stop("Coordinate data must be numeric!")

    group <- as.factor(group)
    coords  <- as.matrix(crd)
    
    k <- nlevels(group)
    if (k < 2)         stop("group must have at least 2 levels!")
    if (nrow(crd) < k) stop("Number of points must be >= number of groups!")

    n_pts   <- nrow(coords)
    grp_int <- as.integer(group)

    # Generate split-position combinations.
    # Note: Instead of (n over k), 
    # the parameterization n * (n-1 over k-1) is used.
    combos <- combn(n_pts - 1L, k - 1L)

    # ---- Center: optimise if NULL, use as given otherwise ----
    if (is.null(cx) || is.null(cy)) {
        x_range <- diff(range(coords[, 1]))
        y_range <- diff(range(coords[, 2]))
        if (x_range == 0) x_range <- 1
        if (y_range == 0) y_range <- 1
        parscale <- c(x_range, y_range)

        starts <- list(c(mean(coords[, 1]), mean(coords[, 2])))
        for (g in seq_len(k)) {
            in_g <- which(grp_int == g)
            if (length(in_g) > 0L) {
                starts <- c(starts,
                            list(c(mean(coords[in_g, 1]),
                                   mean(coords[in_g, 2]))))
            }
        }

        fnToOpt <- function(p) {
            .angular_search(coords, grp_int, p[1], p[2],
                            k, n_pts, combos)$misclass
        }

        best <- NULL
        for (s0 in starts) {
            opt <- optim(par     = s0,
                         fn      = fnToOpt,
                         method  = "Nelder-Mead",
                         control = list(reltol = 1e-8, maxit = 2000,
                                        parscale = parscale))
            if (is.null(best) || opt$value < best$value) best <- opt
        }
        cx <- best$par[1]
        cy <- best$par[2]
    }

    # ---- Final search at chosen center ----
    res         <- .angular_search(coords, grp_int, 
                                   cx, cy, 
                                   k, n_pts, combos)
    best_err    <- res$misclass
    best_cuts   <- res$cuts
    best_margin <- res$margin
    pt_angles   <- res$pt_angles
    
    if (!add) {
        plot(coords, pch = 19, col = "blue",
             xlim = c(-1.5, 1.5), ylim = c(-1, 1),
             asp = 1)
        text(coords, labels = group,
             cex = 0.7, pos = 4)
    }

    # ---- Draw rays from center at the optimal cut angles ----
    usr     <- par("usr")
    ray_len <- 2 * sqrt((usr[2] - usr[1])^2 + (usr[4] - usr[3])^2)

    for (ang in best_cuts) {
        segments(cx, cy,
                 cx + ray_len * cos(ang),
                 cy + ray_len * sin(ang),
                 col = col, lwd = lwd, lty = lty)
    }

    if (!output) return(invisible(NULL))

    # ---- Sector assignment for each point (in original input order) ----
    sc       <- sort(best_cuts %% (2 * pi))
    norm_pts <- pt_angles %% (2 * pi)
    fi       <- findInterval(norm_pts, sc)
    sector   <- ifelse(fi == 0L | fi == k, k, fi)

    # ---- Majority group label in each sector ----
    majority <- sapply(seq_len(k), function(s) {
        pts <- grp_int[sector == s]
        levels(group)[which.max(tabulate(pts, nbins = k))]
    })

    list(
        cuts     = best_cuts,
        margin   = best_margin,
        misclass = best_err,
        sector   = sector,
        majority = majority,
        center   = c(cx, cy)
    )
}

