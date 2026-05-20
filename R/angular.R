# ===================================================================
# ==== Functions for Angular Partitions ====
# ===================================================================


# ---- Internal helpers ----

#' @noRd
.arc_mid <- function(a, b) {
    delta <- (b - a) %% (2 * pi)
    atan2(sin(a + delta / 2), cos(a + delta / 2))
}


#' @noRd
.angular_search <- function(coords, grp_int, cx, cy, k, n_pts, combos) {
    pt_angles <- atan2(coords[, 2] - cy, coords[, 1] - cx)
    pt_radii  <- sqrt((coords[, 1] - cx)^2 + (coords[, 2] - cy)^2)
    ord       <- order(pt_angles)
    s_ang     <- pt_angles[ord]
    s_grp     <- grp_int[ord]
    s_rad     <- pt_radii[ord]

    best_err    <- Inf
    best_margin <- -Inf
    best_cuts   <- NULL

    for (start in 0L:(n_pts - 1L)) {
        ri <- ((seq_len(n_pts) - 1L + start) %% n_pts) + 1L
        rg <- s_grp[ri]
        ra <- s_ang[ri]
        rr <- s_rad[ri]

        for (ci in seq_len(ncol(combos))) {
            sp     <- combos[, ci]
            bounds <- c(0L, sp, n_pts)

            err <- 0L
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


#' K-way angular (wedge) partition of a 2D configuration
#'
#' Finds `k` rays emanating from a center point that partition a 2D
#' configuration into `k` angular sectors, minimising total
#' misclassification (points whose group differs from the majority group
#' in their sector).
#'
#' **Cut search at fixed center.** Points are sorted by their angle from
#' the center. For each circular rotation (n choices, implicitly placing
#' the wrap-around cut) and each of `C(n-1, k-1)` ways to choose the
#' internal split positions, total misclassification is computed.
#' Cut angles are placed at the arc midpoint between adjacent points.
#' Tie-breaker: among k-tuples with the same misclass count, the one with
#' the largest minimum 2D perpendicular distance from a cut ray to the
#' nearest point wins.
#'
#' **Center search.** When `cx` and `cy` are `NULL` (default), the center
#' is optimised by multi-start Nelder-Mead — the brute-force above runs
#' as the inner objective at each candidate center. Starts are the data
#' centroid plus the centroid of each non-empty group. `parscale` is set
#' to the data range. When `cx` and `cy` are supplied, they are used
#' directly (no optimisation).
#'
#' @param crd Numeric matrix or data frame with exactly 2 columns.
#' @param group Factor with `k >= 2` levels.
#' @param cx,cy Center; optimised when either is `NULL`.
#' @param output If `TRUE` (default), return results list.
#' @param col Ray colour (default `"darkorange"`).
#' @param lwd Line width (default `2`).
#' @param lty Line type (default `1`).
#'
#' @return If `output = TRUE`, a list with:
#'   - `cuts` — `numeric[k]`, cut angles in radians
#'   - `margin` — minimum 2D distance from any cut ray to nearest point
#'   - `misclass` — integer
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
                              lty = 1) {

    # ---- Input validation ----
    if (!is.numeric(crd))           stop("Coordinate data must be numeric!")
    if (length(dim(crd)) != 2)      stop("Coordinates must have two dimensions!")
    if (dim(crd)[2] != 2)           stop("Coordinates must be 2-dimensional!")
    if (nrow(crd) != length(group)) stop("nrow(crd) must equal length(group)!")

    group <- as.factor(group)
    k <- nlevels(group)
    if (k < 2)         stop("group must have at least 2 levels!")
    if (nrow(crd) < k) stop("Number of points must be >= number of groups!")

    coords  <- as.matrix(crd)
    n_pts   <- nrow(coords)
    grp_int <- as.integer(group)

    # Pre-generate split-position combinations once.
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
    res         <- .angular_search(coords, grp_int, cx, cy, k, n_pts, combos)
    best_err    <- res$misclass
    best_cuts   <- res$cuts
    best_margin <- res$margin
    pt_angles   <- res$pt_angles

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


#' Three-way angular partition (original fixed-k=3 version)
#'
#' Original implementation specialised to `k = 3` groups. Kept for reference;
#' [angularPartition()] generalises this and is recommended for new code.
#' Algorithm: for each rotation of the angular sort and each internal split
#' pair `(i, j)`, compute total misclassification across the three sectors;
#' keep the minimum.
#'
#' @inheritParams angularPartition
#'
#' @return List with `cuts` (numeric[3] in `(-pi, pi]`), `misclass`,
#'   `sector` (integer with values in 1..3), `majority` (character[3]),
#'   and `center`.
#'
#' @examples
#' \dontrun{
#' angularPartition3(crd, grp)
#' }
#'
#' @export
angularPartition3 <- function(crd,
                               group,
                               cx = NULL,
                               cy = NULL,
                               output = TRUE,
                               col = "darkorange",
                               lwd = 2,
                               lty = 1) {

    # ---- Input validation ----
    if (!is.numeric(crd))           stop("Coordinate data must be numeric!")
    if (length(dim(crd)) != 2)      stop("Coordinates must have two dimensions!")
    if (dim(crd)[2] != 2)           stop("Coordinates must be 2-dimensional!")
    if (nrow(crd) != length(group)) stop("nrow(crd) must equal length(group)!")

    group <- as.factor(group)
    if (nlevels(group) != 3) stop("group must have exactly 3 levels!")
    if (nrow(crd) < 3)       stop("At least 3 points are required!")

    coords <- as.matrix(crd)
    n_pts  <- nrow(coords)

    if (is.null(cx)) cx <- mean(coords[, 1])
    if (is.null(cy)) cy <- mean(coords[, 2])

    grp_int <- as.integer(group)

    pt_angles <- atan2(coords[, 2] - cy, coords[, 1] - cx)
    ord   <- order(pt_angles)
    s_ang <- pt_angles[ord]
    s_grp <- grp_int[ord]

    arc_mid <- function(a, b) {
        delta <- (b - a) %% (2 * pi)
        atan2(sin(a + delta / 2), cos(a + delta / 2))
    }

    best_err  <- Inf
    best_cuts <- NULL

    for (start in 0:(n_pts - 1)) {
        ri <- ((seq_len(n_pts) - 1 + start) %% n_pts) + 1
        rg <- s_grp[ri]
        ra <- s_ang[ri]

        for (i in 1:(n_pts - 2)) {
            for (j in (i + 1):(n_pts - 1)) {
                err <-
                    (i         - max(tabulate(rg[seq_len(i)],        nbins = 3))) +
                    (j - i     - max(tabulate(rg[seq(i + 1, j)],     nbins = 3))) +
                    (n_pts - j - max(tabulate(rg[seq(j + 1, n_pts)], nbins = 3)))

                if (err < best_err) {
                    best_err  <- err
                    best_cuts <- c(
                        arc_mid(ra[i],     ra[i + 1]),
                        arc_mid(ra[j],     ra[j + 1]),
                        arc_mid(ra[n_pts], ra[1])
                    )
                }
            }
        }
    }

    usr     <- par("usr")
    ray_len <- 2 * sqrt((usr[2] - usr[1])^2 + (usr[4] - usr[3])^2)

    for (ang in best_cuts) {
        segments(cx, cy,
                 cx + ray_len * cos(ang),
                 cy + ray_len * sin(ang),
                 col = col, lwd = lwd, lty = lty)
    }

    if (!output) return(invisible(NULL))

    sc       <- sort(best_cuts %% (2 * pi))
    norm_pts <- pt_angles %% (2 * pi)

    sector <- ifelse(norm_pts >= sc[1] & norm_pts < sc[2], 1L,
              ifelse(norm_pts >= sc[2] & norm_pts < sc[3], 2L, 3L))

    majority <- sapply(1:3, function(s) {
        pts <- grp_int[sector == s]
        levels(group)[which.max(tabulate(pts, nbins = 3))]
    })

    list(
        cuts     = best_cuts,
        misclass = best_err,
        sector   = sector,
        majority = majority,
        center   = c(cx, cy)
    )
}
