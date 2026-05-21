# ===================================================================
# ==== Functions for Axial (Parallel-Line) Partitions ====
# ===================================================================


#' Binary LDA separating line for a 2D configuration
#'
#' Fits a linear discriminant analysis to two-group 2D data and draws the
#' classical LDA boundary — perpendicular to LD1 through the midpoint of the
#' class means — on the current plot.
#'
#' When LD1 is along the x-axis (`abs(w[2]) < 1e-10`) the separator is
#' vertical: `slope` is returned as `Inf` and `intercept` carries the
#' x-position of the line.
#'
#' **Not a special case of [axialLines()].** `axialLine()` returns the
#' classical LDA Bayes-rule boundary (closed-form: LD1 direction, midpoint
#' of class means as cut) — optimal under multivariate-normal classes with
#' equal covariances, and the standard reference object in psychometric /
#' Facet Theory pipelines. The return value includes `predicted`, the LDA
#' class assignment per point.
#'
#' [axialLines()] with `k = 2`, in contrast, searches the full
#' `(angle, cut)` space for the empirical minimum-misclass parallel line.
#' Its direction need not be LD1 and its cut need not match the midpoint of
#' means; it can overfit the training sample when classes are non-Gaussian
#' or have unequal covariances. Use it when you want the empirically best
#' linear separator; use `axialLine()` when you want the classical LDA
#' classifier.
#'
#' @param crd Numeric matrix or data frame with exactly 2 columns (x, y).
#' @param group Factor, character, or integer vector with exactly 2 levels.
#' @param fill If `TRUE`, shade the two half-planes (default `FALSE`).
#' @param output If `TRUE` (default), return a list of line parameters.
#' @param col Separator line colour (default `"purple"`).
#' @param cols Length-2 fill colours; used when `fill = TRUE`
#'   (default `c("steelblue", "tomato")`).
#' @param lwd Line width (default `2`).
#' @param lty Line type (default `1`).
#'
#' @return If `output = TRUE`, a list with `slope`, `intercept`, `misclass`
#'   (integer), and `predicted` (LDA-predicted class factor).
#'
#' @seealso [axialLines()] for the empirically-optimal parallel-line
#'   partition over both angle and cut position (any `k >= 2`).
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' crd <- rbind(cbind(rnorm(15, -1), rnorm(15)),
#'              cbind(rnorm(15,  1), rnorm(15)))
#' grp <- factor(c(rep("a", 15), rep("b", 15)))
#' plot(crd, asp = 1)
#' axialLine(crd, grp, fill = TRUE)
#' }
#'
#' @export
axialLine <- function(crd,
                    group,
                    fill = FALSE,
                    output = TRUE,
                    col = "purple",
                    cols = c("steelblue", "tomato"),
                    lwd = 2,
                    lty = 1) {

    # ---- Input validation ----
    if (length(dim(crd)) != 2)       stop("Coordinates must have two dimensions!")
    if (dim(crd)[2] != 2)            stop("Coordinates must be 2-dimensional!")
    if (nrow(crd) != length(group))  stop("nrow(crd) must equal length(group)!")
    if (!is.numeric(as.matrix(crd))) stop("Coordinate data must be numeric!")
    
    group <- as.factor(group)
    if (nlevels(group) != 2) stop("group must have exactly 2 levels!")

    coords  <- as.matrix(crd)
    levels_ <- levels(group)

    # ---- Fit LDA ----
    lda_fit <- lda(coords, grouping = group)

    # ---- Separating line ----
    # Vertical-line guard: when LD1 is along the x-axis (w[2] ~= 0) the
    # separator is vertical. In that case slope = Inf and the intercept
    # field carries the x-position of the line.
    w        <- lda_fit$scaling[, 1]
    mu1      <- colMeans(coords[group == levels_[1], , drop = FALSE])
    mu2      <- colMeans(coords[group == levels_[2], , drop = FALSE])
    midpoint <- (mu1 + mu2) / 2

    vertical <- abs(w[2]) < 1e-10
    if (vertical) {
        slope     <- Inf
        intercept <- midpoint[1]
    } else {
        slope     <- -w[1] / w[2]
        intercept <- midpoint[2] - slope * midpoint[1]
    }

    # ---- Optional half-plane shading ----
    if (fill) {
        usr    <- par("usr")
        xlim_p <- usr[1:2]
        ylim_p <- usr[3:4]

        if (vertical) {
            x_vert    <- intercept
            mu1_left  <- mu1[1] < x_vert
            left_col  <- if (mu1_left) cols[1] else cols[2]
            right_col <- if (mu1_left) cols[2] else cols[1]
            polygon(c(xlim_p[1], x_vert, x_vert, xlim_p[1]),
                    c(ylim_p[1], ylim_p[1], ylim_p[2], ylim_p[2]),
                    col = adjustcolor(left_col, alpha.f = 0.15), border = NA)
            polygon(c(x_vert, xlim_p[2], xlim_p[2], x_vert),
                    c(ylim_p[1], ylim_p[1], ylim_p[2], ylim_p[2]),
                    col = adjustcolor(right_col, alpha.f = 0.15), border = NA)
        } else {
            y_at_xlim <- slope * xlim_p + intercept

            mu1_above <- mu1[2] > slope * mu1[1] + intercept
            above_col <- if (mu1_above) cols[1] else cols[2]
            below_col <- if (mu1_above) cols[2] else cols[1]

            polygon(c(xlim_p[1], xlim_p[2], xlim_p[2], xlim_p[1]),
                    c(y_at_xlim[1], y_at_xlim[2], ylim_p[2], ylim_p[2]),
                    col = adjustcolor(above_col, alpha.f = 0.15), border = NA)
            polygon(c(xlim_p[1], xlim_p[2], xlim_p[2], xlim_p[1]),
                    c(ylim_p[1], ylim_p[1], y_at_xlim[2], y_at_xlim[1]),
                    col = adjustcolor(below_col, alpha.f = 0.15), border = NA)
        }
    }

    # ---- Draw separator line ----
    if (vertical) {
        abline(v = intercept, col = col, lwd = lwd, lty = lty)
    } else {
        abline(a = intercept, b = slope, col = col, lwd = lwd, lty = lty)
    }

    if (!output) return(invisible(NULL))

    # ---- Misclassification ----
    predicted <- predict(lda_fit)$class
    misclass  <- sum(predicted != group)

    list(
        slope     = slope,
        intercept = intercept,
        misclass  = misclass,
        predicted = predicted
    )
}


#' Multi-group parallel-line partition (angle + cuts searched)
#'
#' Partitions a 2D configuration into `k >= 2` groups using `k-1` parallel
#' separating lines. Both the direction of the lines and their positions
#' are searched to minimise empirical misclassification:
#'
#' - **Angle**: grid of `n_angles` equally-spaced angles in `[0, pi)`,
#'   augmented by LDA's LD1 direction. For each angle `theta`, points are
#'   projected onto `(cos theta, sin theta)`.
#' - **Cuts**: for each angle, exact brute-force enumeration over all
#'   `C(n-1, k-1)` ways to split the sorted projections into `k` segments.
#'
#' The `(angle, cuts)` combination with lowest total misclassification wins.
#' Tie-breaker: among configurations with the same misclass count, the one
#' with the largest minimum margin (perpendicular distance from a cut line
#' to the nearest point) is preferred — keeps cut lines visually away from
#' the data points.
#'
#' For `k = 2`, `axialLines()` can differ substantially from [axialLine()]:
#' the latter places the cut at the midpoint of class means on the LD1
#' projection (classical Bayes rule under normality), while `axialLines()`
#' searches the full `(angle, cut)` space for the empirical optimum.
#'
#' Vertical-line guard: when the winning direction has `w[2] ~= 0`, the
#' separators are vertical. `slope` is returned as `Inf` and `intercepts`
#' carry the x-positions of the lines.
#'
#' @param crd Numeric matrix or data frame with exactly 2 columns.
#' @param group Factor, character, or integer vector with `k >= 2` levels.
#' @param fill If `TRUE`, shade the `k` strips between consecutive lines.
#' @param output If `TRUE` (default), return a list of results.
#' @param col Separator line colour (default `"purple"`).
#' @param cols Length-`k` fill colours; auto-generated if `NULL`.
#' @param lwd Line width (default `2`).
#' @param lty Line type (default `1`).
#' @param n_angles Number of angles in `[0, pi)` to scan (default `180L`,
#'   approximately 1-degree resolution); LDA's direction is added as an
#'   extra candidate.
#'
#' @seealso [axialLine()] for the classical LDA Bayes-rule boundary
#'   (closed-form, binary only). Use [axialLine()] when you want the
#'   statistical LDA classifier; use `axialLines()` when you want the
#'   empirically optimal parallel-line partition.
#'
#' @return If `output = TRUE`, a list with:
#'   - `slope` — shared slope of the lines (`Inf` if vertical)
#'   - `intercepts` — `numeric[k-1]`, y-intercepts (or x-positions if vertical)
#'   - `angle` — winning direction in radians
#'   - `margin` — minimum perpendicular distance from the cut line(s) to the
#'     nearest point
#'   - `misclass` — total misclassified points
#'   - `sector` — `integer[n]`, sector `1..k` per point in input order
#'   - `majority` — `character[k]`, majority group per sector
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' crd <- rbind(cbind(rnorm(10, -1), rnorm(10)),
#'              cbind(rnorm(10,  0), rnorm(10)),
#'              cbind(rnorm(10,  1), rnorm(10)))
#' grp <- factor(c(rep("a", 10), rep("b", 10), rep("c", 10)))
#' plot(crd, asp = 1)
#' axialLines(crd, grp, fill = TRUE)
#' }
#'
#' @export
axialLines <- function(crd,
                     group,
                     fill = FALSE,
                     output = TRUE,
                     col = "purple",
                     cols = NULL,
                     lwd = 2,
                     lty = 1,
                     n_angles = 180L) {

    # ---- Input validation ----
    if (length(dim(crd)) != 2) stop("Coordinates must have two dimensions!")
    if (dim(crd)[2] != 2)      stop("Coordinates must be 2-dimensional!")
    if (nrow(crd) != length(group)) stop("nrow(crd) must equal length(group)!")
    if (n_angles < 1L)         stop("n_angles must be >= 1!")
    if (!is.numeric(crd))      stop("Coordinate data must be numeric!")
    
    group <- as.factor(group)
    k     <- nlevels(group)
    if (k < 2) stop("group must have at least 2 levels!")

    coords  <- as.matrix(crd)
    levels_ <- levels(group)
    grp_int <- as.integer(group)
    n_pts   <- nrow(coords)

    # ---- Candidate angles ----
    # Grid in [0, pi) plus LDA's LD1 direction as a high-quality seed.
    lda_fit   <- lda(coords, grouping = group)
    w_lda     <- lda_fit$scaling[, 1]
    theta_lda <- atan2(w_lda[2], w_lda[1]) %% pi
    thetas    <- c(seq(0, pi, length.out = n_angles + 1L)[-(n_angles + 1L)],
                   theta_lda)

    # Pre-compute split-position combinations once (independent of angle).
    if (k >= 3L) combos <- combn(n_pts - 1L, k - 1L)

    # ---- Joint search over (angle, cuts) ----
    # Tie-breaker: among configurations with the same misclassification
    # count, prefer the one with the largest minimum margin. The cut at
    # the midpoint of s_proj[sp] and s_proj[sp+1] has margin
    # (s_proj[sp+1] - s_proj[sp]) / 2; for k cuts, take the min.
    best_err    <- .Machine$integer.max
    best_margin <- -Inf
    best_theta  <- theta_lda
    best_cuts   <- numeric(k - 1L)

    for (theta in thetas) {
        w_t    <- c(cos(theta), sin(theta))
        proj   <- as.numeric(coords %*% w_t)
        ord    <- order(proj)
        s_proj <- proj[ord]
        s_grp  <- grp_int[ord]

        if (k == 2L) {
            for (sp in 1L:(n_pts - 1L)) {
                seg1 <- s_grp[1L:sp]
                seg2 <- s_grp[(sp + 1L):n_pts]
                err  <- (length(seg1) - max(tabulate(seg1, nbins = k))) +
                        (length(seg2) - max(tabulate(seg2, nbins = k)))
                margin <- (s_proj[sp + 1L] - s_proj[sp]) / 2
                if (err < best_err ||
                    (err == best_err && margin > best_margin)) {
                    best_err    <- err
                    best_margin <- margin
                    best_theta  <- theta
                    best_cuts   <- (s_proj[sp] + s_proj[sp + 1L]) / 2
                }
            }
        } else {
            for (ci in seq_len(ncol(combos))) {
                sp     <- combos[, ci]
                bounds <- c(0L, sp, n_pts)
                err    <- 0L
                for (s in seq_len(k)) {
                    seg <- s_grp[(bounds[s] + 1L):bounds[s + 1L]]
                    err <- err + length(seg) - max(tabulate(seg, nbins = k))
                }
                margin <- min(s_proj[sp + 1L] - s_proj[sp]) / 2
                if (err < best_err ||
                    (err == best_err && margin > best_margin)) {
                    best_err    <- err
                    best_margin <- margin
                    best_theta  <- theta
                    best_cuts   <- (s_proj[sp] + s_proj[sp + 1L]) / 2
                }
            }
        }
    }

    # ---- Convert cuts to line geometry ----
    w        <- c(cos(best_theta), sin(best_theta))
    vertical <- abs(w[2]) < 1e-10
    if (vertical) {
        slope      <- Inf
        intercepts <- best_cuts / w[1]
    } else {
        slope      <- -w[1] / w[2]
        intercepts <- best_cuts / w[2]
    }

    # ---- Sector assignment (original point order) ----
    proj_best <- as.numeric(coords %*% w)
    sector    <- findInterval(proj_best, sort(best_cuts)) + 1L

    # ---- Majority labels per sector ----
    majority <- character(k)
    for (s in seq_len(k)) {
        pts_s <- grp_int[sector == s]
        if (length(pts_s) == 0L) {
            majority[s] <- NA_character_
        } else {
            majority[s] <- levels_[which.max(tabulate(pts_s, nbins = k))]
        }
    }

    # ---- Optional strip shading ----
    if (fill) {
        if (is.null(cols)) {
            cols <- hcl.colors(k, palette = "Pastel 1")
        }
        usr   <- par("usr")
        xl    <- usr[1]; xr <- usr[2]
        yb    <- usr[3]; yt <- usr[4]

        cuts_sorted <- sort(best_cuts)
        cuts_ext    <- c(-Inf, cuts_sorted, Inf)

        if (vertical) {
            for (s in seq_len(k)) {
                lo <- cuts_ext[s]
                hi <- cuts_ext[s + 1L]

                x_lo <- if (is.finite(lo)) lo / w[1] else (if (w[1] > 0) xl else xr)
                x_hi <- if (is.finite(hi)) hi / w[1] else (if (w[1] > 0) xr else xl)
                if (x_lo > x_hi) { tmp <- x_lo; x_lo <- x_hi; x_hi <- tmp }
                x_lo <- max(xl, min(xr, x_lo))
                x_hi <- max(xl, min(xr, x_hi))

                polygon(c(x_lo, x_hi, x_hi, x_lo),
                        c(yb,   yb,   yt,   yt),
                        col = adjustcolor(cols[s], alpha.f = 0.15), border = NA)
            }
        } else {
            for (s in seq_len(k)) {
                lo <- cuts_ext[s]
                hi <- cuts_ext[s + 1L]

                y_lo_xl <- if (is.finite(lo)) (lo - w[1] * xl) / w[2] else (if (w[2] > 0) yb else yt)
                y_lo_xr <- if (is.finite(lo)) (lo - w[1] * xr) / w[2] else (if (w[2] > 0) yb else yt)
                y_hi_xl <- if (is.finite(hi)) (hi - w[1] * xl) / w[2] else (if (w[2] > 0) yt else yb)
                y_hi_xr <- if (is.finite(hi)) (hi - w[1] * xr) / w[2] else (if (w[2] > 0) yt else yb)

                y_lo_xl <- max(yb, min(yt, y_lo_xl))
                y_lo_xr <- max(yb, min(yt, y_lo_xr))
                y_hi_xl <- max(yb, min(yt, y_hi_xl))
                y_hi_xr <- max(yb, min(yt, y_hi_xr))

                polygon(c(xl,    xr,    xr,    xl),
                        c(y_lo_xl, y_lo_xr, y_hi_xr, y_hi_xl),
                        col = adjustcolor(cols[s], alpha.f = 0.15), border = NA)
            }
        }
    }

    # ---- Draw k-1 parallel separator lines ----
    for (j in seq_along(intercepts)) {
        if (vertical) {
            abline(v = intercepts[j], col = col, lwd = lwd, lty = lty)
        } else {
            abline(a = intercepts[j], b = slope, col = col, lwd = lwd, lty = lty)
        }
    }

    if (!output) return(invisible(NULL))

    list(
        slope      = slope,
        intercepts = intercepts,
        angle      = best_theta,
        margin     = best_margin,
        misclass   = best_err,
        sector     = sector,
        majority   = majority
    )
}
