# ===================================================================
# ==== Functions for Axial (Parallel-Line) Partitions ====
# ===================================================================


#' Binary LDA separating line for a 2D point configuration
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
#' @param add If `TRUE` (default), add to existing plot; if `FALSE`, call
#'   `plot()` first.
#'
#' @return If `output = TRUE`, a list with `slope`, `intercept`, `misclass`
#'   (integer), `misclass_points` (data frame with columns `x`, `y`, `label`
#'   for each misclassified point), and `predicted` (LDA-predicted class factor).
#'
#' @seealso [axialLines()] for the empirically-optimal parallel-line
#'   partition over both angle and cut position (any `k >= 2`).
#'
#' @examples
#' set.seed(1)
#' crd <- rbind(cbind(rnorm(15, -1), rnorm(15)),
#'              cbind(rnorm(15,  1), rnorm(15)))
#' grp <- factor(c(rep("a", 15), rep("b", 15)))
#' axialLine(crd, grp, fill = TRUE, add = FALSE)
#'
#' @export
axialLine <- function(crd,
                    group,
                    fill = FALSE,
                    output = TRUE,
                    col = "purple",
                    cols = c("steelblue", "tomato"),
                    lwd = 2,
                    lty = 1,
                    add = TRUE) {

    # ---- Input validation ----
    if (any(is.na(crd)))             stop("No NA allowed in crd!")
    if (any(is.na(group)))           stop("No NA allowed in group!")
    if (length(dim(crd)) != 2)       stop("Coordinates must have two dimensions!")
    if (dim(crd)[2] != 2)            stop("Coordinates must be 2-dimensional!")
    if (nrow(crd) != length(group))  stop("nrow(crd) must equal length(group)!")
    if (!is.numeric(as.matrix(crd))) stop("Coordinate data must be numeric!")
    
    group <- as.factor(group)
    if (nlevels(group) != 2) stop("group must have exactly 2 levels!")

    coords  <- as.matrix(crd)
    levels_ <- levels(group)

    # ---- Fit LDA ----
    lda_fit <- MASS::lda(coords, grouping = group)

    # ---- Compute separating line ----
    # Vertical-line guard: when LD1 is along the x-axis (w[2] ~= 0) the
    # separator is vertical. In that case slope = Inf and the intercept
    # field carries the x-position of the line.
    w        <- lda_fit$scaling[, 1]
    mu1      <- lda_fit$means[1, ]
    mu2      <- lda_fit$means[2, ]
    midpoint <- (mu1 + mu2) / 2

    vertical <- abs(w[2]) < 1e-10
    
    if (vertical) {
        slope     <- Inf
        intercept <- midpoint[1]
    } else {
        slope     <- -w[1] / w[2]
        intercept <- midpoint[2] - slope * midpoint[1]
    }

    if (!add) {
        plot(coords, asp = 1)
        graphics::text(coords, labels = group, cex = 0.7, pos = 4)
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
    predicted       <- predict(lda_fit)$class
    misclass_idx    <- which(predicted != group)
    misclass        <- length(misclass_idx)
    misclass_points <- data.frame(
        x     = coords[misclass_idx, 1],
        y     = coords[misclass_idx, 2],
        label = group[misclass_idx]
    )

    list(
        slope           = slope,
        intercept       = intercept,
        misclass        = misclass,
        misclass_points = misclass_points,
        predicted       = predicted
    )
}


#' Multi-group parallel-line partition (angle + cuts searched)
#'
#' Partitions a 2D configuration into `k >= 2` groups using `k-1` parallel
#' separating lines. Both the slopes of the lines and their positions
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
#' For `k = 2`, `axialLines()` can differ from [axialLine()]:
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
#' @param add If `TRUE` (default), add to existing plot; if `FALSE`, call
#'   `plot()` first.
#'
#' @seealso [axialLine()] for the classical LDA Bayes-rule boundary
#'   (closed-form, binary only) if k=2. Use [axialLine()] when you want the
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
#'   - `misclass_points` — data frame (`x`, `y`, `label`) of misclassified points
#'   - `sector` — `integer[n]`, sector `1..k` per point in input order
#'   - `majority` — `character[k]`, majority group per sector
#'
#' @examples
#' set.seed(1)
#' crd <- rbind(cbind(rnorm(10, -1), rnorm(10)),
#'              cbind(rnorm(10,  0), rnorm(10)),
#'              cbind(rnorm(10,  1), rnorm(10)))
#' grp <- factor(c(rep("a", 10), rep("b", 10), rep("c", 10)))
#' axialLines(crd, grp, fill = TRUE, add = FALSE)
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
                     n_angles = 180L,
                     add = TRUE) {

    # ---- Input validation ----
    if (any(is.na(crd)))       stop("No NAs allowed in crd!")
    if (any(is.na(crd)))       stop("No NAs allowed in group!")
    if (length(dim(crd)) != 2) stop("Coordinates must have two dimensions!")
    if (dim(crd)[2] != 2)      stop("Coordinates must be 2-dimensional!")
    if (nrow(crd) != length(group)) stop("nrow(crd) must equal length(group)!")
    if (n_angles < 1L)         stop("n_angles must be >= 1!")
    if (!is.numeric(as.matrix(crd))) stop("Coordinate data must be numeric!")
    
    group <- as.factor(group)
    k     <- nlevels(group)
    if (k < 2) stop("group must have at least 2 levels!")

    # ---- define constants ----
    coords  <- as.matrix(crd)
    levels_ <- levels(group)
    grp_int <- as.integer(group)
    n_pts   <- nrow(coords)
    
    # combos: compute cut-position combinations (independent of angle)
    # each (n over k) selection of k points defines k cut-positions
    if (k >= 3L) combos <- combn(n_pts - 1L, k - 1L)
    
    # ---- Candidate angles to check (0, pi) ----
    # Grid in [0, pi) plus LDA's LD1 direction as a high-quality seed.
    lda_fit   <- lda(coords, grouping = group)
    w_lda     <- lda_fit$scaling[, 1]
    theta_lda <- atan2(w_lda[2], w_lda[1]) %% pi
    thetas    <- c(seq(0, pi, length.out = n_angles + 1L)[-(n_angles + 1L)],
                   theta_lda)

    # ---- Joint search over (angle, cuts) ----
    # For each angle the points are projected and sorted, then cumulative
    # per-group counts (cum[i + 1, g] = number of group-g points among the
    # first i sorted points) let every candidate partition be scored by
    # vector arithmetic with no per-cut re-tabulation. A segment spanning
    # sorted positions a..b has group counts cum[b + 1, ] - cum[a, ], and
    # contributes (b - a + 1) - max_g(count) misclassified points; all
    # C(n - 1, k - 1) partitions are scored together via vectorised pmax.
    #
    # Tie-breaker: among configurations with the same misclassification
    # count, prefer the one with the largest minimum margin. The cut at the
    # midpoint of s_proj[sp] and s_proj[sp + 1] has margin
    # (s_proj[sp + 1] - s_proj[sp]) / 2; for k - 1 cuts, take the min.
    best_err    <- n_pts
    best_margin <- -Inf
    best_theta  <- theta_lda
    best_cuts   <- numeric(k - 1L)

    # loop over all thetas (angles):
    for (theta in thetas) {
        # projection on line: d = x cos theta + y sin theta
        w_t    <- c(cos(theta), sin(theta))
        proj_   <- as.numeric(coords %*% w_t)
        ord    <- order(proj_)
        s_proj <- proj_[ord]
        
        # group sequence along the current theta line:
        s_grp  <- grp_int[ord]

        # cum: table of cumulative per-group counts over the sorted sequence:
        cum <- matrix(0L, nrow = n_pts + 1L, ncol = k)
        for (g in seq_len(k)) cum[-1L, g] <- cumsum(s_grp == g)
        total <- cum[n_pts + 1L, ]

        if (k == 2L) {
            # Single cut after sorted position sp (1 .. n - 1). Inner counts
            # are cum[sp + 1, ], outer counts total - cum[sp + 1, ].
            sp     <- 1L:(n_pts - 1L)
            c_in   <- cum[sp + 1L, , drop = FALSE]
            max_in <- pmax(c_in[, 1L], c_in[, 2L])
            max_ou <- pmax(total[1L] - c_in[, 1L], total[2L] - c_in[, 2L])
            err    <- (sp - max_in) + ((n_pts - sp) - max_ou)
            margin <- (s_proj[sp + 1L] - s_proj[sp]) / 2
            cuts   <- matrix((s_proj[sp] + s_proj[sp + 1L]) / 2, ncol = 1L)
            } 
        else {
            # combos is (k - 1) x M (M is number of possible cuts); 
            # column ci gives the k - 1 cut positions.
            # Accumulate, per segment, the max group count across all M
            # columns; segment s spans positions combos[s - 1] + 1 .. combos[s]
            # (with implicit bounds 0 and n_pts at the ends).
            M       <- ncol(combos)
            sum_max <- integer(M)
            gap     <- matrix(0, nrow = M, ncol = k - 1L)
            
            # loop over the k groups:
            for (s in seq_len(k)) {
                lo    <- if (s == 1L) rep(0L, M)     else combos[s - 1L, ]
                hi    <- if (s == k)  rep(n_pts, M)  else combos[s, ]
                cnt   <- cum[hi + 1L, , drop = FALSE] - 
                          cum[lo + 1L, , drop = FALSE]
                r_max <- cnt[, 1L]
                for (g in 2L:k) r_max <- pmax(r_max, cnt[, g])
                sum_max <- sum_max + r_max
            }
            
            err <- n_pts - sum_max

            # Per-column gaps, margins, and cut midpoints.
            for (s in seq_len(k - 1L))
                gap[, s] <- s_proj[combos[s, ] + 1L] - s_proj[combos[s, ]]
            
            margin <- (if (k - 1L == 1L) gap[, 1L]
                       else do.call(pmin, as.data.frame(gap))) / 2
            
            cuts   <- t((matrix(s_proj[combos],      nrow = k - 1L) +
                         matrix(s_proj[combos + 1L], nrow = k - 1L)) / 2)
        }

        # Best at this angle: lowest err, ties broken by largest margin.
        # which.max returns the first maximum, so the lowest-index (earliest
        # in scan order) column wins ties, matching the original loops.
        m  <- min(err)
        at <- which(err == m)
        j  <- at[which.max(margin[at])]
        if (m < best_err || (m == best_err && margin[j] > best_margin)) {
            best_err    <- m
            best_margin <- margin[j]
            best_theta  <- theta
            best_cuts   <- cuts[j, ]
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

    # ---- Unique group assignment ----
    count_mat <- matrix(0L, nrow = k, ncol = k)
    for (r in seq_len(k)) {
        pts_r <- grp_int[sector == r]
        if (length(pts_r) > 0L)
            count_mat[, r] <- tabulate(pts_r, nbins = k)
    }
    assignment <- .assign_groups(count_mat)
    majority   <- levels_[assignment]

    if (!add) {
        plot(coords, asp = 1)
        graphics::text(coords, labels = group, cex = 0.7, pos = 4)
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

        # Each strip s is the band lo <= w . p <= hi (in projection units),
        # where the same projection w drives the sector assignment above.
        # Build it by clipping the plot rectangle against the two boundary
        # half-planes. This is exact for any slope -- including steep /
        # vertical lines that exit through the top and bottom edges, where
        # corner-clamping a left-edge-to-right-edge polygon would collapse
        # the band and merge strips into a single colour.
        rect_x <- c(xl, xr, xr, xl)
        rect_y <- c(yb, yb, yt, yt)
        for (s in seq_len(k)) {
            lo   <- cuts_ext[s]
            hi   <- cuts_ext[s + 1L]
            poly <- list(x = rect_x, y = rect_y)
            if (is.finite(hi))
                poly <- .clip_halfplane(poly$x, poly$y, w, hi, keep_le = TRUE)
            if (is.finite(lo))
                poly <- .clip_halfplane(poly$x, poly$y, w, lo, keep_le = FALSE)
            if (length(poly$x) >= 3L)
                polygon(poly$x, poly$y,
                        col = adjustcolor(cols[s], alpha.f = 0.15), border = NA)
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

    misclass_idx    <- which(as.character(group) != majority[sector])
    misclass_points <- data.frame(
        x     = coords[misclass_idx, 1],
        y     = coords[misclass_idx, 2],
        label = group[misclass_idx]
    )

    if (!output) return(invisible(NULL))

    list(
        slope           = slope,
        intercepts      = intercepts,
        angle           = best_theta,
        margin          = best_margin,
        misclass        = length(misclass_idx),
        misclass_points = misclass_points,
        sector          = sector,
        majority        = majority
    )
}


#' Clip a polygon against a projection half-plane (Sutherland-Hodgman)
#'
#' Keeps the part of polygon `(px, py)` on one side of the line
#' `w[1] * x + w[2] * y = bound`. With `keep_le = TRUE` the retained side is
#' `w . p <= bound`; with `keep_le = FALSE` it is `w . p >= bound`. Used by
#' [axialLines()] to shade the strip between two parallel cut lines for any
#' slope, including steep / vertical separators.
#'
#' @noRd
.clip_halfplane <- function(px, py, w, bound, keep_le) {
    n <- length(px)
    if (n == 0L) return(list(x = px, y = py))

    # f <= 0 marks the kept side, regardless of orientation.
    sgn <- if (keep_le) 1 else -1
    f   <- sgn * (w[1L] * px + w[2L] * py - bound)

    ox <- numeric(0)
    oy <- numeric(0)
    for (i in seq_len(n)) {
        j  <- if (i == n) 1L else i + 1L
        fi <- f[i]
        fj <- f[j]
        if (fi <= 0) {
            ox <- c(ox, px[i])
            oy <- c(oy, py[i])
        }
        # Edge crosses the boundary: add the intersection point.
        if ((fi < 0 && fj > 0) || (fi > 0 && fj < 0)) {
            t  <- fi / (fi - fj)
            ox <- c(ox, px[i] + t * (px[j] - px[i]))
            oy <- c(oy, py[i] + t * (py[j] - py[i]))
        }
    }
    list(x = ox, y = oy)
}
