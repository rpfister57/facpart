# ===================================================================
# ==== Functions for Ellipses ====
# ===================================================================


#' Fit and draw an ellipse around a 2D configuration
#'
#' Three fitting strategies are available:
#'
#' - `"maxdist"` — covariance-based ellipse, scaled so all points lie inside
#'   with a 5% buffer. Semi-axes are `sqrt(eigenvalue) * k`, where `k` is the
#'   maximum Mahalanobis distance from the center.
#' - `"bestfit"` — minimises the sum of squared normalised elliptic distances
#'   `sqrt((u/a)^2 + (v/b)^2) - 1` via BFGS, starting from the covariance
#'   ellipse.
#' - `"minbound"` — minimum-area enclosing ellipse via [cluster::ellipsoidhull()].
#'
#' The ellipse is drawn on the currently active plot window.
#'
#' @param crd Numeric matrix or data frame with exactly 2 columns (x, y).
#' @param labs Currently unused (kept for backward compatibility).
#' @param mid `"centroid"` (default) places the ellipse center at the data
#'   centroid; `"nullnull"` forces the center to the origin (0, 0).
#'   Only affects `"maxdist"`.
#' @param typ Fitting strategy: `"maxdist"` (default), `"bestfit"`, or
#'   `"minbound"`.
#' @param output If `TRUE` (default), return a list describing the ellipse.
#' @param col Border colour (default `"black"`).
#' @param fill Fill colour, or `NA` (default) for no fill.
#' @param lwd Line width (default `2`).
#' @param lty Line type (default `2`).
#' @param add If `TRUE` (default), add to existing plot; if `FALSE`, call
#'   `plot()` first.
#'
#' @return If `output = TRUE`, a list with elements `cx`, `cy`, `a`, `b`,
#'   `angle` (in radians).
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' crd <- cbind(rnorm(30), 0.5 * rnorm(30))
#' ellipseInConfig(crd, typ = "maxdist", add = FALSE)
#' ellipseInConfig(crd, typ = "bestfit", col = "blue")
#' ellipseInConfig(crd, typ = "minbound", col = "red")
#' }
#'
#' @export
ellipseInConfig <- function(crd,
                            labs = NULL,
                            mid = "centroid",
                            typ = "maxdist",
                            output = TRUE,
                            col = "black",
                            fill = NA,
                            lwd = 2,
                            lty = 2,
                            add = TRUE) {
    available_mid <- c("centroid", "nullnull")
    available_typ <- c("maxdist", "bestfit", "minbound")

    if (!is.numeric(crd)) stop("Coordinate data must be numeric!")
    if (length(dim(crd)) != 2) stop("Coordinates must have 2 dimensions!")
    if (dim(crd)[2] != 2) stop("Coordinates must be 2-dimensional!")
    if (!(mid %in% available_mid)) mid <- "centroid"
    if (!(typ %in% available_typ)) typ <- "maxdist"

    coords <- as.matrix(crd)

    if (mid == "centroid") {
        centroid_x <- mean(coords[ , 1])
        centroid_y <- mean(coords[ , 2])
    }
    else {
        centroid_x <- 0
        centroid_y <- 0
    }

    center <- c(centroid_x, centroid_y)

    if (!add) plot(coords, asp = 1)

    # covariance structure of the point cloud (starting geometry)
    S     <- cov(coords)
    eig   <- eigen(S)   # eigenvalues in decreasing order
    S_inv <- solve(S)   # inverse of S

    # angle of the major axis (eigenvector of largest eigenvalue),
    #   in radians
    eig_angle <- atan2(eig$vectors[2, 1], eig$vectors[1, 1])

    # Mahalanobis distances from center for all points
    centered <- sweep(coords, 2, center)
    mah_dist <- sqrt(rowSums((centered %*% S_inv) * centered))


    # == Maximum distance ellipse "maxdist" ==
    if (typ == "maxdist") {
        eps <- max(mah_dist) * 0.05
        k   <- max(mah_dist) + eps
        a   <- sqrt(eig$values[1]) * k
        b   <- sqrt(eig$values[2]) * k

        draw.ellipse(centroid_x, centroid_y, a = a, b = b,
                     angle = eig_angle * 180 / pi,
                     border = col, col = fill, lwd = lwd, lty = lty)

        if (output) outList <- list(cx = centroid_x, cy = centroid_y,
                                    a = a, b = b, angle = eig_angle)
    }


    # == Best fitting ellipse "bestfit" ==
    if (typ == "bestfit") {

        ellipseResiduals <- function(params, coords) {
            cx  <- params[1]
            cy  <- params[2]
            a   <- abs(params[3])
            b   <- abs(params[4])
            ang <- params[5]
            cos_a <- cos(ang)
            sin_a <- sin(ang)
            u <- (coords[ , 1] - cx) * cos_a + (coords[ , 2] - cy) * sin_a
            v <- -(coords[ , 1] - cx) * sin_a + (coords[ , 2] - cy) * cos_a
            residuals <- sqrt((u / a)^2 + (v / b)^2) - 1
            return(sum(residuals^2))
        }

        k0 <- mean(mah_dist)
        a0 <- sqrt(eig$values[1]) * k0
        b0 <- sqrt(eig$values[2]) * k0

        fit_result <- optim(par = c(centroid_x, centroid_y,
                                    a0, b0, eig_angle),
                            fn  = ellipseResiduals,
                            coords = coords,
                            method = "BFGS")

        fitted_cx  <- fit_result$par[1]
        fitted_cy  <- fit_result$par[2]
        fitted_a   <- abs(fit_result$par[3])
        fitted_b   <- abs(fit_result$par[4])
        fitted_ang <- fit_result$par[5]

        draw.ellipse(fitted_cx, fitted_cy, a = fitted_a, b = fitted_b,
                     angle = fitted_ang * 180 / pi,
                     border = col, col = fill, lwd = lwd, lty = lty)

        if (output) outList <- list(cx = fitted_cx, cy = fitted_cy,
                                    a = fitted_a, b = fitted_b,
                                    angle = fitted_ang)
    }


    # == Minimum bounding ellipse "minbound" ==
    if (typ == "minbound") {
        hull   <- ellipsoidhull(coords)
        min_cx <- hull$loc[1]
        min_cy <- hull$loc[2]
        eig_h  <- eigen(hull$cov)
        ang_h  <- atan2(eig_h$vectors[2, 1], eig_h$vectors[1, 1])
        d      <- sqrt(hull$d2)
        min_a  <- sqrt(eig_h$values[1]) * d
        min_b  <- sqrt(eig_h$values[2]) * d

        draw.ellipse(min_cx, min_cy, a = min_a, b = min_b,
                     angle = ang_h * 180 / pi,
                     border = col, col = fill, lwd = lwd, lty = lty)

        if (output) outList <- list(cx = min_cx, cy = min_cy,
                                    a = min_a, b = min_b, angle = ang_h)
    }

    if (output) return(outList)
}


#' Classify points as inside or outside a given ellipse
#'
#' @param crd Numeric matrix or data frame with exactly 2 columns (x, y).
#' @param cx,cy Center of the ellipse.
#' @param a,b Semi-major and semi-minor axes.
#' @param angle Rotation of the major axis in radians (as returned by
#'   [ellipseInConfig()]).
#'
#' @return A character vector of length `nrow(crd)` with values `"inside"`
#'   or `"outside"`.
#'
#' @examples
#' \dontrun{
#' crd <- cbind(rnorm(30), rnorm(30))
#' inoutEllipse(crd, cx = 0, cy = 0, a = 1, b = 1, angle = 0)
#' }
#'
#' @export
inoutEllipse <- function(crd, cx, cy, a, b, angle) {
    cos_a <- cos(angle)
    sin_a <- sin(angle)
    u <- (crd[ , 1] - cx) * cos_a + (crd[ , 2] - cy) * sin_a
    v <- -(crd[ , 1] - cx) * sin_a + (crd[ , 2] - cy) * cos_a
    return(ifelse((u / a)^2 + (v / b)^2 <= 1, "inside", "outside"))
}
