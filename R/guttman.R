# ===================================================================
# ==== Guttman Monotonicity Coefficient ====
# ===================================================================


#' Guttman's weak monotonicity coefficient mu2
#'
#' Computes Guttman's mu2 coefficient measuring the degree of weak monotone
#' relationship between two numeric vectors. Defined in the HUDAP Manual
#' (2001, p. 219). Pairs where either value is `NA` are excluded (pairwise
#' deletion).
#'
#' @param x Numeric vector.
#' @param y Numeric vector of the same length as `x`.
#'
#' @return A single numeric value in \[-1, 1\]. Returns `NaN` when all pairs
#'   are tied (zero absolute differences).
#'
#' @export
#' @examples
#' mu2(c(1, 2, 3, 4), c(2, 3, 1, 4))
mu2 <- function(x, y) {
    if (!is.numeric(x) || !is.numeric(y)) stop("non-numerical variables!")
    if (length(x) != length(y)) stop("x and y have unequal lengths!")
    valid <- !is.na(x) & !is.na(y)
    x <- x[valid]
    y <- y[valid]
    n <- length(x)
    x_mat <- matrix(rep(x, n), ncol = n)
    y_mat <- matrix(rep(y, n), ncol = n)
    dx <- x_mat - t(x_mat)
    dy <- y_mat - t(y_mat)
    sum(dx * dy) / sum(abs(dx) * abs(dy))
}


#' Matrix of mu2 coefficients for a data frame
#'
#' Computes pairwise Guttman mu2 coefficients for all numeric columns of a
#' data frame, returning either a full symmetric correlation matrix or a
#' lower-triangular distance matrix (`1 - mu2`).
#'
#' @param df A data frame (or coercible object) with numeric columns.
#' @param as_dist If `TRUE`, return a lower-triangular [`dist`] object of
#'   `1 - mu2` values. If `FALSE` (default), return a full symmetric matrix.
#'
#' @return A square numeric matrix with `dimnames` set to the column names of
#'   `df` (when `as_dist = FALSE`), or a [`dist`] object with `Labels` set
#'   to the column names (when `as_dist = TRUE`).
#'
#' @export
#' @examples
#' df <- data.frame(a = 1:5, b = c(2, 3, 1, 5, 4), c = 5:1)
#' mu2df(df)
#' mu2df(df, as_dist = TRUE)
mu2df <- function(df, as_dist = FALSE) {
    df <- as.data.frame(df)
    n_var <- ncol(df)
    for (i in seq_len(n_var))
        if (!is.numeric(df[[i]])) stop("non-numeric variable in data frame!")
    pairs <- combn(n_var, 2)
    n_pairs <- ncol(pairs)
    cor_mat <- matrix(1, nrow = n_var, ncol = n_var,
                      dimnames = list(names(df), names(df)))
    for (i in seq_len(n_pairs)) {
        r <- pairs[1, i]
        s <- pairs[2, i]
        cor_mat[r, s] <- cor_mat[s, r] <- mu2(df[[r]], df[[s]])
    }
    if (as_dist) {
        d <- as.dist(1 - cor_mat, diag = FALSE)
        attr(d, "Labels") <- names(df)
        return(d)
    }
    cor_mat
}
