# Tests for R/axial.R
#
# Both functions draw, and axialLines() reads par("usr") when fill = TRUE, so
# every call is wrapped in pdf(NULL) / dev.off(). Assertions are on the
# returned geometry and misclass fields, never on the plot.

with_null_dev <- function(expr) {
    pdf(NULL)
    on.exit(dev.off(), add = TRUE)
    plot.new()
    plot.window(xlim = c(-10, 10), ylim = c(-10, 10))
    force(expr)
}

# Bijection-constrained misclassification, computed independently of the
# package internals.
bij_err <- function(sector, grp_int, k) {
    perms <- gtools::permutations(k, k)
    best  <- -1L
    for (i in seq_len(nrow(perms))) {
        correct <- sum(grp_int == perms[i, ][sector])
        if (correct > best) best <- correct
    }
    length(grp_int) - best
}

# Slow reference for axialLines(): the same angle grid and the same cut
# enumeration, scored by an explicit k! loop. Shares no code with the package.
ref_axial <- function(coords, grp_int, k, n_angles) {
    n      <- nrow(coords)
    fit    <- MASS::lda(coords, grouping = factor(grp_int))
    w_lda  <- fit$scaling[, 1]
    thetas <- c(seq(0, pi, length.out = n_angles + 1L)[-(n_angles + 1L)],
                atan2(w_lda[2], w_lda[1]) %% pi)
    combos <- combn(n - 1L, k - 1L)
    best   <- Inf
    b_marg <- -Inf
    for (theta in thetas) {
        proj_  <- as.numeric(coords %*% c(cos(theta), sin(theta)))
        ord    <- order(proj_)
        s_proj <- proj_[ord]
        s_grp  <- grp_int[ord]
        for (j in seq_len(ncol(combos))) {
            cp  <- combos[, j]
            bnd <- c(0L, cp, n)
            sec <- integer(n)
            for (s in seq_len(k)) sec[(bnd[s] + 1L):bnd[s + 1L]] <- s
            err <- bij_err(sec, s_grp, k)
            mrg <- min((s_proj[cp + 1L] - s_proj[cp]) / 2)
            if (err < best || (err == best && mrg > b_marg)) {
                best   <- err
                b_marg <- mrg
            }
        }
    }
    list(misclass = best, margin = b_marg)
}


# ---- axialLines(): the search minimises what gets reported ----

test_that("axialLines() reports the error its own lines achieve", {
    set.seed(41)
    for (i in 1:25) {
        k      <- sample(2:4, 1)
        n      <- k * 6
        coords <- cbind(rnorm(n), rnorm(n))
        gi     <- rep(seq_len(k), each = 6)
        grp    <- factor(letters[seq_len(k)][gi])
        res    <- with_null_dev(axialLines(coords, grp, n_angles = 24L))

        # sector must follow from the returned geometry
        proj_ <- as.numeric(coords %*% c(cos(res$angle), sin(res$angle)))
        cuts  <- if (is.infinite(res$slope)) res$intercepts * cos(res$angle)
                 else                        res$intercepts * sin(res$angle)
        expect_equal(res$sector, findInterval(proj_, sort(cuts)) + 1L)

        expect_equal(res$misclass, bij_err(res$sector, gi, k))
        expect_equal(res$misclass, nrow(res$misclass_points))
        expect_length(res$intercepts, k - 1L)
    }
})

test_that("axialLines() attains the brute-force optimum, margin included", {
    set.seed(42)
    for (i in 1:20) {
        k      <- sample(2:3, 1)
        n      <- k * 5
        coords <- cbind(rnorm(n), rnorm(n))
        gi     <- rep(seq_len(k), each = 5)
        grp    <- factor(letters[seq_len(k)][gi])
        got    <- with_null_dev(axialLines(coords, grp, n_angles = 20L))
        want   <- ref_axial(coords, gi, k, 20L)
        expect_equal(got$misclass, want$misclass)
        expect_equal(got$margin,   want$margin)
    }
})

test_that("axialLines() does not score segments by their local majority", {
    # Two segments both dominated by "a" cannot both be labelled "a" once
    # sector/majority are derived, so a score that lets them is unreachable.
    # Whatever the search settles on must be achievable by some bijection.
    set.seed(43)
    for (i in 1:20) {
        k      <- 3L
        n      <- 15L
        coords <- cbind(rnorm(n), rnorm(n))
        # deliberately lopsided groups: one label dominates, which is what
        # makes the unconstrained and bijection-constrained scores diverge
        gi     <- c(rep(1L, 11), 2L, 2L, 3L, 3L)
        grp    <- factor(letters[seq_len(k)][gi])
        res    <- with_null_dev(axialLines(coords, grp, n_angles = 18L))
        expect_equal(res$misclass, bij_err(res$sector, gi, k))
        expect_equal(res$misclass, ref_axial(coords, gi, k, 18L)$misclass)
    }
})

test_that("axialLines() handles k = 2 and perfectly separable data", {
    coords <- rbind(cbind(rnorm(8, -6, 0.2), rnorm(8, 0, 0.2)),
                    cbind(rnorm(8,  6, 0.2), rnorm(8, 0, 0.2)))
    grp    <- factor(rep(c("a", "b"), each = 8))
    res    <- with_null_dev(axialLines(coords, grp))
    expect_equal(res$misclass, 0L)
    expect_length(res$intercepts, 1L)
    expect_equal(nrow(res$misclass_points), 0L)
    # the single cut must fall between the clusters
    expect_gt(res$margin, 1)
})

test_that("axialLines() vertical-line guard returns slope = Inf", {
    # Two clusters separated along x only: the best direction is theta = 0,
    # giving vertical separators.
    coords <- rbind(cbind(rnorm(8, -6, 0.2), rnorm(8, 0, 3)),
                    cbind(rnorm(8,  6, 0.2), rnorm(8, 0, 3)))
    grp    <- factor(rep(c("a", "b"), each = 8))
    res    <- with_null_dev(axialLines(coords, grp, n_angles = 180L))
    expect_equal(res$misclass, 0L)
    expect_true(is.infinite(res$slope) || abs(res$slope) > 50)
    if (is.infinite(res$slope)) {
        # intercepts carry x-positions, between the two clusters
        expect_gt(res$intercepts[1], -6)
        expect_lt(res$intercepts[1],  6)
    }
})


# ---- axialLine(): the LDA classifier ----

test_that("axialLine() returns the LDA boundary and its own error count", {
    set.seed(44)
    coords <- rbind(cbind(rnorm(12, -2), rnorm(12)),
                    cbind(cbind(rnorm(12,  2), rnorm(12))))
    grp    <- factor(rep(c("a", "b"), each = 12))
    res    <- with_null_dev(axialLine(coords, grp))

    # axialLine() is the odd one out: `predicted`, no sector/majority
    expect_true(is.factor(res$predicted))
    expect_null(res$sector)
    expect_equal(res$misclass, sum(res$predicted != grp))
    expect_equal(res$misclass, nrow(res$misclass_points))
    # matches MASS::lda's own predictions
    expect_equal(as.character(res$predicted),
                 as.character(predict(MASS::lda(coords, grouping = grp))$class))
})

test_that("axialLine() requires exactly 2 group levels", {
    coords <- cbind(rnorm(9), rnorm(9))
    grp3   <- factor(rep(c("a", "b", "c"), each = 3))
    expect_error(with_null_dev(axialLine(coords, grp3)), "exactly 2 levels")
})


# ---- Validation ----

test_that("both axial functions reject NAs and malformed input", {
    coords <- cbind(c(1, 2, 3, 9), c(0, 1, 0, 1))
    grp    <- factor(c("a", "a", "b", "b"))

    crd_na <- coords; crd_na[2, 1] <- NA
    grp_na <- grp;    grp_na[2]    <- NA

    for (f in list(axialLine, axialLines)) {
        expect_error(with_null_dev(f(crd_na, grp)),      "NAs? allowed in crd")
        expect_error(with_null_dev(f(coords, grp_na)),   "NAs? allowed in group")
        expect_error(with_null_dev(f(cbind(coords, 1), grp)), "2-dimensional")
        expect_error(with_null_dev(f(coords, grp[-1])),  "must equal")
    }
    expect_error(with_null_dev(axialLines(coords, grp, n_angles = 0L)),
                 "n_angles must be >= 1")
    expect_error(with_null_dev(axialLines(coords, factor(rep("a", 4)))),
                 "at least 2 levels")
})

test_that("output = FALSE returns invisible(NULL)", {
    coords <- cbind(c(1, 2, 3, 9), c(0, 1, 0, 1))
    grp    <- factor(c("a", "a", "b", "b"))
    expect_null(with_null_dev(axialLine(coords,  grp, output = FALSE)))
    expect_null(with_null_dev(axialLines(coords, grp, output = FALSE)))
})
