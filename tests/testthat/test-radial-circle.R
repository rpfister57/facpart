# Tests for R/radial-circle.R
#
# These functions all draw, and several read par("usr"), so every call is
# wrapped in pdf(NULL) / dev.off() to avoid opening the default device and
# leaving a stray Rplots.pdf behind. Assertions are on the returned geometry
# and misclass fields, never on the plot.

# Run `expr` with a null graphics device and an active plot.
with_null_dev <- function(expr) {
    pdf(NULL)
    on.exit(dev.off(), add = TRUE)
    plot.new()
    plot.window(xlim = c(-10, 10), ylim = c(-10, 10))
    force(expr)
}

# Bijection-constrained misclassification of a region assignment, computed
# independently of the package internals.
bij_err <- function(sector, grp_int, k) {
    perms <- gtools::permutations(k, k)
    best  <- -1L
    for (i in seq_len(nrow(perms))) {
        correct <- sum(grp_int == perms[i, ][sector])
        if (correct > best) best <- correct
    }
    length(grp_int) - best
}


# ---- .best_radius(): the nesting constraint must not drop feasible radii ----

test_that(".best_radius() clamps to r_min instead of skipping a feasible split", {
    s_dist <- c(1, 2, 3, 4)
    s_flag <- c(TRUE, TRUE, FALSE, FALSE)

    # Unconstrained: the perfect split sits at the midpoint 2.5.
    expect_equal(.best_radius(s_dist, s_flag, r_min = 0)$misclass, 0L)
    expect_equal(.best_radius(s_dist, s_flag, r_min = 0)$r, 2.5)

    # r_min above that midpoint but below the next distance still admits the
    # perfect split, realised at r = r_min rather than at the midpoint.
    for (r_min in c(2.6, 2.9, 2.999)) {
        res <- .best_radius(s_dist, s_flag, r_min = r_min)
        expect_equal(res$misclass, 0L)
        expect_equal(res$r, r_min)
    }
})

test_that(".best_radius() never returns a radius below r_min", {
    set.seed(11)
    for (i in 1:200) {
        n      <- sample(2:12, 1)
        s_dist <- sort(round(runif(n, 0, 5), 1))
        s_flag <- sample(c(TRUE, FALSE), n, replace = TRUE)
        r_min  <- runif(1, 0, 5)
        expect_gte(.best_radius(s_dist, s_flag, r_min)$r, r_min)
    }
})

test_that(".best_radius() reports the error its radius actually achieves", {
    set.seed(12)
    for (i in 1:300) {
        n      <- sample(2:12, 1)
        # rounding produces plenty of coincident distances
        s_dist <- sort(round(runif(n, 0, 4), 1))
        s_flag <- sample(c(TRUE, FALSE), n, replace = TRUE)
        r_min  <- if (i %% 2L == 0L) 0 else runif(1, 0, 4)
        res    <- .best_radius(s_dist, s_flag, r_min)
        # inside test is `dist <= r`
        actual <- sum(xor(s_dist <= res$r, s_flag))
        expect_equal(res$misclass, actual)
    }
})

test_that(".best_radius() does not split coincident distances", {
    # All four points sit at distance 1, so no radius separates them: any
    # radius >= 1 puts all four inside (2 errors), any smaller radius puts
    # all four outside (2 errors).
    res <- .best_radius(c(1, 1, 1, 1), c(TRUE, FALSE, TRUE, FALSE), r_min = 0)
    expect_equal(res$misclass, 2L)
    expect_equal(sum(xor(c(1, 1, 1, 1) <= res$r, c(TRUE, FALSE, TRUE, FALSE))), 2L)
})


# ---- .cut_radii(): realisable cut positions ----

test_that(".cut_radii() enumerates exactly the realisable positions", {
    # distinct distances: positions 0..n are all realisable
    cp <- .cut_radii(c(1, 2, 3))
    expect_equal(cp$pos, c(0L, 1L, 2L, 3L))
    expect_equal(cp$radius[1:3], c(0.5, 1.5, 2.5))

    # a tie between positions 2 and 3 removes position 2
    cp <- .cut_radii(c(1, 2, 2, 4))
    expect_equal(cp$pos, c(0L, 1L, 3L, 4L))

    # a point at the center removes position 0
    cp <- .cut_radii(c(0, 1, 2))
    expect_false(0L %in% cp$pos)

    # every radius realises the position it is paired with
    set.seed(13)
    for (i in 1:100) {
        s_dist <- sort(round(runif(sample(2:10, 1), 0, 3), 1))
        cp     <- .cut_radii(s_dist)
        expect_equal(vapply(cp$radius, function(r) sum(s_dist <= r), integer(1)),
                     cp$pos)
    }
})


# ---- The search minimises the quantity that is reported ----

test_that("radialCircle() reports the error its own circle achieves", {
    set.seed(21)
    for (i in 1:60) {
        n   <- sample(5:14, 1)
        d   <- round(runif(n, 0.2, 8), 1)
        gi  <- sample(1:2, n, replace = TRUE)
        if (length(unique(gi)) < 2L) next
        crd <- cbind(d, 0)
        grp <- factor(c("a", "b")[gi])
        res <- with_null_dev(radialCircle(crd, grp, cx = 0, cy = 0))

        expect_equal(res$sector, ifelse(d <= res$radius, 1L, 2L))
        expect_equal(res$misclass, bij_err(res$sector, gi, 2L))
        expect_equal(res$misclass, nrow(res$misclass_points))
    }
})

test_that("radialCircle() finds the optimal radius at a fixed center", {
    # Exhaustive reference over every radius that yields a distinct split,
    # including the two degenerate ones.
    best_err <- function(d, gi) {
        cands <- c(0, sort(unique(d)), max(d) + 1)
        min(vapply(cands, function(r) bij_err(ifelse(d <= r, 1L, 2L), gi, 2L),
                   numeric(1)))
    }
    set.seed(22)
    for (i in 1:60) {
        n  <- sample(5:14, 1)
        d  <- round(runif(n, 0.2, 8), 1)
        gi <- sample(1:2, n, replace = TRUE)
        if (length(unique(gi)) < 2L) next
        res <- with_null_dev(
            radialCircle(cbind(d, 0), factor(c("a", "b")[gi]), cx = 0, cy = 0))
        expect_equal(res$misclass, best_err(d, gi))
    }
})

test_that("the search score is the score that gets reported", {
    # The old search scored a split by each region's own local majority,
    # which lets both regions claim the same group and so can report a total
    # no bijection can achieve. Here it scored this configuration 1 while
    # the partition it chose actually costs 2. Whatever the search returns
    # now must be achievable by some region-to-group bijection.
    crd <- cbind(1:8, 0)
    grp <- factor(c("a", "a", "a", "b", "a", "a", "a", "a"))
    gi  <- as.integer(grp)

    d   <- sqrt(rowSums(crd^2))
    ord <- order(d)
    for (k in 2L) {
        res <- .radial_search(d[ord], gi[ord], k)
        sec <- .nested_sector(crd, 0, 0, res$radii, k)
        expect_equal(res$misclass, bij_err(sec, gi, k))
    }
    res <- with_null_dev(radialCircle(crd, grp, cx = 0, cy = 0))
    expect_equal(res$misclass, bij_err(res$sector, gi, 2L))
})

test_that(".refine_radii() never worsens the k-way result and can improve it", {
    # The sequential fit scores each circle in isolation (groups 1..s inside
    # vs the rest outside), which is not the reported criterion, so the
    # refinement pass must be able to improve on it and must never make it
    # worse. Tested on the helper directly: the nesting order the exported
    # function settles on is searched, so reconstructing "the" sequential fit
    # from outside is no longer well defined.
    set.seed(27)
    improved <- 0L
    for (i in 1:60) {
        k   <- sample(3:4, 1)
        n   <- k * 7
        crd <- cbind(rnorm(n), rnorm(n))
        gi  <- rep(seq_len(k), each = 7)

        # a raw level-order sequential fit, centers included
        cxv <- numeric(k - 1L); cyv <- numeric(k - 1L); raw <- numeric(k - 1L)
        prev_cx <- NULL; prev_cy <- NULL; prev_r <- NULL
        for (s in seq_len(k - 1L)) {
            inner <- gi <= s
            circ  <- .optimize_circle(
                crd, inner, prev_cx, prev_cy, prev_r,
                starts = list(c(mean(crd[inner, 1]), mean(crd[inner, 2])),
                              c(mean(crd[, 1]), mean(crd[, 2]))))
            cxv[s] <- circ$cx; cyv[s] <- circ$cy; raw[s] <- circ$r
            prev_cx <- circ$cx; prev_cy <- circ$cy; prev_r <- circ$r
        }

        before <- bij_err(.nested_sector(crd, cxv, cyv, raw, k), gi, k)
        ref    <- .refine_radii(crd, gi, cxv, cyv, raw, k)
        after  <- bij_err(.nested_sector(crd, cxv, cyv, ref, k), gi, k)

        expect_lte(after, before)
        # refinement must keep the circles nested
        for (s in 2:(k - 1L)) {
            need <- sqrt((cxv[s] - cxv[s - 1L])^2 +
                         (cyv[s] - cyv[s - 1L])^2) + ref[s - 1L]
            expect_gte(ref[s], need - 1e-9)
        }
        if (after < before) improved <- improved + 1L
    }
    expect_gt(improved, 0L)
})

test_that("radialCircles() reports the error its own circles achieve", {
    set.seed(23)
    for (i in 1:40) {
        k   <- sample(2:4, 1)
        n   <- k * 6
        crd <- cbind(rnorm(n), rnorm(n))
        gi  <- rep(seq_len(k), each = 6)
        grp <- factor(letters[seq_len(k)][gi])
        for (fixed in c(TRUE, FALSE)) {
            res <- with_null_dev(
                if (fixed) radialCircles(crd, grp, cx = 0, cy = 0)
                else       radialCircles(crd, grp))
            sec <- .nested_sector(crd, res$cx, res$cy, res$radii, k)
            expect_equal(res$sector, sec)
            expect_equal(res$misclass, bij_err(sec, gi, k))
            expect_equal(res$misclass, nrow(res$misclass_points))
        }
    }
})

test_that("concentric radialCircles() is globally optimal over nested radii", {
    # Reference: brute force over all non-decreasing (k-1)-tuples drawn from
    # every radius that yields a distinct split.
    best_err <- function(d, gi, k) {
        cands <- c(0, sort(unique(d)), max(d) + 1)
        grid  <- as.matrix(expand.grid(rep(list(seq_along(cands)), k - 1L)))
        best  <- Inf
        for (i in seq_len(nrow(grid))) {
            rr <- cands[grid[i, ]]
            if (any(diff(rr) < 0)) next
            sec <- rep(k, length(d))
            for (s in (k - 1L):1L) sec[d <= rr[s]] <- s
            best <- min(best, bij_err(sec, gi, k))
        }
        best
    }
    set.seed(24)
    tested <- 0L
    for (i in 1:80) {
        k  <- sample(2:3, 1)
        n  <- sample(5:9, 1)
        d  <- round(runif(n, 0.2, 6), 1)
        gi <- sample(seq_len(k), n, replace = TRUE)
        if (length(unique(gi)) < k) next
        res <- with_null_dev(
            radialCircles(cbind(d, 0), factor(letters[seq_len(k)][gi]),
                          cx = 0, cy = 0))
        expect_equal(res$misclass, best_err(d, gi, k))
        tested <- tested + 1L
    }
    expect_gt(tested, 30L)
})

test_that("radialCircle() and radialCircles(k = 2) agree at the same center", {
    set.seed(25)
    for (i in 1:80) {
        n  <- sample(4:14, 1)
        d  <- round(runif(n, 0.2, 8), 1)
        gi <- sample(1:2, n, replace = TRUE)
        if (length(unique(gi)) < 2L) next
        grp <- factor(c("a", "b")[gi])
        one <- with_null_dev(radialCircle(cbind(d, 0),  grp, cx = 0, cy = 0))
        many <- with_null_dev(radialCircles(cbind(d, 0), grp, cx = 0, cy = 0))
        expect_equal(one$misclass, many$misclass)
        expect_equal(one$radius,   many$radii[1])
        expect_equal(one$sector,   many$sector)
        expect_equal(one$majority, many$majority)
    }
})

test_that("radialCircles() keeps its circles nested", {
    set.seed(26)
    for (i in 1:30) {
        k   <- sample(3:4, 1)
        n   <- k * 7
        crd <- cbind(rnorm(n), rnorm(n))
        grp <- factor(rep(letters[seq_len(k)], each = 7))
        res <- with_null_dev(radialCircles(crd, grp))
        for (s in 2:(k - 1L)) {
            need <- sqrt((res$cx[s] - res$cx[s - 1L])^2 +
                         (res$cy[s] - res$cy[s - 1L])^2) + res$radii[s - 1L]
            expect_gte(res$radii[s], need - 1e-9)
        }
    }
})

test_that("a degenerate partition is returned when it is the true minimum", {
    # 5 points at distances 1..5 labelled a a b a a: no circle separates the
    # single "b", so enclosing everything and calling it "a" (1 error) beats
    # every proper split (2 errors).
    crd <- cbind(1:5, 0)
    grp <- factor(c("a", "a", "b", "a", "a"))
    for (res in list(with_null_dev(radialCircle(crd, grp, cx = 0, cy = 0)),
                     with_null_dev(radialCircles(crd, grp, cx = 0, cy = 0)))) {
        expect_equal(res$misclass, 1L)
        expect_equal(length(unique(res$sector)), 1L)
    }
})

test_that("ties prefer a circle that actually splits the configuration", {
    # Two tight clusters, perfectly separable: the optimum is unique in
    # misclass but many radii achieve it, and the widest gap must win.
    d   <- c(1, 1.1, 1.2, 5, 5.1, 5.2)
    grp <- factor(c("a", "a", "a", "b", "b", "b"))
    res <- with_null_dev(radialCircle(cbind(d, 0), grp, cx = 0, cy = 0))
    expect_equal(res$misclass, 0L)
    expect_gt(res$radius, 1.2)
    expect_lt(res$radius, 5)
    # midpoint of the widest gap
    expect_equal(res$radius, (1.2 + 5) / 2)
})


# ---- Validation parity between the two functions ----

test_that("both functions reject NAs and malformed input", {
    crd <- cbind(c(1, 2, 3, 9), 0)
    grp <- factor(c("a", "a", "b", "b"))

    crd_na <- crd; crd_na[2, 1] <- NA
    grp_na <- grp; grp_na[2]    <- NA

    for (f in list(radialCircle, radialCircles)) {
        expect_error(with_null_dev(f(crd_na, grp, cx = 0, cy = 0)), "NAs.*crd")
        expect_error(with_null_dev(f(crd_na, grp)),                 "NAs.*crd")
        expect_error(with_null_dev(f(crd, grp_na, cx = 0, cy = 0)), "NAs.*group")
        expect_error(with_null_dev(f(crd, grp_na)),                 "NAs.*group")
        expect_error(with_null_dev(f(cbind(crd, 1), grp)),          "2 columns")
        expect_error(with_null_dev(f(crd, grp[-1])),                "must equal")
        expect_error(with_null_dev(f(crd, grp, .method = "BFGS")),  "not available")
    }
})

test_that("output = FALSE returns invisible(NULL)", {
    crd <- cbind(c(1, 2, 3, 9), 0)
    grp <- factor(c("a", "a", "b", "b"))
    expect_null(with_null_dev(radialCircle(crd, grp, cx = 0, cy = 0, output = FALSE)))
    expect_null(with_null_dev(radialCircles(crd, grp, cx = 0, cy = 0, output = FALSE)))
})


# ---- Validation parity across the whole partition family ----

test_that("every partition function rejects NAs with a clear message", {
    crd <- cbind(c(1, 2, 3, 9, 4, 7), c(0, 1, 0, 1, 2, 2))
    g2  <- factor(c("a", "a", "b", "b", "a", "b"))
    g3  <- factor(c("a", "a", "b", "b", "c", "c"))

    crd_na <- crd; crd_na[2, 1] <- NA

    binary <- list(radialCircle, radialEllipse, axialLine)
    k_way  <- list(radialCircles, radialEllipses, axialLines, angularPartition)

    for (f in binary) {
        g_na <- g2; g_na[2] <- NA
        expect_error(with_null_dev(f(crd_na, g2)), "NAs? allowed in crd")
        expect_error(with_null_dev(f(crd, g_na)),  "NAs? allowed in group")
    }
    for (f in k_way) {
        g_na <- g3; g_na[2] <- NA
        expect_error(with_null_dev(f(crd_na, g3)), "NAs? allowed in crd")
        expect_error(with_null_dev(f(crd, g_na)),  "NAs? allowed in group")
    }
    # ellipseInConfig takes no `group`
    expect_error(with_null_dev(ellipseInConfig(crd_na)), "NAs? allowed in crd")
})
