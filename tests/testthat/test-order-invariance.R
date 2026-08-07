# The partition functions must not depend on the factor level order.
#
# radialCircles() and radialEllipses() fit their boundaries sequentially, which
# needs an inside-to-outside ordering of the groups. Taking that from the
# factor level order made the result depend on how the groups happened to be
# named -- with R's default alphabetical levels, a configuration nested
# inner < mid < outer but labelled so that "mid" sorts first was fitted as the
# wrong nesting, and .assign_groups() cannot repair a nesting the geometry
# search never looked for. The ordering is now searched instead
# (.nesting_orders() in utils.R), so relabelling or reordering the levels
# cannot change the answer. These tests pin that down.

with_null_dev <- function(expr) {
    pdf(NULL)
    on.exit(dev.off(), add = TRUE)
    plot.new()
    plot.window(xlim = c(-10, 10), ylim = c(-10, 10))
    force(expr)
}

# Everything about a result that must be invariant: the geometry, the induced
# partition of the points into regions, the misclass count, and which group
# owns each region. Deliberately compares the partition as sets of point
# indices, not the raw `sector` codes, plus `majority` -- together those pin
# down the region-to-group mapping without assuming a labelling.
invariant_key <- function(res) {
    # every returned field that is geometry rather than labelling
    drop <- c("misclass", "misclass_points", "sector", "majority", "predicted")
    geom <- unlist(res[setdiff(names(res), drop)], use.names = TRUE)
    list(misclass = if (is.list(res$misclass)) res$misclass$n else res$misclass,
         geom     = if (length(geom)) signif(as.numeric(geom), 10) else numeric(0),
         part     = split(seq_along(res$sector), res$sector),
         majority = res$majority)
}

# Call `f` once per permutation of the level order and return TRUE when every
# result carries the same invariant key.
same_under_all_level_orders <- function(f, crd, lab, k, extra = list()) {
    perms <- gtools::permutations(k, k)
    uniq  <- sort(unique(lab))
    ref   <- NULL
    for (i in seq_len(nrow(perms))) {
        grp <- factor(lab, levels = uniq[perms[i, ]])
        key <- invariant_key(
            with_null_dev(do.call(f, c(list(crd, grp), extra))))
        if (is.null(ref)) ref <- key
        else if (!isTRUE(all.equal(ref, key))) return(FALSE)
    }
    TRUE
}

# k concentric rings with `sd` noise; group g sits on ring g.
rings <- function(k, per = 7L, seed = 1L, sd = 0.15) {
    set.seed(seed)
    n    <- k * per
    ring <- rep(seq_len(k), each = per)
    th   <- runif(n, 0, 2 * pi)
    list(crd = cbind(ring * cos(th), ring * sin(th)) +
               matrix(rnorm(2 * n, 0, sd), n),
         lab = letters[seq_len(k)][ring])
}


# ---- The motivating case ----

test_that("a nesting whose alphabetical levels are in the wrong order is found", {
    # 2 inner / 7 middle / 3 outer, labelled so that the *middle* group sorts
    # first alphabetically -- the shape of the FacetsGutt91 configuration.
    set.seed(8)
    inner <- cbind(rnorm(2, 0, 0.05), rnorm(2, 0, 0.05))
    th_m  <- seq(0, 2 * pi, length.out = 8)[-8]
    mid   <- cbind(cos(th_m), sin(th_m))
    th_o  <- seq(0, 2 * pi, length.out = 4)[-4]
    outer <- cbind(2.5 * cos(th_o), 2.5 * sin(th_o))
    crd   <- rbind(inner, mid, outer)
    lab   <- c(rep("infer", 2), rep("applFact", 7), rep("applSTM", 3))

    alpha <- factor(lab)
    expect_equal(levels(alpha), c("applFact", "applSTM", "infer"))

    circ <- with_null_dev(radialCircles(crd, alpha))
    expect_equal(circ$misclass, 0L)
    expect_equal(circ$majority, c("infer", "applFact", "applSTM"))

    ell <- with_null_dev(radialEllipses(crd, alpha))
    expect_equal(ell$misclass$n, 0L)
    expect_equal(ell$majority, c("infer", "applFact", "applSTM"))

    # and the alphabetical order gives exactly what the hand-set order gives
    fixed <- factor(lab, levels = c("infer", "applFact", "applSTM"))
    expect_equal(invariant_key(with_null_dev(radialCircles(crd, fixed))),
                 invariant_key(circ))
})


# ---- Invariance across every level order ----

test_that("radialCircles() is invariant under all k! level orders", {
    for (k in 2:4) {
        r <- rings(k, seed = 10L + k)
        expect_true(same_under_all_level_orders(radialCircles, r$crd, r$lab, k))
    }
})

test_that("radialEllipses() is invariant under all k! level orders", {
    for (k in 2:4) {
        r <- rings(k, seed = 20L + k)
        expect_true(same_under_all_level_orders(radialEllipses, r$crd, r$lab, k))
    }
})

test_that("concentric radialCircles() is invariant under all k! level orders", {
    for (k in 2:4) {
        r <- rings(k, seed = 30L + k)
        expect_true(same_under_all_level_orders(
            radialCircles, r$crd, r$lab, k, list(cx = 0, cy = 0)))
    }
})

test_that("radialEllipses(ellipse = ) is invariant under all k! level orders", {
    for (k in 3:4) {
        r <- rings(k, seed = 40L + k)
        expect_true(same_under_all_level_orders(
            radialEllipses, r$crd, r$lab, k,
            list(ellipse = c(0, 0, 1.4, 1.0, 0.3))))
    }
})

test_that("invariance holds with no radial structure at all", {
    # The searched ordering must be well defined even when no ordering is
    # "right": these configurations are pure noise, so ties abound and any
    # index-based tie-break would leak the level order.
    set.seed(50)
    for (i in 1:12) {
        k   <- sample(2:4, 1)
        n   <- k * 6
        crd <- cbind(rnorm(n), rnorm(n))
        lab <- letters[seq_len(k)][rep(seq_len(k), each = 6)]
        expect_true(same_under_all_level_orders(radialCircles,  crd, lab, k))
        expect_true(same_under_all_level_orders(radialEllipses, crd, lab, k))
    }
})

test_that("axialLines() and angularPartition() are invariant too", {
    # Neither uses level order to fit, but both report `majority` through
    # .assign_groups(), whose tie-break used to be index-based.
    set.seed(60)
    for (k in 2:3) {
        n   <- k * 7
        crd <- cbind(rnorm(n), rnorm(n))
        lab <- letters[seq_len(k)][rep(seq_len(k), each = 7)]
        expect_true(same_under_all_level_orders(
            axialLines, crd, lab, k, list(n_angles = 24L)))
        expect_true(same_under_all_level_orders(angularPartition, crd, lab, k))
    }
})


# ---- The pieces the invariance rests on ----

test_that(".nesting_orders() is invariant and puts the data-derived order first", {
    r <- rings(3L, seed = 70L)
    uniq <- sort(unique(r$lab))
    ref  <- NULL
    for (p in seq_len(6)) {
        perm <- gtools::permutations(3, 3)[p, ]
        grp  <- factor(r$lab, levels = uniq[perm])
        cand <- .nesting_orders(r$crd, grp)
        # translate index-based orders into level *names* to compare across
        # different level orders
        named <- lapply(cand, function(o) levels(grp)[o])
        if (is.null(ref)) ref <- named else expect_equal(named, ref)
        # candidate 1 is the mean-radial order, which here is a < b < c
        expect_equal(named[[1L]], c("a", "b", "c"))
        expect_length(cand, 6L)
        expect_true(attr(cand, "exhaustive"))
    }
})

test_that(".nesting_orders() falls back past max_perm and says so", {
    r    <- rings(3L, seed = 71L)
    grp  <- factor(r$lab, levels = sort(unique(r$lab)))
    cand <- .nesting_orders(r$crd, grp, max_perm = 2L)
    expect_length(cand, 1L)
    expect_false(attr(cand, "exhaustive"))
    # under the cap, nothing is dropped
    expect_true(attr(.nesting_orders(r$crd, grp), "exhaustive"))
})

test_that("radialCircles() warns rather than silently weakening past the cap", {
    # k = 6 is past max_perm = 120, so only the data-derived order is fitted.
    r   <- rings(6L, per = 5L, seed = 72L)
    grp <- factor(r$lab, levels = sort(unique(r$lab)))
    expect_warning(with_null_dev(radialCircles(r$crd, grp, output = FALSE)),
                   "too many group orderings")
    expect_warning(with_null_dev(radialEllipses(r$crd, grp, output = FALSE)),
                   "too many group orderings")
    # k = 5 is still exhaustive, so no warning
    r5   <- rings(5L, per = 5L, seed = 73L)
    grp5 <- factor(r5$lab, levels = sort(unique(r5$lab)))
    expect_silent(with_null_dev(radialCircles(r5$crd, grp5, output = FALSE)))
})

test_that(".assign_groups() breaks ties by level name, not by group index", {
    # A fully tied table: every bijection scores the same, so the tie-break
    # decides. Passing the names must give the same answer whichever order
    # the names arrive in.
    count_mat <- matrix(1L, 3L, 3L)
    a <- .assign_groups(count_mat, c("x", "y", "z"))
    b <- .assign_groups(count_mat, c("z", "y", "x"))
    expect_equal(c("x", "y", "z")[a], c("z", "y", "x")[b])

    # without names the historical first-in-permutation-order winner is kept
    expect_equal(.assign_groups(count_mat), c(1L, 2L, 3L))

    # a clear winner is unaffected by the tie-break
    cm <- matrix(c(9L, 0L, 0L,
                   0L, 9L, 0L,
                   0L, 0L, 9L), 3L, 3L, byrow = TRUE)
    expect_equal(.assign_groups(cm), c(1L, 2L, 3L))
    expect_equal(.assign_groups(cm, c("p", "q", "r")), c(1L, 2L, 3L))
})
