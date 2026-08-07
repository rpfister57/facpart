# test angular.R functions

# ---- Helpers --------------------------------------------------------------

# angularPartition() always draws: it reads par("usr") and calls segments(),
# so it needs an active device with an existing plot. Running it bare under
# Rscript silently opens the default device and leaves a stray Rplots.pdf in
# the working directory, so every call that gets past input validation goes
# through here. `code` is a promise and is forced only after plot().
with_plot <- function(crd, code) {
    pdf(NULL)
    on.exit(dev.off(), add = TRUE)
    plot(crd, asp = 1)
    force(code)
}


# k disjoint angular wedges with clear gaps between them: perfectly
# separable, so the optimal partition has zero misclassification.
wedge_config <- function(k = 3L, per = 6L) {
    centers <- 2 * pi * (seq_len(k) - 1L) / k
    half    <- pi / (2 * k)
    ang     <- unlist(lapply(centers,
                             function(c0) seq(c0 - half, c0 + half,
                                              length.out = per)))
    rad     <- rep(seq(0.6, 1.4, length.out = per), k)
    list(crd = cbind(rad * cos(ang), rad * sin(ang)),
         grp = factor(rep(letters[seq_len(k)], each = per)))
}


# Regular n-gon on the unit circle; labels recycled over the vertices.
ngon_config <- function(n = 12L, labels = c("a", "b")) {
    ang <- 2 * pi * (seq_len(n) - 1L) / n
    list(crd = cbind(cos(ang), sin(ang)),
         grp = factor(rep(labels, length.out = n)))
}


# Independent reference implementation of the search criterion: enumerate
# every choice of k cut-gaps and every bijective arc-to-group assignment,
# keep the best. Deliberately the slow, obvious version -- angularPartition()
# gets the same answer via a bound-and-prune scan over a vectorised candidate
# axis, and this is what pins that optimisation to the right quantity.
ref_angular_misclass <- function(crd, group, cx, cy) {
    group  <- factor(group)
    k      <- nlevels(group)
    n      <- nrow(crd)
    g      <- as.integer(group)[order(atan2(crd[, 2] - cy, crd[, 1] - cx))]
    perms  <- gtools::permutations(k, k)
    combos <- combn(n, k)
    best   <- -1L
    for (j in seq_len(ncol(combos))) {
        cuts <- combos[, j]
        arcs <- vector("list", k)
        for (s in seq_len(k)) {
            # interior arc s spans sorted positions cuts[s]+1 .. cuts[s+1];
            # arc k wraps cuts[k]+1 .. n, 1 .. cuts[1]
            if (s < k) {
                idx <- if (cuts[s] + 1L <= cuts[s + 1L])
                    seq(cuts[s] + 1L, cuts[s + 1L]) else integer(0)
            } else {
                idx <- c(if (cuts[k] + 1L <= n) seq(cuts[k] + 1L, n) else integer(0),
                         seq_len(cuts[1L]))
            }
            arcs[[s]] <- tabulate(g[idx], nbins = k)
        }
        for (p in seq_len(nrow(perms))) {
            correct <- 0L
            for (s in seq_len(k)) correct <- correct + arcs[[s]][perms[p, s]]
            if (correct > best) best <- correct
        }
    }
    n - best
}


# ---- Result shape ---------------------------------------------------------

test_that("angularPartition returns the documented result shape", {
    cfg <- wedge_config()
    res <- with_plot(cfg$crd, angularPartition(cfg$crd, cfg$grp))

    expect_named(res, c("cuts", "margin", "misclass", "misclass_points",
                        "sector", "majority", "center", "pt_angles"))
    expect_length(res$cuts, 3L)              # k cut angles
    expect_length(res$margin, 1L)
    expect_length(res$center, 2L)
    expect_length(res$sector, nrow(cfg$crd))
    expect_length(res$pt_angles, nrow(cfg$crd))
    expect_length(res$majority, 3L)

    expect_type(res$misclass, "integer")     # integer, not the list that
    expect_length(res$misclass, 1L)          # radialEllipse() returns
    expect_type(res$majority, "character")
    expect_true(is.numeric(res$cuts))
    expect_true(res$margin > 0)

    # angles are atan2 output, i.e. in (-pi, pi]
    expect_true(all(res$pt_angles > -pi & res$pt_angles <= pi))
    expect_true(all(res$cuts > -pi & res$cuts <= pi))
})


test_that("misclass_points honours the highlightMisclass() contract", {
    cfg <- ngon_config(12L, c("a", "b"))
    res <- with_plot(cfg$crd, angularPartition(cfg$crd, cfg$grp))

    expect_s3_class(res$misclass_points, "data.frame")
    expect_named(res$misclass_points, c("x", "y", "label"))
    expect_identical(res$misclass, nrow(res$misclass_points))
    expect_gt(res$misclass, 0L)              # alternating labels: not separable

    # the listed points are exactly the ones whose group differs from the
    # majority group of their sector, with their original coordinates
    bad <- which(as.character(cfg$grp) != res$majority[res$sector])
    expect_identical(res$misclass, length(bad))
    expect_equal(res$misclass_points$x, cfg$crd[bad, 1])
    expect_equal(res$misclass_points$y, cfg$crd[bad, 2])
    expect_identical(as.character(res$misclass_points$label),
                     as.character(cfg$grp[bad]))

    expect_null(with_plot(cfg$crd, highlightMisclass(res)))
})


# ---- Correctness of the search -------------------------------------------

test_that("clean wedges are separated perfectly", {
    cfg <- wedge_config(k = 3L, per = 6L)
    res <- with_plot(cfg$crd, angularPartition(cfg$crd, cfg$grp))

    expect_identical(res$misclass, 0L)
    expect_equal(nrow(res$misclass_points), 0L)

    # the wedge gaps are centred on -pi/3, pi/3 and pi, and the nearest
    # point sits at radius 0.6 across a pi/3 gap: 0.6 * sin(pi/6) = 0.3
    expect_equal(sort(res$cuts), c(-pi / 3, pi / 3, pi))
    expect_equal(res$margin, 0.3)
    expect_equal(res$center, c(0, 0), tolerance = 1e-6)

    # each group owns exactly one sector, and vice versa
    tab <- table(res$sector, cfg$grp)
    expect_identical(sum(tab > 0), 3L)
    expect_true(all(rowSums(tab > 0) == 1L))
    expect_true(all(colSums(tab > 0) == 1L))
    expect_setequal(res$majority, levels(cfg$grp))
})


test_that("misclass matches an independent brute-force optimum", {
    # fixed centre, so both sides solve exactly the same problem and the
    # centre optimisation is out of the picture
    cases <- list(
        alt2 = ngon_config(12L, c("a", "b")),
        alt3 = ngon_config(12L, c("a", "b", "c")),
        alt4 = ngon_config(12L, c("a", "b", "c", "d")),
        blk3 = list(crd = ngon_config(12L)$crd,
                    grp = factor(rep(c("a", "b", "c"), each = 4L)))
    )
    for (nm in names(cases)) {
        cfg <- cases[[nm]]
        for (ctr in list(c(0, 0), c(0.2, -0.1))) {
            res <- with_plot(cfg$crd,
                             angularPartition(cfg$crd, cfg$grp,
                                              cx = ctr[1], cy = ctr[2]))
            expect_identical(res$misclass,
                             ref_angular_misclass(cfg$crd, cfg$grp,
                                                  ctr[1], ctr[2]),
                             info = paste(nm, "at centre",
                                          paste(ctr, collapse = ",")))
        }
    }
})


test_that("sector assignment follows from cuts and pt_angles", {
    cfg <- ngon_config(14L, c("a", "b", "c"))
    res <- with_plot(cfg$crd, angularPartition(cfg$crd, cfg$grp))
    k   <- nlevels(cfg$grp)

    expect_true(all(res$sector %in% seq_len(k)))

    # recompute the sector of every point from the returned geometry alone
    sc  <- sort(res$cuts %% (2 * pi))
    fi  <- findInterval(res$pt_angles %% (2 * pi), sc)
    expect_identical(as.integer(ifelse(fi == 0L | fi == k, k, fi)),
                     as.integer(res$sector))

    # pt_angles are measured from the returned centre
    expect_equal(res$pt_angles,
                 atan2(cfg$crd[, 2] - res$center[2],
                       cfg$crd[, 1] - res$center[1]))
})


test_that("results are reported in original row order", {
    cfg <- wedge_config(k = 3L, per = 6L)
    prm <- c(13:18, 1:6, 7:12)

    res <- with_plot(cfg$crd, angularPartition(cfg$crd, cfg$grp))
    shf <- with_plot(cfg$crd[prm, ],
                     angularPartition(cfg$crd[prm, ], cfg$grp[prm]))

    expect_identical(shf$sector, res$sector[prm])
    expect_identical(shf$misclass, res$misclass)
    expect_equal(shf$margin, res$margin)
})


test_that("repeated calls are deterministic", {
    cfg <- wedge_config(k = 3L, per = 5L)
    a <- with_plot(cfg$crd, angularPartition(cfg$crd, cfg$grp))
    b <- with_plot(cfg$crd, angularPartition(cfg$crd, cfg$grp))
    expect_identical(a, b)
})


test_that("golden values for a deterministic non-separable case", {
    # regression guard: alternating labels on a regular 12-gon at a fixed
    # centre. Any change to the search must reproduce these bit-for-bit.
    cfg <- ngon_config(12L, c("a", "b"))
    res <- with_plot(cfg$crd,
                     angularPartition(cfg$crd, cfg$grp, cx = 0, cy = 0))

    expect_identical(res$misclass, 5L)
    expect_equal(res$cuts, c(-1.83259571459404658, -0.26179938779914969))
    expect_equal(res$margin, sin(pi / 12))   # unit radius, pi/6 gap
    expect_identical(res$sector, c(2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L,
                                   1L, 1L, 1L))
    expect_identical(res$majority, c("b", "a"))

    cfg3 <- ngon_config(12L, c("a", "b", "c"))
    res3 <- with_plot(cfg3$crd,
                      angularPartition(cfg3$crd, cfg3$grp, cx = 0, cy = 0))

    expect_identical(res3$misclass, 6L)
    expect_equal(res3$cuts, c(-1.8325957145940466, -1.3089969389957472,
                              2.3561944901923448))
    expect_identical(res3$sector, c(3L, 3L, 3L, 3L, 3L, 1L, 1L, 1L, 1L,
                                    2L, 3L, 3L))
    expect_identical(res3$majority, c("c", "a", "b"))
})


# ---- Arguments ------------------------------------------------------------

test_that("supplied cx/cy are used verbatim", {
    cfg <- wedge_config()
    res <- with_plot(cfg$crd,
                     angularPartition(cfg$crd, cfg$grp, cx = 0.05, cy = -0.05))
    expect_identical(res$center, c(0.05, -0.05))

    # only one of the two is enough to trigger the centre optimisation
    half <- with_plot(cfg$crd, angularPartition(cfg$crd, cfg$grp, cx = 0.05))
    expect_false(isTRUE(all.equal(half$center, c(0.05, -0.05))))
})


test_that("output = FALSE returns invisible NULL", {
    cfg <- wedge_config()
    expect_null(with_plot(cfg$crd,
                          angularPartition(cfg$crd, cfg$grp, output = FALSE)))
})


test_that("group may be a character vector", {
    cfg <- wedge_config()
    keys <- c("cuts", "margin", "misclass", "sector", "majority", "center")
    fac <- with_plot(cfg$crd, angularPartition(cfg$crd, cfg$grp))
    chr <- with_plot(cfg$crd,
                     angularPartition(cfg$crd, as.character(cfg$grp)))
    expect_identical(chr[keys], fac[keys])
})


test_that("fill and add do not change the numeric result", {
    cfg <- wedge_config()
    keys <- c("cuts", "margin", "misclass", "sector", "majority", "center")
    plain <- with_plot(cfg$crd, angularPartition(cfg$crd, cfg$grp))

    filled <- with_plot(cfg$crd,
                        angularPartition(cfg$crd, cfg$grp, fill = TRUE))
    expect_identical(filled[keys], plain[keys])

    # add = FALSE draws its own plot, so it needs a device but no plot
    pdf(NULL)
    on.exit(dev.off(), add = TRUE)
    fresh <- angularPartition(cfg$crd, cfg$grp, add = FALSE)
    expect_identical(fresh[keys], plain[keys])
})


test_that("k = 2 and a data frame crd are handled", {
    cfg <- wedge_config(k = 2L, per = 8L)
    res <- with_plot(cfg$crd, angularPartition(cfg$crd, cfg$grp))
    expect_identical(res$misclass, 0L)
    expect_length(res$cuts, 2L)

    df  <- as.data.frame(cfg$crd)
    dfr <- with_plot(df, angularPartition(df, cfg$grp))
    expect_identical(dfr$misclass, res$misclass)
    expect_identical(dfr$sector, res$sector)
})


# ---- Input validation ----------------------------------------------------

test_that("angularPartition rejects invalid input", {
    crd <- cbind(c(1, 0, -1, 0, 0.7, -0.7), c(0, 1, 0, -1, 0.7, -0.7))
    grp <- factor(c("a", "b", "a", "b", "a", "b"))

    # validation runs before anything is drawn, so no device is needed here
    expect_error(angularPartition(cbind(c("a", "b"), c("c", "d")),
                                  factor(c("x", "y"))),
                 "must be numeric", fixed = TRUE)
    expect_error(angularPartition(rbind(crd, c(NA, 0)),
                                  factor(c(as.character(grp), "a"))),
                 "No NAs allowed in crd!", fixed = TRUE)
    expect_error(angularPartition(crd,
                                  factor(c("a", "b", "a", "b", "a", NA))),
                 "No NAs allowed in group!", fixed = TRUE)
    expect_error(angularPartition(1:6, grp),
                 "must have two dimensions", fixed = TRUE)
    expect_error(angularPartition(matrix(1:6, ncol = 1L), grp),
                 "Coordinates must have 2 columns!", fixed = TRUE)
    expect_error(angularPartition(cbind(1:6, 1:6, 1:6), grp),
                 "Coordinates must have 2 columns!", fixed = TRUE)
    expect_error(angularPartition(crd, factor(c("a", "b", "a"))),
                 "nrow(crd) must equal length(group)!", fixed = TRUE)
    expect_error(angularPartition(crd, factor(rep("a", 6L))),
                 "group must have at least 2 levels!", fixed = TRUE)
})
