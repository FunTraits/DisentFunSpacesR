# The faster implementations must give the same results as the original ones.

test_that(".procuste_rtest reproduces ade4::procuste.rtest exactly", {
  set.seed(11)
  X <- matrix(rnorm(300), 100); Y <- matrix(rnorm(200), 100)  # unrelated: p not minimal
  rownames(X) <- rownames(Y) <- paste0("s", 1:100)
  set.seed(5); a <- ade4::procuste.rtest(as.data.frame(X), as.data.frame(Y), nrepet = 199)
  set.seed(5); b <- DisentFunSpaces:::.procuste_rtest(X, Y, nrepet = 199)
  expect_equal(b$obs, a$obs, tolerance = 1e-12)
  expect_equal(b$sim, a$sim, tolerance = 1e-12)
  expect_identical(b$pvalue, a$pvalue)
  expect_gt(a$pvalue, 0.005)
  skip_on_os("windows")
  op <- options(DisentFunSpaces.parallel_min_axes = 1); on.exit(options(op))
  set.seed(5); c2 <- DisentFunSpaces:::.procuste_rtest(X, Y, nrepet = 199, cores = 2)
  expect_equal(c2$sim, a$sim, tolerance = 1e-12)
  expect_identical(c2$pvalue, a$pvalue)
})

test_that("compare_spaces gives the same results with several cores", {
  skip_on_os("windows")
  sp <- build_spaces(ex(), envfit_perm = 0)
  a <- compare_spaces(sp, nrepet = 49, n_axes = 3)
  b <- compare_spaces(sp, nrepet = 49, n_axes = 3, cores = 2)
  expect_identical(a$table$p_value, b$table$p_value)
  expect_equal(a$table$r, b$table$r)
})

test_that("vectorised CorrelType equals one cor()/lm() per trait and axis", {
  ref <- function(pcoa, type) {
    vec <- pcoa$vectors; type <- as.data.frame(type)
    A <- matrix(NA_real_, ncol(type), ncol(vec))
    for (i in seq_len(ncol(type))) for (j in seq_len(ncol(vec))) {
      v <- type[[i]]; y <- vec[, j]
      A[i, j] <- if (is.factor(v) || is.character(v)) {
        v <- factor(v); ok <- !is.na(v)
        if (nlevels(droplevels(v[ok])) < 2) NA else sqrt(summary(stats::lm(y[ok] ~ v[ok]))$r.squared)
      } else if (isTRUE(stats::sd(v, na.rm = TRUE) > 0)) stats::cor(y, v, use = "pairwise.complete.obs") else NA
    }
    A
  }
  ts <- sim_tr(missing = 0.1)
  sp <- suppressWarnings(build_space(ts, envfit_perm = 0))
  for (tab in sp$blocks$x) {
    expect_equal(unname(CorrelType(sp$PCoA, tab)), ref(sp$PCoA, tab), tolerance = 1e-10)
  }
})

test_that("partial PCoA gives the same leading axes as the full PCoA", {
  skip_if_not_installed("RSpectra")
  tr <- ex()
  full <- build_space(tr, max_axes = 6, envfit_perm = 0)
  part <- build_space(tr, max_axes = 6, envfit_perm = 0, eigen = "partial")
  expect_true(isTRUE(part$PCoA$partial))
  expect_equal(abs(unname(part$coords)), abs(unname(full$coords)), tolerance = 1e-8)
  expect_equal(part$eigen$relative, full$eigen$relative, tolerance = 1e-10)
  expect_true(is.na(part$negative_eigen$n))
  # default: full decomposition kept
  expect_null(full$PCoA$partial)
  # all positive axes (max_axes = NULL): same axes as the full decomposition
  fa <- build_space(tr, c("Diet", "Locomotion"), envfit_perm = 0)
  pa <- build_space(tr, c("Diet", "Locomotion"), envfit_perm = 0, eigen = "partial")
  expect_true(isTRUE(pa$PCoA$partial))
  expect_equal(ncol(pa$coords), ncol(fa$coords))
  expect_equal(abs(unname(pa$coords)), abs(unname(fa$coords)), tolerance = 1e-8)
  expect_equal(pa$eigen$relative, fa$eigen$relative, tolerance = 1e-10)
  # more positive axes than the cap: falls back to the full decomposition
  op <- options(DisentFunSpaces.partial_k_max = 3); on.exit(options(op))
  expect_message(fb <- build_space(tr, c("Diet", "Locomotion"), envfit_perm = 0, eigen = "partial"),
                 "full decomposition")
  expect_null(fb$PCoA$partial)
  expect_equal(fb$coords, fa$coords)
})

test_that("streamed accumulation of distances handles NA and weights", {
  set.seed(3)
  ds <- list(c(0.1, NA, 0.3, 0.5), c(0.2, 0.4, NA, 0.1), c(NA, NA, 0.2, 0.9))
  w <- c(1, 0.5, 2)
  num <- 0; den <- 0
  for (i in 1:3) { ok <- !is.na(ds[[i]]); d <- ds[[i]]; d[!ok] <- 0
    num <- num + w[i] * d * d; den <- den + w[i] * ok }
  expect_identical(DisentFunSpaces:::.wmean_sq_stream(1:3, function(i) ds[[i]], w), num / den)
  ds2 <- list(c(0.1, 0.2), c(0.3, 0.4))
  expect_identical(DisentFunSpaces:::.wmean_sq_stream(1:2, function(i) ds2[[i]], c(1, 1)),
                   (ds2[[1]]^2 + ds2[[2]]^2) / 2)
})

test_that("assemblage metrics and hull tests match the matrix-based versions", {
  set.seed(4)
  X <- matrix(rnorm(400), 100)
  X <- rbind(X, X[1:5, ])                       # duplicated positions
  D <- as.matrix(dist(X)); diag(D) <- NA
  m <- DisentFunSpaces:::.assemblage_metrics(X, 3, 1)
  expect_equal(unname(m["FOri"]), mean(apply(D, 1, min, na.rm = TRUE)), tolerance = 1e-12)
  expect_equal(unname(m["MPD"]), mean(D, na.rm = TRUE), tolerance = 1e-12)
  expect_equal(DisentFunSpaces:::.mean_pairwise_dist(X, max_n = 10, block = 17),
               mean(D, na.rm = TRUE), tolerance = 1e-10)
  # hull: half-space test agrees with the volume criterion
  Xn <- X[1:80, 1:3]; P <- rbind(X[81:105, 1:3], Xn[chull(Xn[, 1:2])[1:3], ])
  vol <- vapply(seq_len(nrow(P)), function(i) {
    v0 <- DisentFunSpaces:::.hull_vol(Xn); v1 <- DisentFunSpaces:::.hull_vol(rbind(Xn, P[i, , drop = FALSE]))
    v1 <= v0 * (1 + 1e-7)
  }, logical(1))
  expect_identical(DisentFunSpaces:::.in_hull(Xn, P), vol)
})

test_that("permute null model: cached distances give the same results", {
  sp <- build_space(ex(), envfit_perm = 0)
  st <- simulate_introductions(sp, prop = 0.15, bias = "peripheral", seed = 2)
  a <- nonnative_assemblage(sp, st, n_axes = 6, hull_axes = 3, null = "permute", nrep = 29, seed = 3)
  op <- options(DisentFunSpaces.max_n_dist_matrix = 0); on.exit(options(op))
  b <- nonnative_assemblage(sp, st, n_axes = 6, hull_axes = 3, null = "permute", nrep = 29, seed = 3)
  expect_equal(a, b, tolerance = 1e-10)
  # the cached replicate metrics equal the direct computation
  X <- sp$coords[, 1:6]; pc <- DisentFunSpaces:::.permute_cache(X)
  R <- c(3, 10, 11, 50, 200)
  expect_equal(DisentFunSpaces:::.permute_metrics(pc, R, X[-R, ], 3, 1),
               DisentFunSpaces:::.assemblage_metrics(X[-R, ], 3, 1), tolerance = 1e-12)
})
