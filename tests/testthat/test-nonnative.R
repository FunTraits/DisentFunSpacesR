nn_setup <- function() {
  sp <- build_spaces(ex(), spaces = list(Diet = "Diet", Loco = "Locomotion",
                                          Combined = c("Diet", "Locomotion")))
  st <- simulate_introductions(sp$Combined, prop = 0.15, bias = "peripheral", seed = 1)
  list(sp = sp, st = st)
}

test_that("simulate_introductions assigns the requested proportion", {
  s <- nn_setup()
  expect_equal(sum(s$st == "non-native"), round(0.15 * length(s$st)))
  expect_setequal(names(s$st), rownames(s$sp$Combined$coords))
})

test_that("non-native species metrics follow their definitions", {
  s <- nn_setup()
  imp <- stats::setNames(rep(c("MC", "MV"), length.out = sum(s$st != "native")),
                         names(s$st)[s$st != "native"])
  r <- nonnative_species(s$sp$Combined, s$st, k = 3, impact = imp)
  expect_s3_class(r, "fs_nonnative")
  X <- s$sp$Combined$coords
  nat <- names(s$st)[s$st == "native"]; nn <- names(s$st)[s$st != "native"]
  D <- as.matrix(stats::dist(X))
  i <- r$nonnative$species[1]
  expect_equal(r$nonnative$FNov_raw[1], min(D[i, nat]), tolerance = 1e-8)
  expect_equal(r$nonnative$mean_dist_native[1], mean(D[i, nat]), tolerance = 1e-8)
  expect_equal(r$nonnative$FSp_native_raw[1], sqrt(sum((X[i, ] - colMeans(X[nat, ]))^2)), tolerance = 1e-8)
  j <- r$native$species[2]
  expect_equal(r$native$FUn_before_raw[2], sum(sort(D[j, setdiff(nat, j)])[1:3]), tolerance = 1e-8)
  expect_equal(r$native$FUn_after_raw[2], sum(sort(D[j, setdiff(rownames(X), j)])[1:3]), tolerance = 1e-8)
  expect_true(all(r$native$dFUn_rel <= 1e-12))
  expect_equal(r$nonnative$FUSI, log(1 + r$nonnative$W * r$nonnative$FNov +
                                       r$nonnative$W * r$nonnative$FSp_native))
  expect_equal(sort(unique(r$nonnative$W)), c(0.2, 1))
  ## a non-native species identical to a native is never outside the native hull
  d <- nonnative_species(build_spaces(ex(), spaces = list(D = "Diet"))$D, s$st, hull_axes = 3)
  expect_true(all(d$nonnative$inside_native_hull[d$nonnative$FNov_raw == 0] %in% c(TRUE, NA)))
  ## logical status gives the same result
  r2 <- nonnative_species(s$sp$Combined, s$st != "native", k = 3)
  expect_equal(r2$nonnative$FNov_raw, r$nonnative$FNov_raw)
})

test_that("assemblage-level change and null models", {
  s <- nn_setup()
  a <- nonnative_assemblage(s$sp$Combined, s$st, n_axes = 4, hull_axes = 3,
                            null = "permute", nrep = 19, seed = 1)
  expect_equal(nrow(a), 1)
  expect_equal(a$S_cur - a$S_hist, sum(s$st != "native"))
  expect_true(a$FRic_cur >= a$FRic_hist)
  expect_true(all(c("ses_FRic", "p_FRic", "ses_FDis") %in% names(a)))
  expect_warning(nonnative_assemblage(s$sp$Combined, s$st, n_axes = 3, hull_axes = 3,
                                      null = "pool", nrep = 3), "uninformative")
  ## several assemblages, extirpations and per-site status
  sp <- rownames(s$sp$Combined$coords)
  set.seed(2)
  comm <- matrix(stats::rbinom(5 * length(sp), 1, 0.4), 5, dimnames = list(paste0("s", 1:5), sp))
  ext <- names(s$st)[s$st == "native"][1:10]
  b <- nonnative_assemblage(s$sp$Combined, s$st, comm = comm, extirpated = ext,
                            n_axes = 4, hull_axes = 3, null = "pool", nrep = 9, seed = 1)
  expect_equal(nrow(b), 5)
  expect_equal(b$S_cur, b$n_native - b$n_extirpated + b$n_nonnative)
  stm <- matrix("native", 5, length(sp), dimnames = dimnames(comm))
  stm[1, 1:20] <- "introduced"
  c3 <- nonnative_assemblage(s$sp$Combined, stm, comm = comm, n_axes = 4, hull_axes = 3)
  expect_equal(c3$n_nonnative[1], sum(comm[1, 1:20]))
  expect_equal(c3$n_nonnative[2], 0)
})

test_that("nonnative_all, summaries and plots", {
  s <- nn_setup()
  x <- nonnative_all(s$sp, s$st, n_axes = 4, hull_axes = 3, null = "permute", nrep = 9, seed = 1)
  expect_s3_class(x, "fs_nonnative_all")
  expect_setequal(unique(x$nonnative$space), names(s$sp))
  sm <- summary_nonnative(x)
  expect_equal(nrow(sm), 3)
  expect_true(all(sm$p_outside_hull >= 0 & sm$p_outside_hull <= 1))
  expect_equal(nrow(compare_nonnative(x, "FNov")$pairs), 3)
  expect_equal(nrow(compare_nonnative(x, "dFUn_rel")$pairs), 3)
  expect_s3_class(plot_nonnative(s$sp$Combined, s$st, label = 3), "ggplot")
  expect_s3_class(plot_nonnative_change(x, "ses_FRic"), "ggplot")
})
