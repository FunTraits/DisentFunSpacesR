grp_setup <- function() {
  sp <- build_spaces(ex(), spaces = list(Diet = "Diet", Loco = "Locomotion",
                                          Combined = c("Diet", "Locomotion")))
  s <- rownames(sp$Combined$coords)
  ord <- stats::setNames(birds_traits$order, birds_traits$species)[s]
  iucn <- stats::setNames(birds_traits$iucn, birds_traits$species)[s]
  list(sp = sp, ord = ord, iucn = iucn,
       rec = list(threatened = c("VU", "EN", "CR"), other = c("LC", "NT")))
}

test_that("group labels are parsed: exclusive, multi-label, matrix, recoding", {
  sp <- c("a", "b", "c", "d")
  g <- DisentFunSpaces:::.as_groups(c(a = "x", b = "y", c = NA, d = "x"), sp)
  expect_true(g$exclusive); expect_equal(colSums(g$M), c(x = 2, y = 1))
  g <- DisentFunSpaces:::.as_groups(c(a = "food; aquarium", b = "food", c = "", d = "pet"), sp, sep = ";")
  expect_false(g$exclusive); expect_equal(colSums(g$M)[["food"]], 2); expect_false(any(g$M["c", ]))
  m <- matrix(c(1, 0, 1, 0, 1, 1, 0, 0), 4, dimnames = list(sp, c("u1", "u2")))
  g <- DisentFunSpaces:::.as_groups(m, sp)
  expect_equal(colSums(g$M), c(u1 = 2, u2 = 2)); expect_false(g$exclusive)
  g <- DisentFunSpaces:::.as_groups(c(a = "LC", b = "EN", c = "CR", d = "DD"), sp,
                                    recode = list(threatened = c("EN", "CR")), levels = c("threatened", "LC"))
  expect_equal(colSums(g$M), c(threatened = 2, LC = 1))
  expect_equal(DisentFunSpaces:::.parse_recode("thr = VU, EN\nother = LC;NT\n bad line"),
               list(thr = c("VU", "EN"), other = c("LC", "NT")))
  expect_null(DisentFunSpaces:::.parse_recode(""))
})

test_that("group_overlap follows its definitions and the global tests match vegan", {
  s <- grp_setup(); lv <- c("Anseriformes", "Accipitriformes", "Strigiformes")
  g <- group_overlap(s$sp$Combined, s$ord, levels = lv, hull_axes = 3, nperm = 99, seed = 1)
  expect_s3_class(g, "fs_groups")
  expect_equal(nrow(g$pairs), 3); expect_equal(g$groups$n, as.vector(table(s$ord)[lv]))
  X <- s$sp$Combined$coords; i <- s$ord %in% lv
  ## PERMANOVA pseudo-F and R2 = vegan::adonis2 with Euclidean distances
  a <- vegan::adonis2(stats::dist(X[i, ]) ~ s$ord[i], permutations = 0)
  expect_equal(g$test$permanova_R2, a$R2[1], tolerance = 1e-8)
  expect_equal(g$test$permanova_F, a$F[1], tolerance = 1e-8)
  bd <- stats::anova(vegan::betadisper(stats::dist(X[i, ]), factor(s$ord[i]), type = "centroid"))
  expect_equal(g$test$permdisp_F, bd$`F value`[1], tolerance = 1e-6)
  ## overlap of a group with itself is 1; Jaccard components add up
  p <- g$pairs
  expect_equal(p$beta_jac, p$beta_jtu + p$beta_jne)
  expect_true(all(p$overlap >= 0 & p$overlap <= 1))
  A <- X[s$ord == "Anseriformes", 1:3]
  expect_equal(DisentFunSpaces:::.hull_inter(A, A), DisentFunSpaces:::.hull_vol(A), tolerance = 1e-6)
  ## centroid distance
  B <- X[s$ord == "Accipitriformes", ]
  expect_equal(p$centroid_dist[p$group1 == "Anseriformes" & p$group2 == "Accipitriformes"],
               sqrt(sum((colMeans(X[s$ord == "Anseriformes", ]) - colMeans(B))^2)))
  ## null model adds SES; non-exclusive groups have no global test
  g2 <- group_overlap(s$sp$Diet, s$ord, levels = lv[1:2], hull_axes = 3, null = "permute", nrep = 9,
                      nperm = 0, seed = 1)
  expect_true(all(c("ses_overlap", "p_overlap", "ses_centroid_dist") %in% names(g2$pairs)))
  multi <- stats::setNames(ifelse(s$ord %in% lv[1:2], paste(s$ord, "big", sep = ";"), "small"), names(s$ord))
  g3 <- group_overlap(s$sp$Diet, multi, sep = ";", hull_axes = 3, nperm = 0)
  expect_false(g3$settings$exclusive); expect_true(is.na(g3$test$permanova_R2))
})

test_that("group_species generalises nonnative_species", {
  s <- grp_setup()
  st <- simulate_introductions(s$sp$Combined, prop = 0.15, seed = 1)
  a <- nonnative_species(s$sp$Combined, st, k = 3)
  b <- group_species(s$sp$Combined, st, reference = "native", k = 3)
  expect_equal(b$focal$FNov_raw, a$nonnative$FNov_raw)
  expect_equal(b$focal$inside_ref_hull, a$nonnative$inside_native_hull)
  expect_equal(b$reference$dFUn_rel, a$native$dFUn_rel)
  ## several focal groups, each against the reference; pooled
  r <- group_species(s$sp$Combined, s$iucn, recode = s$rec, reference = "other")
  expect_equal(unique(r$focal$group), c("DD", "threatened")[c("DD", "threatened") %in% r$focal$group])
  r2 <- group_species(s$sp$Combined, s$iucn, recode = s$rec, levels = c("threatened", "other"),
                      reference = "other")
  expect_equal(unique(r2$focal$group), "threatened")
  expect_equal(nrow(r2$focal), sum(s$iucn %in% c("VU", "EN", "CR")))
  ## weights: FUSW built as FUSI
  w <- stats::setNames(rep(0.5, nrow(r2$focal)), r2$focal$species)
  r3 <- group_species(s$sp$Combined, s$iucn, recode = s$rec, levels = c("threatened", "other"),
                      reference = "other", weight = w)
  expect_equal(r3$focal$FUSW, log(1 + 0.5 * r3$focal$FNov + 0.5 * r3$focal$FSp_ref))
  cc <- group_contribution(s$sp$Combined, s$iucn, recode = s$rec, levels = c("threatened", "other"),
                           reference = "other", hull_axes = 3)
  expect_equal(cc$S_all - cc$S_ref, cc$n_focal)
})

test_that("groups_all, summaries and plots", {
  s <- grp_setup()
  x <- groups_all(s$sp, s$iucn, recode = s$rec, levels = c("threatened", "other"), reference = "other",
                  hull_axes = 3, nperm = 19, null = "permute", nrep = 9, seed = 1)
  expect_s3_class(x, "fs_groups_all")
  expect_setequal(unique(x$pairs$space), names(s$sp))
  sg <- summary_groups(x); expect_equal(nrow(sg), 3)
  sf <- summary_groups(x, "focal"); expect_true(all(c("FNov_rel", "FRic_contribution_rel") %in% names(sf)))
  expect_s3_class(plot_groups(s$sp$Combined, s$ord, levels = c("Anseriformes", "Strigiformes")), "ggplot")
  expect_s3_class(plot_groups(s$sp$Combined, s$iucn, recode = s$rec, levels = c("threatened", "other"),
                              reference = "other", hull_axes = 3), "ggplot")
  expect_s3_class(plot_group_overlap(x, "ses_overlap"), "ggplot")
  y <- groups_all(s$sp, s$ord, levels = c("Anseriformes", "Strigiformes"), hull_axes = 3, nperm = 0)
  expect_null(y$focal); expect_error(summary_groups(y, "focal"), "reference")
})

test_that("kernel densities: overlap of identical sets is 1, disjoint sets ~0", {
  s <- grp_setup(); X <- s$sp$Combined$coords
  ks <- DisentFunSpaces:::.kde_setup(X, 3)
  p1 <- DisentFunSpaces:::.kde_density(X, ks)
  expect_equal(sum(p1), 1); expect_equal(DisentFunSpaces:::.kde_overlap(p1, p1), 1)
  lo <- X[, 1] < stats::quantile(X[, 1], 0.1); hi <- X[, 1] > stats::quantile(X[, 1], 0.9)
  ks2 <- DisentFunSpaces:::.kde_setup(X, 2, bw = 0.2)
  expect_lt(DisentFunSpaces:::.kde_overlap(DisentFunSpaces:::.kde_density(X[lo, ], ks2),
                                           DisentFunSpaces:::.kde_density(X[hi, ], ks2)), 0.05)
  g <- group_overlap(s$sp$Combined, s$ord, levels = c("Anseriformes", "Accipitriformes"), hull_axes = 3,
                     density = TRUE, null = "permute", nrep = 9, nperm = 0, seed = 1)
  expect_true(all(c("overlap_kde", "ses_overlap_kde") %in% names(g$pairs)))
  expect_true("FRic_kde" %in% names(g$groups))
  ## abundance weights change the density
  ab <- stats::setNames(seq_len(nrow(X)), rownames(X))
  g2 <- group_overlap(s$sp$Combined, s$ord, levels = c("Anseriformes", "Accipitriformes"), hull_axes = 3,
                      density = TRUE, abundance = ab, nperm = 0)
  expect_false(isTRUE(all.equal(g$pairs$overlap_kde, g2$pairs$overlap_kde)))
})

test_that("group_assemblage equals nonnative_assemblage with per-site status", {
  s <- grp_setup(); sp <- rownames(s$sp$Combined$coords)
  set.seed(2)
  long <- do.call(rbind, lapply(1:6, function(i) {
    x <- sample(sp, 60)
    data.frame(site = paste0("s", i), species = x,
               status = sample(c("native", "introduced"), 60, TRUE, prob = c(0.8, 0.2)))
  }))
  cm <- community_from_long(long, group_col = "status")
  expect_equal(dim(cm$comm), dim(cm$groups))
  a <- group_assemblage(s$sp$Combined, cm, reference = "native", hull_axes = 3)
  st <- cm$groups; st[is.na(st)] <- "native"
  b <- nonnative_assemblage(s$sp$Combined, st, comm = (cm$comm > 0) * 1, hull_axes = 3)
  expect_equal(a$d_FRic, b$d_FRic); expect_equal(a$FDis_all, b$FDis_cur)
  x <- groups_all(s$sp, NULL, reference = "native", comm = cm, hull_axes = 3, nperm = 0)
  expect_equal(nrow(x$sites), 6 * length(s$sp))
  expect_s3_class(plot_group_sites(x), "ggplot")
})

test_that("covariate_space: R2 matches lm, autocorrelation detects structure", {
  s <- grp_setup(); X <- s$sp$Combined$coords
  x <- stats::setNames(X[, 1] + stats::rnorm(nrow(X), 0, 0.01), rownames(X))
  cv <- covariate_space(s$sp$Combined, x, fit_axes = 3, nperm = 49, seed = 1)
  expect_equal(cv$summary$R2, summary(stats::lm(x ~ X[, 1:3]))$r.squared)
  expect_gt(cv$summary$knn_cor, 0.5); expect_lt(cv$summary$p_fit, 0.05)
  noise <- stats::setNames(stats::rnorm(nrow(X)), rownames(X))
  cv2 <- covariate_space(s$sp, noise, nperm = 49, seed = 1)
  expect_equal(nrow(cv2$summary), 3)
  expect_s3_class(plot_covariate(s$sp$Combined, x), "ggplot")
  expect_error(covariate_space(s$sp$Combined, -abs(x), transform = "log10"), "positive")
})

test_that("extinction scenarios follow their definitions", {
  s <- grp_setup()
  e <- extinction_scenarios(s$sp, s$iucn, hull_axes = 3, nrep = 9, steps = 4, seed = 1)
  expect_s3_class(e, "fs_extinction")
  nthr <- sum(s$iucn %in% c("VU", "EN", "CR"))
  expect_equal(unique(e$summary$n_removed), nthr)
  c0 <- e$curve[e$curve$n_removed == 0, ]
  expect_true(all(abs(c0$mean - 1) < 1e-12))
  expect_true(all(e$curve$mean[e$curve$metric == "FRic"] <= 1 + 1e-12))
  expect_equal(nrow(e$stochastic), 4 * length(s$sp))
  expect_s3_class(plot_extinction(e, "FDis"), "ggplot")
  ## numeric probabilities and a user priority
  P <- stats::setNames(iucn_probabilities()[s$iucn], names(s$iucn))
  e2 <- extinction_scenarios(s$sp$Diet, P, priority = -P, max_n = 10, steps = 2, nrep = 5,
                             metrics = "FRic", stochastic = FALSE)
  expect_null(e2$stochastic); expect_equal(max(e2$curve$n_removed), 10)
})

test_that("groups_from_columns turns indicator columns into non-exclusive groups", {
  d <- data.frame(species = c("a", "b", "c", "d", "e"), Fisheries = c(1, 0, 1, 0, NA),
                  Aquarium = c("no", "no", "yes", "no", NA), check.names = FALSE)
  d$`All uses` <- c(1, 1, 1, 0, NA)
  g <- groups_from_columns(d, c("Fisheries", "Aquarium"), "species")
  expect_equal(as.vector(unname(g)), c("Fisheries", "none", "Fisheries;Aquarium", "none", NA))
  expect_equal(attr(g, "sep"), ";")
  gr <- DisentFunSpaces:::.as_groups(g, d$species)   # sep taken from the attribute
  expect_false(gr$exclusive)
  expect_equal(colSums(gr$M)[c("none", "Fisheries", "Aquarium")], c(none = 2, Fisheries = 2, Aquarium = 1))
  expect_true(is.na(groups_from_columns(d, "Fisheries", "species", none = "")[["b"]]))
  expect_error(groups_from_columns(d, "species", "species"), "not a 0 / 1")
  expect_error(groups_from_columns(d, "Fisheries", "species", none = "Fisheries"), "also the name")
  expect_true(DisentFunSpaces:::.is_indicator(d$Aquarium))
  expect_false(DisentFunSpaces:::.is_indicator(c(0.2, 1.5)))
  d$Aquarium <- c(0, 1, 1, 0, NA)
  expect_equal(names(DisentFunSpaces:::.union_columns(d, c("Fisheries", "Aquarium", "All uses"))), "All uses")
  # analysis with a reference group of unused species
  s <- grp_setup()
  sp <- names(s$ord)
  dd <- data.frame(species = sp, u1 = as.numeric(s$ord == "Anseriformes"),
                   u2 = as.numeric(s$ord %in% c("Anseriformes", "Accipitriformes")))
  gg <- groups_from_columns(dd, c("u1", "u2"), "species")
  r <- groups_all(s$sp["Diet"], gg, reference = "none", hull_axes = 3, nperm = 0)
  expect_setequal(unique(r$focal$group), c("u1", "u2"))
  expect_true(is.na(r$test$permanova_R2))
})

test_that("group centroids are the mean coordinates and are drawn by plot_groups", {
  s <- grp_setup()
  ce <- group_centroids(s$sp$Combined, s$iucn, recode = s$rec, n_axes = 3)
  expect_named(ce, c("space", "group", "n", "PCoA1", "PCoA2", "PCoA3"))
  X <- s$sp$Combined$coords
  thr <- names(s$iucn)[s$iucn %in% c("VU", "EN", "CR")]
  expect_equal(unlist(ce[ce$group == "threatened", c("PCoA1", "PCoA2", "PCoA3")]),
               colMeans(X[thr, 1:3]), ignore_attr = TRUE)
  expect_equal(ce$n[ce$group == "threatened"], length(thr))
  all_sp <- group_centroids(s$sp, s$iucn, recode = s$rec)
  expect_setequal(unique(all_sp$space), names(s$sp))
  n_layers <- function(p) length(p$layers)
  expect_equal(n_layers(plot_groups(s$sp$Combined, s$iucn, recode = s$rec)),
               n_layers(plot_groups(s$sp$Combined, s$iucn, recode = s$rec, centroids = FALSE)) + 1)
  expect_equal(n_layers(plot_groups(s$sp$Combined, s$iucn, recode = s$rec, reference = "other")),
               n_layers(plot_groups(s$sp$Combined, s$iucn, recode = s$rec, reference = "other",
                                    centroids = FALSE)) + 2)
})
