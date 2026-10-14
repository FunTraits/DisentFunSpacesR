test_that("spaces are built and compared", {
  tr <- ex()
  sp <- build_spaces(tr)
  expect_s3_class(sp, "fs_spaces")
  expect_named(sp, c("Locomotion", "Reproduction", "Diet", "Combined"))
  expect_equal(nrow(sp$Diet$coords), length(tr$species))
  expect_true(all(sp$Combined$eigen$relative > 0))
  cmp <- compare_spaces(sp, nrepet = 19)
  expect_equal(nrow(cmp$table), choose(4, 2))
  expect_true(all(cmp$table$r >= 0 & cmp$table$r <= 1))
  expect_equal(unname(diag(cmp$matrix)), rep(1, 4))
  cv <- compare_spaces(sp[c("Diet", "Combined")], method = "vegan", nrepet = 19, mantel = TRUE, n_axes = 3)
  expect_true(all(c("m2", "mantel_r") %in% names(cv$table)))
})

test_that("a space compared with itself has Procrustes r = 1", {
  sp <- build_spaces(ex(), spaces = list(A = "Diet", B = "Diet"))
  expect_equal(compare_spaces(sp, nrepet = 9)$table$r, 1, tolerance = 1e-8)
  expect_equal(max(procrustes_residuals(sp$A, sp$B)$residual), 0, tolerance = 1e-6)
})

test_that("species metrics follow the published definitions", {
  sp <- build_spaces(ex(), spaces = list(Diet = "Diet", Loco = "Locomotion"))
  m <- species_metrics(sp$Diet, k = 5)
  expect_equal(range(m$FUn), c(0, 1)); expect_equal(range(m$FSp), c(0, 1))
  X <- sp$Diet$coords
  expect_equal(m$FSp_raw[1], sqrt(sum((X[1, ] - colMeans(X))^2)))
  iu <- stats::setNames(birds_traits$iucn, birds_traits$species)
  all <- species_metrics_all(sp, iucn = iu)
  expect_true("FUSE" %in% names(all))
  i <- 1
  expect_equal(all$FUSE[i], log(1 + all$P[i] * all$FSp[i] + all$P[i] * all$FUn[i]))
  cm <- compare_metrics(all, "FUn", top = 0.1)
  expect_equal(nrow(cm$pairs), 1)
})

test_that("plots are ggplot objects", {
  sp <- build_spaces(ex(), spaces = list(Diet = "Diet", Loco = "Locomotion"))
  expect_s3_class(plot_space(sp$Loco), "ggplot")
  expect_s3_class(plot_comparison(compare_spaces(sp, nrepet = 9)), "ggplot")
  expect_s3_class(plot_procrustes(procrustes_residuals(sp$Diet, sp$Loco)), "ggplot")
  expect_s3_class(plot_quality(space_quality(sp$Loco, 4)), "ggplot")
  expect_s3_class(plot_eigen(sp), "ggplot")
})

test_that("imputation fills all gaps and keeps observed values", {
  skip_if_not_installed("missForest")
  tr <- ex(FALSE)
  set.seed(1); idx <- sample(length(tr$species), 20)
  tr$data$clutch_size[idx] <- NA
  im <- impute_traits(tr, ntree = 20, maxiter = 3)
  expect_false(anyNA(im$data))
  obs <- !is.na(tr$data$incubation_days)
  expect_equal(im$data$incubation_days[obs], tr$data$incubation_days[obs])
  diet <- im$data[, grep("^diet_", names(im$data))]
  expect_equal(unname(rowSums(diet)), rep(1, length(tr$species)))
})

test_that("generated script parses", {
  p <- list(data_file = "x.csv", trait_info_file = "trait_info.csv", species_col = "species",
            iucn_col = "iucn", sep = ",", missing = "impute", phylo_file = "tree.nwk",
            n_phylo_axes = 10, ntree = 100, seed = 1,
            spaces = list(Diet = "Diet", Combined = c("Diet", "Locomotion")),
            weighting = "variable", correction = "none", method = "ade4", nrepet = 999,
            n_axes = NULL, mantel = FALSE, k = 5)
  expect_silent(parse(text = analysis_script(p)))
})

test_that("species lacking a whole dimension are excluded from that space only", {
  tr <- ex(FALSE)
  no_diet <- tr$species[1:5]
  tr$data[1:5, grep("^diet_", names(tr$data))] <- NA
  expect_message(sp <- build_spaces(tr, spaces = list(Diet = "Diet", Loco = "Locomotion")),
                 "excluded")
  expect_setequal(sp$Diet$excluded, no_diet)
  expect_equal(nrow(sp$Loco$coords), 333)
  cmp <- compare_spaces(sp, nrepet = 9)
  expect_equal(cmp$table$n_species, 333 - length(no_diet))
})

test_that("build_space reproduces the published PCoAGraph call", {
  tr <- ex()
  ti <- tr$trait_info; dat <- tr$data
  loco <- dat[, ti$trait[ti$dimension == "Locomotion"]]
  diet <- ade4::prep.fuzzy(dat[, ti$trait[ti$block == "diet"]], 10, labels = "diet")
  # published-style call: list of tables + groups, dist.ktab + ape::pcoa
  ref <- PCoAGraph(x = list(loco, diet), datax = cbind.data.frame(loco, diet),
                   namesx = tr$species, groups = c("Q", "F"), envfit_perm = 0)
  sp <- build_space(tr, dimensions = c("Locomotion", "Diet"), envfit_perm = 0)
  expect_equal(sp$blocks$groups, c("Q", "F"))
  expect_equal(as.vector(sp$dist), as.vector(ref$dist), tolerance = 1e-12)
  expect_equal(sp$PCoA$values$Eigenvalues, ref$PCoA$values$Eigenvalues, tolerance = 1e-10)
  expect_equal(abs(unname(sp$coords)), abs(unname(ref$PCoA$vectors)), tolerance = 1e-8)
  expect_equal(sort(rownames(sp$PCoACor)), sort(rownames(ref$PCoACor)))
  expect_true(all(c("Axis.1", "Axis.2") %in% names(sp$varPlot)))
  # same as a direct dist.ktab call
  dk <- ade4::dist.ktab(ade4::ktab.list.df(list(loco, diet)), c("Q", "F"),
                        scann = FALSE, option = "scaledBYrange")
  expect_equal(as.vector(sp$dist), as.vector(dk), tolerance = 1e-12)
})

test_that("CorrelType matches cor() for numeric traits", {
  sp <- build_space(ex(), "Locomotion", envfit_perm = 0)
  v <- sp$PCoA$vectors
  expect_equal(unname(sp$PCoACor["wing_length", 1:3]),
               unname(cor(v[, 1:3], ex()$data$wing_length)[, 1]), tolerance = 1e-12)
  ts <- sim_tr()
  nom <- ts$trait_info$trait[ts$trait_info$type == "N"][1]
  sf <- build_space(ts, ts$trait_info$dimension[ts$trait_info$trait == nom], envfit_perm = 0)
  expect_true(all(sf$PCoACor[nom, ] >= 0))             # correlation ratio
})

test_that("max_axes truncates coordinates but keeps the full PCoA", {
  tr <- ex()
  full <- build_space(tr, "Locomotion", envfit_perm = 0)
  part <- build_space(tr, "Locomotion", max_axes = 3, envfit_perm = 0)
  expect_equal(ncol(part$coords), 3)
  expect_equal(unname(part$coords), unname(full$coords[, 1:3]))
  expect_equal(ncol(part$PCoA$vectors), ncol(full$PCoA$vectors))
})

test_that("dimension weighting averages per-dimension dist.ktab distances", {
  tr <- ex()
  dims <- c("Diet", "Reproduction")
  sp <- build_space(tr, dims, weighting = "dimension", envfit_perm = 0)
  dd <- lapply(dims, function(dm) as.vector(trait_distance(tr, dm)))
  expect_equal(as.vector(sp$dist), sqrt((dd[[1]]^2 + dd[[2]]^2) / 2), tolerance = 1e-12)
})

test_that("unknown dimensions give a clear error", {
  expect_error(build_space(ex(), "Nope"), "Unknown dimension")
})

test_that("columns of the data table can be used as imputation predictors", {
  skip_if_not_installed("missForest")
  tr <- ex(FALSE)
  pred <- birds_traits[, paste0("Eigen.", 1:5)]
  rownames(pred) <- birds_traits$species
  im <- impute_traits(tr, covariates = pred, ntree = 20, maxiter = 3)
  expect_false(anyNA(im$data))
  expect_equal(im$imputation$covariates, paste0("Eigen.", 1:5))
  expect_false(any(grepl("Eigen", names(im$data))))          # predictors are not traits
  im2 <- impute_traits(tr, phylo = birds_tree, covariates = pred[, 1:2], n_phylo_axes = 3,
                       ntree = 20, maxiter = 3)
  expect_false(anyNA(im2$data))
  bad <- pred; rownames(bad) <- NULL
  expect_error(impute_traits(tr, covariates = bad), "row names")
  both <- birds_traits[, c("body_mass", "Eigen.1")]; rownames(both) <- birds_traits$species
  expect_error(impute_traits(tr, covariates = both), "both as trait and as predictor")
})

test_that("undefined distances stop before the PCoA with a clear message", {
  df <- data.frame(sp = paste0("s", 1:6),
                   a1 = c(1, 0, 1, NA, NA, 1), a2 = c(0, 1, 1, NA, NA, 0),
                   b1 = c(1, 1, 0, 1, 0, NA), b2 = c(0, 0, 1, 1, 1, NA))
  ti <- data.frame(trait = c("a1", "a2", "b1", "b2"), type = "B", dimension = "H",
                   block = c("a", "a", "b", "b"))
  tr <- prepare_traits(df, ti, "sp")
  expect_error(build_space(tr, envfit_perm = 0), "share no non-missing trait")
})

test_that("imputed multichoice blocks keep at least one modality", {
  skip_if_not_installed("missForest")
  tr <- sim_tr()                                  # no multichoice trait in the bird data
  b1 <- tr$trait_info$block[tr$trait_info$type == "B"][1]
  hab <- tr$trait_info$trait[tr$trait_info$block == b1]
  set.seed(3); idx <- sample(nrow(tr$data), 30)
  tr$data[idx, hab] <- NA
  im <- impute_traits(tr, ntree = 20, maxiter = 3)
  m <- as.matrix(im$data[, hab])
  expect_true(all(m %in% c(0, 1)))
  expect_true(all(rowSums(m) >= 1))
  obs <- setdiff(seq_len(nrow(m)), idx)
  expect_equal(m[obs, ], as.matrix(tr$data[obs, hab]))
})

test_that("engine = 'fast' gives the same space as dist.ktab", {
  tr <- ex(FALSE)
  a <- build_space(tr, c("Locomotion", "Diet", "Reproduction"), envfit_perm = 0)
  b <- build_space(tr, c("Locomotion", "Diet", "Reproduction"), envfit_perm = 0, engine = "fast")
  expect_equal(as.vector(a$dist), as.vector(b$dist), tolerance = 1e-10)
  expect_equal(a$eigen, b$eigen, tolerance = 1e-8)
})

test_that("identical profiles get exactly FUn = 0 and are counted", {
  df <- data.frame(sp = paste0("s", 1:12), d = c(rep("a", 8), "b", "b", "c", "a"),
                   e = c(rep("x", 8), "y", "x", "y", "y"))
  ti <- data.frame(trait = c("d", "e"), type = "N", dimension = "A")
  sp <- build_space(prepare_traits(df, ti, "sp"), envfit_perm = 0)
  m <- species_metrics(sp, k = 3)
  expect_true(all(m$FUn_raw[1:8] == 0))
  expect_equal(m$n_identical[1:8], rep(7L, 8))
  expect_equal(length(unique(m$FSp_raw[1:8])), 1)
  mm <- rbind(cbind(space = "A", m), cbind(space = "B", m))
  cm <- suppressWarnings(compare_metrics(mm, "FUn", top = 0.5))
  expect_true(all(c("n_top", "informative") %in% names(cm$top_info)))
  # 8 species tied at FUn = 0 cannot all be in a top-10 % set: flagged
  expect_warning(compare_metrics(mm, "FSp", top = 0.1), NA)
  expect_warning(compare_metrics(transform(mm, FUn = 0), "FUn", top = 0.1), "uninformative")
})

test_that("the memory-lean PCoA gives the same result as ape::pcoa", {
  for (tr in list(ex(), ex(FALSE))) {
    d <- trait_distance(tr)
    a <- ape::pcoa(d)
    b <- DisentFunSpaces:::.pcoa_lean(d)
    expect_equal(a$values, b$values, tolerance = 1e-8)
    expect_equal(abs(unname(a$vectors)), abs(unname(b$vectors)), tolerance = 1e-6)
    expect_equal(a$trace, b$trace, tolerance = 1e-10)
  }
  old <- options(DisentFunSpaces.lean_pcoa_n = 10)
  on.exit(options(old))
  sp1 <- build_space(ex(), "Locomotion", envfit_perm = 0)
  options(DisentFunSpaces.lean_pcoa_n = 1e6)
  sp2 <- build_space(ex(), "Locomotion", envfit_perm = 0)
  expect_equal(abs(unname(sp1$coords)), abs(unname(sp2$coords)), tolerance = 1e-6)
})
