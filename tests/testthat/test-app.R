test_that("app server runs the full workflow on the example data", {
  skip_if_not_installed("shiny")
  shiny::testServer(DisentFunSpaces:::app_server, {
    session$setInputs(src = "example", species_col = "species", iucn_col = "iucn", sep = ",")
    expect_equal(nrow(raw()), 2000)
    expect_equal(nrow(rv$ti), 26)
    session$setInputs(validate = 1)
    expect_s3_class(rv$traits, "fs_traits")
    session$setInputs(missing = "drop", run_missing = 1)
    expect_lt(length(rv$traits_final$species), 2000)
    session$setInputs(spaces_sel = c("Diet", "Locomotion", "Combined"), weighting = "variable",
                      correction = "none", run_spaces = 1)
    expect_named(rv$spaces, c("Diet", "Locomotion", "Combined"))
    session$setInputs(pr_method = "ade4", nrepet = 19, axes_mode = "all", mantel = FALSE, run_compare = 1)
    expect_equal(nrow(rv$cmp$table), 3)
    session$setInputs(k_nn = 5, metric = "FUSE", top = 0.1, run_metrics = 1)
    expect_true("FUSE" %in% names(rv$metrics))
    session$setInputs(n_phylo = 10, ntree = 100, seed = 1)
    expect_true(grepl("build_spaces", analysis_script(script_params())))
    session$setInputs(nn_src = "simulate", nn_prop = 0.1, nn_bias = "peripheral",
                      nn_sim_space = "Combined", nn_hull_axes = 3, nn_null = "permute",
                      nn_nrep = 9, run_nonnative = 1)
    expect_s3_class(rv$nonnative, "fs_nonnative_all")
    expect_setequal(unique(rv$nonnative$nonnative$space), names(rv$spaces))
  })
})

test_that("app handles an uploaded semicolon CSV with guessed trait types", {
  f <- tempfile(fileext = ".csv")
  d <- birds_traits[seq(1, 2000, by = 8), setdiff(names(birds_traits), c("order", "family"))]
  utils::write.table(d, f, sep = ";", row.names = FALSE)
  shiny::testServer(DisentFunSpaces:::app_server, {
    session$setInputs(src = "upload", sep = ";", file = list(datapath = f, name = "my_traits.csv"))
    session$setInputs(species_col = "species", iucn_col = "iucn")
    expect_equal(ncol(raw()), ncol(d))
    expect_false("iucn" %in% rv$ti$trait)
    expect_equal(rv$ti$type[rv$ti$trait == "diet_inv"], "F")
    # assign two dimensions with the bulk tool
    session$setInputs(bulk_dim = "Diet", bulk_rows = grep("^diet_", rv$ti$trait, value = TRUE), bulk_apply = 1)
    session$setInputs(validate = 1)
    expect_setequal(unique(rv$traits$trait_info$dimension), c("Diet", "All traits"))
    session$setInputs(missing = "keep", run_missing = 1)
    session$setInputs(spaces_sel = c("Diet", "All traits", "Combined"), weighting = "variable",
                      correction = "none", max_axes = 5, run_spaces = 1)
    expect_equal(ncol(rv$spaces$Combined$coords), 5)
  })
})

test_that("app imputes with predictor columns chosen from the data table", {
  skip_if_not_installed("missForest")
  shiny::testServer(DisentFunSpaces:::app_server, {
    session$setInputs(src = "example", species_col = "species", iucn_col = "iucn", sep = ",")
    session$setInputs(validate = 1)
    session$setInputs(missing = "impute", cov_cols = c("Eigen.1", "Eigen.2", "Eigen.3"),
                      n_phylo = 10, ntree = 20, seed = 1, run_missing = 1)
    expect_false(anyNA(rv$traits_final$data))
    expect_equal(rv$traits_final$imputation$covariates, c("Eigen.1", "Eigen.2", "Eigen.3"))
    expect_match(rv$missing_msg, "Eigen.1")
    session$setInputs(spaces_sel = c("Diet", "Combined"), weighting = "variable",
                      correction = "none", max_axes = NA, run_spaces = 1)
    session$setInputs(pr_method = "ade4", nrepet = 9, axes_mode = "all", mantel = FALSE, k_nn = 5)
    scr <- analysis_script(script_params())
    expect_match(scr, "covariates = predictors", fixed = TRUE)
    expect_silent(parse(text = scr))
  })
})

test_that("app loads the full 8,413-species example with its trait table", {
  shiny::testServer(DisentFunSpaces:::app_server, {
    session$setInputs(src = "example_all", species_col = "species", iucn_col = "iucn", sep = ",")
    expect_equal(nrow(raw()), 8413)
    expect_equal(nrow(rv$ti), 26)
    expect_equal(data_name(), "birds_traits_all.csv")
  })
})

test_that("app uses the large-data settings (fast engine, partial PCoA, cores)", {
  skip_if_not_installed("RSpectra")
  op <- options(DisentFunSpaces.app_large_n = 100); on.exit(options(op))
  shiny::testServer(DisentFunSpaces:::app_server, {
    session$setInputs(src = "example", species_col = "species", iucn_col = "iucn", sep = ",")
    session$setInputs(validate = 1)
    session$setInputs(missing = "drop", run_missing = 1)
    session$setInputs(spaces_sel = c("Diet", "Locomotion", "Combined"), weighting = "variable",
                      correction = "none", max_axes = 6, engine = "auto", eigen_mode = "partial",
                      run_spaces = 1)
    expect_true(isTRUE(rv$spaces$Combined$PCoA$partial))
    expect_equal(ncol(rv$spaces$Combined$coords), 6)
    expect_equal(rv$build_settings$engine, "ade4")       # < 1000 species
    # all positive axes with the partial decomposition = full decomposition
    session$setInputs(max_axes = NA, eigen_mode = "partial", run_spaces = 2)
    pa <- rv$spaces$Locomotion
    session$setInputs(eigen_mode = "full", run_spaces = 3)
    fu <- rv$spaces$Locomotion
    expect_equal(ncol(pa$coords), ncol(fu$coords))
    expect_equal(abs(unname(pa$coords)), abs(unname(fu$coords)), tolerance = 1e-8)
    session$setInputs(eigen_mode = "partial", max_axes = 6, run_spaces = 4)
    session$setInputs(pr_method = "ade4", nrepet = 19, axes_mode = "all", mantel = FALSE,
                      cores = 2, run_compare = 1)
    expect_equal(nrow(rv$cmp$table), 3)
    session$setInputs(k_nn = 5, n_phylo = 10, ntree = 100, seed = 1)
    scr <- analysis_script(script_params())
    expect_match(scr, 'eigen = "partial"', fixed = TRUE)
    expect_silent(parse(text = scr))
  })
})
