test_that("app server runs the full workflow on the example data", {
  skip_if_not_installed("shiny")
  shiny::testServer(DisentFunSpaces:::app_server, {
    session$setInputs(src = "example", species_col = "species", iucn_col = "iucn", sep = ",")
    expect_equal(nrow(raw()), 1000)
    expect_equal(nrow(rv$ti), 26)
    session$setInputs(validate = 1)
    expect_s3_class(rv$traits, "fs_traits")
    session$setInputs(missing = "drop", run_missing = 1)
    expect_lt(length(rv$traits_final$species), 1000)
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
  d <- birds_traits[seq(1, 1000, by = 4), setdiff(names(birds_traits), c("order", "family"))]
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
