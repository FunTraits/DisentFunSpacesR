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
    session$setInputs(grp_src = "simulate", nn_prop = 0.1, nn_bias = "peripheral",
                      nn_sim_space = "Combined", grp_ref = "native", grp_hull_axes = 3,
                      grp_null = "permute", grp_nrep = 9, run_groups = 1)
    expect_s3_class(rv$groups, "fs_groups_all")
    expect_setequal(unique(rv$groups$focal$space), names(rv$spaces))
    expect_equal(unique(rv$groups$focal$group), "non-native")
    # any grouping column: taxonomic order, groups only then with a reference
    session$setInputs(group_cols = c("order", "family"), grp_src = "column", grp_col = "order",
                      grp_sep = "", grp_recode = "")
    expect_true("Passeriformes" %in% names(grp_counts()$n))
    session$setInputs(grp_levels = c("Anseriformes", "Accipitriformes", "Strigiformes"),
                      grp_ref = character(0), grp_null = "none", run_groups = 2)
    expect_equal(nrow(rv$groups$pairs), 3 * length(rv$spaces))
    expect_null(rv$groups$focal)
    expect_true(all(is.finite(rv$groups$test$permanova_R2)))
    # threat status merged with the recoding box, reference = not threatened
    session$setInputs(grp_col = "iucn", grp_recode = "threatened = VU, EN, CR\nnot threatened = LC, NT")
    expect_setequal(names(grp_counts()$n), c("threatened", "not threatened"))
    session$setInputs(grp_levels = c("threatened", "not threatened"), grp_ref = "not threatened",
                      run_groups = 3)
    expect_equal(unique(rv$groups$focal$group), "threatened")
    expect_false(is.null(rv$groups$contribution))
    expect_equal(prj_params()$groups$recode$threatened, c("VU", "EN", "CR"))
    expect_equal(prj_params()$groups$reference, "not threatened")
    # densities
    session$setInputs(grp_density = TRUE, grp_density_axes = 2, grp_bw = 1, grp_abund = "", run_groups = 4)
    expect_true("overlap_kde" %in% names(rv$groups$pairs))
    expect_true(prj_params()$groups$density)
    # continuous variable
    session$setInputs(group_cols = c("order", "Eigen.1"), grp_src = "continuous", grp_col = "Eigen.1",
                      cov_transform = "none", cov_fit_axes = 3, cov_nperm = 19, run_groups = 5)
    expect_s3_class(rv$covariate, "fs_covariate")
    expect_equal(prj_params()$covariate$col, "Eigen.1")
    # extinction scenarios (tab 6)
    session$setInputs(ext_steps = 3, ext_nrep = 5, ext_metric = "FRic", run_ext = 1)
    expect_s3_class(rv$extinction, "fs_extinction")
    expect_equal(prj_params()$extinction$nrep, 5)
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

test_that("app home page and export tab follow the session", {
  shiny::testServer(DisentFunSpaces:::app_server, {
    session$setInputs(src = "example", species_col = "species", iucn_col = "iucn", sep = ",")
    expect_equal(which(steps_done()), 1L)
    expect_equal(next_step(), 2L)
    session$setInputs(validate = 1)
    session$setInputs(missing = "drop", run_missing = 1)
    session$setInputs(spaces_sel = c("Diet", "Locomotion", "Combined"), weighting = "variable",
                      correction = "none", run_spaces = 1)
    session$setInputs(pr_method = "ade4", nrepet = 19, axes_mode = "all", mantel = FALSE, run_compare = 1)
    expect_equal(which(steps_done()), 1:5)
    expect_equal(prj_steps(), "compare")
    session$setInputs(k_nn = 5, n_phylo = 10, ntree = 100, seed = 1, top = 0.1,
                      prj_title = "My analysis", prj_author = "Me", prj_name = "my project!")
    expect_equal(prj_name(), "my_project_")
    p <- prj_params()
    expect_equal(p$nrepet, 19)
    pr <- write_analysis_project(raw(), current_ti(), file.path(tempfile(), prj_name()), params = p,
                                 steps = prj_steps(), title = "My analysis", run = FALSE,
                                 zip = FALSE, quiet = TRUE)
    s <- list.files(file.path(pr$dir, "scripts"))
    expect_true("04_ANALYSIS_compare_spaces.R" %in% s)
    expect_false("05_ANALYSIS_species_metrics.R" %in% s)
    e <- new.env(); sys.source(file.path(pr$dir, "utils/parameters.R"), envir = e)
    expect_named(e$PARAMS$spaces, c("Diet", "Locomotion", "Combined"))
    expect_equal(e$PARAMS$missing, "drop")
  })
})

test_that("separator and decimal mark of text files are detected", {
  f <- tempfile(fileext = ".csv")
  writeLines(c("Species;Family;REs;VEp", "A a;Cyprinidae;0.138;0.21", "B b;Cichlidae;0.131;0.22"), f)
  d <- DisentFunSpaces:::.read_table_fast(f, "auto")
  expect_equal(ncol(d), 4); expect_equal(attr(d, "sep"), ";"); expect_true(is.numeric(d$REs))
  # wrong separator chosen by the user: read again with the detected one
  d <- DisentFunSpaces:::.read_table_fast(f, ",")
  expect_equal(ncol(d), 4)
  # decimal comma
  writeLines(c("Species;REs;VEp", "A a;0,138;0,21", "B b;1,5;2,25"), f)
  d <- DisentFunSpaces:::.read_table_fast(f, "auto")
  expect_equal(attr(d, "dec"), ","); expect_equal(d$REs, c(0.138, 1.5))
  writeLines(c("Species\tREs", "A a\t0.1"), f)
  expect_equal(attr(DisentFunSpaces:::.read_table_fast(f, "auto"), "sep"), "\t")
})
