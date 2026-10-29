test_that("write_analysis_project writes a complete, parseable project", {
  dir <- file.path(withr_tmp <- tempfile("prj"), "birds_spaces")
  on.exit(unlink(withr_tmp, recursive = TRUE))
  res <- write_analysis_project(
    birds_traits[1:300, ], birds_trait_info, dir,
    params = list(species_col = "species", iucn_col = "iucn", missing = "drop", nrepet = 19,
                  nonnative = list(source = "simulate", prop = 0.1, bias = "peripheral")),
    steps = c("compare", "metrics", "nonnative"), title = "Test", author = "Me",
    run = FALSE, zip = FALSE, quiet = TRUE)
  expect_true(is.na(res$ok))
  f <- list.files(res$dir, recursive = TRUE, all.files = TRUE)
  expect_true(all(c("run_all.R", "README.md", "birds_spaces.Rproj", ".gitignore", "session_info.txt",
                    "data/raw/traits.csv", "data/raw/trait_info.csv", "data/metadata.csv",
                    "utils/parameters.R", "utils/library.R", "utils/functions.R",
                    "manuscript/Methods_Results.Rmd",
                    "scripts/01_DATA_prepare_traits.R", "scripts/04_ANALYSIS_compare_spaces.R",
                    "scripts/06_ANALYSIS_nonnative.R") %in% f))
  for (r in f[grepl("\\.R$", f)]) expect_silent(parse(file.path(res$dir, r)))
  # every placeholder filled
  txt <- unlist(lapply(file.path(res$dir, f[grepl("\\.(R|md|Rmd)$", f)]), readLines, warn = FALSE))
  expect_false(any(grepl("{{", txt, fixed = TRUE)))
  # the parameters are read back as given
  e <- new.env(); sys.source(file.path(res$dir, "utils/parameters.R"), envir = e)
  expect_equal(e$PARAMS$missing, "drop")
  expect_equal(e$PARAMS$nrepet, 19L)
  expect_equal(e$PARAMS$nonnative$bias, "peripheral")
  expect_named(e$PARAMS$spaces, c("Locomotion", "Reproduction", "Diet", "Combined"))
  # the archived table gives the same trait object as the original one
  d <- utils::read.csv(file.path(res$dir, "data/raw/traits.csv"), check.names = FALSE,
                       stringsAsFactors = FALSE, na.strings = c("NA", ""))
  ti <- utils::read.csv(file.path(res$dir, "data/raw/trait_info.csv"), stringsAsFactors = FALSE)
  a <- prepare_traits(d, ti, species_col = "species")
  b <- prepare_traits(birds_traits[1:300, ], birds_trait_info, species_col = "species")
  expect_equal(a$data, b$data, ignore_attr = TRUE)
})

test_that("optional steps are left out and the project refuses a non-empty folder", {
  dir <- tempfile("prj2")
  on.exit(unlink(dir, recursive = TRUE))
  res <- write_analysis_project(birds_traits[1:200, ], birds_trait_info, dir,
                                params = list(missing = "drop"), steps = character(0),
                                run = FALSE, zip = FALSE, quiet = TRUE)
  s <- list.files(file.path(res$dir, "scripts"))
  expect_false(any(grepl("^0[456]_", s)))
  expect_error(write_analysis_project(birds_traits, birds_trait_info, dir, run = FALSE, quiet = TRUE),
               "not empty")
  expect_error(write_analysis_project(birds_traits, birds_trait_info, tempfile(), steps = "nonnative",
                                      run = FALSE, quiet = TRUE), "non-native")
})

test_that("the project runs on its own and zips (full pipeline)", {
  skip_on_cran()
  skip_if_not(nzchar(Sys.which("Rscript")) || file.exists(file.path(R.home("bin"), "Rscript")))
  base <- tempfile("prj3"); dir.create(base)
  on.exit(unlink(base, recursive = TRUE))
  res <- write_analysis_project(birds_traits[1:250, ], birds_trait_info, file.path(base, "run"),
                                params = list(iucn_col = "iucn", missing = "drop", nrepet = 9),
                                steps = c("compare", "metrics"), quiet = TRUE)
  expect_true(res$ok, info = res$error)
  expect_true(file.exists(res$zip))
  tabs <- list.files(file.path(res$dir, "results", "tables"))
  expect_true(all(c("Table1_traits.csv", "Table2_spaces.csv", "Table3_procrustes.csv",
                    "Table4_metric_congruence.csv", "TableS5_species_metrics.csv") %in% tabs))
  expect_true(file.exists(file.path(res$dir, "results", "figures", "Fig2_procrustes.png")))
  cmp <- readRDS(file.path(res$dir, "data", "processed", "comparison.rds"))
  expect_equal(nrow(cmp$table), 6)
})

test_that("the project tree is drawn", {
  d <- tempfile("tree"); dir.create(file.path(d, "a", "b"), recursive = TRUE)
  file.create(file.path(d, "a", "x.R"))
  on.exit(unlink(d, recursive = TRUE))
  tr <- DisentFunSpaces:::.project_tree(d)
  expect_match(tr, "x.R")
  expect_equal(length(strsplit(tr, "\n")[[1]]), 4)
})
