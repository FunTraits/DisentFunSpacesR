#' Write (and run) a reproducible analysis project
#'
#' Writes a self-contained, reproducible project for a functional-space
#' analysis, organised as a research repository, then runs its own scripts
#' and zips it. This is what the "Export" tab of the Shiny application does:
#' it does not save what is on screen, it writes the project with the
#' settings of the tabs and RUNS it, so that the archive is exactly what
#' anyone re-running it obtains.
#'
#' The project contains:
#' \preformatted{
#' <project>/
#' |-- scripts/      00_START_GeneralScript.R, 01_DATA_prepare_traits.R,
#' |                 02_DATA_missing_values.R, 03_ANALYSIS_functional_spaces.R,
#' |                 04_ANALYSIS_compare_spaces.R, 05_ANALYSIS_species_metrics.R,
#' |                 06_ANALYSIS_groups.R (the steps requested), setup_environment.R
#' |-- utils/        library.R, parameters.R (every setting), functions.R
#' |-- data/raw/     traits.csv, trait_info.csv (+ phylogeny.tre)
#' |-- data/processed/  objects passed between scripts (.rds)
#' |-- data/metadata.csv  definition of every variable
#' |-- results/figures/, results/tables/
#' |-- manuscript/   Methods_Results.Rmd (+ .docx once rendered)
#' |-- run_all.R, <project>.Rproj, README.md, session_info.txt, .gitignore
#' }
#'
#' @param data Trait table (data frame, one row per species), as given to
#'   [prepare_traits()]. It may hold further columns (IUCN category, grouping variables,
#'   predictors of missing values); they are archived with it.
#' @param trait_info Trait-description table (see [check_trait_info()]).
#' @param dir Folder of the project (created; must not exist or be empty).
#'   Its name is the project name.
#' @param params Named list of settings, with the names used by
#'   [analysis_script()]: `species_col`, `iucn_col`,
#'   `missing` (`"keep"`, `"drop"`, `"impute"`), `cov_cols`,
#'   `n_phylo_axes`, `ntree`, `seed`, `spaces` (named list of dimensions),
#'   `weighting`, `correction`, `max_axes`, `engine`, `eigen`, `keep_dist`,
#'   `method`, `nrepet`, `n_axes`, `mantel`, `cores`, `k`, and also `top`
#'   (proportion of species of the top sets, default 0.1), `coords_axes`
#'   (axes in the coordinate table, default 10) and `groups` (list:
#'   `source` = `"column"` or `"simulate"`, `col` (grouping column), `sep`,
#'   `recode`, `levels`, `reference`, `pool_focal`, `prop`, `bias`,
#'   `sim_space`, `hull_axes`, `null`, `nrep`, `nperm`, `density`,
#'   `density_axes`, `bw`, `abundance_col`, `sites`; see [groups_all()]),
#'   `covariate` (continuous variable that is not a trait: list `col`,
#'   `transform`, `fit_axes`, `nperm`; see [covariate_space()]) and
#'   `extinction` (with the metrics step and an IUCN column: list `steps`,
#'   `nrep`, `hull_axes`; see [extinction_scenarios()]).
#'   The settings of version 0.7.0 (`status_col` and `nonnative` = list(`source`,
#'   `native`, ...)) are still accepted and translated (reference = native
#'   label). Missing settings take the defaults of the published analysis.
#' @param assemblages Optional long table of assemblages (one row per
#'   species in an assemblage; see [community_from_long()]), used by the
#'   species-groups step when `params$groups$sites` gives its columns
#'   (`site_col`, `species_col`, `group_col`); archived as
#'   `data/raw/assemblages.csv`.
#' @param phylo Optional phylogeny (`phylo` object or Newick file) used to
#'   impute missing values; archived as `data/raw/phylogeny.tre`.
#' @param steps Optional analyses included after the spaces: any of
#'   `"compare"` (Procrustes tests), `"metrics"` (FUn, FSp, FUSE) and
#'   `"groups"` (species groups, [groups_all()]; `"nonnative"` is an alias)
#'   (none: `character(0)`).
#' @param title,author Title and author(s), used in the README and report.
#' @param data_source Name of the original data file (for the README).
#' @param run Run the project's scripts (`run_all.R`, in a separate R
#'   process) once written.
#' @param zip Zip the project (`<dir>.zip`, next to it).
#' @param quiet No progress messages.
#' @return A list: `dir`, `zip` (path or `NULL`), `log` (path of the run
#'   log), `ok` (`TRUE` when the run completed, `NA` when not run) and
#'   `error` (message, when it did not).
#' @seealso [run_app()], [analysis_script()]
#' @export
#' @examples
#' \donttest{
#' data(birds_traits); data(birds_trait_info)
#' res <- write_analysis_project(
#'   birds_traits[1:300, ], birds_trait_info,
#'   dir = file.path(tempdir(), "birds_spaces"),
#'   params = list(species_col = "species", iucn_col = "iucn",
#'                 missing = "drop", nrepet = 99),
#'   steps = c("compare", "metrics"), title = "Functional spaces of birds",
#'   run = FALSE, zip = FALSE)
#' list.files(res$dir, recursive = TRUE)
#' }
write_analysis_project <- function(data, trait_info, dir, params = list(),
                                   phylo = NULL, assemblages = NULL, steps = c("compare", "metrics"),
                                   title = "Functional spaces analysis", author = "",
                                   data_source = NULL, run = TRUE, zip = TRUE,
                                   quiet = FALSE) {
  say <- function(...) if (!quiet) message(...)
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  if (dir.exists(dir) && length(list.files(dir, all.files = TRUE, no.. = TRUE)))
    stop("`dir` exists and is not empty: ", dir, call. = FALSE)
  steps <- unique(sub("^nonnative$", "groups", as.character(steps)))
  bad <- setdiff(steps, c("compare", "metrics", "groups"))
  if (length(bad)) stop("Unknown step(s): ", paste(bad, collapse = ", "),
                        " (use \"compare\", \"metrics\", \"groups\").", call. = FALSE)
  tpl <- system.file("templates", "analysis_project", package = "DisentFunSpaces")
  if (!nzchar(tpl)) stop("Templates not found in the DisentFunSpaces installation.", call. = FALSE)

  # --- settings -----------------------------------------------------------------
  ti <- check_trait_info(trait_info, data)
  p <- .project_params(params, ti, data, has_phylo = !is.null(phylo))
  if (!p$species_col %in% names(data))
    stop("Species column '", p$species_col, "' not found in `data`.", call. = FALSE)
  if ("groups" %in% steps && is.null(p$groups) && is.null(p$covariate))
    stop("The species-groups step needs `params$groups$col` (a grouping column; for non-native ",
         "species, the status column) or `params$groups$source = \"simulate\"`.", call. = FALSE)
  if ("groups" %in% steps && identical(p$groups$source, "column") && !all(p$groups$col %in% names(data)))
    stop("Grouping column(s) '", paste(setdiff(p$groups$col, names(data)), collapse = "', '"),
         "' not found in `data`.", call. = FALSE)
  if (!is.null(p$covariate) && !p$covariate$col %in% names(data))
    stop("Column '", p$covariate$col, "' not found in `data`.", call. = FALSE)
  if (!is.null(p$groups$sites) && is.null(assemblages)) {
    warning("`params$groups$sites` given without `assemblages`: assemblages left out.", call. = FALSE)
    p$groups$sites <- NULL
  }
  if (!"groups" %in% steps) { p$groups <- NULL; p$covariate <- NULL }
  if (!"metrics" %in% steps || is.null(p$iucn_col)) p$extinction <- NULL
  n_spaces <- length(p$spaces)
  if ("compare" %in% steps && n_spaces < 2) steps <- setdiff(steps, "compare")

  project <- basename(normalizePath(dir, mustWork = FALSE))
  for (d in c("scripts", "utils", "data/raw", "data/processed", "results/figures",
              "results/tables", "manuscript"))
    dir.create(file.path(dir, d), recursive = TRUE, showWarnings = FALSE)

  # --- data/raw -----------------------------------------------------------------
  utils::write.csv(data, file.path(dir, "data/raw/traits.csv"), row.names = FALSE, na = "",
                   fileEncoding = "UTF-8")
  ti_out <- as.data.frame(trait_info, stringsAsFactors = FALSE)
  utils::write.csv(ti_out, file.path(dir, "data/raw/trait_info.csv"), row.names = FALSE, na = "")
  if (!is.null(phylo)) {
    if (is.character(phylo)) phylo <- ape::read.tree(phylo)
    ape::write.tree(phylo, file.path(dir, "data/raw/phylogeny.tre"))
  }
  if (!is.null(p$groups$sites))
    utils::write.csv(as.data.frame(assemblages), file.path(dir, "data/raw/assemblages.csv"),
                     row.names = FALSE, na = "", fileEncoding = "UTF-8")
  utils::write.csv(.project_metadata(), file.path(dir, "data/metadata.csv"), row.names = FALSE)

  # --- scripts ------------------------------------------------------------------
  scripts <- c("01_DATA_prepare_traits.R", "02_DATA_missing_values.R",
               "03_ANALYSIS_functional_spaces.R",
               if ("compare" %in% steps) "04_ANALYSIS_compare_spaces.R",
               if ("metrics" %in% steps) "05_ANALYSIS_species_metrics.R",
               if ("groups" %in% steps) "06_ANALYSIS_groups.R")
  cran <- c("ggplot2", "ade4", "ape", "vegan",
            if (identical(p$missing, "impute")) "missForest",
            if (identical(p$eigen, "partial")) "RSpectra")
  dims <- unique(ti$dimension)
  extra <- c(if (!is.null(p$iucn_col)) sprintf("`%s` (IUCN category, for FUSE)", p$iucn_col),
             if (!is.null(p$groups$col)) sprintf("%s (species groups, not a trait)",
               paste0("`", p$groups$col, "`", collapse = ", ")),
             if (!is.null(p$covariate$col)) sprintf("`%s` (continuous variable, not a trait)", p$covariate$col),
             if (length(p$cov_cols)) sprintf("%s (predictors of missing values, not traits)",
                                             paste0("`", p$cov_cols, "`", collapse = ", ")))
  vals <- list(
    TITLE = title, AUTHOR = if (nzchar(author)) author else "[author]",
    DATE = format(Sys.Date(), "%Y-%m-%d"), PROJECT = project,
    DFS_VERSION = as.character(utils::packageVersion("DisentFunSpaces")),
    CRAN_PACKAGES = paste0('"', cran, '"', collapse = ", "),
    N_SPECIES = format(nrow(data), big.mark = ","), N_TRAITS = length(unique(ti$block)),
    N_DIMS = length(dims), DIMENSIONS = paste(dims, collapse = ", "),
    SPECIES_COL = p$species_col,
    DATA_SOURCE = if (is.null(data_source)) "data frame given in R" else data_source,
    EXTRA_COLS = if (length(extra)) paste0(" Other columns: ", paste(extra, collapse = "; "), ".") else "",
    PHYLO_TREE = if (!is.null(phylo)) " (+ phylogeny.tre)" else "",
    PHYLO_SENTENCE = if (!is.null(phylo))
      "\n\n`data/raw/phylogeny.tre` is the phylogeny whose eigenvectors predict missing values." else "",
    SCRIPT_TREE = paste(sprintf("\u2502   \u251c\u2500\u2500 %-34s # %s", scripts, .PROJECT_SCRIPTS[scripts]),
                        collapse = "\n"),
    STEP_TABLE = paste(sprintf("| %d | `%s` | %s |", seq_along(scripts), scripts,
                               .PROJECT_OUTPUTS[scripts]), collapse = "\n"),
    SETTINGS = .project_settings_md(p, steps))
  cp <- function(from, to = from) {
    txt <- readLines(file.path(tpl, from), warn = FALSE, encoding = "UTF-8")
    writeLines(enc2utf8(.fill_template(txt, vals)), file.path(dir, to), useBytes = TRUE)
  }
  for (s in c("00_START_GeneralScript.R", scripts, "setup_environment.R"))
    cp(file.path("scripts", s))
  cp("utils/library.R"); cp("utils/functions.R")
  writeLines(enc2utf8(.project_parameters_r(p, title)), file.path(dir, "utils/parameters.R"),
             useBytes = TRUE)
  cp("run_all.R"); cp("README.md")
  cp("manuscript/Methods_Results.Rmd")
  # named rproj.txt / gitignore.txt in the package: .Rbuildignore drops *.Rproj
  file.copy(file.path(tpl, "rproj.txt"), file.path(dir, paste0(project, ".Rproj")))
  file.copy(file.path(tpl, "gitignore.txt"), file.path(dir, ".gitignore"))

  pk <- c("DisentFunSpaces", "ade4", "ape", "vegan", "ggplot2", "missForest", "RSpectra",
          "rmarkdown", "knitr")
  writeLines(c(paste("Project:", title),
               paste("Written:", format(Sys.time(), "%Y-%m-%d %H:%M %Z")),
               "(replaced by sessionInfo() at the end of run_all.R)", "",
               paste("R", paste(R.version$major, R.version$minor, sep = "."), "--", R.version$platform),
               vapply(pk, function(x) sprintf("%-16s %s", x,
                 tryCatch(as.character(utils::packageVersion(x)), error = function(e) "not installed")),
                 character(1))),
             file.path(dir, "session_info.txt"))

  out <- list(dir = normalizePath(dir), zip = NULL, log = NULL, ok = NA, error = NULL)

  # --- run the project's own scripts ----------------------------------------------
  if (run) {
    log <- file.path(out$dir, "results", "run_all.log")
    say("Running the project's scripts (run_all.R)...")
    rs <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
    old <- Sys.getenv("R_LIBS", unset = NA)
    Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
    on.exit(if (is.na(old)) Sys.unsetenv("R_LIBS") else Sys.setenv(R_LIBS = old), add = TRUE)
    owd <- setwd(out$dir)
    on.exit(setwd(owd), add = TRUE)
    status <- tryCatch(suppressWarnings(system2(rs, c("--no-save", "--no-restore", "run_all.R"),
                                                stdout = log, stderr = log)),
                       error = function(e) { cat(conditionMessage(e), file = log, append = TRUE); 1L })
    setwd(owd)
    out$log <- log
    lines <- tryCatch(readLines(log, warn = FALSE), error = function(e) character(0))
    out$ok <- identical(as.integer(status), 0L) && any(grepl("run successfully", lines))
    if (!out$ok) {
      out$error <- paste(utils::tail(lines, 15), collapse = "\n")
      if (!quiet) warning("The project did not run to the end:\n", out$error, call. = FALSE)
    }
  }

  # --- zip --------------------------------------------------------------------------
  if (zip) {
    zf <- paste0(out$dir, ".zip")
    if (file.exists(zf)) unlink(zf)
    owd2 <- setwd(dirname(out$dir))
    on.exit(setwd(owd2), add = TRUE)
    files <- list.files(project, recursive = TRUE, all.files = TRUE, full.names = TRUE)
    if (requireNamespace("zip", quietly = TRUE)) zip::zip(zf, files) else
      utils::zip(zf, files, flags = "-r9Xq")
    setwd(owd2)
    out$zip <- zf
    say("Project archive: ", zf)
  }
  out
}

# --- what each script does and writes (README) --------------------------------
.PROJECT_SCRIPTS <- c(
  "01_DATA_prepare_traits.R" = "Reads and checks the traits",
  "02_DATA_missing_values.R" = "Missing values (keep / drop / impute)",
  "03_ANALYSIS_functional_spaces.R" = "Gower distance + PCoA per space",
  "04_ANALYSIS_compare_spaces.R" = "Procrustes tests between spaces",
  "05_ANALYSIS_species_metrics.R" = "FUn, FSp, FUSE per space",
  "06_ANALYSIS_groups.R" = "Species groups (overlap, focal vs reference)")
.PROJECT_OUTPUTS <- c(
  "01_DATA_prepare_traits.R" = "Table1_traits",
  "02_DATA_missing_values.R" = "TableS1_missing_values (+ TableS2_imputation_error)",
  "03_ANALYSIS_functional_spaces.R" = "Table2_spaces, TableS3/S4, Fig1_space_*, FigS1_inertia, FigS2_quality",
  "04_ANALYSIS_compare_spaces.R" = "Table3_procrustes, Fig2_procrustes",
  "05_ANALYSIS_species_metrics.R" = "Table4_metric_congruence, TableS5_species_metrics, Fig3_metric_congruence (+ Table8, TableS10-S11, Fig6_extinction_*)",
  "06_ANALYSIS_groups.R" = "Table5_groups_summary, Table6_group_pairs (+ Table7, Table9, TableS6-S13), Fig4_groups_*, Fig5_group_overlap (+ Fig7, Fig8)")

# Settings completed with the defaults of the published analysis.
.project_params <- function(p, ti, data, has_phylo = FALSE) {
  `%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a
  nz <- function(x) if (is.null(x) || !length(x) || is.na(x[1]) || !nzchar(x[1])) NULL else x
  dims <- unique(ti$dimension)
  default_sp <- stats::setNames(as.list(dims), dims)
  if (length(dims) > 1) default_sp$Combined <- dims
  n <- nrow(data)
  out <- list(
    data_file = "data/raw/traits.csv",
    trait_info_file = "data/raw/trait_info.csv",
    phylo_file = if (has_phylo) "data/raw/phylogeny.tre" else NULL,
    species_col = p$species_col %||% "species",
    iucn_col = nz(p$iucn_col),
    seed = p$seed %||% 20251029,
    missing = p$missing %||% "keep",
    cov_cols = if (identical(p$missing, "impute")) nz(p$cov_cols) else NULL,
    n_phylo_axes = as.integer(p$n_phylo_axes %||% 10), ntree = as.integer(p$ntree %||% 100),
    spaces = p$spaces %||% default_sp,
    weighting = p$weighting %||% "variable", correction = p$correction %||% "none",
    max_axes = if (is.null(p$max_axes) || is.na(p$max_axes)) NULL else as.integer(p$max_axes),
    engine = p$engine %||% "ade4", eigen = p$eigen %||% "full",
    keep_dist = p$keep_dist %||% (n < 5000),
    method = p$method %||% "ade4", nrepet = as.integer(p$nrepet %||% 999),
    n_axes = if (is.null(p$n_axes) || is.na(p$n_axes)) NULL else as.integer(p$n_axes),
    mantel = isTRUE(p$mantel), cores = as.integer(p$cores %||% 1),
    k = as.integer(p$k %||% 5), top = p$top %||% 0.1, coords_axes = as.integer(p$coords_axes %||% 10))
  gp <- p$groups
  # settings of version 0.7.0: status column + nonnative list
  if (is.null(gp) && (!is.null(p$nonnative) || !is.null(nz(p$status_col)))) {
    nn <- p$nonnative %||% list()
    gp <- list(source = nn$source %||% "column", col = nz(p$status_col),
               reference = nn$native %||% "native", prop = nn$prop, bias = nn$bias,
               sim_space = nn$sim_space, hull_axes = nn$hull_axes, null = nn$null, nrep = nn$nrep)
    if (identical(gp$null, "pool")) gp$null <- "permute"
  }
  if (is.null(gp) && !is.null(nz(p$group_col))) gp <- list(col = p$group_col)
  if (!is.null(gp)) {
    src <- gp$source %||% "column"
    col <- nz(gp$col %||% p$group_col)
    if (src == "column" && is.null(col)) gp <- NULL else
      gp <- list(source = src, col = if (src == "column") col else NULL,
                 none = if (src == "column" && (length(col) > 1 || !is.null(gp$none))) (gp$none %||% "") else NULL,
                 sep = nz(gp$sep), recode = if (length(gp$recode)) gp$recode else NULL,
                 levels = if (src == "column" && length(gp$levels)) as.character(gp$levels) else NULL,
                 reference = if (src == "simulate") (nz(gp$reference) %||% "native") else nz(gp$reference),
                 pool_focal = isTRUE(gp$pool_focal),
                 prop = gp$prop %||% 0.1, bias = gp$bias %||% "random",
                 sim_space = gp$sim_space %||% utils::tail(names(out$spaces), 1),
                 hull_axes = as.integer(gp$hull_axes %||% 4), null = gp$null %||% "none",
                 nrep = as.integer(gp$nrep %||% 99), nperm = as.integer(gp$nperm %||% 999),
                 density = isTRUE(gp$density), density_axes = as.integer(gp$density_axes %||% 2),
                 bw = gp$bw %||% 1, abundance_col = nz(gp$abundance_col),
                 sites = if (length(gp$sites) && !is.null(nz(gp$sites$site_col)))
                   list(site_col = gp$sites$site_col, species_col = gp$sites$species_col %||% "species",
                        group_col = nz(gp$sites$group_col)) else NULL)
  }
  out$groups <- gp
  cv <- p$covariate
  out$covariate <- if (!is.null(nz(cv$col))) list(col = cv$col, transform = cv$transform %||% "none",
    fit_axes = as.integer(cv$fit_axes %||% 4), nperm = as.integer(cv$nperm %||% 999)) else NULL
  ex <- p$extinction
  out$extinction <- if (!is.null(ex)) list(steps = as.integer(ex$steps %||% 10),
    nrep = as.integer(ex$nrep %||% 99), hull_axes = as.integer(ex$hull_axes %||% 4)) else NULL
  out
}

.fill_template <- function(txt, vals) {
  txt <- enc2utf8(txt)
  for (k in names(vals))
    txt <- gsub(paste0("{{", k, "}}"), enc2utf8(as.character(vals[[k]])), txt, fixed = TRUE,
                useBytes = TRUE)
  Encoding(txt) <- "UTF-8"
  txt
}

.dep <- function(x) paste(deparse(x, width.cutoff = 500L,
                                  control = c("keepNA", "keepInteger", "niceNames")), collapse = "")

# utils/parameters.R: every setting in one commented list.
.project_parameters_r <- function(p, title) {
  line <- function(name, comment, last = FALSE)
    sprintf("  %-15s = %s%s  # %s", name, .dep(p[[name]]), if (last) "" else ",", comment)
  sp <- paste0("list(\n", paste0("    `", names(p$spaces), "` = ", vapply(p$spaces, .dep, ""),
                                 collapse = ",\n"), "\n  )")
  lst <- function(x) if (is.null(x)) "NULL" else paste0("list(\n",
    paste0("    ", names(x), " = ", vapply(x, .dep, ""), collapse = ",\n"), "\n  )")
  nn <- lst(p$groups)
  c("# ============================================================================",
    sprintf("# Parameters of \"%s\"", title),
    sprintf("# Written by DisentFunSpaces %s on %s. Every setting of the analysis is here:",
            utils::packageVersion("DisentFunSpaces"), format(Sys.Date())),
    "# change it here (not in the scripts), then source(\"run_all.R\").",
    "# ============================================================================",
    "PARAMS <- list(",
    "  # --- data ---",
    line("data_file", "trait table (one row per species)"),
    line("trait_info_file", "type, dimension, transform, weight of each trait"),
    line("phylo_file", "phylogeny used to impute missing values (or NULL)"),
    line("species_col", "species column"),
    line("iucn_col", "IUCN category column, for FUSE (or NULL)"),
    line("seed", "random seed (imputation, permutations, simulations)"),
    "  # --- missing values ---",
    line("missing", "\"keep\" (pairwise Gower), \"drop\" or \"impute\" (missForest)"),
    line("cov_cols", "predictor columns of the trait table (imputation)"),
    line("n_phylo_axes", "phylogenetic eigenvectors used as predictors"),
    line("ntree", "trees per random forest"),
    "  # --- functional spaces ---",
    paste0("  spaces          = ", sp, ",  # dimensions of each space"),
    line("weighting", "\"variable\" (each trait equally), \"dimension\", \"dimension_scaled\""),
    line("correction", "negative eigenvalues: \"none\", \"cailliez\", \"lingoes\""),
    line("max_axes", "PCoA axes kept (NULL = all positive axes)"),
    line("engine", "distance: \"ade4\" (dist.ktab) or \"fast\" (same values)"),
    line("eigen", "PCoA: \"full\" or \"partial\" (RSpectra, same axes)"),
    line("keep_dist", "keep distances (quality, Mantel tests)"),
    line("coords_axes", "axes written in TableS3_pcoa_coordinates"),
    "  # --- comparison of spaces ---",
    line("method", "Procrustes test: \"ade4\" (procuste.rtest) or \"vegan\" (protest)"),
    line("nrepet", "permutations"),
    line("n_axes", "axes compared and used by the metrics (NULL = all)"),
    line("mantel", "add Mantel tests"),
    line("cores", "CPU cores for the permutations (same results)"),
    "  # --- species metrics ---",
    line("k", "nearest neighbours (FUn)"),
    line("top", "proportion of species of the top sets (congruence)"),
    "  # --- species groups (any categorical attribute that is not a trait) ---",
    "  # col: grouping column; sep: separator of multiple labels; recode: merged labels;",
    "  # levels: groups analysed; reference: reference group(s) of the species-level analysis",
    paste0("  groups          = ", nn, ",  # NULL: not analysed"),
    "  # continuous variable that is not a trait (col, transform, fit_axes, nperm)",
    paste0("  covariate       = ", lst(p$covariate), ",  # NULL: not analysed"),
    "  # extinction scenarios (metrics step, IUCN column): steps, nrep, hull_axes",
    paste0("  extinction      = ", lst(p$extinction), "  # NULL: not analysed"),
    ")")
}

.project_settings_md <- function(p, steps) {
  b <- function(k, v) sprintf("- **%s:** %s", k, v)
  sp <- paste(sprintf("%s (%s)", names(p$spaces), vapply(p$spaces, paste, "", collapse = " + ")),
              collapse = "; ")
  paste(c(
    b("Missing values", switch(p$missing, keep = "kept (pairwise Gower distance)",
                               drop = "incomplete species dropped",
                               impute = sprintf("imputed with missForest (%d trees%s%s)", p$ntree,
                                 if (length(p$cov_cols)) paste0(", predictors: ", paste(p$cov_cols, collapse = ", ")) else "",
                                 if (!is.null(p$phylo_file)) sprintf(", %d phylogenetic eigenvectors", p$n_phylo_axes) else ""))),
    b("Spaces", sp),
    b("Distance", sprintf("generalised Gower (%s engine), %s weighting", p$engine, p$weighting)),
    b("PCoA", sprintf("correction %s, %s, %s eigendecomposition", p$correction,
                      if (is.null(p$max_axes)) "all positive axes" else paste(p$max_axes, "axes"), p$eigen)),
    if ("compare" %in% steps) b("Procrustes tests", sprintf("%s, %d permutations, %s%s", p$method, p$nrepet,
                                if (is.null(p$n_axes)) "all axes" else paste(p$n_axes, "axes"),
                                if (p$mantel) ", + Mantel tests" else "")),
    if ("metrics" %in% steps) b("Species metrics", sprintf("k = %d nearest neighbours, top %g%%%s", p$k,
                                100 * p$top, if (!is.null(p$iucn_col)) ", FUSE from IUCN categories" else "")),
    if (!is.null(p$groups)) b("Species groups", if (p$groups$source == "simulate")
      sprintf("simulated native / non-native status (%g%%, %s bias, %s space), reference: native",
              100 * p$groups$prop, p$groups$bias, p$groups$sim_space) else
      sprintf("%s '%s'%s%s, reference: %s, null model: %s",
              if (length(p$groups$col) > 1) "indicator columns" else "column",
              paste(p$groups$col, collapse = "', '"),
              if (length(p$groups$levels)) sprintf(" (%s)", paste(p$groups$levels, collapse = ", ")) else "",
              if (!is.null(p$groups$sep)) sprintf(", labels separated by '%s'", p$groups$sep) else "",
              if (is.null(p$groups$reference)) "none" else paste(p$groups$reference, collapse = " + "),
              p$groups$null)),
    if (!is.null(p$groups) && isTRUE(p$groups$density)) b("Densities", sprintf(
      "kernel densities on %d axes (bandwidth x %g%s)", p$groups$density_axes, p$groups$bw,
      if (!is.null(p$groups$abundance_col)) sprintf(", weighted by '%s'", p$groups$abundance_col) else "")),
    if (!is.null(p$groups$sites)) b("Assemblages", sprintf("data/raw/assemblages.csv (columns '%s', '%s'%s)",
      p$groups$sites$site_col, p$groups$sites$species_col,
      if (!is.null(p$groups$sites$group_col)) sprintf(", groups per assemblage in '%s'", p$groups$sites$group_col) else "")),
    if (!is.null(p$covariate)) b("Continuous variable", sprintf("'%s' (transform: %s; fit on %d axes, %d permutations)",
      p$covariate$col, p$covariate$transform, p$covariate$fit_axes, p$covariate$nperm)),
    if (!is.null(p$extinction)) b("Extinction scenarios", sprintf("%d steps, %d replicates, hulls on %d axes",
      p$extinction$steps, p$extinction$nrep, p$extinction$hull_axes)),
    b("Random seed", p$seed)), collapse = "\n")
}

# data/metadata.csv: definition of the variables of the data and result tables.
.project_metadata <- function() {
  v <- function(file, variable, description)
    data.frame(file = file, variable = variable, description = description, stringsAsFactors = FALSE)
  rbind(
    v("data/raw/trait_info.csv", c("trait", "type", "dimension", "block", "transform", "weight",
                                   "levels", "circ_min", "circ_max", "use"),
      c("column of the trait table", "Q continuous, O ordinal, N nominal, D binary, B multichoice, F fuzzy-coded, C circular",
        "functional dimension", "columns of one fuzzy-coded or multichoice trait share a block",
        "transformation of continuous traits", "relative weight in the distance",
        "ordered levels of an ordinal trait (a;b;c)", "lower limit of a circular trait",
        "upper limit of a circular trait", "trait used")),
    v("results/tables/Table2_spaces.csv", c("space", "dimensions", "n_species", "n_variables",
                                            "PCoA1_pct ... PCoA4_pct", "cum4_pct"),
      c("functional space", "dimensions combined", "species in the space", "trait variables",
        "inertia of each of the first PCoA axes (%)", "cumulative inertia of the first four axes (%)")),
    v("results/tables/TableS3_pcoa_coordinates.csv", c("space", "species", "PCoA1 ..."),
      c("functional space", "species", "coordinates on the PCoA axes")),
    v("results/tables/Table3_procrustes.csv", c("space1, space2", "n_species", "n_axes", "r", "m2",
                                                "p_value", "signif"),
      c("spaces compared", "species shared", "axes used in each space", "Procrustes correlation",
        "Procrustes sum of squares (vegan only)", "permutation p-value", "significance code")),
    v("results/tables/TableS5_species_metrics.csv", c("FUn", "FSp", "FUn_raw", "FSp_raw",
                                                      "n_identical", "iucn", "P", "FUSE"),
      c("functional uniqueness: sum of distances to the k nearest neighbours, rescaled to [0, 1] within the space",
        "functional specialisation: distance to the centroid, rescaled to [0, 1]",
        "FUn before rescaling", "FSp before rescaling", "other species with the same position",
        "IUCN category", "extinction probability of the category (EDGE2)",
        "log(1 + P x FSp + P x FUn)")),
    v("results/tables/Table4_metric_congruence.csv", c("metric", "spearman", "jaccard_top", "n_shared_top"),
      c("species metric compared", "Spearman correlation between the rankings in two spaces",
        "Jaccard overlap of the top sets", "species in both top sets")),
    v("results/tables/Table5_groups_summary.csv", c("permanova_R2, permanova_p", "permdisp_p",
                                                    "overlap_median, overlap_min", "segregation_median"),
      c("share of variance explained by the groups and permutation p-value (PERMANOVA)",
        "test of homogeneity of dispersions (PERMDISP)",
        "median and minimum hull overlap (intersection / union) of the pairs of groups",
        "median ratio of observed / expected same-group neighbours")),
    v("results/tables/TableS6_groups.csv", c("n", "FRic", "FDis", "FOri", "MPD", "p_outside_others",
                                             "prop_V_exclusive", "knn_same", "segregation"),
      c("species of the group", "hull volume relative to all species", "mean distance to the group centroid",
        "mean nearest-neighbour distance", "mean pairwise distance",
        "proportion of species outside the hull of the other species",
        "share of the group hull not shared with the other species",
        "mean proportion of the k nearest neighbours in the same group",
        "knn_same / its expectation under random labels")),
    v("results/tables/Table6_group_pairs.csv", c("V1, V2, V_inter", "overlap", "beta_jac, beta_jtu, beta_jne",
                                                 "p1_in_hull2, p2_in_hull1", "centroid_dist", "ses_*, p_*"),
      c("hull volumes of the two groups and of their intersection (relative to all species)",
        "intersection / union", "functional Jaccard dissimilarity, turnover and nestedness components",
        "proportion of the species of one group inside the hull of the other",
        "distance between group centroids", "standardised effect size and p-value (label permutations)")),
    v("results/tables/TableS7_focal_species.csv", c("FNov", "FSp_ref", "inside_ref_hull", "nearest_ref"),
      c("novelty: distance to the nearest reference species (rescaled)",
        "distance to the reference centroid (rescaled)", "inside the convex hull of reference species",
        "nearest reference species")),
    v("results/tables/TableS8_reference_species.csv", c("dFUn_rel", "n_focal_knn"),
      c("relative change in FUn of reference species once focal species are added",
        "focal species among the k nearest neighbours")),
    v("results/tables/TableS9_group_contribution.csv", c("FRic_*", "FDis_*", "FOri_*", "MPD_*",
                                                         "d_*", "ses_*, p_*"),
      c("functional richness of the reference group alone (ref) and with the focal group (all)",
        "functional dispersion", "functional originality", "mean pairwise distance",
        "contribution of the focal group (all - ref)", "standardised effect size and p-value against the null model")),
    v("results/tables/TableS13_group_assemblages.csv", c("assemblage, group", "*_ref, *_all, d_*",
                                                         "overlap", "p_focal_outside"),
      c("assemblage and focal group", "metrics of the reference species alone, with the focal species, and contribution",
        "hull overlap of focal and reference species in the assemblage",
        "proportion of focal species outside the reference hull")),
    v("results/tables/Table8_extinction_scenarios.csv", c("loss_risk", "loss_random", "ses, p"),
      c("relative loss after removing the threatened species by decreasing extinction probability",
        "same number of species removed at random", "standardised effect size and one-sided p-value")),
    v("results/tables/TableS11_stochastic_extinctions.csv", c("n_lost_mean", "loss", "loss_random"),
      c("mean number of species lost (EDGE2 probabilities)", "relative loss", "loss for as many random species")),
    v("results/tables/Table9_continuous_variable.csv", c("R2, R2_adj, p_fit", "knn_cor, p_knn",
                                                         "rho_FUn, rho_FSp"),
      c("linear fit of the variable on the first axes, permutation p-value",
        "correlation with the mean of the k nearest neighbours (functional autocorrelation)",
        "Spearman correlations with functional uniqueness and specialisation")))
}

# An ASCII tree of a folder (for the Shiny app).
.project_tree <- function(dir, max_files = 400) {
  root <- basename(dir)
  utf8 <- isTRUE(l10n_info()[["UTF-8"]])
  g <- if (utf8) c(mid = "\u251c\u2500\u2500 ", end = "\u2514\u2500\u2500 ", bar = "\u2502   ") else
    c(mid = "|-- ", end = "`-- ", bar = "|   ")
  walk <- function(d, prefix) {
    items <- sort(list.files(d, all.files = TRUE, no.. = TRUE))
    isdir <- dir.exists(file.path(d, items))
    items <- c(items[isdir], items[!isdir]); isdir <- c(isdir[isdir], isdir[!isdir])
    out <- character(0)
    for (i in seq_along(items)) {
      last <- i == length(items)
      out <- c(out, paste0(prefix, if (last) g[["end"]] else g[["mid"]],
                           items[i], if (isdir[i]) "/" else ""))
      if (isdir[i]) out <- c(out, walk(file.path(d, items[i]),
                                       paste0(prefix, if (last) "    " else g[["bar"]])))
    }
    out
  }
  lines <- c(paste0(root, "/"), walk(dir, ""))
  if (length(lines) > max_files) lines <- c(utils::head(lines, max_files), "...")
  paste(lines, collapse = "\n")
}
