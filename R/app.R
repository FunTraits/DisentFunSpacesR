#' Launch the DisentFunSpaces Shiny application
#'
#' Interactive version of the workflow: load a trait table, declare the type
#' and functional dimension of each trait, handle missing values, build the
#' functional spaces, compare them with Procrustes tests and compute
#' species-level metrics. Results, the trait-information table and an R
#' script reproducing the analysis can be downloaded.
#'
#' @param launch.browser Passed to [shiny::runApp()].
#' @param max_upload_mb Maximum upload size (MB).
#' @details Large datasets: from 1000 species the app uses, by default
#'   ("Performance" panel of tab 4), the fast distance engine (same values as
#'   `ade4::dist.ktab()`, much less memory) and a partial eigendecomposition
#'   (package 'RSpectra') that computes only the axes kept - by default all
#'   positive-eigenvalue axes, as in the published analysis. From 3000
#'   species (`options(DisentFunSpaces.app_large_n)`) the app shows memory
#'   estimates and advice. Text files are read with
#'   `data.table::fread()` when installed. Procrustes permutations can use
#'   several cores (macOS / Linux) without changing the results. Computations
#'   on matrices (PCoA, Procrustes) are much faster with an optimised BLAS
#'   (e.g. Accelerate on macOS, OpenBLAS on Linux).
#' @param ... Further arguments to [shiny::runApp()].
#' @return Called for its side effect.
#' @export
#' @examples
#' if (interactive()) run_app()
run_app <- function(launch.browser = TRUE, max_upload_mb = 1000, ...) {
  old <- options(shiny.maxRequestSize = max_upload_mb * 1024^2)
  on.exit(options(old), add = TRUE)
  shiny::runApp(fs_app(), launch.browser = launch.browser, ...)
}

# Large-data settings used by the app ----------------------------------------
# number of species from which the app treats a dataset as "large" (fast
# distance engine, partial PCoA with a limited number of axes)
.app_large_n <- function() getOption("DisentFunSpaces.app_large_n", 3000)
.app_cores <- function() {
  if (.Platform$OS.type != "unix") return(1L)
  n <- tryCatch(parallel::detectCores(logical = FALSE), error = function(e) NA)
  if (is.na(n)) n <- 1L
  as.integer(max(1L, min(8L, n - 1L)))
}
# reference (unoptimised) BLAS in use? PCoA and Procrustes tests are then
# several times slower than with an optimised BLAS (Accelerate, OpenBLAS, MKL)
.reference_blas <- function() {
  b <- tryCatch(utils::sessionInfo()$BLAS, error = function(e) "")
  length(b) == 1 && nzchar(b) && grepl("libRblas", b, fixed = TRUE)
}
# guess the field separator (and decimal mark) of a delimited text file from
# its first lines: the candidate (";", tab, ",", "|") that splits every line
# into the same number (> 1) of fields, the most frequent one winning
.detect_sep <- function(path, n = 20L) {
  l <- tryCatch(readLines(path, n = n, warn = FALSE, encoding = "UTF-8"),
                error = function(e) character(0))
  l <- sub("^\ufeff", "", l)
  l <- l[nzchar(trimws(l))]
  if (!length(l)) return(list(sep = ",", dec = "."))
  strip <- function(x) gsub('"[^"]*"', "", x)          # ignore quoted fields
  l2 <- strip(l)
  cand <- c(";", "\t", ",", "|")
  cnt <- vapply(cand, function(ch) {
    k <- lengths(regmatches(l2, gregexpr(ch, l2, fixed = TRUE)))
    if (k[1] > 0 && length(unique(k)) == 1L) k[1] else if (k[1] > 0) k[1] / 2 else 0
  }, numeric(1))
  sep <- if (max(cnt) > 0) cand[which.max(cnt)] else ","
  # decimal comma (French / European CSV): numbers like 0,25 between separators
  dec <- "."
  if (sep != ",") {
    body <- l2[-1]
    if (length(body) &&
        any(grepl("(^|[^0-9.])[0-9]+,[0-9]+($|[^0-9])", body)) &&
        !any(grepl("(^|[^0-9,])[0-9]+\\.[0-9]+($|[^0-9])", body))) dec <- ","
  }
  list(sep = sep, dec = dec)
}

# fast reader for large text files ('data.table' if installed). `sep = "auto"`
# guesses the separator; a table read as a single column with the separator
# chosen by the user is read again with the guessed one.
.read_table_fast <- function(path, sep = "auto", dec = NULL) {
  guess <- .detect_sep(path)
  if (identical(sep, "auto") || is.null(sep) || !nzchar(sep)) sep <- guess$sep
  if (is.null(dec)) dec <- if (identical(sep, guess$sep)) guess$dec else "."
  rd <- function(sep, dec) {
    if (requireNamespace("data.table", quietly = TRUE)) {
      d <- data.table::fread(path, sep = sep, dec = dec, header = TRUE,
                             na.strings = c("NA", ""), check.names = FALSE,
                             data.table = FALSE, showProgress = FALSE,
                             integer64 = "double", quote = "\"", encoding = "UTF-8")
    } else {
      d <- utils::read.table(path, header = TRUE, sep = sep, dec = dec, quote = "\"",
                             stringsAsFactors = FALSE, check.names = FALSE,
                             na.strings = c("NA", ""), comment.char = "",
                             fileEncoding = "UTF-8")
    }
    names(d) <- sub("^\ufeff", "", names(d))
    d
  }
  d <- rd(sep, dec)
  if (ncol(d) == 1L && !identical(sep, guess$sep)) {
    sep <- guess$sep; dec <- guess$dec
    d <- rd(sep, dec)
  }
  attr(d, "sep") <- sep; attr(d, "dec") <- dec
  d
}

#' @rdname run_app
#' @return `fs_app()` returns a `shiny.appobj`.
#' @export
fs_app <- function() shiny::shinyApp(ui = app_ui(), server = app_server)

# Small UI helpers ------------------------------------------------------------
.card_x <- function(title, ...) bslib::card(fill = FALSE, bslib::card_header(title),
                                            bslib::card_body(..., fillable = FALSE))
.CITATION_TEXT <- paste(
  "Toussaint A., Tedesco P., Grenouillet G., Kasari-Toussaint L., Brosse S. (2026)",
  "Disentangling functional spaces: a multidimensional framework for ecology and",
  "conservation. Ecography e08578. doi:10.1002/ecog.08578")
# column names proposed as grouping variables (not traits)
.GROUP_COL_GUESS <- c("status", "origin", "native_status", "nativeness", "alien", "introduced",
                      "order", "family", "use", "uses", "human_use", "migration", "migratory",
                      "activity", "diel", "diurnal", "group", "groups",
                      "fisheries", "fishery", "aquaculture", "aquarium", "game fish", "gamefish",
                      "game_fish", "sport fish", "bait", "commercial", "all uses", "all_uses",
                      "any use", "used", "harvested", "threatened")
# the steps of the workflow: tab value, label, what it does
.APP_STEPS <- data.frame(
  tab = c("1. Data", "2. Traits", "3. Missing values", "4. Functional spaces",
          "5. Compare spaces", "6. Species metrics", "7. Species groups", "8. Export"),
  label = c("Data", "Traits", "Missing values", "Functional spaces", "Compare spaces",
            "Species metrics", "Species groups", "Export"),
  what = c("A trait table: one row per species, one column per trait (or per modality).",
           "The type and the functional dimension of each trait.",
           "Keep them (pairwise Gower distance), drop incomplete species, or impute them.",
           "One space per dimension and per combination: generalised Gower distance + PCoA.",
           "Procrustes tests: do the dimensions place species alike?",
           "Functional uniqueness, specialisation and FUSE in every space.",
           "Groups of species (origin, uses, migration...): overlap, novelty, contribution (optional).",
           "A reproducible project: data, scripts, results and report in one .zip."),
  optional = c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE),
  stringsAsFactors = FALSE)

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------
app_ui <- function() {
  sh <- shiny::tags
  type_choices <- stats::setNames(trait_types()$code,
                                  paste0(trait_types()$code, " - ", trait_types()$label))
  bslib::page_navbar(
    title = sh$span(.brand_mark(), "DisentFunSpaces"),
    id = "nav",
    navbar_options = bslib::navbar_options(bg = "#ffffff"),
    selected = "home",
    theme = .app_theme(),
    fillable = FALSE,
    window_title = "DisentFunSpaces",
    header = .app_deps(),

    # HOME ----------------------------------------------------------------------
    bslib::nav_panel("Home", value = "home",
      sh$div(class = "home-hero",
        sh$div(
          sh$div(class = "home-kicker", "Toussaint et al. 2026 \u00b7 Ecography"),
          sh$h1("Disentangling functional spaces"),
          sh$p(class = "home-lead",
            "Build one functional space per functional dimension of your traits (locomotion,",
            " diet, reproduction...) and for their combinations, compare them with Procrustes",
            " tests, and find the species that are functionally unique, specialised or",
            " endangered in each. Any mixed-type trait table: continuous, ordinal, nominal,",
            " binary, multichoice, fuzzy-coded or circular traits."),
          shiny::uiOutput("home_state"),
          sh$div(class = "home-actions",
            shiny::actionButton("home_example", "Start with the example: 2,000 bird species",
                                class = "btn-primary"),
            shiny::actionButton("home_upload", "Upload a trait table"),
            shiny::uiOutput("home_next", inline = TRUE))),
        sh$div(.home_art_svg())),
      bslib::layout_columns(col_widths = c(7, 5), fill = FALSE, style = "align-items: start;",
        sh$div(class = "home-col",
          .card_x("Workflow", shiny::uiOutput("home_steps"))),
        sh$div(class = "home-col",
          .card_x("This session", shiny::uiOutput("home_stats")),
          .card_x("From the R console",
            sh$div(class = "home-code", sh$pre(paste(
              "library(DisentFunSpaces)",
              "tr <- prepare_traits(birds_traits, birds_trait_info, \"species\")",
              "spaces <- build_spaces(drop_incomplete(tr))",
              "compare_spaces(spaces)",
              "",
              "# the whole analysis as a reproducible project",
              "write_analysis_project(birds_traits, birds_trait_info,",
              "  dir = \"~/birds_spaces\", params = list(missing = \"drop\"))",
              sep = "\n"))),
            sh$p(class = "hint", "Every tab calls these package functions; the Export tab writes",
                 " them into a project that runs on its own.")),
          .card_x("Citation",
            sh$p(class = "home-cite", .CITATION_TEXT),
            sh$p(class = "hint", sh$a("doi:10.1002/ecog.08578", href = "https://doi.org/10.1002/ecog.08578",
                                      target = "_blank"), " \u00b7 ",
                 sh$a("Code and data of the paper", href = "https://github.com/FunTraits/DisentFunSpaces",
                      target = "_blank"))))
      )
    ),

    # 1. DATA -----------------------------------------------------------------
    bslib::nav_panel(.step_title(1, "Data"), value = "1. Data",
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(width = 320,
          shiny::radioButtons("src", "Data source",
            c("Example: 2000 bird species (fast)" = "example",
              "Example: all 8,413 bird species of Toussaint et al. 2026" = "example_all",
              "Upload a file" = "upload")),
          shiny::conditionalPanel("input.src == 'example_all'",
            shiny::tags$p(class = "step-help",
              "Full dataset of the Ecography paper (22 % missing values, impute them in tab 3).",
              " With the default large-data settings of tab 4 (fast distance, partial PCoA),",
              " each space is built in less than a minute.")),
          shiny::conditionalPanel("input.src == 'upload'",
            shiny::fileInput("file", "Trait table (.csv, .tsv, .txt, .xlsx; up to 1 GB)",
                             accept = c(".csv", ".tsv", ".txt", ".xlsx", ".xls")),
            shiny::selectInput("sep", "Separator (text files)",
                               c("Auto-detect" = "auto", "Comma" = ",", "Semicolon" = ";",
                                 "Tab" = "\t", "Pipe" = "|"), selected = "auto")),
          shiny::selectInput("species_col", "Species column (a species may be repeated: one row per record)",
                             choices = NULL),
          shiny::selectInput("iucn_col", "IUCN category column (optional, for FUSE)", choices = NULL),
          shiny::selectizeInput("group_cols", "Other variables, not traits (optional: groups, continuous attributes)", choices = NULL,
                                multiple = TRUE,
                                options = list(placeholder = "e.g. status, use, migration")),
          sh$p(class = "step-help",
               "One row per species, one column per trait (or per modality for fuzzy-coded",
               " and multichoice traits). Intraspecific records must be aggregated beforehand.",
               " Declare here every column that is not a trait (e.g. 0 / 1 columns of human uses):",
               " otherwise it is used as a binary trait to build the spaces.")
        ),
        bslib::card(bslib::card_header("Preview"), shiny::uiOutput("data_info"),
                    shiny::uiOutput("perf_info"),
                    DT::DTOutput("data_preview"))
      )
    ),

    # 2. TRAITS ---------------------------------------------------------------
    bslib::nav_panel(.step_title(2, "Traits"), value = "2. Traits",
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(width = 320,
          shiny::actionButton("guess", "Guess trait types", class = "btn-outline-primary w-100"),
          shiny::fileInput("ti_file", "...or load a trait-information CSV", accept = ".csv"),
          shiny::hr(),
          shiny::textInput("bulk_dim", "Assign dimension to selected rows", placeholder = "e.g. Diet"),
          shiny::selectizeInput("bulk_rows", "Rows (traits)", choices = NULL, multiple = TRUE),
          shiny::actionButton("bulk_apply", "Apply", class = "btn-sm btn-outline-secondary"),
          shiny::hr(),
          shiny::actionButton("validate", "Validate traits", class = "btn-primary w-100"),
          shiny::br(), shiny::br(),
          shiny::downloadButton("dl_ti", "Download trait-information CSV", class = "btn-sm w-100")
        ),
        bslib::navset_card_tab(
          bslib::nav_panel("Trait description",
            sh$p(class = "step-help",
              "For each column: tick 'use', choose the type, and give the functional dimension",
              " (e.g. Locomotion, Diet, Reproduction). Columns of a fuzzy-coded (F) or multichoice (B)",
              " trait must share the same 'block' name. 'levels' gives the order of an ordinal trait",
              " stored as text (a;b;c). Circular traits need the cycle limits (e.g. 0 and 24 for hours)."),
            sh$div(style = "max-height: 62vh; overflow-y: auto;", shiny::uiOutput("ti_editor"))),
          bslib::nav_panel("Validation", shiny::uiOutput("validation_msg"),
                           DT::DTOutput("trait_summary")),
          bslib::nav_panel("Trait types", DT::DTOutput("types_table"))
        )
      )
    ),

    # 3. MISSING --------------------------------------------------------------
    bslib::nav_panel(.step_title(3, "Missing"), value = "3. Missing values",
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(width = 320,
          shiny::radioButtons("missing", "Handling of missing values", c(
            "Keep (pairwise Gower distance)" = "keep",
            "Drop incomplete species" = "drop",
            "Impute with missForest" = "impute"), selected = "keep"),
          shiny::conditionalPanel("input.missing == 'impute'",
            shiny::selectizeInput("cov_cols", "Predictor columns from the data table (optional)",
              choices = NULL, multiple = TRUE,
              options = list(placeholder = "e.g. Eigen.1, Eigen.2 ...")),
            shiny::tags$p(class = "step-help",
              "Columns of the initial table used only to predict missing values (e.g. phylogenetic",
              " eigenvectors, as in the published analysis, or taxonomy): they are not traits and do",
              " not enter functional distances. Columns matching 'eigen' are proposed by default."),
            shiny::fileInput("phylo", "...and/or a phylogeny (Newick, optional)",
                             accept = c(".tre", ".nwk", ".newick", ".txt")),
            shiny::numericInput("n_phylo", "Eigenvectors computed from the phylogeny", 10, 1, 50),
            shiny::numericInput("ntree", "Trees per forest", 100, 10, 1000),
            shiny::numericInput("seed", "Random seed", 20251029)),
          shiny::actionButton("run_missing", "Apply", class = "btn-primary w-100")
        ),
        bslib::card(bslib::card_header("Missing values per trait"), shiny::uiOutput("missing_msg"),
                    DT::DTOutput("missing_table"))
      )
    ),

    # 4. SPACES ---------------------------------------------------------------
    bslib::nav_panel(.step_title(4, "Spaces"), value = "4. Functional spaces",
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(width = 330,
          shiny::checkboxGroupInput("spaces_sel", "Spaces to build", choices = NULL),
          bslib::accordion(open = FALSE,
            bslib::accordion_panel("Add a custom combination",
              shiny::selectizeInput("custom_dims", "Dimensions", choices = NULL, multiple = TRUE),
              shiny::textInput("custom_name", "Name", placeholder = "e.g. Diet+Locomotion"),
              shiny::actionButton("add_custom", "Add", class = "btn-sm btn-outline-secondary"))),
          shiny::radioButtons("weighting", "Weighting in combined spaces", c(
            "Each trait equally (published analysis)" = "variable",
            "Each dimension equally" = "dimension",
            "Each dimension equally, rescaled (types balanced)" = "dimension_scaled")),
          shiny::selectInput("correction", "Negative eigenvalues correction",
                             c("None (published analysis)" = "none", "Cailliez" = "cailliez", "Lingoes" = "lingoes")),
          shiny::numericInput("max_axes", "Max. PCoA axes kept (empty = all, as published)",
                              value = NA, min = 2, step = 1),
          bslib::accordion(open = FALSE,
            bslib::accordion_panel("Performance (large datasets)",
              shiny::selectInput("engine", "Distance computation", c(
                "Auto (fast from 1000 species)" = "auto",
                "ade4::dist.ktab (published analysis)" = "ade4",
                "Fast (same values, low memory)" = "fast")),
              shiny::selectInput("eigen_mode", "PCoA eigendecomposition", c(
                "Auto (partial from 1000 species)" = "auto",
                "Full (published analysis)" = "full",
                "Partial: axes kept only (RSpectra)" = "partial")),
              shiny::tags$p(class = "step-help",
                "Fast distances and the partial PCoA give the same distances, coordinates and",
                " inertia as the published analysis (signs of axes are arbitrary). The partial PCoA",
                " computes all positive axes (or the max. number of axes set above) without the full",
                " eigendecomposition; it does not report negative eigenvalues. Requires 'RSpectra'."))),
          shiny::tags$p(class = "step-help",
            "Distances: ade4::dist.ktab (one table per dimension x type); PCoA: ape::pcoa (PCoAGraph)."),
          shiny::actionButton("run_spaces", "Build spaces", class = "btn-primary w-100")
        ),
        shiny::uiOutput("spaces_msg"),
        DT::DTOutput("spaces_table"),
        bslib::layout_columns(col_widths = c(8, 4),
          bslib::card(full_screen = TRUE,
            bslib::card_header(bslib::layout_columns(col_widths = c(4, 2, 2, 4),
              shiny::selectInput("plot_space", NULL, choices = NULL),
              shiny::numericInput("ax1", "x", 1, 1, 20), shiny::numericInput("ax2", "y", 2, 1, 20),
              shiny::sliderInput("min_cor", "Arrows |r| >=", 0, 1, 0.3, 0.05))),
            shiny::checkboxInput("color_iucn", "Colour species by IUCN category", FALSE),
            .fig_output("space_plot", height = "520px"),
            shiny::downloadButton("dl_space_plot", "PDF", class = "btn-sm")),
          bslib::navset_card_tab(
            bslib::nav_panel("Inertia", .fig_output("eigen_plot", height = "300px")),
            bslib::nav_panel("Quality", .fig_output("quality_plot", height = "300px"),
                             sh$p(class = "step-help", "mAD: Maire et al. (2015) Global Ecol. Biogeogr.")),
            bslib::nav_panel("Trait correlations", DT::DTOutput("cor_table"))
          )
        )
      )
    ),

    # 5. PROCRUSTES -----------------------------------------------------------
    bslib::nav_panel(.step_title(5, "Compare"), value = "5. Compare spaces",
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(width = 320,
          shiny::radioButtons("pr_method", "Procrustes test", c(
            "ade4::procuste.rtest (published analysis)" = "ade4",
            "vegan::protest (symmetric)" = "vegan")),
          shiny::numericInput("nrepet", "Permutations", 999, 9, 99999),
          shiny::radioButtons("axes_mode", "PCoA axes used",
                              c("All positive axes (published)" = "all", "First k axes" = "k")),
          shiny::conditionalPanel("input.axes_mode == 'k'",
            shiny::numericInput("n_axes", "k", 4, 1, 50)),
          shiny::checkboxInput("mantel", "Add Mantel tests on distances", FALSE),
          if (.Platform$OS.type == "unix")
            shiny::numericInput("cores", "CPU cores for permutations (same results)", .app_cores(), 1, 64),
          shiny::actionButton("run_compare", "Run comparisons", class = "btn-primary w-100")
        ),
        bslib::layout_columns(col_widths = c(6, 6),
          bslib::card(full_screen = TRUE, bslib::card_header("Procrustes correlations"),
                      .fig_output("cmp_plot", height = "420px"),
                      shiny::downloadButton("dl_cmp_plot", "PDF", class = "btn-sm")),
          bslib::card(bslib::card_header("Tests"), DT::DTOutput("cmp_table"))
        ),
        bslib::card(full_screen = TRUE,
          bslib::card_header("Species-level Procrustes residuals"),
          bslib::layout_columns(col_widths = c(3, 3, -6),
            shiny::selectInput("res_a", "Space A", choices = NULL),
            shiny::selectInput("res_b", "Space B", choices = NULL)),
          bslib::layout_columns(col_widths = c(7, 5),
            .fig_output("res_plot", height = "460px"),
            DT::DTOutput("res_table")))
      )
    ),

    # 6. METRICS --------------------------------------------------------------
    bslib::nav_panel(.step_title(6, "Metrics"), value = "6. Species metrics",
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(width = 320,
          shiny::numericInput("k_nn", "k nearest neighbours (FUn)", 5, 1, 50),
          shiny::selectInput("metric", "Metric compared across spaces", c("FUn", "FSp")),
          shiny::sliderInput("top", "Top proportion of species", 0.01, 0.5, 0.1, 0.01),
          shiny::actionButton("run_metrics", "Compute", class = "btn-primary w-100"),
          sh$p(class = "step-help",
            "FUn: sum of distances to the k nearest neighbours; FSp: distance to the centroid;",
            " both rescaled to [0, 1] within each space. FUSE = log(1 + P x FSp + P x FUn) with",
            " P the EDGE2 extinction probability of the IUCN category.")
        ),
        bslib::layout_columns(col_widths = c(6, 6),
          bslib::card(bslib::card_header("Congruence of species rankings between spaces"),
                      DT::DTOutput("metric_pairs")),
          bslib::card(full_screen = TRUE,
            bslib::card_header(bslib::layout_columns(col_widths = c(6, 6),
              shiny::selectInput("m_a", NULL, choices = NULL),
              shiny::selectInput("m_b", NULL, choices = NULL))),
            .fig_output("metric_scatter", height = "380px"))
        ),
        bslib::card(bslib::card_header("Species metrics (all spaces)"), DT::DTOutput("metric_table")),
        bslib::card(full_screen = TRUE,
          bslib::card_header(bslib::layout_columns(col_widths = c(4, 2, 2, 2, 2),
            sh$span("Extinction scenarios (needs the IUCN column)"),
            shiny::numericInput("ext_steps", "Steps", 10, 2, 50),
            shiny::numericInput("ext_nrep", "Replicates", 49, 9, 9999),
            shiny::selectInput("ext_metric", "Metric", c("FRic", "FDis", "FOri", "MPD")),
            shiny::actionButton("run_ext", "Run", class = "btn-primary btn-sm w-100", style = "margin-top: 1.9rem;"))),
          sh$p(class = "step-help",
            "Species are removed by decreasing extinction probability (CR, EN, VU; ties at random)",
            " or at random, up to the number of threatened species; the stochastic scenario draws",
            " extinctions with the EDGE2 probabilities. Metrics are relative to all species",
            " (hulls on the axes set in tab 7)."),
          .fig_output("ext_plot", height = "330px"),
          bslib::layout_columns(col_widths = c(6, 6),
            sh$div(sh$h6("Loss after removing the threatened species"), DT::DTOutput("ext_summary")),
            sh$div(sh$h6("Stochastic extinctions (EDGE2 probabilities)"), DT::DTOutput("ext_stoch"))))
      )
    ),

    # 7. SPECIES GROUPS -------------------------------------------------------
    bslib::nav_panel(.step_title(7, "Groups"), value = "7. Species groups",
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(width = 330,
          shiny::radioButtons("grp_src", "Non-trait variable",
            c("Categories (groups of species)" = "column",
              "Indicator columns (0 / 1, one group per column)" = "binary",
              "Continuous variable" = "continuous",
              "Simulate native / non-native status (exploration)" = "simulate")),
          shiny::conditionalPanel("input.grp_src == 'binary'",
            shiny::selectizeInput("grp_bin_cols", "Indicator columns (one group each)", choices = NULL,
                                  multiple = TRUE, options = list(placeholder = "declare the 0 / 1 columns in tab 1")),
            shiny::textInput("grp_none", "Group of species set to 0 in all columns", "none",
                             placeholder = "empty: leave them out"),
            sh$p(class = "step-help",
                 "E.g. one column per human use (Fisheries, Aquaculture, Aquarium...): a species belongs",
                 " to every group set to 1 (groups are not exclusive) and to the group above if set to 0",
                 " everywhere (natural reference: unused species). A column that is the union of the",
                 " others ('All uses') is redundant with them: analyse it alone (used vs none) rather",
                 " than with the other columns."),
            shiny::uiOutput("grp_bin_info")),
          shiny::conditionalPanel("input.grp_src == 'column' || input.grp_src == 'continuous'",
            shiny::selectInput("grp_col", "Variable", choices = NULL),
            sh$p(class = "step-help",
                 "Any attribute that is not a trait: origin, human use, migration, activity period,",
                 " threat status, taxonomy, range size, abundance... Declare the columns in tab 1",
                 " (the IUCN column is also offered).")),
          shiny::conditionalPanel("input.grp_src == 'column'",
            shiny::textInput("grp_sep", "Separator of multiple labels (optional)", "",
                             placeholder = "e.g. ;  for 'food;aquarium'")),
          shiny::conditionalPanel("input.grp_src == 'column' || input.grp_src == 'binary'",
            shiny::textAreaInput("grp_recode", "Merge labels (optional, one group per line)", "",
                                 rows = 2, placeholder = "threatened = VU, EN, CR\nother = LC, NT"),
            shiny::uiOutput("grp_levels_info"),
            shiny::selectizeInput("grp_levels", "Groups analysed", choices = NULL, multiple = TRUE)),
          shiny::conditionalPanel("input.grp_src == 'continuous'",
            shiny::selectInput("cov_transform", "Transformation",
                               c("None" = "none", "log10" = "log10", "log10(1 + x)" = "log1p", "Ranks" = "rank")),
            shiny::numericInput("cov_fit_axes", "Axes of the linear fit", 4, 1, 20),
            shiny::numericInput("cov_nperm", "Permutations", 999, 9, 99999)),
          shiny::conditionalPanel("input.grp_src == 'simulate'",
            shiny::sliderInput("nn_prop", "Proportion of non-native species", 0.02, 0.5, 0.1, 0.01),
            shiny::selectInput("nn_bias", "Position of introduced species",
                               c("random", "peripheral", "central")),
            shiny::selectInput("nn_sim_space", "Space used for the bias", choices = NULL)),
          shiny::conditionalPanel("input.grp_src != 'continuous'",
            shiny::selectizeInput("grp_ref", "Reference group (species level; optional)",
                                  choices = NULL, multiple = TRUE,
                                  options = list(placeholder = "none: compare groups only")),
            shiny::conditionalPanel("input.grp_ref && input.grp_ref.length > 0",
              shiny::checkboxInput("grp_pool", "Pool the other groups into one focal group", FALSE)),
            shiny::numericInput("grp_hull_axes", "Axes for convex hulls", 4, 2, 8),
            shiny::checkboxInput("grp_density", "Also overlap of kernel densities (robust to outliers)", FALSE),
            shiny::conditionalPanel("input.grp_density",
              bslib::layout_columns(col_widths = c(6, 6),
                shiny::numericInput("grp_density_axes", "Axes (2-4)", 2, 2, 4),
                shiny::numericInput("grp_bw", "Bandwidth x", 1, 0.1, 5, step = 0.1)),
              shiny::selectInput("grp_abund", "Weight species by (optional)", choices = NULL)),
            shiny::selectInput("grp_null", "Null model",
              c("None" = "none", "Permute group labels" = "permute")),
            shiny::numericInput("grp_nrep", "Null replicates", 99, 9, 9999),
            sh$details(sh$summary("Assemblages (optional)"),
              shiny::fileInput("grp_sites_file", "Long table: one row per species in an assemblage",
                               accept = c(".csv", ".tsv", ".txt")),
              shiny::uiOutput("grp_sites_cols"),
              sh$p(class = "step-help",
                   "Columns: assemblage (site, basin...), species and, optionally, the group of the species",
                   " in that assemblage (e.g. native here, introduced there). Needs a reference group."))),
          shiny::actionButton("run_groups", "Compute", class = "btn-primary w-100"),
          sh$p(class = "step-help",
            "Groups: convex hull volume (FRic), dispersion, pairwise overlap of hulls",
            " (intersection / union, functional Jaccard turnover and nestedness) and of kernel",
            " densities, PERMANOVA / PERMDISP on the coordinates. With a reference group (e.g.",
            " natives, residents, unexploited species): novelty of each focal species (FNov =",
            " distance to the nearest reference species), position in the reference hull, loss of",
            " uniqueness of reference species (dFUn_rel), contribution of the focal group to FRic,",
            " FDis, FOri and MPD, overall and per assemblage. Continuous variable: linear fit on the",
            " axes, functional autocorrelation (k nearest neighbours) and correlations with FUn and",
            " FSp (uses n axes / k of tabs 5-6).")
        ),
        # --- continuous variable
        shiny::conditionalPanel("input.grp_src == 'continuous'",
          bslib::card(bslib::card_header("Continuous variable in each space"),
                      DT::DTOutput("cov_summary")),
          bslib::card(full_screen = TRUE,
            bslib::card_header(bslib::layout_columns(col_widths = c(6, 3, 3),
              shiny::selectInput("cov_space", NULL, choices = NULL),
              shiny::numericInput("cov_ax1", NULL, 1, 1, 20),
              shiny::numericInput("cov_ax2", NULL, 2, 1, 20))),
            .fig_output("cov_plot", height = "440px")),
          bslib::card(bslib::card_header("Species"), DT::DTOutput("cov_species"))),
        # --- groups
        shiny::conditionalPanel("input.grp_src != 'continuous'",
        bslib::card(bslib::card_header("Separation of the groups in each space"),
                    DT::DTOutput("grp_summary"),
                    shiny::uiOutput("grp_focal_summary_ui")),
        bslib::layout_columns(col_widths = c(7, 5),
          bslib::card(full_screen = TRUE,
            bslib::card_header(bslib::layout_columns(col_widths = c(4, 2, 2, 2, 2),
              shiny::selectInput("grp_space", NULL, choices = NULL),
              shiny::numericInput("grp_ax1", NULL, 1, 1, 20),
              shiny::numericInput("grp_ax2", NULL, 2, 1, 20),
              shiny::checkboxInput("grp_contours", "Densities", FALSE),
              shiny::checkboxInput("grp_centroids", "Centroids", TRUE))),
            .fig_output("grp_plot", height = "440px")),
          bslib::card(full_screen = TRUE,
            bslib::card_header(bslib::layout_columns(col_widths = c(5, 7),
              sh$span("Pairwise overlap across spaces"),
              shiny::selectInput("grp_pair_metric", NULL,
                c("Hull overlap (intersection / union)" = "overlap",
                  "SES of the overlap (null model)" = "ses_overlap",
                  "Overlap of kernel densities" = "overlap_kde",
                  "SES of the density overlap" = "ses_overlap_kde",
                  "Jaccard dissimilarity" = "beta_jac",
                  "Turnover component" = "beta_jtu",
                  "Nestedness component" = "beta_jne",
                  "Distance between centroids" = "centroid_dist",
                  "SES of the centroid distance" = "ses_centroid_dist")))),
            .fig_output("grp_overlap_plot", height = "440px"))
        ),
        bslib::navset_card_tab(
          bslib::nav_panel("Groups", DT::DTOutput("grp_groups_table")),
          bslib::nav_panel("Pairs of groups", DT::DTOutput("grp_pairs_table")),
          bslib::nav_panel("Centroids", DT::DTOutput("grp_centroids_table")),
          bslib::nav_panel("Species mixing", DT::DTOutput("grp_species_table")),
          bslib::nav_panel("Focal species", DT::DTOutput("grp_focal_table")),
          bslib::nav_panel("Reference species", DT::DTOutput("grp_ref_table")),
          bslib::nav_panel("Contribution of focal groups", DT::DTOutput("grp_contrib_table")),
          bslib::nav_panel("Assemblages",
            bslib::layout_columns(col_widths = c(4, -8),
              shiny::selectInput("grp_site_metric", NULL,
                c("Relative contribution to FRic" = "FRic_change_rel", "SES of the FRic contribution" = "ses_FRic",
                  "Contribution to FDis" = "d_FDis", "SES of the FDis contribution" = "ses_FDis",
                  "Contribution to FOri" = "d_FOri", "Contribution to MPD" = "d_MPD",
                  "Hull overlap focal / reference" = "overlap",
                  "Focal species outside the reference hull" = "p_focal_outside"))),
            .fig_output("grp_sites_plot", height = "340px"),
            DT::DTOutput("grp_sites_table"))))
      )
    ),

    # 8. EXPORT ---------------------------------------------------------------
    bslib::nav_panel(.step_title(8, "Export"), value = "8. Export",
      bslib::layout_columns(col_widths = c(4, 8), fill = FALSE, style = "align-items: start;",
        sh$div(class = "home-col",
          .card_x("Reproducible project (.zip)",
            shiny::textInput("prj_title", "Title", "Functional spaces analysis", width = "100%"),
            shiny::textInput("prj_author", "Author(s)", "", width = "100%"),
            shiny::textInput("prj_name", "Project (folder) name",
                             paste0("DisentFunSpaces_", format(Sys.Date(), "%Y%m%d")), width = "100%"),
            sh$label(class = "form-label", "Analyses included (those run in the tabs)"),
            shiny::uiOutput("prj_steps"),
            shiny::checkboxInput("prj_run", "Run the project's scripts before zipping (recommended)", TRUE),
            shiny::uiOutput("prj_run_hint"),
            shiny::textInput("prj_save_dir", "Also save the .zip in this folder (optional)", "",
                             placeholder = "/path/to/folder", width = "100%"),
            sh$p(class = "hint",
              "The project is written with the settings of the tabs, then its own scripts are RUN",
              " (run_all.R, in a separate R session): tables, figures and the Methods & Results",
              " report come from them, so the archive reproduces itself."),
            shiny::actionButton("prj_build", "Build the project", class = "btn-primary w-100"),
            shiny::uiOutput("prj_status"),
            shiny::uiOutput("prj_dl_ui")),
          .card_x("Individual files",
            sh$p(class = "hint", "Quick downloads of what is on screen; the project above holds",
                 " all of them, with the code that produced them."),
            sh$div(class = "d-grid gap-2 export-files",
              shiny::downloadButton("dl_script", "Single R script (.R)", class = "btn-sm"),
              shiny::downloadButton("dl_ti2", "Trait information (.csv)", class = "btn-sm"),
              shiny::downloadButton("dl_coords", "PCoA coordinates (.csv)", class = "btn-sm"),
              shiny::downloadButton("dl_cmp", "Procrustes tests (.csv)", class = "btn-sm"),
              shiny::downloadButton("dl_metrics", "Species metrics (.csv)", class = "btn-sm"),
              shiny::downloadButton("dl_nonnative", "Species groups (.zip of .csv)", class = "btn-sm"),
              shiny::downloadButton("dl_rds", "All R objects (.rds)", class = "btn-sm")))),
        .card_x("Contents of the project",
          bslib::navset_underline(id = "prj_tabs",
            bslib::nav_panel("Architecture", shiny::uiOutput("prj_tree")),
            bslib::nav_panel("README", shiny::uiOutput("prj_readme")),
            bslib::nav_panel("Methods & Results", shiny::uiOutput("prj_report")),
            bslib::nav_panel("parameters.R", shiny::uiOutput("prj_params")),
            bslib::nav_panel("Single script", shiny::verbatimTextOutput("script_preview")),
            bslib::nav_panel("Run log", shiny::uiOutput("prj_log"))))
      )
    ),
    bslib::nav_spacer(),
    bslib::nav_menu(title = "Help", align = "right",
      bslib::nav_item(sh$a("Package documentation (README)", target = "_blank",
                           href = "https://github.com/FunTraits/DisentFunSpacesR#readme")),
      bslib::nav_item(sh$a("The paper (Ecography 2026)", target = "_blank",
                           href = "https://doi.org/10.1002/ecog.08578")),
      bslib::nav_item(sh$a("Report a problem", target = "_blank",
                           href = "https://github.com/FunTraits/DisentFunSpacesR/issues"))),
    bslib::nav_item(sh$span(class = "version", paste0("v", utils::packageVersion("DisentFunSpaces"))))
  )
}

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------
app_server <- function(input, output, session) {
  rv <- shiny::reactiveValues(ti = NULL, traits = NULL, traits_final = NULL, spaces = NULL,
                              cmp = NULL, metrics = NULL, groups = NULL, grp_set = NULL,
                              extinction = NULL, covariate = NULL,
                              custom = list(), ti_version = 0, spaces_version = 0,
                              valid_msg = NULL, missing_msg = NULL, phylo_name = NULL)
  notify_error <- function(e) shiny::showNotification(conditionMessage(e), type = "error", duration = 10)
  sig3 <- function(t) { num <- vapply(t, is.numeric, logical(1)); t[num] <- lapply(t[num], signif, 3); t }
  dt_page <- function(t) DT::datatable(sig3(t), rownames = FALSE,
                                       options = list(pageLength = 10, scrollX = TRUE))

  # ---- 1. data --------------------------------------------------------------
  raw_file <- shiny::reactive({
    if (input$src %in% c("example", "example_all")) {
      nm <- if (input$src == "example") "birds_traits" else "birds_traits_all"
      e <- new.env(); utils::data(list = nm, package = "DisentFunSpaces", envir = e)
      return(e[[nm]])
    }
    shiny::req(input$file)
    path <- input$file$datapath
    tryCatch({
      if (grepl("\\.xlsx?$", input$file$name, ignore.case = TRUE)) {
        if (!requireNamespace("readxl", quietly = TRUE)) stop("Install 'readxl' to read Excel files.")
        as.data.frame(readxl::read_excel(path))
      } else {
        .read_table_fast(path, input$sep)
      }
    }, error = function(e) { notify_error(e); NULL })
  })
  # the table used by the analyses: species repeated in several rows (records:
  # individuals, populations...) get one identifier per record (record_ids());
  # the original names are kept in attribute "species"
  raw <- shiny::reactive({
    d <- raw_file(); shiny::req(d)
    sc <- input$species_col
    if (!is.null(sc) && sc %in% names(d)) {
      sp <- as.character(d[[sc]])
      if (anyDuplicated(sp[!is.na(sp)])) {
        d[[sc]] <- record_ids(sp)
        attr(d, "species") <- sp
      }
    }
    d
  })
  # species of each record (named by record), NULL without repeated species
  rec_map <- shiny::reactive({
    d <- raw(); sp <- attr(d, "species")
    if (is.null(sp)) NULL else stats::setNames(sp, d[[input$species_col]])
  })
  data_name <- shiny::reactive(switch(input$src, example = "birds_traits.csv",
                                     example_all = "birds_traits_all.csv", input$file$name))

  shiny::observeEvent(raw_file(), {
    d <- raw_file(); shiny::req(d)
    cols <- names(d)
    guess_sp <- cols[which(tolower(cols) %in% c("species", "sp", "taxon", "scientificname",
                                                "scientificnamestd", "binomial"))[1]]
    if (is.na(guess_sp)) guess_sp <- cols[1]
    guess_iucn <- cols[which(tolower(cols) %in% c("iucn", "category", "iucn_category", "redlist"))[1]]
    shiny::updateSelectInput(session, "species_col", choices = cols, selected = guess_sp)
    shiny::updateSelectInput(session, "iucn_col", choices = c("(none)" = "", cols),
                             selected = if (is.na(guess_iucn)) "" else guess_iucn)
    guess_st <- cols[tolower(cols) %in% .GROUP_COL_GUESS]
    shiny::updateSelectizeInput(session, "group_cols", choices = cols, selected = guess_st)
    # the full bird dataset has species pairs sharing no reproduction trait:
    # propose imputation (as in the published analysis)
    if (identical(input$src, "example_all")) shiny::updateRadioButtons(session, "missing", selected = "impute")
    rv$ti <- NULL; rv$traits <- NULL; rv$traits_final <- NULL
    rv$spaces <- NULL; rv$cmp <- NULL; rv$metrics <- NULL; rv$groups <- NULL
    rv$extinction <- NULL; rv$covariate <- NULL
  })

  read_sep <- shiny::reactive({
    d <- raw_file()
    s <- attr(d, "sep"); if (is.null(s)) s <- if (identical(input$sep, "auto")) "," else input$sep
    dc <- attr(d, "dec"); if (is.null(dc)) dc <- "."
    list(sep = s, dec = dc)
  })
  output$data_info <- shiny::renderUI({
    d <- raw_file(); shiny::req(d)
    rm <- rec_map()
    rs <- read_sep()
    sep_txt <- if (!is.null(attr(d, "sep")))
      sprintf(" - separator: %s, decimal mark: '%s'",
              c(";" = "semicolon", "\t" = "tab", "," = "comma", "|" = "pipe")[[rs$sep]], rs$dec)
    else ""
    shiny::tagList(
      shiny::tags$p(sprintf("%d %s x %d columns - %.1f%% missing cells%s",
                            nrow(d), if (is.null(rm)) "species" else "rows", ncol(d),
                            100 * mean(is.na(d)), sep_txt)),
      if (!is.null(rm)) {
        n_rep <- sum(table(rm) > 1)
        shiny::tags$div(class = "alert alert-info", style = "font-size: 0.9rem;",
          shiny::tags$b("Repeated species. "),
          sprintf(paste0("%d rows for %d species (%d species in several rows). Each row is a record",
                         " (individual, population...) placed in the spaces with its own identifier",
                         " (e.g. '%s'). Groups and IUCN categories are read per row; assemblages",
                         " list species, which are given all their records."),
                  length(rm), length(unique(rm)), n_rep, names(rm)[duplicated(rm)][1]))
      },
      if (ncol(d) == 1L) shiny::tags$div(class = "alert alert-warning",
        "Only one column was read: check the separator of the file."))
  })
  output$perf_info <- shiny::renderUI({
    d <- raw_file(); shiny::req(d)
    n <- nrow(d)
    if (n < .app_large_n()) return(NULL)
    gb <- 8 * n^2 / 1024^3
    shiny::tags$div(class = "alert alert-info", style = "font-size: 0.9rem;",
      shiny::tags$b("Large dataset. "),
      sprintf(paste0("One %d x %d distance matrix takes %.1f GB of RAM (plan about %.0f GB free).",
                     " By default (tab 4, 'Performance') the app computes distances with the fast",
                     " engine and only the positive PCoA axes, without a full eigendecomposition."),
              n, n, gb, max(2, ceiling(4 * gb))),
      if (!requireNamespace("RSpectra", quietly = TRUE)) shiny::tags$span(
        shiny::tags$b(" Install the package 'RSpectra'"),
        " (install.packages('RSpectra')): without it the full eigendecomposition is used,",
        " which takes about ten minutes per space for 8,000 species."),
      if (.reference_blas()) shiny::tags$span(
        " R uses its reference BLAS: an optimised BLAS (Accelerate on macOS, OpenBLAS on",
        " Linux) makes PCoA and Procrustes tests several times faster."))
  })
  output$data_preview <- DT::renderDT(DT::datatable(raw_file(), options = list(scrollX = TRUE, pageLength = 10),
                                                    rownames = FALSE), server = TRUE)

  # ---- 2. trait information --------------------------------------------------
  trait_cols <- shiny::reactive({
    shiny::req(raw(), input$species_col)
    setdiff(names(raw()), c(input$species_col, input$iucn_col, input$group_cols))
  })

  do_guess <- function() {
    d <- raw(); shiny::req(d, input$species_col)
    g <- guess_trait_info(d[, c(input$species_col, trait_cols()), drop = FALSE],
                          species_col = input$species_col)
    rv$ti <- g; rv$ti_version <- rv$ti_version + 1
  }
  shiny::observeEvent(input$guess, do_guess())
  shiny::observeEvent(list(raw(), input$species_col, input$iucn_col, input$group_cols), {
    sc <- input$species_col
    if (is.null(sc) || !nzchar(sc)) return()
    if (isTRUE(input$src %in% c("example", "example_all")) && identical(sc, "species")) {
      e <- new.env(); utils::data("birds_trait_info", package = "DisentFunSpaces", envir = e)
      rv$ti <- e$birds_trait_info; rv$ti_version <- rv$ti_version + 1
    } else if (!is.null(raw()) && sc %in% names(raw())) {
      do_guess()
    }
  }, ignoreInit = FALSE)

  shiny::observeEvent(input$ti_file, {
    ti <- tryCatch(utils::read.csv(input$ti_file$datapath, stringsAsFactors = FALSE),
                   error = function(e) { notify_error(e); NULL })
    shiny::req(ti)
    base <- guess_trait_info(raw()[, c(input$species_col, trait_cols()), drop = FALSE],
                             species_col = input$species_col)
    base$use <- FALSE
    m <- match(base$trait, ti$trait)
    for (col in intersect(names(ti), names(base))) {
      if (col == "trait") next
      base[[col]][!is.na(m)] <- ti[[col]][m[!is.na(m)]]
    }
    if (!"use" %in% names(ti)) base$use[!is.na(m)] <- TRUE
    rv$ti <- base; rv$ti_version <- rv$ti_version + 1
    shiny::showNotification(sprintf("%d of %d traits matched.", sum(!is.na(m)), nrow(base)))
  })

  # read the editor inputs back into a table
  current_ti <- shiny::reactive({
    ti <- rv$ti; shiny::req(ti)
    v <- rv$ti_version
    for (i in seq_len(nrow(ti))) {
      id <- function(f) paste0("ti_", f, "_", v, "_", i)
      if (is.null(input[[id("type")]])) next
      ti$use[i] <- isTRUE(input[[id("use")]])
      ti$type[i] <- input[[id("type")]]
      ti$dimension[i] <- input[[id("dim")]]
      ti$block[i] <- input[[id("block")]]
      ti$transform[i] <- input[[id("tr")]]
      ti$weight[i] <- input[[id("w")]]
      ti$levels[i] <- input[[id("lev")]]
      ti$circ_min[i] <- input[[id("cmin")]]
      ti$circ_max[i] <- input[[id("cmax")]]
    }
    ti
  })

  output$ti_editor <- shiny::renderUI({
    ti <- rv$ti; shiny::req(ti)
    v <- rv$ti_version
    tcodes <- trait_types()$code
    sh <- shiny::tags
    rows <- lapply(seq_len(nrow(ti)), function(i) {
      id <- function(f) paste0("ti_", f, "_", v, "_", i)
      sh$tr(
        sh$td(shiny::checkboxInput(id("use"), NULL, isTRUE(ti$use[i]), width = "30px")),
        sh$td(sh$code(ti$trait[i])),
        sh$td(shiny::selectInput(id("type"), NULL, tcodes, ti$type[i], selectize = FALSE, width = "70px")),
        sh$td(shiny::textInput(id("dim"), NULL, ti$dimension[i], width = "150px")),
        sh$td(shiny::textInput(id("block"), NULL, ti$block[i], width = "130px")),
        sh$td(shiny::selectInput(id("tr"), NULL, .transforms, ti$transform[i], selectize = FALSE, width = "90px")),
        sh$td(shiny::numericInput(id("w"), NULL, ti$weight[i], min = 0, step = 0.1, width = "70px")),
        sh$td(shiny::textInput(id("lev"), NULL, ifelse(is.na(ti$levels[i]), "", ti$levels[i]), width = "150px")),
        sh$td(shiny::numericInput(id("cmin"), NULL, ti$circ_min[i], width = "70px")),
        sh$td(shiny::numericInput(id("cmax"), NULL, ti$circ_max[i], width = "70px"))
      )
    })
    sh$table(class = "table table-sm ti-table",
      sh$thead(sh$tr(lapply(c("use", "trait", "type", "dimension", "block", "transform",
                              "weight", "levels (ordinal)", "circ. min", "circ. max"), sh$th))),
      sh$tbody(rows))
  })

  shiny::observe({
    ti <- rv$ti; shiny::req(ti)
    shiny::updateSelectizeInput(session, "bulk_rows", choices = ti$trait)
  })
  shiny::observeEvent(input$bulk_apply, {
    shiny::req(input$bulk_rows, nzchar(input$bulk_dim))
    ti <- current_ti()
    ti$dimension[ti$trait %in% input$bulk_rows] <- input$bulk_dim
    rv$ti <- ti; rv$ti_version <- rv$ti_version + 1
  })

  output$types_table <- DT::renderDT(DT::datatable(trait_types(), rownames = FALSE,
                                                   options = list(dom = "t", pageLength = 10)))

  shiny::observeEvent(input$validate, {
    ti <- current_ti()
    res <- tryCatch({
      w <- NULL
      tr <- withCallingHandlers(
        prepare_traits(raw(), ti, species_col = input$species_col, taxon = attr(raw(), "species")),
        warning = function(wn) { w <<- c(w, conditionMessage(wn)); invokeRestart("muffleWarning") })
      list(tr = tr, w = w)
    }, error = function(e) { notify_error(e); NULL })
    if (is.null(res)) { rv$valid_msg <- NULL; return() }
    rv$ti <- ti
    rv$traits <- res$tr; rv$traits_final <- NULL; rv$spaces <- NULL; rv$cmp <- NULL; rv$metrics <- NULL
    rv$valid_msg <- res$w
    dims <- unique(res$tr$trait_info$dimension)
    sp <- default_spaces(res$tr)
    shiny::updateCheckboxGroupInput(session, "spaces_sel", choices = names(sp), selected = names(sp))
    shiny::updateSelectizeInput(session, "custom_dims", choices = dims)
    rv$custom <- list()
    if (length(res$tr$species) >= .app_large_n())
      shiny::showNotification(sprintf(paste0("%d species: the spaces will be built with fast distances and a",
                                             " partial PCoA (tab 4, 'Performance')."),
                                      length(res$tr$species)), type = "message", duration = 10)
    shiny::showNotification(sprintf("%d species, %d traits in %d dimension(s).",
                                    length(res$tr$species), nrow(trait_summary(res$tr)), length(dims)),
                            type = "message")
    bslib::nav_select("nav", "3. Missing values")
  })
  output$validation_msg <- shiny::renderUI({
    if (is.null(rv$traits)) return(shiny::tags$p(class = "step-help", "Not validated yet."))
    shiny::tagList(
      shiny::tags$div(class = "alert alert-success", "Traits validated."),
      if (length(rv$valid_msg)) shiny::tags$div(class = "alert alert-warning",
                                                shiny::HTML(paste(rv$valid_msg, collapse = "<br>"))))
  })
  output$trait_summary <- DT::renderDT({
    shiny::req(rv$traits)
    DT::datatable(trait_summary(rv$traits), rownames = FALSE, options = list(dom = "t", pageLength = 100))
  })
  ti_download <- function() shiny::downloadHandler(
    filename = function() "trait_info.csv",
    content = function(file) utils::write.csv(current_ti(), file, row.names = FALSE))
  output$dl_ti <- ti_download()
  output$dl_ti2 <- ti_download()

  # ---- 3. missing values -----------------------------------------------------
  output$missing_table <- DT::renderDT({
    shiny::req(rv$traits)
    tr <- if (!is.null(rv$traits_final)) rv$traits_final else rv$traits
    DT::datatable(trait_summary(tr), rownames = FALSE, options = list(dom = "t", pageLength = 100))
  })
  output$missing_msg <- shiny::renderUI({
    if (is.null(rv$traits)) return(shiny::tags$p(class = "step-help", "Validate the traits first (tab 2)."))
    shiny::tags$p(rv$missing_msg)
  })
  shiny::observeEvent(input$run_missing, {
    shiny::req(rv$traits)
    tr <- rv$traits
    out <- tryCatch(switch(input$missing,
      keep = { rv$missing_msg <- "Missing values kept: distances use the traits available for each pair of species."; tr },
      drop = {
        r <- drop_incomplete(tr)
        rv$missing_msg <- sprintf("%d of %d species kept.", length(r$species), length(tr$species)); r
      },
      impute = {
        tree <- NULL
        if (!is.null(input$phylo)) tree <- ape::read.tree(input$phylo$datapath)
        rv$phylo_name <- if (!is.null(input$phylo)) input$phylo$name else NULL
        r <- shiny::withProgress(message = "missForest imputation...", value = 0.3,
          impute_traits(tr, phylo = tree, covariates = cov_data(),
                        n_phylo_axes = input$n_phylo, seed = input$seed,
                        ntree = input$ntree, verbose = FALSE))
        oob <- r$imputation$oob
        preds <- c(r$imputation$covariates,
                   if (!is.null(tree)) sprintf("%d eigenvectors of %s", input$n_phylo, input$phylo$name))
        rv$missing_msg <- sprintf("Imputation done (predictors: %s). Out-of-bag error: %s.",
                                  if (length(preds)) paste(preds, collapse = ", ") else "traits only",
                                  paste(names(oob), signif(oob, 3), collapse = ", ")); r
      }), error = function(e) { notify_error(e); NULL })
    shiny::req(out)
    rv$traits_final <- out; rv$spaces <- NULL; rv$cmp <- NULL; rv$metrics <- NULL
    bslib::nav_select("nav", "4. Functional spaces")
  })
  traits_used <- shiny::reactive(if (!is.null(rv$traits_final)) rv$traits_final else rv$traits)

  # candidate predictor columns: every column of the table that is not a used trait
  shiny::observe({
    d <- raw(); shiny::req(d, input$species_col)
    used <- if (!is.null(rv$traits)) rv$traits$trait_info$trait else
      if (!is.null(rv$ti)) rv$ti$trait[rv$ti$use] else character(0)
    cand <- setdiff(names(d), c(input$species_col, used))
    prev <- shiny::isolate(input$cov_cols)
    sel <- if (length(prev)) intersect(prev, cand) else cand[grepl("eigen", cand, ignore.case = TRUE)]
    shiny::updateSelectizeInput(session, "cov_cols", choices = cand, selected = sel)
  })
  cov_data <- function() {
    cols <- input$cov_cols
    if (!length(cols)) return(NULL)
    d <- raw(); shiny::req(d, input$species_col %in% names(d))
    cols <- intersect(setdiff(cols, rv$traits$trait_info$trait), names(d))
    if (!length(cols)) return(NULL)
    out <- d[, cols, drop = FALSE]
    rownames(out) <- as.character(d[[input$species_col]])
    out
  }

  # ---- 4. spaces --------------------------------------------------------------
  shiny::observeEvent(input$add_custom, {
    shiny::req(length(input$custom_dims) >= 1)
    nm <- if (nzchar(input$custom_name)) input$custom_name else paste(input$custom_dims, collapse = "+")
    rv$custom[[nm]] <- input$custom_dims
    sp <- c(default_spaces(rv$traits), rv$custom)
    shiny::updateCheckboxGroupInput(session, "spaces_sel", choices = names(sp),
                                    selected = c(input$spaces_sel, nm))
  })
  shiny::observeEvent(input$run_spaces, {
    tr <- traits_used(); shiny::req(tr, length(input$spaces_sel) >= 1)
    all_sp <- c(default_spaces(tr), rv$custom)[input$spaces_sel]
    res <- tryCatch(shiny::withProgress(message = "Building spaces", value = 0, {
      mx <- if (is.null(input$max_axes) || is.na(input$max_axes)) NULL else as.integer(input$max_axes)
      n <- length(tr$species)
      eng <- if (is.null(input$engine)) "auto" else input$engine
      if (eng == "auto") eng <- if (n >= 1000) "fast" else "ade4"
      em <- if (is.null(input$eigen_mode)) "auto" else input$eigen_mode
      eig <- if (em == "auto") { if (n >= 1000) "partial" else "full" } else em
      if (input$correction != "none") eig <- "full"
      rv$build_settings <- list(engine = eng, eigen = eig, keep_dist = n < 5000)
      build_spaces(tr, spaces = all_sp, weighting = input$weighting, correction = input$correction,
                   max_axes = mx, envfit_perm = 0, engine = eng, eigen = eig,
                   keep_dist = n < 5000,
                   progress = function(i, n, nm) shiny::setProgress((i - 1) / n, detail = nm))
    }), error = function(e) { notify_error(e); NULL })
    shiny::req(res)
    rv$spaces <- res; rv$cmp <- NULL; rv$metrics <- NULL; rv$groups <- NULL
    rv$extinction <- NULL; rv$covariate <- NULL
    rv$spaces_version <- rv$spaces_version + 1
    nms <- names(res)
    shiny::updateSelectInput(session, "plot_space", choices = nms, selected = nms[length(nms)])
    shiny::updateSelectInput(session, "res_a", choices = nms, selected = nms[1])
    shiny::updateSelectInput(session, "res_b", choices = nms, selected = nms[length(nms)])
    shiny::updateSelectInput(session, "m_a", choices = nms, selected = nms[1])
    shiny::updateSelectInput(session, "m_b", choices = nms, selected = nms[length(nms)])
  })
  output$spaces_msg <- shiny::renderUI({
    if (is.null(traits_used())) return(shiny::tags$p(class = "step-help", "Validate the traits first (tab 2)."))
    if (is.null(rv$spaces)) return(shiny::tags$p(class = "step-help", "Choose the spaces and click 'Build spaces'."))
    NULL
  })
  output$spaces_table <- DT::renderDT({
    shiny::req(rv$spaces)
    DT::datatable(space_summary(rv$spaces), rownames = FALSE, options = list(dom = "t"))
  })
  sel_space <- shiny::reactive({ shiny::req(rv$spaces, input$plot_space %in% names(rv$spaces)); rv$spaces[[input$plot_space]] })
  iucn_vec <- shiny::reactive({
    if (is.null(input$iucn_col) || !nzchar(input$iucn_col)) return(NULL)
    d <- raw()
    # after a new import, input$iucn_col / input$species_col may still hold the previous
    # dataset's column names until updateSelectInput() has round-tripped to the client
    shiny::req(d, input$iucn_col %in% names(d), input$species_col %in% names(d))
    stats::setNames(toupper(as.character(d[[input$iucn_col]])), as.character(d[[input$species_col]]))
  })
  min_cor_d <- shiny::debounce(shiny::reactive(input$min_cor), 400)
  space_plot <- shiny::reactive({
    s <- sel_space()
    ax <- c(input$ax1, input$ax2)
    shiny::validate(shiny::need(max(ax) <= ncol(s$coords), sprintf("This space has %d axes.", ncol(s$coords))),
                    shiny::need(ax[1] != ax[2], "Choose two different axes."))
    col <- NULL
    if (isTRUE(input$color_iucn) && !is.null(iucn_vec())) {
      col <- factor(iucn_vec(), levels = c("LC", "NT", "VU", "EN", "CR", "DD"))
      names(col) <- names(iucn_vec())
    }
    p <- plot_space(s, axes = ax, min_cor = min_cor_d(), color_by = col, arrows = is.null(col))
    if (!is.null(col)) p <- p + ggplot2::scale_colour_manual(
      values = c(LC = "#3B9AB2", NT = "#78B7C5", VU = "#EBCC2A", EN = "#E1AF00", CR = "#F21A00", DD = "grey60"),
      drop = TRUE, na.value = "grey80")
    p
  })
  output$space_plot <- shiny::bindCache(.render_fig(space_plot()),
    rv$spaces_version, input$plot_space, input$ax1, input$ax2, min_cor_d(), input$color_iucn,
    input$iucn_col)
  output$dl_space_plot <- shiny::downloadHandler(
    filename = function() paste0("space_", make.names(input$plot_space), ".pdf"),
    content = function(file) ggplot2::ggsave(file, space_plot(), width = 7, height = 6))
  output$eigen_plot <- shiny::bindCache(
    .render_fig({ shiny::req(rv$spaces); plot_eigen(rv$spaces) }), rv$spaces_version)
  output$quality_plot <- shiny::bindCache(.render_fig({
    s <- sel_space()
    shiny::validate(shiny::need(!is.null(s$dist), "Distance not stored (> 5000 species)."))
    plot_quality(stats::setNames(list(space_quality(s)), s$name))
  }), rv$spaces_version, input$plot_space)
  output$cor_table <- DT::renderDT({
    s <- sel_space()
    cc <- s$cor[, c("trait", "dimension", "type", intersect(paste0("PCoA", 1:4), names(s$cor)))]
    num <- vapply(cc, is.numeric, logical(1)); cc[num] <- lapply(cc[num], round, 2)
    DT::datatable(cc, rownames = FALSE, options = list(dom = "tp", pageLength = 15))
  })

  # ---- 5. comparison ----------------------------------------------------------
  n_axes_val <- shiny::reactive(if (input$axes_mode == "all") NULL else input$n_axes)
  shiny::observeEvent(input$run_compare, {
    shiny::req(rv$spaces)
    if (length(rv$spaces) < 2) { shiny::showNotification("Build at least two spaces.", type = "warning"); return() }
    p_max <- max(vapply(rv$spaces, function(s) ncol(s$coords), integer(1)))
    if (!is.null(n_axes_val())) p_max <- min(p_max, n_axes_val())
    n_sp <- nrow(rv$spaces[[1]]$coords)
    if (p_max > 200 && n_sp >= 2000)
      shiny::showNotification(sprintf(paste0("Spaces with up to %d axes: each permutation computes the SVD of a",
                                             " %d x %d matrix, which can take a long time. Use 'First k axes',",
                                             " or fewer axes in tab 4, for a quick exploration."),
                                      p_max, p_max, p_max), type = "warning", duration = 15)
    cores <- if (is.null(input$cores) || is.na(input$cores)) 1L else as.integer(input$cores)
    res <- tryCatch(shiny::withProgress(message = "Procrustes tests", value = 0, {
      compare_spaces(rv$spaces, n_axes = n_axes_val(), method = input$pr_method,
                     nrepet = input$nrepet, mantel = input$mantel, cores = cores,
                     progress = function(i, n, nm) shiny::setProgress((i - 1) / n, detail = nm))
    }), error = function(e) { notify_error(e); NULL })
    rv$cmp <- res
  })
  output$cmp_plot <- .render_fig({
    shiny::validate(shiny::need(rv$cmp, "Click 'Run comparisons'."))
    plot_comparison(rv$cmp)
  })
  output$dl_cmp_plot <- shiny::downloadHandler(filename = "procrustes.pdf",
    content = function(file) ggplot2::ggsave(file, plot_comparison(rv$cmp), width = 6, height = 5))
  output$cmp_table <- DT::renderDT({
    shiny::req(rv$cmp)
    t <- rv$cmp$table
    num <- vapply(t, is.numeric, logical(1)); t[num] <- lapply(t[num], signif, 3)
    DT::datatable(t, rownames = FALSE, options = list(dom = "tp", scrollX = TRUE, pageLength = 15))
  })
  residuals_r <- shiny::reactive({
    shiny::req(rv$spaces, input$res_a %in% names(rv$spaces), input$res_b %in% names(rv$spaces))
    shiny::validate(shiny::need(input$res_a != input$res_b, "Choose two different spaces."))
    procrustes_residuals(rv$spaces[[input$res_a]], rv$spaces[[input$res_b]], n_axes = n_axes_val())
  })
  output$res_plot <- shiny::bindCache(.render_fig(plot_procrustes(residuals_r())),
    rv$spaces_version, input$res_a, input$res_b, n_axes_val())
  output$res_table <- DT::renderDT({
    r <- residuals_r(); r$residual <- signif(r$residual, 3)
    DT::datatable(r, rownames = FALSE, options = list(pageLength = 12, dom = "ftp"))
  })

  # ---- 6. metrics --------------------------------------------------------------
  shiny::observe({
    ch <- c("FUn", "FSp", if (!is.null(iucn_vec())) "FUSE")
    shiny::updateSelectInput(session, "metric", choices = ch, selected = shiny::isolate(input$metric))
  })
  shiny::observeEvent(input$run_metrics, {
    shiny::req(rv$spaces)
    rv$metrics <- tryCatch(shiny::withProgress(message = "Species metrics", value = 0.3,
                             species_metrics_all(rv$spaces, n_axes = n_axes_val(), k = input$k_nn,
                                                 iucn = iucn_vec())),
                           error = function(e) { notify_error(e); NULL })
  })
  shiny::observeEvent(input$run_ext, {
    shiny::req(rv$spaces)
    iu <- iucn_vec()
    if (is.null(iu)) { shiny::showNotification("Choose the IUCN column in tab 1.", type = "warning"); return() }
    rv$extinction <- tryCatch(shiny::withProgress(message = "Extinction scenarios", value = 0.3,
      extinction_scenarios(rv$spaces, iu, steps = input$ext_steps, n_axes = n_axes_val(),
                           hull_axes = if (is.null(input$grp_hull_axes)) 4 else input$grp_hull_axes,
                           nrep = input$ext_nrep, seed = input$seed)),
      error = function(e) { notify_error(e); NULL })
  })
  output$ext_plot <- .render_fig({
    shiny::validate(shiny::need(rv$extinction, "Click 'Run'."))
    plot_extinction(rv$extinction, input$ext_metric)
  })
  output$ext_summary <- DT::renderDT({
    shiny::req(rv$extinction)
    t <- rv$extinction$summary; t <- t[t$metric == input$ext_metric, ]
    DT::datatable(sig3(t), rownames = FALSE, options = list(dom = "t", scrollX = TRUE))
  })
  output$ext_stoch <- DT::renderDT({
    shiny::req(rv$extinction$stochastic)
    t <- rv$extinction$stochastic; t <- t[t$metric == input$ext_metric, ]
    DT::datatable(sig3(t), rownames = FALSE, options = list(dom = "t", scrollX = TRUE))
  })
  metric_cmp <- shiny::reactive({
    shiny::req(rv$metrics, input$metric %in% names(rv$metrics))
    compare_metrics(rv$metrics, input$metric, top = input$top)
  })
  output$metric_pairs <- DT::renderDT({
    shiny::validate(shiny::need(rv$metrics, "Click 'Compute'."))
    t <- metric_cmp()$pairs
    t$spearman <- round(t$spearman, 3); t$jaccard_top <- round(t$jaccard_top, 3)
    DT::datatable(t, rownames = FALSE, options = list(dom = "tp", pageLength = 15))
  })
  output$metric_scatter <- .render_fig({
    shiny::req(rv$metrics, input$m_a, input$m_b)
    m <- rv$metrics; met <- input$metric
    a <- m[m$space == input$m_a, c("species", met)]; b <- m[m$space == input$m_b, c("species", met)]
    d <- merge(a, b, by = "species"); names(d) <- c("species", "x", "y")
    top <- metric_cmp()$top_sets
    d$group <- ifelse(d$species %in% top[[input$m_a]] & d$species %in% top[[input$m_b]], "Top in both",
               ifelse(d$species %in% top[[input$m_a]], paste("Top in", input$m_a),
               ifelse(d$species %in% top[[input$m_b]], paste("Top in", input$m_b), "Other")))
    d$hover <- sprintf("%s<br>%s = %.3g / %.3g", d$species, met, d$x, d$y)
    ggplot2::ggplot(d, ggplot2::aes(x = x, y = y, colour = group)) +
      .geom_point_hover(ggplot2::aes(text = hover), alpha = 0.7) +
      ggplot2::scale_colour_manual(values = stats::setNames(
        c("grey70", "#C62828", "#1565C0", "#6A1B9A"),
        c("Other", "Top in both", paste("Top in", input$m_a), paste("Top in", input$m_b)))) +
      ggplot2::labs(x = paste(met, "-", input$m_a), y = paste(met, "-", input$m_b), colour = NULL) +
      ggplot2::theme_classic(base_size = 11) + ggplot2::theme(legend.position = "bottom")
  })
  output$metric_table <- DT::renderDT({
    shiny::req(rv$metrics)
    m <- rv$metrics
    num <- vapply(m, is.numeric, logical(1)); m[num] <- lapply(m[num], signif, 3)
    DT::datatable(m, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE))
  })

  # ---- 7. species groups -----------------------------------------------------
  shiny::observe({
    gc <- unique(c(input$group_cols, if (isTRUE(nzchar(input$iucn_col))) input$iucn_col))
    shiny::updateSelectInput(session, "grp_col", choices = if (length(gc)) gc else c("(declare columns in tab 1)" = ""),
      selected = if (isTRUE(shiny::isolate(input$grp_col) %in% gc)) shiny::isolate(input$grp_col) else gc[1])
  })
  grp_species <- shiny::reactive({
    if (!is.null(rv$spaces)) rownames(rv$spaces[[1]]$coords) else {
      d <- raw(); shiny::req(d, input$species_col %in% names(d)); as.character(d[[input$species_col]])
    }
  })
  # indicator (0 / 1) columns among the declared non-trait columns
  grp_bin_avail <- shiny::reactive({
    d <- raw(); shiny::req(d)
    gc <- intersect(input$group_cols, names(d))
    gc[vapply(gc, function(cl) .is_indicator(d[[cl]]), logical(1))]
  })
  shiny::observe({
    av <- tryCatch(grp_bin_avail(), error = function(e) character(0))
    old <- shiny::isolate(input$grp_bin_cols)
    sel <- if (length(old) && all(old %in% av)) old else {
      un <- tryCatch(names(.union_columns(raw(), av)), error = function(e) character(0))
      if (length(setdiff(av, un))) setdiff(av, un) else av
    }
    shiny::updateSelectizeInput(session, "grp_bin_cols", choices = av, selected = sel)
  })
  output$grp_bin_info <- shiny::renderUI({
    av <- tryCatch(grp_bin_avail(), error = function(e) character(0))
    if (!length(av)) return(shiny::tags$div(class = "alert alert-warning", style = "font-size: 0.85rem;",
      "No 0 / 1 column among the non-trait variables: declare them in tab 1",
      " ('Other variables, not traits'), then rebuild the spaces."))
    un <- tryCatch(.union_columns(raw(), intersect(input$grp_bin_cols, av)), error = function(e) integer(0))
    if (length(un)) shiny::tags$div(class = "alert alert-warning", style = "font-size: 0.85rem;",
      sprintf(paste("'%s' is the union of the other columns selected%s: it duplicates them as a",
                    "group; remove it, or analyse it alone (used vs none)."),
              names(un)[1], if (un[1] > 0) sprintf(" (except for %d species: check them)", un[1]) else ""))
  })
  grp_raw <- shiny::reactive({
    if (identical(input$grp_src, "binary")) {
      cols <- input$grp_bin_cols
      if (!length(cols)) return(NULL)
      d <- raw(); shiny::req(d, all(cols %in% names(d)), input$species_col %in% names(d))
      return(tryCatch(groups_from_columns(d, cols, input$species_col, none = trimws(input$grp_none)),
                      error = function(e) { notify_error(e); NULL }))
    }
    if (!identical(input$grp_src, "column")) return(NULL)
    gcol <- input$grp_col
    if (is.null(gcol) || !nzchar(gcol)) return(NULL)
    d <- raw(); shiny::req(d, gcol %in% names(d), input$species_col %in% names(d))
    stats::setNames(as.character(d[[gcol]]), as.character(d[[input$species_col]]))
  })
  cov_vec <- shiny::reactive({
    gcol <- input$grp_col
    if (is.null(gcol) || !nzchar(gcol)) return(NULL)
    d <- raw(); shiny::req(d, gcol %in% names(d), input$species_col %in% names(d))
    v <- suppressWarnings(as.numeric(as.character(d[[gcol]])))
    if (all(is.na(v))) return(NULL)
    stats::setNames(v, as.character(d[[input$species_col]]))
  })
  shiny::observe({
    d <- raw(); shiny::req(d)
    num <- names(d)[vapply(d, is.numeric, logical(1))]
    num <- setdiff(num, c(input$species_col, if (!is.null(rv$ti)) rv$ti$trait[rv$ti$use %in% TRUE]))
    shiny::updateSelectInput(session, "grp_abund", choices = c("(none)" = "", num), selected = "")
  })
  # assemblages: long table (assemblage, species, optional group)
  grp_sites_raw <- shiny::reactive({
    f <- input$grp_sites_file; if (is.null(f)) return(NULL)
    tryCatch(.read_table_fast(f$datapath, "auto"), error = function(e) { notify_error(e); NULL })
  })
  output$grp_sites_cols <- shiny::renderUI({
    d <- grp_sites_raw(); if (is.null(d)) return(NULL)
    cols <- names(d); lc <- tolower(cols)
    pick <- function(cand, def) { i <- which(lc %in% cand)[1]; if (is.na(i)) def else cols[i] }
    shiny::tagList(
      shiny::tags$p(class = "step-help", sprintf("%d rows, %d columns.", nrow(d), ncol(d))),
      shiny::selectInput("grp_site_col", "Assemblage column", cols,
                         pick(c("site", "assemblage", "basin", "plot", "locality", "country"), cols[1])),
      shiny::selectInput("grp_site_sp_col", "Species column", cols,
                         pick(c("species", "sp", "taxon", "binomial"), cols[min(2, length(cols))])),
      shiny::selectInput("grp_site_grp_col", "Group in the assemblage (optional)", c("(same as the variable)" = "", cols),
                         pick(c("status", "origin", "group"), "")))
  })
  grp_comm <- shiny::reactive({
    d <- grp_sites_raw(); if (is.null(d)) return(NULL)
    shiny::req(input$grp_site_col %in% names(d), input$grp_site_sp_col %in% names(d))
    gc <- if (isTRUE(nzchar(input$grp_site_grp_col))) input$grp_site_grp_col else NULL
    tryCatch(community_from_long(d, input$grp_site_col, input$grp_site_sp_col, gc, records = rec_map()),
             error = function(e) { notify_error(e); NULL })
  })
  grp_sep <- shiny::reactive({
    if (identical(input$grp_src, "binary")) return(";")
    if (is.null(input$grp_sep) || !nzchar(trimws(input$grp_sep))) NULL else trimws(input$grp_sep)
  })
  grp_recode <- shiny::reactive(.parse_recode(input$grp_recode))
  grp_counts <- shiny::reactive({
    g <- grp_raw(); shiny::req(g)
    gr <- tryCatch(.as_groups(g, grp_species(), sep = grp_sep(), recode = grp_recode()),
                   error = function(e) NULL)
    shiny::req(gr)
    list(n = sort(colSums(gr$M), decreasing = TRUE), exclusive = gr$exclusive,
         unassigned = sum(rowSums(gr$M) == 0))
  })
  shiny::observe({
    if (identical(input$grp_src, "simulate")) {
      shiny::updateSelectizeInput(session, "grp_ref", choices = c("native", "non-native"), selected = "native")
      return()
    }
    if (identical(input$grp_src, "continuous")) return()
    cn <- tryCatch(grp_counts(), error = function(e) NULL)
    if (is.null(cn)) {
      shiny::updateSelectizeInput(session, "grp_levels", choices = character(0))
      shiny::updateSelectizeInput(session, "grp_ref", choices = character(0)); return()
    }
    lv <- names(cn$n)
    lab <- stats::setNames(lv, sprintf("%s (%d)", lv, cn$n))
    sel <- lv[cn$n > 2][seq_len(min(8, sum(cn$n > 2)))]
    shiny::updateSelectizeInput(session, "grp_levels", choices = lab, selected = sel)
    ref <- lv[tolower(lv) %in% c("native", "resident", "none", "no", "not threatened")]
    shiny::updateSelectizeInput(session, "grp_ref", choices = lab, selected = ref[seq_len(min(1, length(ref)))])
  })
  output$grp_levels_info <- shiny::renderUI({
    cn <- tryCatch(grp_counts(), error = function(e) NULL)
    if (is.null(cn)) return(NULL)
    shiny::tags$p(class = "step-help", sprintf("%d groups (%s), %d species without a group.",
      length(cn$n), if (cn$exclusive) "one label per species" else "several labels per species",
      cn$unassigned))
  })
  shiny::observe({
    shiny::req(rv$spaces)
    nm <- names(rv$spaces)
    shiny::updateSelectInput(session, "cov_space", choices = nm, selected = utils::tail(nm, 1))
    shiny::updateSelectInput(session, "grp_space", choices = nm,
      selected = if (isTRUE(shiny::isolate(input$grp_space) %in% nm)) shiny::isolate(input$grp_space) else utils::tail(nm, 1))
    shiny::updateSelectInput(session, "nn_sim_space", choices = nm, selected = utils::tail(nm, 1))
  })
  shiny::observeEvent(input$run_groups, {
    shiny::req(rv$spaces)
    sp <- rownames(rv$spaces[[1]]$coords)
    if (identical(input$grp_src, "continuous")) {
      v <- cov_vec()
      if (is.null(v)) {
        shiny::showNotification("Choose a numeric variable (declared in tab 1).", type = "warning"); return()
      }
      rv$cov_set <- list(x = v, col = input$grp_col, transform = input$cov_transform)
      rv$covariate <- tryCatch(shiny::withProgress(message = "Continuous variable", value = 0.3,
        covariate_space(rv$spaces, v, transform = input$cov_transform, n_axes = n_axes_val(),
                        k = if (is.null(input$k_nn)) 5 else input$k_nn, fit_axes = input$cov_fit_axes,
                        nperm = input$cov_nperm, seed = input$seed)),
        error = function(e) { notify_error(e); NULL })
      return()
    }
    if (identical(input$grp_src, "simulate")) {
      ss <- if (isTRUE(input$nn_sim_space %in% names(rv$spaces))) input$nn_sim_space else length(rv$spaces)
      g <- simulate_introductions(rv$spaces[[ss]], prop = input$nn_prop, bias = input$nn_bias,
                                  seed = input$seed)
      set <- list(source = "simulate", col = NULL, groups = g, sep = NULL, recode = NULL, levels = NULL,
                  reference = if (length(input$grp_ref)) input$grp_ref else "native")
    } else {
      g <- grp_raw()
      if (is.null(g)) {
        shiny::showNotification(if (identical(input$grp_src, "binary"))
          "Choose the indicator (0 / 1) columns, declared in tab 1 as variables that are not traits." else
          "Declare a grouping column in tab 1, or simulate a native / non-native status.",
                                type = "warning"); return()
      }
      lv <- input$grp_levels
      if (length(lv) < 1) { shiny::showNotification("Choose the groups analysed.", type = "warning"); return() }
      ref <- intersect(input$grp_ref, lv)
      if (length(input$grp_ref) && !length(ref))
        shiny::showNotification("The reference group is not among the groups analysed: species-level analysis skipped.",
                                type = "warning")
      if (length(ref) && !length(setdiff(lv, ref))) {
        shiny::showNotification("Add at least one group other than the reference.", type = "warning"); return()
      }
      bin <- identical(input$grp_src, "binary")
      set <- list(source = "column", col = if (bin) input$grp_bin_cols else input$grp_col, groups = g,
                  sep = grp_sep(), recode = grp_recode(), levels = lv, reference = if (length(ref)) ref else NULL,
                  none = if (bin) trimws(input$grp_none))
    }
    if (is.null(set$reference) && length(set$levels) == 1) {
      shiny::showNotification("One group only: choose a reference group or more groups.", type = "warning"); return()
    }
    set$density <- isTRUE(input$grp_density)
    set$density_axes <- if (is.null(input$grp_density_axes)) 2 else input$grp_density_axes
    set$bw <- if (is.null(input$grp_bw)) 1 else input$grp_bw
    set$abundance_col <- if (set$density && isTRUE(nzchar(input$grp_abund)) &&
                             isTRUE(input$grp_abund %in% names(d)) &&
                             isTRUE(input$species_col %in% names(d))) input$grp_abund else NULL
    ab <- if (!is.null(set$abundance_col)) {
      d <- raw(); stats::setNames(as.numeric(d[[set$abundance_col]]), as.character(d[[input$species_col]]))
    }
    cm <- if (!is.null(set$reference)) grp_comm() else NULL
    if (!is.null(grp_sites_raw()) && is.null(set$reference))
      shiny::showNotification("Assemblages need a reference group: skipped.", type = "warning")
    set$sites <- if (!is.null(cm)) list(site_col = input$grp_site_col, species_col = input$grp_site_sp_col,
                                         group_col = if (isTRUE(nzchar(input$grp_site_grp_col))) input$grp_site_grp_col,
                                         data = grp_sites_raw())
    rv$grp_set <- set
    rv$groups <- tryCatch(shiny::withProgress(message = "Species groups", value = 0.3, {
      groups_all(rv$spaces, set$groups, reference = set$reference, pool_focal = isTRUE(input$grp_pool),
                 sep = set$sep, recode = set$recode, levels = set$levels, n_axes = n_axes_val(),
                 k = if (is.null(input$k_nn)) 5 else input$k_nn,
                 hull_axes = input$grp_hull_axes, null = input$grp_null, nrep = input$grp_nrep,
                 seed = input$seed, density = set$density, density_axes = set$density_axes,
                 bw = set$bw, abundance = ab, comm = cm)
    }), error = function(e) { notify_error(e); NULL })
  })
  # ---- continuous variable outputs
  output$cov_summary <- DT::renderDT({
    shiny::validate(shiny::need(rv$covariate, "Click 'Compute'."))
    DT::datatable(sig3(rv$covariate$summary), rownames = FALSE, options = list(dom = "t", scrollX = TRUE))
  })
  output$cov_plot <- .render_fig({
    shiny::req(rv$covariate, input$cov_space %in% names(rv$spaces))
    s <- rv$spaces[[input$cov_space]]; ax <- c(input$cov_ax1, input$cov_ax2)
    shiny::validate(shiny::need(max(ax) <= ncol(s$coords), sprintf("This space has %d axes.", ncol(s$coords))),
                    shiny::need(ax[1] != ax[2], "Choose two different axes."))
    plot_covariate(s, rv$cov_set$x, axes = ax, transform = rv$cov_set$transform, name = rv$cov_set$col)
  })
  output$cov_species <- DT::renderDT({ shiny::req(rv$covariate); dt_page(rv$covariate$species) })
  # ---- assemblages
  output$grp_sites_plot <- .render_fig({
    shiny::validate(shiny::need(rv$groups$sites, "Upload assemblages (sidebar) and choose a reference group."),
                    shiny::need(input$grp_site_metric %in% names(rv$groups$sites),
                                "Not computed: choose the 'Permute group labels' null model."))
    plot_group_sites(rv$groups, input$grp_site_metric)
  })
  output$grp_sites_table <- DT::renderDT({ shiny::req(rv$groups$sites); dt_page(rv$groups$sites) })
  output$grp_summary <- DT::renderDT({
    shiny::validate(shiny::need(rv$groups, "Click 'Compute'."))
    DT::datatable(sig3(summary_groups(rv$groups)), rownames = FALSE,
                  options = list(dom = "t", scrollX = TRUE))
  })
  output$grp_focal_summary_ui <- shiny::renderUI({
    if (is.null(rv$groups$focal)) return(NULL)
    shiny::tagList(shiny::tags$h6(sprintf("Focal groups relative to the reference (%s)",
                                          rv$groups$settings$reference), style = "margin-top: 1rem;"),
                   DT::DTOutput("grp_focal_summary"))
  })
  output$grp_focal_summary <- DT::renderDT({
    shiny::req(rv$groups$focal)
    DT::datatable(sig3(summary_groups(rv$groups, "focal")), rownames = FALSE,
                  options = list(dom = "tp", pageLength = 12, scrollX = TRUE))
  })
  grp_plot_r <- shiny::reactive({
    shiny::req(rv$groups, input$grp_space %in% names(rv$spaces))
    s <- rv$spaces[[input$grp_space]]; set <- rv$grp_set
    ax <- c(input$grp_ax1, input$grp_ax2)
    shiny::validate(shiny::need(max(ax) <= ncol(s$coords), sprintf("This space has %d axes.", ncol(s$coords))),
                    shiny::need(ax[1] != ax[2], "Choose two different axes."))
    x <- rv$groups
    if (!is.null(set$reference) && !is.null(x$focal)) {
      r <- list(focal = x$focal[x$focal$space == input$grp_space, -1],
                settings = list(hull_axes = min(input$grp_hull_axes, ncol(s$coords))))
      plot_groups(s, set$groups, axes = ax, sep = set$sep, recode = set$recode, levels = set$levels,
                  reference = set$reference, pool_focal = isTRUE(input$grp_pool), x = r, label = 5,
                  centroids = isTRUE(input$grp_centroids))
    } else plot_groups(s, set$groups, axes = ax, sep = set$sep, recode = set$recode, levels = set$levels,
                       contours = isTRUE(input$grp_contours), centroids = isTRUE(input$grp_centroids))
  })
  output$grp_plot <- .render_fig(grp_plot_r())
  output$grp_overlap_plot <- .render_fig({
    shiny::req(rv$groups)
    p <- rv$groups$pairs
    shiny::validate(shiny::need(!is.null(p) && nrow(p), "Needs at least two groups."),
                    shiny::need(input$grp_pair_metric %in% names(p),
                                "Not computed: choose the 'Permute group labels' null model."))
    plot_group_overlap(rv$groups, input$grp_pair_metric)
  })
  output$grp_groups_table <- DT::renderDT({ shiny::req(rv$groups); dt_page(rv$groups$groups) })
  output$grp_pairs_table <- DT::renderDT({ shiny::req(rv$groups$pairs); dt_page(rv$groups$pairs) })
  # centroids of the groups in every space (axes retained for the metrics)
  grp_centroids_r <- shiny::reactive({
    shiny::req(rv$groups, rv$spaces); set <- rv$grp_set
    tryCatch(group_centroids(rv$spaces, set$groups, sep = set$sep, recode = set$recode,
                             levels = set$levels, n_axes = n_axes_val()),
             error = function(e) NULL)
  })
  output$grp_centroids_table <- DT::renderDT({ shiny::req(grp_centroids_r()); dt_page(grp_centroids_r()) })
  output$grp_species_table <- DT::renderDT({ shiny::req(rv$groups); dt_page(rv$groups$species) })
  output$grp_focal_table <- DT::renderDT({
    shiny::validate(shiny::need(rv$groups$focal, "Choose a reference group to analyse focal species."))
    dt_page(rv$groups$focal)
  })
  output$grp_ref_table <- DT::renderDT({
    shiny::validate(shiny::need(rv$groups$reference, "Choose a reference group."))
    dt_page(rv$groups$reference)
  })
  output$grp_contrib_table <- DT::renderDT({
    shiny::validate(shiny::need(rv$groups$contribution, "Choose a reference group."))
    dt_page(rv$groups$contribution)
  })

  # ---- 8. export -------------------------------------------------------------
  script_params <- shiny::reactive({
    shiny::req(rv$spaces)
    sp_all <- c(default_spaces(rv$traits), rv$custom)
    list(data_file = data_name(), trait_info_file = "trait_info.csv",
         species_col = input$species_col, records = !is.null(rec_map()),
         iucn_col = if (isTRUE(nzchar(input$iucn_col))) input$iucn_col else NULL,
         group_col = if (identical(input$grp_src, "binary") && length(input$grp_bin_cols)) input$grp_bin_cols else
           if (isTRUE(nzchar(input$grp_col)) && !identical(input$grp_src, "simulate")) input$grp_col else NULL,
         group_none = if (identical(input$grp_src, "binary")) trimws(input$grp_none),
         group_cols = input$group_cols,
         group_reference = rv$grp_set$reference, group_sep = rv$grp_set$sep,
         group_recode = rv$grp_set$recode, group_levels = rv$grp_set$levels,
         hull_axes = if (is.null(input$grp_hull_axes)) 4 else input$grp_hull_axes,
         sep = read_sep()$sep, dec = read_sep()$dec, missing = input$missing,
         phylo_file = rv$phylo_name, n_phylo_axes = input$n_phylo, ntree = input$ntree,
         cov_cols = if (identical(input$missing, "impute")) names(cov_data()) else NULL,
         seed = input$seed, spaces = lapply(rv$spaces, `[[`, "dimensions"),
         weighting = input$weighting, correction = input$correction,
         max_axes = if (is.null(input$max_axes) || is.na(input$max_axes)) NULL else as.integer(input$max_axes),
         engine = rv$build_settings$engine, eigen = rv$build_settings$eigen,
         keep_dist = rv$build_settings$keep_dist,
         cores = if (is.null(input$cores) || is.na(input$cores)) 1L else as.integer(input$cores),
         method = input$pr_method, nrepet = input$nrepet, n_axes = n_axes_val(),
         mantel = input$mantel, k = input$k_nn)
  })
  output$script_preview <- shiny::renderText({
    shiny::validate(shiny::need(rv$spaces, "Build the spaces to generate the script."))
    analysis_script(script_params())
  })
  output$dl_script <- shiny::downloadHandler(filename = "functional_spaces_analysis.R",
    content = function(file) writeLines(analysis_script(script_params()), file))
  output$dl_coords <- shiny::downloadHandler(filename = "pcoa_coordinates.csv",
    content = function(file) {
      shiny::req(rv$spaces)
      utils::write.csv(coords_table(rv$spaces), file, row.names = FALSE)
    })
  output$dl_cmp <- shiny::downloadHandler(filename = "procrustes_tests.csv",
    content = function(file) { shiny::req(rv$cmp); utils::write.csv(rv$cmp$table, file, row.names = FALSE) })
  output$dl_metrics <- shiny::downloadHandler(filename = "species_metrics.csv",
    content = function(file) { shiny::req(rv$metrics); utils::write.csv(rv$metrics, file, row.names = FALSE) })
  output$dl_nonnative <- shiny::downloadHandler(filename = "species_groups.zip",
    content = function(file) {
      shiny::req(!is.null(rv$groups) || !is.null(rv$covariate))
      d <- tempfile("groups_"); dir.create(d)
      x <- rv$groups
      if (!is.null(rv$covariate)) {
        utils::write.csv(rv$covariate$summary, file.path(d, "continuous_summary.csv"), row.names = FALSE)
        utils::write.csv(rv$covariate$species, file.path(d, "continuous_species.csv"), row.names = FALSE)
      }
      tabs <- if (is.null(x)) list() else list(groups_summary = summary_groups(x), groups = x$groups, pairs = x$pairs,
                   centroids = tryCatch(grp_centroids_r(), error = function(e) NULL),
                   tests = x$test, species_mixing = x$species, focal_summary =
                     if (!is.null(x$focal)) summary_groups(x, "focal"),
                   focal_species = x$focal, reference_species = x$reference,
                   contribution = x$contribution, assemblages = x$sites)
      tabs <- tabs[!vapply(tabs, is.null, logical(1))]
      f <- file.path(d, paste0(names(tabs), ".csv"))
      for (i in seq_along(tabs)) utils::write.csv(tabs[[i]], f[i], row.names = FALSE)
      f <- list.files(d, full.names = TRUE)
      if (requireNamespace("zip", quietly = TRUE)) zip::zip(file, basename(f), root = d) else
        utils::zip(file, f, flags = "-j9X")
    })
  output$dl_rds <- shiny::downloadHandler(filename = "functional_spaces_results.rds",
    content = function(file) saveRDS(list(trait_info = current_ti(), traits = traits_used(),
                                          spaces = rv$spaces, comparison = rv$cmp,
                                          metrics = rv$metrics, extinction = rv$extinction,
                                          groups = rv$groups, continuous = rv$covariate), file))
  # ---- home page and progress through the steps ----------------------------
  shiny::observeEvent(input$go_tab, bslib::nav_select("nav", input$go_tab))
  shiny::observeEvent(input$home_example, {
    shiny::updateRadioButtons(session, "src", selected = "example")
    bslib::nav_select("nav", "1. Data")
  })
  shiny::observeEvent(input$home_upload, {
    shiny::updateRadioButtons(session, "src", selected = "upload")
    bslib::nav_select("nav", "1. Data")
  })
  steps_done <- shiny::reactive({
    d <- tryCatch(raw(), error = function(e) NULL)
    c(!is.null(d), !is.null(rv$traits), !is.null(rv$traits_final), !is.null(rv$spaces),
      !is.null(rv$cmp), !is.null(rv$metrics), !is.null(rv$groups) || !is.null(rv$covariate),
      isTRUE(!is.null(rv$prj_result$zip)))
  })
  next_step <- shiny::reactive({
    done <- steps_done()
    core <- which(!done & !.APP_STEPS$optional)
    if (length(core)) core[1] else NA_integer_
  })
  shiny::observe(session$sendCustomMessage("dfs_steps", list(done = as.list(which(steps_done())))))
  output$home_state <- shiny::renderUI({
    done <- steps_done(); d <- tryCatch(raw(), error = function(e) NULL)
    n_done <- sum(done[!.APP_STEPS$optional])
    txt <- if (is.null(d)) "No data loaded yet." else
      sprintf("Data: %s (%s species). %s %s", data_name(), format(nrow(d), big.mark = ","),
              if (!is.null(rv$spaces)) sprintf("%d functional space%s built.", length(rv$spaces),
                                               if (length(rv$spaces) > 1) "s" else "")
              else sprintf("Next: %s.", tolower(.APP_STEPS$label[next_step()])),
              if (is.na(next_step())) "Project exported." else "")
    shiny::tagList(
      shiny::tags$p(class = "home-state", txt),
      shiny::tags$div(class = "home-progress", title = sprintf("%d of 5 core steps done", n_done),
                      shiny::tags$div(class = "done", style = sprintf("width:%.0f%%;", 100 * n_done / 5))))
  })
  output$home_next <- shiny::renderUI({
    k <- next_step()
    if (is.na(k) || k <= 1) return(NULL)
    shiny::tags$a(class = "btn btn-link", href = "#", `data-go` = .APP_STEPS$tab[k],
                  paste0("Continue: ", .APP_STEPS$label[k], " \u2192"))
  })
  output$home_steps <- shiny::renderUI({
    done <- steps_done(); k <- next_step(); if (is.na(k)) k <- 0L
    shiny::tags$ol(class = "home-steps", lapply(seq_len(nrow(.APP_STEPS)), function(i) {
      cls <- if (done[i]) "done" else if (i == k) "next" else ""
      st <- if (done[i]) "done" else if (i == k) "next" else if (.APP_STEPS$optional[i]) "optional" else ""
      shiny::tags$li(class = cls, `data-go` = .APP_STEPS$tab[i], tabindex = 0, role = "link",
        shiny::tags$span(class = "n", i),
        shiny::tags$div(shiny::tags$div(class = "t", .APP_STEPS$label[i]),
                        shiny::tags$p(.APP_STEPS$what[i])),
        shiny::tags$span(class = "s", st))
    }))
  })
  output$home_stats <- shiny::renderUI({
    d <- tryCatch(raw(), error = function(e) NULL)
    tr <- traits_used()
    stat <- function(v, l) shiny::tags$div(class = "home-stat",
      shiny::tags$div(class = "v", v), shiny::tags$div(class = "l", l))
    shiny::tags$div(class = "home-stats",
      stat(if (is.null(d)) "\u2013" else format(if (is.null(tr)) nrow(d) else length(tr$species),
                                                big.mark = ","), "species"),
      stat(if (is.null(tr)) "\u2013" else nrow(trait_summary(tr)), "traits used"),
      stat(if (is.null(tr)) "\u2013" else length(unique(tr$trait_info$dimension)), "dimensions"),
      stat(if (is.null(rv$spaces)) "\u2013" else length(rv$spaces), "spaces built"))
  })

  # ---- 8. reproducible project ----------------------------------------------
  prj_steps <- shiny::reactive(c(if (!is.null(rv$cmp)) "compare", if (!is.null(rv$metrics)) "metrics",
                                 if (!is.null(rv$groups) || !is.null(rv$covariate)) "groups"))
  prj_params <- shiny::reactive({
    p <- script_params()
    p$top <- if (is.null(input$top)) 0.1 else input$top
    p$groups <- if (!is.null(rv$groups)) {
      set <- rv$grp_set
      list(source = set$source, col = set$col, none = set$none, sep = set$sep, recode = set$recode, levels = set$levels,
           reference = set$reference, pool_focal = isTRUE(input$grp_pool),
           prop = input$nn_prop, bias = input$nn_bias, sim_space = input$nn_sim_space,
           hull_axes = input$grp_hull_axes, null = input$grp_null, nrep = input$grp_nrep,
           density = set$density, density_axes = set$density_axes, bw = set$bw,
           abundance_col = set$abundance_col,
           sites = if (!is.null(set$sites)) set$sites[c("site_col", "species_col", "group_col")])
    }
    p$covariate <- if (!is.null(rv$covariate)) list(col = rv$cov_set$col, transform = rv$cov_set$transform,
      fit_axes = input$cov_fit_axes, nperm = input$cov_nperm)
    p$extinction <- if (!is.null(rv$extinction)) list(steps = input$ext_steps, nrep = input$ext_nrep,
      hull_axes = if (is.null(input$grp_hull_axes)) 4 else input$grp_hull_axes)
    p
  })
  prj_name <- shiny::reactive({
    nm <- gsub("[^A-Za-z0-9_.-]+", "_", trimws(if (is.null(input$prj_name)) "" else input$prj_name))
    if (!nzchar(nm)) "DisentFunSpaces_analysis" else nm
  })
  prj_phylo <- shiny::reactive(if (identical(input$missing, "impute") && !is.null(input$phylo))
    input$phylo$datapath else NULL)
  # a preview of the project (written without running it) for the tabs
  prev_preview <- NULL
  prj_preview <- shiny::reactive({
    shiny::req(rv$spaces)
    p <- prj_params(); st <- prj_steps(); nm <- prj_name()
    ti <- current_ti(); d <- raw_file()
    if (!is.null(prev_preview)) unlink(prev_preview, recursive = TRUE)
    prev_preview <<- tempfile("prj_preview_")
    dir <- file.path(prev_preview, nm)
    tryCatch(write_analysis_project(d, ti, dir, params = p, phylo = prj_phylo(), steps = st,
                                    assemblages = if (!is.null(rv$groups)) rv$grp_set$sites$data,
                                    title = input$prj_title, author = input$prj_author,
                                    data_source = data_name(), run = FALSE, zip = FALSE, quiet = TRUE),
             error = function(e) { notify_error(e); NULL })
  }) |> shiny::debounce(800)
  prj_shown <- shiny::reactive({
    r <- rv$prj_result
    if (!is.null(r) && isTRUE(file.exists(r$dir))) r else prj_preview()
  })

  output$prj_steps <- shiny::renderUI({
    st <- prj_steps()
    item <- function(on, label) shiny::tags$li(class = if (on) "on" else "off",
                                               shiny::tags$span(class = "dot"), label)
    shiny::tagList(shiny::tags$ul(class = "export-steps",
      item(!is.null(rv$traits), "Traits and missing values (tabs 2-3)"),
      item(!is.null(rv$spaces), sprintf("Functional spaces (tab 4)%s",
                                        if (!is.null(rv$spaces)) sprintf(": %s", paste(names(rv$spaces), collapse = ", ")) else "")),
      item("compare" %in% st, "Procrustes tests (tab 5)"),
      item("metrics" %in% st, "Species metrics (tab 6)"),
      item("groups" %in% st, "Species groups (tab 7)")),
      if (is.null(rv$spaces)) shiny::tags$p(class = "hint", "Build the spaces (tab 4) first."))
  })
  output$prj_run_hint <- shiny::renderUI({
    if (!isTRUE(input$prj_run)) return(shiny::tags$p(class = "hint",
      "Not run: the archive holds the data and scripts only; run run_all.R to get the results."))
    n <- if (!is.null(rv$spaces)) nrow(rv$spaces[[1]]$coords) else 0
    if (n >= .app_large_n()) shiny::tags$p(class = "hint", sprintf(paste0(
      "%s species: the run repeats every step of the analysis and can take as long as in the",
      " tabs (imputation, PCoA, permutations). The app stays usable meanwhile."), format(n, big.mark = ",")))
  })

  shiny::observeEvent(input$prj_build, {
    if (is.null(rv$spaces)) {
      shiny::showNotification("Build the functional spaces (tab 4) first.", type = "warning"); return()
    }
    if (!is.null(rv$prj_job) && rv$prj_job$is_alive()) return()
    base <- tempfile("prj_"); dir.create(base)
    dir <- file.path(base, prj_name())
    args <- list(data = raw_file(), trait_info = current_ti(), dir = dir, params = prj_params(),
                 phylo = if (!is.null(prj_phylo())) ape::read.tree(prj_phylo()), steps = prj_steps(),
                 assemblages = if (!is.null(rv$groups)) rv$grp_set$sites$data,
                 title = input$prj_title, author = input$prj_author, data_source = data_name(),
                 run = isTRUE(input$prj_run), zip = TRUE, quiet = TRUE)
    rv$prj_result <- NULL; rv$prj_dir <- dir
    if (requireNamespace("callr", quietly = TRUE)) {
      f <- file.path(base, "args.rds"); saveRDS(args, f)
      rv$prj_job <- callr::r_bg(function(f, libs) {
        .libPaths(libs)
        do.call(DisentFunSpaces::write_analysis_project, readRDS(f))
      }, args = list(f = f, libs = .libPaths()), supervise = TRUE)
    } else {
      rv$prj_result <- tryCatch(shiny::withProgress(message = "Building the project...", value = 0.5,
                                                    do.call(write_analysis_project, args)),
                                error = function(e) list(ok = FALSE, error = conditionMessage(e)))
      save_copy(rv$prj_result)
    }
  })
  save_copy <- function(res) {
    sd <- path.expand(trimws(if (is.null(input$prj_save_dir)) "" else input$prj_save_dir))
    if (nzchar(sd) && !is.null(res$zip) && dir.exists(sd)) {
      ok <- file.copy(res$zip, file.path(sd, basename(res$zip)), overwrite = TRUE)
      if (ok) shiny::showNotification(paste("Saved:", file.path(sd, basename(res$zip))), duration = 10)
    } else if (nzchar(sd) && !dir.exists(sd)) {
      shiny::showNotification(paste("Folder not found, .zip not copied:", sd), type = "warning")
    }
  }
  shiny::observe({
    j <- rv$prj_job; shiny::req(j)
    if (j$is_alive()) { shiny::invalidateLater(1000); return() }
    res <- tryCatch(j$get_result(), error = function(e) list(ok = FALSE, error = conditionMessage(e)))
    shiny::isolate({ rv$prj_job <- NULL; rv$prj_result <- res; save_copy(res) })
  })
  output$prj_status <- shiny::renderUI({
    if (!is.null(rv$prj_job)) {
      shiny::invalidateLater(1500)
      log <- file.path(if (is.null(rv$prj_dir)) "" else rv$prj_dir, "results", "run_all.log")
      last <- if (file.exists(log)) utils::tail(grep("Loading|Done|Report", readLines(log, warn = FALSE),
                                                     value = TRUE), 1) else "Writing the project..."
      return(shiny::tags$div(class = "alert alert-info mt-3", shiny::tags$b("Building\u2026 "),
                             if (length(last)) last else "Running run_all.R..."))
    }
    r <- rv$prj_result; if (is.null(r)) return(NULL)
    if (isTRUE(r$ok)) shiny::tags$div(class = "alert alert-success mt-3",
      "Project built and run: data, scripts, tables, figures and report are in the archive.")
    else if (is.na(if (is.null(r$ok)) FALSE else r$ok)) shiny::tags$div(class = "alert alert-secondary mt-3",
      "Project written (not run): run run_all.R to produce the results.")
    else shiny::tags$div(class = "alert alert-danger mt-3", shiny::tags$b("The project did not run to the end."),
      shiny::tags$pre(style = "white-space:pre-wrap;font-size:11px;margin-top:6px;", if (is.null(r$error)) "" else r$error),
      "The archive still holds the data and scripts: see the Run log tab, or run run_all.R.")
  })
  output$prj_dl_ui <- shiny::renderUI({
    r <- rv$prj_result; shiny::req(r, r$zip, file.exists(r$zip))
    shiny::tagList(
      shiny::downloadButton("prj_dl", paste("Download", basename(r$zip)), class = "btn-success w-100 mt-2"),
      shiny::tags$p(class = "hint mt-1", sprintf("%.1f MB", file.size(r$zip) / 1024^2)))
  })
  output$prj_dl <- shiny::downloadHandler(filename = function() basename(rv$prj_result$zip),
                                          content = function(f) file.copy(rv$prj_result$zip, f))

  hint <- function(...) shiny::tags$p(class = "hint", ...)
  output$prj_tree <- shiny::renderUI({
    if (is.null(rv$spaces)) return(hint("Build the spaces (tab 4) to see the project."))
    r <- prj_shown(); shiny::req(r)
    shiny::tagList(
      hint(if (identical(r, rv$prj_result)) "The project as built (with its results)." else
        "Preview: the project as it will be written; results/ and data/processed/ are filled when it runs."),
      shiny::tags$div(class = "tree", .project_tree(r$dir)))
  })
  md_view <- function(file) {
    if (!file.exists(file)) return(hint("Not written yet."))
    html <- tryCatch(markdown_html(file), error = function(e) NULL)
    if (is.null(html)) return(shiny::tags$pre(paste(readLines(file, warn = FALSE, encoding = "UTF-8"), collapse = "\n")))
    shiny::tags$div(class = "report", shiny::HTML(html))
  }
  markdown_html <- function(file) {
    if (!requireNamespace("rmarkdown", quietly = TRUE) || !rmarkdown::pandoc_available()) return(NULL)
    out <- tempfile(fileext = ".html")
    rmarkdown::pandoc_convert(file, from = "gfm", to = "html", output = out)
    paste(readLines(out, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  }
  output$prj_readme <- shiny::renderUI({
    if (is.null(rv$spaces)) return(hint("Build the spaces (tab 4) to see the README."))
    r <- prj_shown(); shiny::req(r); md_view(file.path(r$dir, "README.md"))
  })
  output$prj_params <- shiny::renderUI({
    if (is.null(rv$spaces)) return(hint("Build the spaces (tab 4) to see the parameters."))
    r <- prj_shown(); shiny::req(r)
    shiny::tagList(hint("utils/parameters.R: every setting of the analysis, in one place."),
      shiny::tags$pre(style = "max-height:560px;overflow:auto;font-size:12px;",
                      paste(readLines(file.path(r$dir, "utils", "parameters.R"), warn = FALSE,
                                      encoding = "UTF-8"), collapse = "\n")))
  })
  output$prj_report <- shiny::renderUI({
    r <- rv$prj_result
    if (is.null(r)) return(hint("Build the project to read the Methods & Results it writes."))
    docx <- file.path(r$dir, "report", "Methods_Results.docx")
    if (!file.exists(docx)) return(hint(if (!isTRUE(r$ok))
      "No report: the project was not run, or stopped before the report (see the Run log tab)." else
      "No rendered report (pandoc not found): knit report/Methods_Results.Rmd in RStudio."))
    html <- tempfile(fileext = ".html")
    ok <- tryCatch({ rmarkdown::pandoc_convert(docx, to = "html", output = html,
                                               options = c("--embed-resources", "--standalone")); TRUE },
                   error = function(e) tryCatch({ rmarkdown::pandoc_convert(docx, to = "html", output = html,
                                                  options = "--self-contained"); TRUE }, error = function(e) FALSE))
    if (!ok) return(hint("The report is in report/Methods_Results.docx."))
    body <- paste(readLines(html, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
    body <- sub("(?s).*<body[^>]*>", "", body, perl = TRUE); body <- sub("(?s)</body>.*", "", body, perl = TRUE)
    shiny::tags$div(class = "report", shiny::HTML(body))
  })
  output$prj_log <- shiny::renderUI({
    r <- rv$prj_result
    if (!is.null(rv$prj_job)) shiny::invalidateLater(2000)
    log <- if (!is.null(r$log)) r$log else file.path(if (is.null(rv$prj_dir)) "" else rv$prj_dir,
                                                       "results", "run_all.log")
    if (!file.exists(log)) return(hint("The log of run_all.R appears here once the project runs."))
    shiny::tags$pre(class = "runlog", paste(readLines(log, warn = FALSE), collapse = "\n"))
  })
}
