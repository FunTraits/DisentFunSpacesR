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
# fast reader for large text files ('data.table' if installed)
.read_table_fast <- function(path, sep) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    d <- data.table::fread(path, sep = sep, header = TRUE, na.strings = c("NA", ""),
                           check.names = FALSE, data.table = FALSE, showProgress = FALSE,
                           integer64 = "double", quote = "\"")
    return(d)
  }
  utils::read.table(path, header = TRUE, sep = sep, quote = "\"",
                    stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("NA", ""),
                    comment.char = "")
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
# the steps of the workflow: tab value, label, what it does
.APP_STEPS <- data.frame(
  tab = c("1. Data", "2. Traits", "3. Missing values", "4. Functional spaces",
          "5. Compare spaces", "6. Species metrics", "7. Non-native species", "8. Export"),
  label = c("Data", "Traits", "Missing values", "Functional spaces", "Compare spaces",
            "Species metrics", "Non-native species", "Export"),
  what = c("A trait table: one row per species, one column per trait (or per modality).",
           "The type and the functional dimension of each trait.",
           "Keep them (pairwise Gower distance), drop incomplete species, or impute them.",
           "One space per dimension and per combination: generalised Gower distance + PCoA.",
           "Procrustes tests: do the dimensions place species alike?",
           "Functional uniqueness, specialisation and FUSE in every space.",
           "Novelty of non-native species and their effect on the native fauna (optional).",
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
    bg = "#ffffff",
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
                               c("Comma" = ",", "Semicolon" = ";", "Tab" = "\t"))),
          shiny::selectInput("species_col", "Species column", choices = NULL),
          shiny::selectInput("iucn_col", "IUCN category column (optional, for FUSE)", choices = NULL),
          shiny::selectInput("status_col", "Native / non-native status column (optional)", choices = NULL),
          sh$p(class = "step-help",
               "One row per species, one column per trait (or per modality for fuzzy-coded",
               " and multichoice traits). Intraspecific records must be aggregated beforehand.")
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
            shiny::plotOutput("space_plot", height = "520px"),
            shiny::downloadButton("dl_space_plot", "PDF", class = "btn-sm")),
          bslib::navset_card_tab(
            bslib::nav_panel("Inertia", shiny::plotOutput("eigen_plot", height = "300px")),
            bslib::nav_panel("Quality", shiny::plotOutput("quality_plot", height = "300px"),
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
                      shiny::plotOutput("cmp_plot", height = "420px"),
                      shiny::downloadButton("dl_cmp_plot", "PDF", class = "btn-sm")),
          bslib::card(bslib::card_header("Tests"), DT::DTOutput("cmp_table"))
        ),
        bslib::card(full_screen = TRUE,
          bslib::card_header("Species-level Procrustes residuals"),
          bslib::layout_columns(col_widths = c(3, 3, 6),
            shiny::selectInput("res_a", "Space A", choices = NULL),
            shiny::selectInput("res_b", "Space B", choices = NULL), NULL),
          bslib::layout_columns(col_widths = c(7, 5),
            shiny::plotOutput("res_plot", height = "460px"),
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
            shiny::plotOutput("metric_scatter", height = "380px"))
        ),
        bslib::card(bslib::card_header("Species metrics (all spaces)"), DT::DTOutput("metric_table"))
      )
    ),

    # 7. NON-NATIVE SPECIES ----------------------------------------------------
    bslib::nav_panel(.step_title(7, "Non-natives"), value = "7. Non-native species",
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(width = 320,
          shiny::radioButtons("nn_src", "Species status",
            c("Status column of the data" = "column",
              "Simulate introductions (exploration)" = "simulate")),
          shiny::conditionalPanel("input.nn_src == 'column'",
            shiny::textInput("nn_native", "Label of native species", "native"),
            sh$p(class = "step-help", "Any other non-missing value is treated as non-native.")),
          shiny::conditionalPanel("input.nn_src == 'simulate'",
            shiny::sliderInput("nn_prop", "Proportion of non-native species", 0.02, 0.5, 0.1, 0.01),
            shiny::selectInput("nn_bias", "Position of introduced species",
                               c("random", "peripheral", "central")),
            shiny::selectInput("nn_sim_space", "Space used for the bias", choices = NULL)),
          shiny::numericInput("nn_hull_axes", "Axes for convex hulls", 4, 2, 8),
          shiny::selectInput("nn_null", "Null model (assemblage change)",
            c("None" = "none", "Permute status" = "permute",
              "Random colonists (needs a wider pool)" = "pool")),
          shiny::numericInput("nn_nrep", "Null replicates", 99, 9, 9999),
          shiny::actionButton("run_nonnative", "Compute", class = "btn-primary w-100"),
          sh$p(class = "step-help",
            "Non-native species: FNov = distance to the nearest native species; FSp_native =",
            " distance to the native centroid; position inside/outside the native convex hull.",
            " Natives: relative change in FUn once non-native species are added (dFUn_rel).",
            " Whole fauna: change in FRic, FDis, FOri and MPD (uses n axes / k of tabs 5-6).")
        ),
        bslib::card(bslib::card_header("Effect of non-native species in each space"),
                    DT::DTOutput("nn_summary")),
        bslib::layout_columns(col_widths = c(7, 5),
          bslib::card(full_screen = TRUE,
            bslib::card_header(bslib::layout_columns(col_widths = c(6, 3, 3),
              shiny::selectInput("nn_space", NULL, choices = NULL),
              shiny::numericInput("nn_ax1", NULL, 1, 1, 20),
              shiny::numericInput("nn_ax2", NULL, 2, 1, 20))),
            shiny::plotOutput("nn_plot", height = "420px")),
          bslib::card(bslib::card_header("Assemblage-level change"),
            DT::DTOutput("nn_assemblage"))
        ),
        bslib::navset_card_tab(
          bslib::nav_panel("Non-native species", DT::DTOutput("nn_species_table")),
          bslib::nav_panel("Native species", DT::DTOutput("nn_native_table")))
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
              shiny::downloadButton("dl_nonnative", "Non-native analysis (.csv)", class = "btn-sm"),
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
                              cmp = NULL, metrics = NULL, nonnative = NULL,
                              nn_status = NULL, custom = list(), ti_version = 0, spaces_version = 0,
                              valid_msg = NULL, missing_msg = NULL, phylo_name = NULL)
  notify_error <- function(e) shiny::showNotification(conditionMessage(e), type = "error", duration = 10)

  # ---- 1. data --------------------------------------------------------------
  raw <- shiny::reactive({
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
  data_name <- shiny::reactive(switch(input$src, example = "birds_traits.csv",
                                     example_all = "birds_traits_all.csv", input$file$name))

  shiny::observeEvent(raw(), {
    d <- raw(); shiny::req(d)
    cols <- names(d)
    guess_sp <- cols[which(tolower(cols) %in% c("species", "sp", "taxon", "scientificname",
                                                "scientificnamestd", "binomial"))[1]]
    if (is.na(guess_sp)) guess_sp <- cols[1]
    guess_iucn <- cols[which(tolower(cols) %in% c("iucn", "category", "iucn_category", "redlist"))[1]]
    shiny::updateSelectInput(session, "species_col", choices = cols, selected = guess_sp)
    shiny::updateSelectInput(session, "iucn_col", choices = c("(none)" = "", cols),
                             selected = if (is.na(guess_iucn)) "" else guess_iucn)
    guess_st <- cols[which(tolower(cols) %in% c("status", "origin", "native_status",
                                                "nativeness", "alien", "introduced"))[1]]
    shiny::updateSelectInput(session, "status_col", choices = c("(none)" = "", cols),
                             selected = if (is.na(guess_st)) "" else guess_st)
    # the full bird dataset has species pairs sharing no reproduction trait:
    # propose imputation (as in the published analysis)
    if (identical(input$src, "example_all")) shiny::updateRadioButtons(session, "missing", selected = "impute")
    rv$ti <- NULL; rv$traits <- NULL; rv$traits_final <- NULL
    rv$spaces <- NULL; rv$cmp <- NULL; rv$metrics <- NULL; rv$nonnative <- NULL
  })

  output$data_info <- shiny::renderUI({
    d <- raw(); shiny::req(d)
    shiny::tags$p(sprintf("%d species x %d columns - %.1f%% missing cells",
                          nrow(d), ncol(d), 100 * mean(is.na(d))))
  })
  output$perf_info <- shiny::renderUI({
    d <- raw(); shiny::req(d)
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
  output$data_preview <- DT::renderDT(DT::datatable(raw(), options = list(scrollX = TRUE, pageLength = 10),
                                                    rownames = FALSE), server = TRUE)

  # ---- 2. trait information --------------------------------------------------
  trait_cols <- shiny::reactive({
    shiny::req(raw(), input$species_col)
    setdiff(names(raw()), c(input$species_col, input$iucn_col, input$status_col))
  })

  do_guess <- function() {
    d <- raw(); shiny::req(d, input$species_col)
    g <- guess_trait_info(d[, c(input$species_col, trait_cols()), drop = FALSE],
                          species_col = input$species_col)
    rv$ti <- g; rv$ti_version <- rv$ti_version + 1
  }
  shiny::observeEvent(input$guess, do_guess())
  shiny::observeEvent(list(raw(), input$species_col, input$iucn_col, input$status_col), {
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
        prepare_traits(raw(), ti, species_col = input$species_col),
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
    d <- raw()
    cols <- setdiff(cols, rv$traits$trait_info$trait)
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
    rv$spaces <- res; rv$cmp <- NULL; rv$metrics <- NULL; rv$nonnative <- NULL
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
    d <- raw(); stats::setNames(toupper(as.character(d[[input$iucn_col]])), as.character(d[[input$species_col]]))
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
  output$space_plot <- shiny::bindCache(shiny::renderPlot(space_plot(), res = 96),
    rv$spaces_version, input$plot_space, input$ax1, input$ax2, min_cor_d(), input$color_iucn,
    input$iucn_col)
  output$dl_space_plot <- shiny::downloadHandler(
    filename = function() paste0("space_", make.names(input$plot_space), ".pdf"),
    content = function(file) ggplot2::ggsave(file, space_plot(), width = 7, height = 6))
  output$eigen_plot <- shiny::bindCache(
    shiny::renderPlot({ shiny::req(rv$spaces); plot_eigen(rv$spaces) }, res = 96), rv$spaces_version)
  output$quality_plot <- shiny::bindCache(shiny::renderPlot({
    s <- sel_space()
    shiny::validate(shiny::need(!is.null(s$dist), "Distance not stored (> 5000 species)."))
    plot_quality(stats::setNames(list(space_quality(s)), s$name))
  }, res = 96), rv$spaces_version, input$plot_space)
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
  output$cmp_plot <- shiny::renderPlot({
    shiny::validate(shiny::need(rv$cmp, "Click 'Run comparisons'."))
    plot_comparison(rv$cmp)
  }, res = 96)
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
  output$res_plot <- shiny::bindCache(shiny::renderPlot(plot_procrustes(residuals_r()), res = 96),
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
  output$metric_scatter <- shiny::renderPlot({
    shiny::req(rv$metrics, input$m_a, input$m_b)
    m <- rv$metrics; met <- input$metric
    a <- m[m$space == input$m_a, c("species", met)]; b <- m[m$space == input$m_b, c("species", met)]
    d <- merge(a, b, by = "species"); names(d) <- c("species", "x", "y")
    top <- metric_cmp()$top_sets
    d$group <- ifelse(d$species %in% top[[input$m_a]] & d$species %in% top[[input$m_b]], "Top in both",
               ifelse(d$species %in% top[[input$m_a]], paste("Top in", input$m_a),
               ifelse(d$species %in% top[[input$m_b]], paste("Top in", input$m_b), "Other")))
    ggplot2::ggplot(d, ggplot2::aes(x = x, y = y, colour = group)) +
      ggplot2::geom_point(alpha = 0.7) +
      ggplot2::scale_colour_manual(values = stats::setNames(
        c("grey70", "#C62828", "#1565C0", "#6A1B9A"),
        c("Other", "Top in both", paste("Top in", input$m_a), paste("Top in", input$m_b)))) +
      ggplot2::labs(x = paste(met, "-", input$m_a), y = paste(met, "-", input$m_b), colour = NULL) +
      ggplot2::theme_classic(base_size = 11) + ggplot2::theme(legend.position = "bottom")
  }, res = 96)
  output$metric_table <- DT::renderDT({
    shiny::req(rv$metrics)
    m <- rv$metrics
    num <- vapply(m, is.numeric, logical(1)); m[num] <- lapply(m[num], signif, 3)
    DT::datatable(m, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE))
  })

  # ---- 7. non-native species ---------------------------------------------------
  status_vec <- shiny::reactive({
    if (is.null(input$status_col) || !nzchar(input$status_col)) return(NULL)
    d <- raw(); shiny::req(input$status_col %in% names(d))
    stats::setNames(as.character(d[[input$status_col]]), as.character(d[[input$species_col]]))
  })
  shiny::observe({
    shiny::req(rv$spaces)
    nm <- names(rv$spaces)
    shiny::updateSelectInput(session, "nn_space", choices = nm,
      selected = if (isTRUE(shiny::isolate(input$nn_space) %in% nm)) shiny::isolate(input$nn_space) else utils::tail(nm, 1))
    shiny::updateSelectInput(session, "nn_sim_space", choices = nm, selected = utils::tail(nm, 1))
  })
  shiny::observeEvent(input$run_nonnative, {
    shiny::req(rv$spaces)
    sp <- rownames(rv$spaces[[1]]$coords)
    st <- if (identical(input$nn_src, "simulate")) {
      ss <- if (isTRUE(input$nn_sim_space %in% names(rv$spaces))) input$nn_sim_space else length(rv$spaces)
      simulate_introductions(rv$spaces[[ss]], prop = input$nn_prop, bias = input$nn_bias,
                             seed = input$seed)
    } else {
      v <- status_vec()
      if (is.null(v)) {
        shiny::showNotification("Choose a status column in tab 1, or simulate introductions.",
                                type = "warning"); return()
      }
      stats::setNames(v[sp], sp)
    }
    nat <- if (identical(input$nn_src, "simulate")) "native" else input$nn_native
    rv$nn_status <- st; rv$nn_native_label <- nat
    rv$nonnative <- tryCatch(shiny::withProgress(message = "Non-native species", value = 0.3, {
      nonnative_all(rv$spaces, st, native = nat, n_axes = n_axes_val(),
                    k = if (is.null(input$k_nn)) 5 else input$k_nn,
                    hull_axes = input$nn_hull_axes, null = input$nn_null,
                    nrep = input$nn_nrep, seed = input$seed)
    }), error = function(e) { notify_error(e); NULL })
  })
  sig3 <- function(t) { num <- vapply(t, is.numeric, logical(1)); t[num] <- lapply(t[num], signif, 3); t }
  output$nn_summary <- DT::renderDT({
    shiny::validate(shiny::need(rv$nonnative, "Click 'Compute'."))
    DT::datatable(sig3(summary_nonnative(rv$nonnative)), rownames = FALSE,
                  options = list(dom = "t", scrollX = TRUE))
  })
  nn_plot_r <- shiny::reactive({
    shiny::req(rv$nonnative, input$nn_space %in% names(rv$spaces))
    s <- rv$spaces[[input$nn_space]]
    ax <- c(input$nn_ax1, input$nn_ax2)
    shiny::validate(shiny::need(max(ax) <= ncol(s$coords), sprintf("This space has %d axes.", ncol(s$coords))),
                    shiny::need(ax[1] != ax[2], "Choose two different axes."))
    x <- rv$nonnative
    r <- list(nonnative = x$nonnative[x$nonnative$space == input$nn_space, -1],
              settings = list(hull_axes = min(input$nn_hull_axes, ncol(s$coords))))
    plot_nonnative(s, rv$nn_status, axes = ax, x = r, native = rv$nn_native_label, label = 5)
  })
  output$nn_plot <- shiny::renderPlot(nn_plot_r(), res = 96)
  output$nn_assemblage <- DT::renderDT({
    shiny::req(rv$nonnative$assemblage)
    a <- rv$nonnative$assemblage
    keep <- intersect(c("space", "n_native", "n_nonnative", "FRic_hist", "FRic_cur", "FRic_change_rel",
                        "d_FDis", "d_FOri", "d_MPD", "p_nn_outside", "ses_FRic", "p_FRic",
                        "ses_FDis", "p_FDis"), names(a))
    DT::datatable(sig3(a[, keep]), rownames = FALSE, options = list(dom = "t", scrollX = TRUE))
  })
  output$nn_species_table <- DT::renderDT({
    shiny::req(rv$nonnative)
    DT::datatable(sig3(rv$nonnative$nonnative), rownames = FALSE,
                  options = list(pageLength = 10, scrollX = TRUE))
  })
  output$nn_native_table <- DT::renderDT({
    shiny::req(rv$nonnative)
    DT::datatable(sig3(rv$nonnative$native), rownames = FALSE,
                  options = list(pageLength = 10, scrollX = TRUE))
  })

  # ---- 8. export -------------------------------------------------------------
  script_params <- shiny::reactive({
    shiny::req(rv$spaces)
    sp_all <- c(default_spaces(rv$traits), rv$custom)
    list(data_file = data_name(), trait_info_file = "trait_info.csv",
         species_col = input$species_col,
         iucn_col = if (isTRUE(nzchar(input$iucn_col))) input$iucn_col else NULL,
         status_col = if (isTRUE(nzchar(input$status_col))) input$status_col else NULL,
         native_label = if (is.null(input$nn_native)) "native" else input$nn_native,
         hull_axes = if (is.null(input$nn_hull_axes)) 4 else input$nn_hull_axes,
         sep = input$sep, missing = input$missing,
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
  output$dl_nonnative <- shiny::downloadHandler(filename = "nonnative_species.csv",
    content = function(file) {
      shiny::req(rv$nonnative)
      a <- rv$nonnative$nonnative; a$group <- "non-native"
      b <- rv$nonnative$native; b$group <- "native"
      cols <- union(names(a), names(b))
      a[setdiff(cols, names(a))] <- NA; b[setdiff(cols, names(b))] <- NA
      utils::write.csv(rbind(a[cols], b[cols]), file, row.names = FALSE)
    })
  output$dl_rds <- shiny::downloadHandler(filename = "functional_spaces_results.rds",
    content = function(file) saveRDS(list(trait_info = current_ti(), traits = traits_used(),
                                          spaces = rv$spaces, comparison = rv$cmp,
                                          metrics = rv$metrics, nonnative = rv$nonnative), file))
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
      !is.null(rv$cmp), !is.null(rv$metrics), !is.null(rv$nonnative),
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
                                 if (!is.null(rv$nonnative)) "nonnative"))
  prj_params <- shiny::reactive({
    p <- script_params()
    p$top <- if (is.null(input$top)) 0.1 else input$top
    p$nonnative <- if (!is.null(rv$nonnative)) list(
      source = if (identical(input$nn_src, "simulate")) "simulate" else "column",
      native = rv$nn_native_label, prop = input$nn_prop, bias = input$nn_bias,
      sim_space = input$nn_sim_space, hull_axes = input$nn_hull_axes,
      null = input$nn_null, nrep = input$nn_nrep)
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
    ti <- current_ti(); d <- raw()
    if (!is.null(prev_preview)) unlink(prev_preview, recursive = TRUE)
    prev_preview <<- tempfile("prj_preview_")
    dir <- file.path(prev_preview, nm)
    tryCatch(write_analysis_project(d, ti, dir, params = p, phylo = prj_phylo(), steps = st,
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
      item("nonnative" %in% st, "Non-native species (tab 7)")),
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
    args <- list(data = raw(), trait_info = current_ti(), dir = dir, params = prj_params(),
                 phylo = if (!is.null(prj_phylo())) ape::read.tree(prj_phylo()), steps = prj_steps(),
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
    docx <- file.path(r$dir, "manuscript", "Methods_Results.docx")
    if (!file.exists(docx)) return(hint(if (!isTRUE(r$ok))
      "No report: the project was not run, or stopped before the report (see the Run log tab)." else
      "No rendered report (pandoc not found): knit manuscript/Methods_Results.Rmd in RStudio."))
    html <- tempfile(fileext = ".html")
    ok <- tryCatch({ rmarkdown::pandoc_convert(docx, to = "html", output = html,
                                               options = c("--embed-resources", "--standalone")); TRUE },
                   error = function(e) tryCatch({ rmarkdown::pandoc_convert(docx, to = "html", output = html,
                                                  options = "--self-contained"); TRUE }, error = function(e) FALSE))
    if (!ok) return(hint("The report is in manuscript/Methods_Results.docx."))
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
