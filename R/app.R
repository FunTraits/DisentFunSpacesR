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
#' @param ... Further arguments to [shiny::runApp()].
#' @return Called for its side effect.
#' @export
#' @examples
#' if (interactive()) run_app()
run_app <- function(launch.browser = TRUE, max_upload_mb = 200, ...) {
  options(shiny.maxRequestSize = max_upload_mb * 1024^2)
  shiny::runApp(fs_app(), launch.browser = launch.browser, ...)
}

#' @rdname run_app
#' @return `fs_app()` returns a `shiny.appobj`.
#' @export
fs_app <- function() shiny::shinyApp(ui = app_ui(), server = app_server)

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------
app_ui <- function() {
  sh <- shiny::tags
  type_choices <- stats::setNames(trait_types()$code,
                                  paste0(trait_types()$code, " - ", trait_types()$label))
  bslib::page_navbar(
    title = "DisentFunSpaces",
    id = "nav",
    bg = "#1F4E79",
    inverse = TRUE,
    theme = bslib::bs_theme(version = 5, primary = "#1F4E79",
                            base_font = bslib::font_collection("system-ui", "-apple-system",
                                                               "Segoe UI", "sans-serif")),
    header = sh$style(shiny::HTML("
      .ti-table td, .ti-table th {padding: 2px 6px; vertical-align: middle; font-size: 0.85rem;}
      .ti-table .form-group {margin-bottom: 0;}
      .ti-table select, .ti-table input {font-size: 0.85rem; padding: 2px 4px; height: 30px;}
      .ti-table th {position: sticky; top: 0; background: #f4f6f8; z-index: 2;}
      .step-help {color: #555; font-size: 0.9rem;}
      .shiny-notification {font-size: 0.9rem;}")),

    # 1. DATA -----------------------------------------------------------------
    bslib::nav_panel("1. Data",
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(width = 320,
          shiny::radioButtons("src", "Data source",
            c("Example (2000 bird species, real data)" = "example", "Upload a file" = "upload")),
          shiny::conditionalPanel("input.src == 'upload'",
            shiny::fileInput("file", "Trait table (.csv, .tsv, .txt, .xlsx)",
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
                    DT::DTOutput("data_preview"))
      )
    ),

    # 2. TRAITS ---------------------------------------------------------------
    bslib::nav_panel("2. Traits",
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
    bslib::nav_panel("3. Missing values",
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
    bslib::nav_panel("4. Functional spaces",
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
    bslib::nav_panel("5. Compare spaces",
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
    bslib::nav_panel("6. Species metrics",
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
    bslib::nav_panel("7. Non-native species",
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
    bslib::nav_panel("8. Export",
      bslib::card(
        bslib::card_header("Download results"),
        sh$p("Every table can also be copied from its own tab. The R script reproduces the",
             " analysis with the package functions, using the trait-information CSV."),
        bslib::layout_column_wrap(width = 1 / 3,
          shiny::downloadButton("dl_script", "R script (.R)"),
          shiny::downloadButton("dl_ti2", "Trait information (.csv)"),
          shiny::downloadButton("dl_coords", "PCoA coordinates (.csv)"),
          shiny::downloadButton("dl_cmp", "Procrustes table (.csv)"),
          shiny::downloadButton("dl_metrics", "Species metrics (.csv)"),
          shiny::downloadButton("dl_nonnative", "Non-native analysis (.csv)"),
          shiny::downloadButton("dl_rds", "All R objects (.rds)")),
        shiny::hr(),
        sh$p(class = "step-help",
          "Please cite: Toussaint A., Tedesco P., Grenouillet G., Kasari-Toussaint L., Brosse S. (2026)",
          " Disentangling functional spaces: a multidimensional framework for ecology and conservation.",
          " Ecography e08578. doi:10.1002/ecog.08578"),
        shiny::verbatimTextOutput("script_preview")
      )
    ),
    bslib::nav_spacer(),
    bslib::nav_item(sh$a("Paper", href = "https://doi.org/10.1002/ecog.08578", target = "_blank"))
  )
}

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------
app_server <- function(input, output, session) {
  rv <- shiny::reactiveValues(ti = NULL, traits = NULL, traits_final = NULL, spaces = NULL,
                              cmp = NULL, metrics = NULL, nonnative = NULL,
                              nn_status = NULL, custom = list(), ti_version = 0,
                              valid_msg = NULL, missing_msg = NULL, phylo_name = NULL)
  notify_error <- function(e) shiny::showNotification(conditionMessage(e), type = "error", duration = 10)

  # ---- 1. data --------------------------------------------------------------
  raw <- shiny::reactive({
    if (input$src == "example") {
      e <- new.env(); utils::data("birds_traits", package = "DisentFunSpaces", envir = e)
      return(e$birds_traits)
    }
    shiny::req(input$file)
    path <- input$file$datapath
    tryCatch({
      if (grepl("\\.xlsx?$", input$file$name, ignore.case = TRUE)) {
        if (!requireNamespace("readxl", quietly = TRUE)) stop("Install 'readxl' to read Excel files.")
        as.data.frame(readxl::read_excel(path))
      } else {
        utils::read.table(path, header = TRUE, sep = input$sep, quote = "\"",
                          stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("NA", ""),
                          comment.char = "")
      }
    }, error = function(e) { notify_error(e); NULL })
  })
  data_name <- shiny::reactive(if (input$src == "example") "birds_traits.csv" else input$file$name)

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
    rv$ti <- NULL; rv$traits <- NULL; rv$traits_final <- NULL
    rv$spaces <- NULL; rv$cmp <- NULL; rv$metrics <- NULL; rv$nonnative <- NULL
  })

  output$data_info <- shiny::renderUI({
    d <- raw(); shiny::req(d)
    shiny::tags$p(sprintf("%d species x %d columns - %.1f%% missing cells",
                          nrow(d), ncol(d), 100 * mean(is.na(d))))
  })
  output$data_preview <- DT::renderDT(DT::datatable(raw(), options = list(scrollX = TRUE, pageLength = 10),
                                                    rownames = FALSE))

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
    if (identical(input$src, "example") && identical(sc, "species")) {
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
      build_spaces(tr, spaces = all_sp, weighting = input$weighting, correction = input$correction,
                   max_axes = mx, envfit_perm = 0,
                   progress = function(i, n, nm) shiny::setProgress(i / n, detail = nm))
    }), error = function(e) { notify_error(e); NULL })
    shiny::req(res)
    rv$spaces <- res; rv$cmp <- NULL; rv$metrics <- NULL
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
    p <- plot_space(s, axes = ax, min_cor = input$min_cor, color_by = col, arrows = is.null(col))
    if (!is.null(col)) p <- p + ggplot2::scale_colour_manual(
      values = c(LC = "#3B9AB2", NT = "#78B7C5", VU = "#EBCC2A", EN = "#E1AF00", CR = "#F21A00", DD = "grey60"),
      drop = TRUE, na.value = "grey80")
    p
  })
  output$space_plot <- shiny::renderPlot(space_plot(), res = 96)
  output$dl_space_plot <- shiny::downloadHandler(
    filename = function() paste0("space_", make.names(input$plot_space), ".pdf"),
    content = function(file) ggplot2::ggsave(file, space_plot(), width = 7, height = 6))
  output$eigen_plot <- shiny::renderPlot({ shiny::req(rv$spaces); plot_eigen(rv$spaces) }, res = 96)
  output$quality_plot <- shiny::renderPlot({
    s <- sel_space()
    shiny::validate(shiny::need(!is.null(s$dist), "Distance not stored (> 5000 species)."))
    plot_quality(stats::setNames(list(space_quality(s)), s$name))
  }, res = 96)
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
    res <- tryCatch(shiny::withProgress(message = "Procrustes tests", value = 0, {
      compare_spaces(rv$spaces, n_axes = n_axes_val(), method = input$pr_method,
                     nrepet = input$nrepet, mantel = input$mantel,
                     progress = function(i, n, nm) shiny::setProgress(i / n, detail = nm))
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
  output$res_plot <- shiny::renderPlot(plot_procrustes(residuals_r()), res = 96)
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
    rv$metrics <- tryCatch(species_metrics_all(rv$spaces, n_axes = n_axes_val(), k = input$k_nn,
                                               iucn = iucn_vec()),
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
}
