# ============================================================================
# Packages used by the pipeline
# ============================================================================
# Nothing is installed here: missing packages stop the pipeline with the list
# of what to install (see scripts/setup_environment.R). Versions that
# produced the results: session_info.txt.

cran <- c({{CRAN_PACKAGES}}, "DisentFunSpaces")
missing_pkgs <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs))
  stop("Missing packages: ", paste(missing_pkgs, collapse = ", "), ". ",
       "Install them from CRAN (see scripts/setup_environment.R), then run again.",
       call. = FALSE)

suppressPackageStartupMessages({
  library(DisentFunSpaces)
  library(ggplot2)
})

if (as.character(utils::packageVersion("DisentFunSpaces")) != "{{DFS_VERSION}}")
  message("Note: this project was produced with DisentFunSpaces {{DFS_VERSION}}; ",
          "version ", utils::packageVersion("DisentFunSpaces"), " is installed.")
