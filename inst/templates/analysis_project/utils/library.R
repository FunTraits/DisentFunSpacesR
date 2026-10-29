# ============================================================================
# Packages used by the pipeline
# ============================================================================
# CRAN packages that are missing are installed; DisentFunSpaces itself is
# installed from GitHub when missing. Versions that produced the results:
# session_info.txt.

cran <- c({{CRAN_PACKAGES}})
missing_pkgs <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) install.packages(missing_pkgs)

if (!requireNamespace("DisentFunSpaces", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  remotes::install_github("FunTraits/DisentFunSpacesR")
}

suppressPackageStartupMessages({
  library(DisentFunSpaces)
  library(ggplot2)
})

if (as.character(utils::packageVersion("DisentFunSpaces")) != "{{DFS_VERSION}}")
  message("Note: this project was produced with DisentFunSpaces {{DFS_VERSION}}; ",
          "version ", utils::packageVersion("DisentFunSpaces"), " is installed.")
