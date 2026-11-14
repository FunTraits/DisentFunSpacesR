# ============================================================================
# setup_environment.R -- checks which packages of the pipeline are missing
# (not run by run_all.R). This script installs nothing itself: it lists the
# missing packages, to be installed by the user from CRAN (e.g. with the
# 'Packages' pane of RStudio). Optionally, freeze
# the versions with renv afterwards: renv::init(); renv::snapshot()
# ============================================================================
pkgs <- c({{CRAN_PACKAGES}}, "DisentFunSpaces", "rmarkdown", "knitr")
missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) {
  message("Packages to install from CRAN:\n  ",
          paste0("\"", missing_pkgs, "\"", collapse = ", "))
} else {
  message("All the packages of the pipeline are installed.")
}
