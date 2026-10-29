# ============================================================================
# setup_environment.R -- installs the packages of the pipeline once (not run
# by run_all.R). Optionally, freeze the versions with renv afterwards:
#   renv::init(); renv::snapshot()
# ============================================================================
cran <- c({{CRAN_PACKAGES}}, "remotes", "rmarkdown", "knitr", "here")
install.packages(setdiff(cran, rownames(installed.packages())))
remotes::install_github("FunTraits/DisentFunSpacesR")
