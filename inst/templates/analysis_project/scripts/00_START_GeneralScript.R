# ============================================================================
# 00_START_GeneralScript.R -- loads the packages, the parameters and the
# helpers. Sourced at the top of every script of the pipeline.
# ============================================================================
source("utils/library.R")
source("utils/parameters.R")
source("utils/functions.R")
set.seed(PARAMS$seed)
