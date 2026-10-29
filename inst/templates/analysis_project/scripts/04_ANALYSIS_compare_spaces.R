# ============================================================================
# 04_ANALYSIS_compare_spaces.R -- pairwise comparison of the functional spaces
# with Procrustes tests (and optional Mantel tests on the distances).
#   Inputs : data/processed/spaces.rds
#   Outputs: data/processed/comparison.rds
#            results/tables/Table3_procrustes.csv
#            results/figures/Fig2_procrustes
# ============================================================================
source("scripts/00_START_GeneralScript.R")

spaces <- read_object("spaces")
cmp <- compare_spaces(spaces, n_axes = PARAMS$n_axes, method = PARAMS$method,
                      nrepet = PARAMS$nrepet, mantel = PARAMS$mantel,
                      cores = PARAMS$cores)
print(cmp)
save_object(cmp, "comparison")
save_table(cmp$table, "Table3_procrustes")
save_figure(plot_comparison(cmp), "Fig2_procrustes", width = 6, height = 5)
