# ============================================================================
# 05_ANALYSIS_species_metrics.R -- species-level functional uniqueness (FUn),
# specialisation (FSp) and, with IUCN categories, FUSE in every space; do the
# spaces identify the same distinct species?
#   Inputs : data/processed/spaces.rds, data/raw/traits.csv (IUCN column)
#   Outputs: data/processed/metrics.rds
#            results/tables/Table4_metric_congruence.csv, TableS5_species_metrics.csv
#            results/figures/Fig3_metric_congruence
# ============================================================================
source("scripts/00_START_GeneralScript.R")

spaces <- read_object("spaces")
iucn <- species_vector(read_traits_table(), PARAMS$iucn_col)
metrics <- species_metrics_all(spaces, n_axes = PARAMS$n_axes, k = PARAMS$k, iucn = iucn)
save_object(metrics, "metrics")
save_table(metrics, "TableS5_species_metrics")

if (length(spaces) > 1) {
  mets <- intersect(c("FUn", "FSp", "FUSE"), names(metrics))
  cong <- do.call(rbind, lapply(mets, function(m)
    cbind(metric = m, compare_metrics(metrics, m, top = PARAMS$top)$pairs)))
  save_table(cong, "Table4_metric_congruence")

  cong$pair <- paste(cong$space1, "vs", cong$space2)
  p <- ggplot(cong, aes(x = metric, y = pair, fill = spearman)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.2f", spearman)), size = 3) +
    scale_fill_gradient2(limits = c(-1, 1), low = "#b2182b", mid = "white", high = "#2166ac") +
    labs(x = NULL, y = NULL, fill = "Spearman rho",
         title = "Congruence of species rankings between spaces") +
    theme_minimal(base_size = 11)
  save_figure(p, "Fig3_metric_congruence", width = 6.5, height = 0.4 * nrow(cong) / length(mets) + 2.5)
}
