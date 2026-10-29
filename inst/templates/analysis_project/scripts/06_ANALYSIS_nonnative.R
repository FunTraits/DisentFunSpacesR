# ============================================================================
# 06_ANALYSIS_nonnative.R -- non-native species relative to the native fauna
# in every space: functional novelty (FNov), position in the native convex
# hull, loss of uniqueness of native species, assemblage-level change.
#   Inputs : data/processed/spaces.rds, data/raw/traits.csv (status column)
#   Outputs: data/processed/nonnative.rds
#            results/tables/Table5_nonnative_summary.csv, Table6_assemblage_change.csv,
#            TableS6_nonnative_species.csv, TableS7_native_species.csv
#            results/figures/Fig4_nonnative_<space>
# ============================================================================
source("scripts/00_START_GeneralScript.R")

spaces <- read_object("spaces")
nn_par <- PARAMS$nonnative
sp <- rownames(spaces[[1]]$coords)

if (identical(nn_par$source, "simulate")) {
  # exploration: introductions simulated in one space, as in the application
  status <- simulate_introductions(spaces[[nn_par$sim_space]], prop = nn_par$prop,
                                   bias = nn_par$bias, seed = PARAMS$seed)
  utils::write.csv(data.frame(species = names(status), status = unname(status)),
                   file.path("data", "processed", "species_status_simulated.csv"), row.names = FALSE)
} else {
  status <- species_vector(read_traits_table(), PARAMS$status_col)
  status <- stats::setNames(status[sp], sp)
}

nn <- nonnative_all(spaces, status, native = nn_par$native, n_axes = PARAMS$n_axes,
                    k = PARAMS$k, hull_axes = nn_par$hull_axes, null = nn_par$null,
                    nrep = nn_par$nrep, seed = PARAMS$seed)
print(nn)
save_object(nn, "nonnative")
save_object(status, "species_status")

save_table(summary_nonnative(nn), "Table5_nonnative_summary")
if (!is.null(nn$assemblage)) save_table(nn$assemblage, "Table6_assemblage_change")
save_table(nn$nonnative, "TableS6_nonnative_species")
save_table(nn$native, "TableS7_native_species")

for (nm in names(spaces)) {
  s <- spaces[[nm]]
  r <- list(nonnative = nn$nonnative[nn$nonnative$space == nm, -1],
            settings = list(hull_axes = min(nn_par$hull_axes, ncol(s$coords))))
  save_figure(plot_nonnative(s, status, x = r, native = nn_par$native, label = 5),
              paste0("Fig4_nonnative_", safe_name(nm)))
}
