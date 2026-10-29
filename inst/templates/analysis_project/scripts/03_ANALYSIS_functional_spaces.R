# ============================================================================
# 03_ANALYSIS_functional_spaces.R -- one functional space per dimension (and
# per combination of dimensions): generalised Gower distance between species
# (ade4::dist.ktab) and Principal Coordinates Analysis (ape::pcoa).
#   Inputs : data/processed/traits_final.rds
#   Outputs: data/processed/spaces.rds
#            results/tables/Table2_spaces.csv, TableS3_pcoa_coordinates.csv,
#            TableS4_trait_axis_correlations.csv
#            results/figures/Fig1_space_<space>, FigS1_inertia, FigS2_quality
# ============================================================================
source("scripts/00_START_GeneralScript.R")

traits <- read_object("traits_final")
spaces <- build_spaces(traits, spaces = PARAMS$spaces,
                       weighting = PARAMS$weighting, correction = PARAMS$correction,
                       max_axes = PARAMS$max_axes, engine = PARAMS$engine,
                       eigen = PARAMS$eigen, keep_dist = PARAMS$keep_dist,
                       envfit_perm = 0)
print(spaces)
save_object(spaces, "spaces")

save_table(space_summary(spaces), "Table2_spaces")
save_table(coords_table(spaces, n_axes = PARAMS$coords_axes), "TableS3_pcoa_coordinates")
save_table(bind_fill(lapply(names(spaces), function(nm)
  cbind(space = nm, spaces[[nm]]$cor, stringsAsFactors = FALSE))),
  "TableS4_trait_axis_correlations")

for (nm in names(spaces))
  save_figure(plot_space(spaces[[nm]]), paste0("Fig1_space_", safe_name(nm)))
save_figure(plot_eigen(spaces), "FigS1_inertia", width = 7, height = 4.5)

with_dist <- Filter(function(s) !is.null(s$dist), spaces)
if (length(with_dist))
  save_figure(plot_quality(lapply(with_dist, space_quality)), "FigS2_quality",
              width = 7, height = 4.5)
