# ============================================================================
# 06_ANALYSIS_groups.R -- species groups in every space. Groups are defined by
# any categorical attribute that is not a trait (origin, human use, migration,
# activity period, threat status, taxonomy...).
#   Groups   : hull volume, dispersion, pairwise overlap of convex hulls
#              (functional Jaccard, turnover / nestedness), PERMANOVA / PERMDISP.
#   Reference: when a reference group is set (e.g. natives), novelty (FNov) of
#              each focal species, position in the reference hull, loss of
#              uniqueness of reference species, contribution of the focal group.
#   Inputs : data/processed/spaces.rds, data/raw/traits.csv (grouping column)
#   Outputs: data/processed/groups.rds
#            results/tables/Table5_groups_summary.csv, Table6_group_pairs.csv,
#            TableS6_groups.csv (+ Table7_focal_summary.csv, TableS7_focal_species.csv,
#            TableS8_reference_species.csv, TableS9_group_contribution.csv)
#            results/figures/Fig4_groups_<space>, Fig5_group_overlap
#            (+ continuous variable: Table9, TableS12, Fig7_continuous_<space>;
#               assemblages: TableS13_group_assemblages, Fig8_group_assemblages)
# ============================================================================
source("scripts/00_START_GeneralScript.R")

spaces <- read_object("spaces")
gp <- PARAMS$groups
sp <- rownames(spaces[[1]]$coords)
traits_raw <- read_traits_table()

# --- continuous variable that is not a trait ---------------------------------
cv <- PARAMS$covariate
if (!is.null(cv)) {
  x <- suppressWarnings(as.numeric(species_vector(traits_raw, cv$col)))
  names(x) <- as.character(traits_raw[[PARAMS$species_col]])
  cov <- covariate_space(spaces, x, transform = cv$transform, n_axes = PARAMS$n_axes, k = PARAMS$k,
                         fit_axes = cv$fit_axes, nperm = cv$nperm, seed = PARAMS$seed)
  print(cov)
  save_object(cov, "continuous")
  save_table(cov$summary, "Table9_continuous_variable")
  save_table(cov$species, "TableS12_continuous_species")
  for (nm in names(spaces))
    save_figure(plot_covariate(spaces[[nm]], x, transform = cv$transform, name = cv$col),
                paste0("Fig7_continuous_", safe_name(nm)))
}
# --- species groups -------------------------------------------------------------
if (!is.null(gp)) {

  if (identical(gp$source, "simulate")) {
    # exploration: native / non-native status simulated in one space, as in the application
    groups <- simulate_introductions(spaces[[gp$sim_space]], prop = gp$prop,
                                     bias = gp$bias, seed = PARAMS$seed)
    utils::write.csv(data.frame(species = names(groups), group = unname(groups)),
                     file.path("data", "processed", "species_groups_simulated.csv"), row.names = FALSE)
  } else {
    groups <- species_vector(traits_raw, gp$col)
    groups <- stats::setNames(groups[sp], sp)
  }
  # abundance weights of the densities, and assemblages (long table)
  abundance <- if (isTRUE(gp$density) && !is.null(gp$abundance_col))
    stats::setNames(as.numeric(traits_raw[[gp$abundance_col]]), as.character(traits_raw[[PARAMS$species_col]]))
  comm <- if (!is.null(gp$sites)) {
    a <- utils::read.csv(file.path("data", "raw", "assemblages.csv"), check.names = FALSE,
                         stringsAsFactors = FALSE)
    community_from_long(a, gp$sites$site_col, gp$sites$species_col, gp$sites$group_col)
  }

  grp <- groups_all(spaces, groups, reference = gp$reference, pool_focal = gp$pool_focal,
                    sep = gp$sep, recode = gp$recode, levels = gp$levels,
                    n_axes = PARAMS$n_axes, k = PARAMS$k, hull_axes = gp$hull_axes,
                    null = gp$null, nrep = gp$nrep, nperm = gp$nperm, seed = PARAMS$seed,
                    density = isTRUE(gp$density), density_axes = gp$density_axes, bw = gp$bw,
                    abundance = abundance, comm = comm)
  print(grp)
  save_object(grp, "groups")
  save_object(groups, "species_groups")

  save_table(summary_groups(grp), "Table5_groups_summary")
  save_table(grp$groups, "TableS6_groups")
  if (!is.null(grp$pairs) && nrow(grp$pairs)) {
    save_table(grp$pairs, "Table6_group_pairs")
    save_figure(plot_group_overlap(grp, "overlap"), "Fig5_group_overlap")
  }
  if (!is.null(grp$focal)) {
    save_table(summary_groups(grp, "focal"), "Table7_focal_summary")
    save_table(grp$focal, "TableS7_focal_species")
    save_table(grp$reference, "TableS8_reference_species")
    if (!is.null(grp$contribution)) save_table(grp$contribution, "TableS9_group_contribution")
  }
  if (!is.null(grp$sites)) {
    save_table(grp$sites, "TableS13_group_assemblages")
    save_figure(plot_group_sites(grp, "FRic_change_rel"), "Fig8_group_assemblages", width = 7, height = 4)
  }

  for (nm in names(spaces)) {
    s <- spaces[[nm]]
    fig <- if (!is.null(grp$focal)) {
      r <- list(focal = grp$focal[grp$focal$space == nm, -1],
                settings = list(hull_axes = min(gp$hull_axes, ncol(s$coords))))
      plot_groups(s, groups, sep = gp$sep, recode = gp$recode, levels = gp$levels,
                  reference = gp$reference, pool_focal = gp$pool_focal, x = r, label = 5)
    } else plot_groups(s, groups, sep = gp$sep, recode = gp$recode, levels = gp$levels)
    save_figure(fig, paste0("Fig4_groups_", safe_name(nm)))
  }
}
