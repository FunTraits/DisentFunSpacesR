# ============================================================================
# 02_DATA_missing_values.R -- handling of missing trait values:
#   "keep"   the generalised Gower distance uses, for each pair of species,
#            the traits known for both;
#   "drop"   species with a missing value are removed;
#   "impute" missing values are predicted with missForest, from the other
#            traits, optional predictor columns and phylogenetic eigenvectors.
#   Inputs : data/processed/traits_prepared.rds (+ data/raw/phylogeny.tre)
#   Outputs: data/processed/traits_final.rds
#            results/tables/TableS1_missing_values.csv
#            results/tables/TableS2_imputation_error.csv (imputation only)
# ============================================================================
source("scripts/00_START_GeneralScript.R")

traits <- read_object("traits_prepared")
missing_before <- trait_summary(traits)[, c("trait", "dimension", "type", "pct_missing")]

if (PARAMS$missing == "drop") {
  traits <- drop_incomplete(traits)

} else if (PARAMS$missing == "impute") {
  tree <- if (!is.null(PARAMS$phylo_file)) ape::read.tree(PARAMS$phylo_file)
  predictors <- NULL
  if (length(PARAMS$cov_cols)) {
    traits_raw <- read_traits_table()
    predictors <- traits_raw[, PARAMS$cov_cols, drop = FALSE]
    rownames(predictors) <- as.character(traits_raw[[PARAMS$species_col]])
  }
  traits <- impute_traits(traits, phylo = tree, covariates = predictors,
                          n_phylo_axes = PARAMS$n_phylo_axes, ntree = PARAMS$ntree,
                          seed = PARAMS$seed, verbose = FALSE)
  oob <- traits$imputation$oob
  save_table(data.frame(error = names(oob), value = unname(oob)), "TableS2_imputation_error")
}

missing_after <- trait_summary(traits)[, c("trait", "pct_missing")]
names(missing_before)[4] <- "pct_missing_raw"
names(missing_after)[2] <- "pct_missing_analysed"
tab <- merge(missing_before, missing_after, by = "trait", sort = FALSE)
tab$n_species_analysed <- length(traits$species)
save_table(tab, "TableS1_missing_values")
save_object(traits, "traits_final")
