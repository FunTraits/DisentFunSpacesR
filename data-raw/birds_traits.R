## Build the example datasets `birds_traits`, `birds_trait_info`, `birds_tree`
## ---------------------------------------------------------------------------
## Real data: 1000 bird species selected among the species analysed in
## Toussaint et al. (2026, Ecography, doi:10.1002/ecog.08578), described by
## the traits used in that paper, with their OBSERVED (non-imputed,
## untransformed) values:
##   * Locomotion (AVONET; Tobias et al. 2022, CC BY 4.0): tarsus length, wing
##     length, Kipp's distance, tail length, body mass;
##   * Reproduction (Amniote database; Myhrvold et al. 2015, CC0): clutch
##     size, clutches per year, adult body mass, incubation time, longevity,
##     fledging age, egg mass, body length;
##   * Diet (EltonTraits 1.0; Wilman et al. 2014, CC0): 10 diet categories
##     (fuzzy-coded, % summing to 100);
##   * Beak (AVONET; validation space of the paper, `use = FALSE` by
##     default): culmen length, beak width, beak depth.
## Also: IUCN Red List category (API v4, rredlist), the first 10 phylogenetic
## eigenvectors (PCoA of the square-rooted cophenetic distances of the full
## tree, as used for imputation in the paper) and the phylogeny pruned to
## the 1000 species.
##
## Species selection: species of the analysis present in the phylogeny and in
## the eigenvector table; all species with no missing value (281) plus the
## species with the fewest missing values, up to 1000 (ties at random, fixed
## seed). The sample is therefore biased towards well-studied species.
##
## Run from the package root:  source("data-raw/birds_traits.R")
## ---------------------------------------------------------------------------

FUSE <- path.expand(Sys.getenv(
  "FUSE_BIRD_DIR",
  "~/Library/CloudStorage/OneDrive-Personnel/iCloud Drive/00_Papier_2_Submitted/FUSE_bird"))
N_SPECIES <- 1000
SEED <- 2026

ready <- utils::read.csv(file.path(FUSE, "data/processed/phenoBirdsImputedREADY.csv"),
                         check.names = FALSE)                   # species of the analysis
raw   <- utils::read.csv(file.path(FUSE, "data/processed/phenoBirds.csv"),
                         check.names = FALSE)                   # observed values
names(raw)[1] <- "species"
tree  <- ape::read.tree(file.path(FUSE, "data/processed/BirdPhylogeny.tre"))
eig   <- utils::read.delim(file.path(FUSE, "data/processed/pcoaPhylogenyAves.txt"))
avonet <- as.data.frame(readxl::read_excel(file.path(FUSE, "data/raw/AVONET.xlsx"),
                                           sheet = "AVONET1_BirdLife"))   # taxonomy

## ---- species pool and random draw ------------------------------------------------
d <- raw[match(ready$scientificNameStd, raw$species), ]
stopifnot(identical(d$species, ready$scientificNameStd))
d_all <- d; ready_all <- ready                       # all 8,413 species (birds_traits_all)
## Selection: all species with no missing value for the 13 continuous traits
## and the diet (281 species), completed to 1000 with the species having the
## fewest missing values (ties at the cut-off broken at random, fixed seed).
loco_cols  <- c("Tarsus.Length", "Wing.Length", "Kipps.Distance", "Tail.Length", "Mass")
repro_cols <- c("litter_or_clutch_size_n", "litters_or_clutches_per_y", "adult_body_mass_g",
                "incubation_d", "longevity_y", "fledging_age_d", "egg_mass_g", "adult_svl_cm")
diet_cols  <- grep("^Diet\\.", names(d), value = TRUE)
pool <- d$species %in% tree$tip.label & d$species %in% eig$scientificNameStd &
  stats::complete.cases(d[, diet_cols])
d <- d[pool, ]; ready <- ready[pool, ]
n_na <- rowSums(is.na(d[, c(loco_cols, repro_cols)]))
message(nrow(d), " species in the pool; ", sum(n_na == 0), " with no missing value")
set.seed(SEED)
ord <- order(n_na, stats::runif(length(n_na)))      # fewest NA first, random ties
pick <- sort(d$species[ord[seq_len(N_SPECIES)]])
message("maximum number of missing traits among selected species: ",
        max(n_na[match(pick, d$species)]))
i <- match(pick, d$species)
d <- d[i, ]; ready <- ready[i, ]
av <- avonet[match(gsub("_", " ", pick), avonet$Species1), ]   # BirdLife order and family
stopifnot(!anyNA(av$Species1))

## ---- assemble the trait tables ------------------------------------------------------
## Phylogenetic eigenvectors are those of the published analysis (Eigen.1-10 of
## phenoBirdsImputedREADY.csv; identical to pcoaPhylogenyAves.txt for the
## species present in the tree).
make_table <- function(d, ready, av) {
  out <- data.frame(
  species = ready$scientificNameStd,
  order   = av$Order1,
  family  = av$Family1,
  # Locomotion (AVONET; mm, g)
  tarsus_length  = d$Tarsus.Length,
  wing_length    = d$Wing.Length,
  kipps_distance = d$Kipps.Distance,
  tail_length    = d$Tail.Length,
  body_mass      = d$Mass,
  # Reproduction (Amniote database; n, n per year, g, days, years, days, g, cm)
  clutch_size      = d$litter_or_clutch_size_n,
  clutches_per_year = d$litters_or_clutches_per_y,
  adult_body_mass  = d$adult_body_mass_g,
  incubation_days  = d$incubation_d,
  longevity_years  = d$longevity_y,
  fledging_days    = d$fledging_age_d,
  egg_mass         = d$egg_mass_g,
  body_length      = d$adult_svl_cm,
  # Diet (EltonTraits; % of diet, sums to 100)
  diet_inv    = d$Diet.Inv,   diet_vend  = d$Diet.Vend,  diet_vect  = d$Diet.Vect,
  diet_vfish  = d$Diet.Vfish, diet_vunk  = d$Diet.Vunk,  diet_scav  = d$Diet.Scav,
  diet_fruit  = d$Diet.Fruit, diet_nect  = d$Diet.Nect,  diet_seed  = d$Diet.Seed,
  diet_planto = d$Diet.PlantO,
  # Beak (AVONET; mm) - validation space of the paper
  beak_culmen = d$Beak.Length_Culmen,
  beak_width  = d$Beak.Width,
  beak_depth  = d$Beak.Depth,
  # threat status and imputation predictors (not traits)
  iucn = ready$iucn_category,
  ready[, paste0("Eigen.", 1:10)],
  stringsAsFactors = FALSE, row.names = NULL
)
  num <- vapply(out, is.numeric, logical(1))
  out[num] <- lapply(out[num], function(x) signif(x, 6))
  out
}
birds_traits <- make_table(d, ready, av)
stopifnot(identical(birds_traits$species, pick))
av_all <- avonet[match(gsub("_", " ", ready_all$scientificNameStd), avonet$Species1), ]
stopifnot(!anyNA(av_all$Species1))
birds_traits_all <- make_table(d_all, ready_all, av_all)


## ---- trait-information table ----------------------------------------------------
tr <- function(trait, type, dimension, block = trait, transform = "none", use = TRUE)
  data.frame(trait, type, dimension, block, transform, weight = 1, levels = NA_character_,
             circ_min = NA_real_, circ_max = NA_real_, use = use, stringsAsFactors = FALSE)
diet <- grep("^diet_", names(birds_traits), value = TRUE)
birds_trait_info <- rbind(
  tr(c("tarsus_length", "wing_length", "kipps_distance", "tail_length", "body_mass"),
     "Q", "Locomotion", transform = "log10"),
  tr(c("clutch_size", "clutches_per_year", "adult_body_mass", "incubation_days",
       "longevity_years", "fledging_days", "egg_mass", "body_length"),
     "Q", "Reproduction", transform = "log10"),
  tr(diet, "F", "Diet", block = "diet"),
  tr(c("beak_culmen", "beak_width", "beak_depth"), "Q", "Beak", transform = "log10", use = FALSE)
)
rownames(birds_trait_info) <- NULL

## ---- phylogeny -----------------------------------------------------------------------
birds_tree <- ape::keep.tip(tree, pick)
birds_tree$node.label <- NULL
stopifnot(setequal(birds_tree$tip.label, birds_traits$species))

## ---- checks and save -------------------------------------------------------------
DisentFunSpaces::check_trait_info(birds_trait_info, birds_traits)
print(table(birds_traits$iucn))
print(round(100 * colMeans(is.na(birds_traits[birds_trait_info$trait])), 1))
print(mean(stats::complete.cases(birds_traits[birds_trait_info$trait[birds_trait_info$use]])))
save(birds_traits, file = "data/birds_traits.rda", compress = "xz")
save(birds_trait_info, file = "data/birds_trait_info.rda", compress = "xz")
save(birds_tree, file = "data/birds_tree.rda", compress = "xz")
DisentFunSpaces::check_trait_info(birds_trait_info, birds_traits_all)
save(birds_traits_all, file = "data/birds_traits_all.rda", compress = "xz")
