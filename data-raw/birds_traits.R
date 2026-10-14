## Build the example datasets `birds_traits`, `birds_trait_info`, `birds_tree`
## ---------------------------------------------------------------------------
## Real data: 2000 bird species drawn at random from the species analysed in
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
## the 2000 species.
##
## Species pool: the 8,413 species of the analysis present in the phylogeny
## and in the eigenvector table (7,619 species). 2000 species are drawn by
## simple random sampling (fixed seed).
##
## Run from the package root:  source("data-raw/birds_traits.R")
## ---------------------------------------------------------------------------

FUSE <- path.expand(Sys.getenv(
  "FUSE_BIRD_DIR",
  "~/Library/CloudStorage/OneDrive-Personnel/iCloud Drive/00_Papier_2_Submitted/FUSE_bird"))
N_SPECIES <- 2000
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
pool <- d$species %in% tree$tip.label & d$species %in% eig$scientificNameStd
d <- d[pool, ]; ready <- ready[pool, ]
message(nrow(d), " species in the pool")
set.seed(SEED)
pick <- sort(sample(d$species, N_SPECIES))
i <- match(pick, d$species)
d <- d[i, ]; ready <- ready[i, ]
av <- avonet[match(gsub("_", " ", pick), avonet$Species1), ]   # BirdLife order and family
stopifnot(!anyNA(av$Species1))

## ---- assemble the trait table -----------------------------------------------------
birds_traits <- data.frame(
  species = pick,
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
  eig[match(pick, eig$scientificNameStd), paste0("Eigen.", 1:10)],
  stringsAsFactors = FALSE, row.names = NULL
)
num <- vapply(birds_traits, is.numeric, logical(1))
birds_traits[num] <- lapply(birds_traits[num], function(x) signif(x, 6))

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
