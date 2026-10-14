#' DisentFunSpaces: build and compare functional spaces from mixed-type traits
#'
#' Workflow (Toussaint et al. 2026, Ecography, doi:10.1002/ecog.08578):
#' 1. describe traits: [guess_trait_info()] / [check_trait_info()]
#'    (type and functional dimension of every trait);
#' 2. [prepare_traits()], then optionally [impute_traits()] or
#'    [drop_incomplete()];
#' 3. [build_spaces()]: one PCoA space per dimension and for any combination
#'    (generalised Gower distance, [trait_distance()]);
#' 4. [compare_spaces()]: pairwise Procrustes tests;
#'    [procrustes_residuals()] for species-level discrepancies;
#' 5. [species_metrics_all()] and [compare_metrics()]: FUn, FSp, FUSE.
#'
#' [run_app()] launches the Shiny application.
#' @keywords internal
"_PACKAGE"

#' Example dataset: traits of 2000 bird species
#'
#' Observed (non-imputed, untransformed) traits of 2000 bird species selected
#' among the species analysed by Toussaint et al. (2026, *Ecography*,
#' \doi{10.1002/ecog.08578}), described by the traits used in that paper:
#' locomotion (five continuous traits, AVONET), reproduction (eight
#' continuous traits, Amniote database) and diet (one fuzzy-coded trait with
#' ten categories, EltonTraits). Three beak traits, used in the paper to
#' build a validation space, are included with `use = FALSE` in
#' [birds_trait_info]. The table also holds the IUCN Red List category of
#' each species and the first ten phylogenetic eigenvectors, used as
#' imputation predictors in the paper. Script: `data-raw/birds_traits.R`
#' (source repository).
#'
#' **Selection.** All 281 species of the published analysis present in the
#' phylogeny that have no missing value for the 23 traits, plus 1719 species
#' drawn at random (fixed seed) among those with at least one trait informed
#' in each dimension (locomotion, reproduction, diet), for a total of 2000
#' species. BirdLife taxonomy; 31 orders and 182 families (1065
#' Passeriformes).
#'
#' **Missing values** (19 % of trait values) concern reproduction traits
#' only: `adult_body_mass` (2 %), `clutch_size` (26 %), `egg_mass` (45 %),
#' `incubation_days` (69 %), `fledging_days` (73 %), `clutches_per_year`
#' (74 %), `body_length` (75 %) and `longevity_years` (75 %). As in the
#' published analysis, they should be imputed (see [impute_traits()] and the
#' vignette).
#'
#' @format A data frame with 2000 rows and 40 columns:
#' \describe{
#'   \item{species}{Binomial name, BirdLife taxonomy (`Genus_species`).}
#'   \item{order, family}{BirdLife taxonomy from AVONET (not traits).}
#'   \item{tarsus_length, wing_length, kipps_distance, tail_length}{Locomotion:
#'     lengths (mm), AVONET.}
#'   \item{body_mass}{Locomotion: body mass (g), AVONET.}
#'   \item{clutch_size}{Reproduction: number of eggs per clutch.}
#'   \item{clutches_per_year}{Reproduction: number of clutches per year.}
#'   \item{adult_body_mass}{Reproduction: adult body mass (g).}
#'   \item{incubation_days}{Reproduction: incubation time (days).}
#'   \item{longevity_years}{Reproduction: longevity (years).}
#'   \item{fledging_days}{Reproduction: fledging age (days).}
#'   \item{egg_mass}{Reproduction: egg mass (g).}
#'   \item{body_length}{Reproduction: body length, tip of beak to cloaca
#'     (cm). All reproduction traits from the Amniote database.}
#'   \item{diet_inv, diet_vend, diet_vect, diet_vfish, diet_vunk, diet_scav,
#'     diet_fruit, diet_nect, diet_seed, diet_planto}{Diet: percentage of
#'     invertebrates, endothermic vertebrates, ectothermic vertebrates, fish,
#'     unknown vertebrates, carrion, fruit, nectar, seeds and other plant
#'     material (fuzzy-coded, sums to 100), EltonTraits.}
#'   \item{beak_culmen, beak_width, beak_depth}{Beak: culmen length, beak
#'     width and depth (mm), AVONET; validation space, not used by
#'     default.}
#'   \item{iucn}{IUCN Red List category (`LC`, `NT`, `VU`, `EN`, `CR`),
#'     retrieved from the IUCN Red List API v4 (package rredlist) in 2026
#'     for BirdLife names; not a trait.}
#'   \item{Eigen.1, Eigen.2, Eigen.3, Eigen.4, Eigen.5, Eigen.6, Eigen.7,
#'     Eigen.8, Eigen.9, Eigen.10}{First ten axes of a principal coordinates
#'     analysis (`stats::cmdscale()`) of the square-rooted cophenetic
#'     distances of the full bird phylogeny (10,985 tips), as used for
#'     imputation in the paper (Penone et al. 2014); predictors, not traits.}
#' }
#'
#' @source
#' * Tobias, J. A. et al. (2022) AVONET: morphological, ecological and
#'   geographical data for all birds. *Ecology Letters* 25: 581-597,
#'   \doi{10.1111/ele.13898} (CC BY 4.0).
#' * Wilman, H. et al. (2014) EltonTraits 1.0: species-level foraging
#'   attributes of the world's birds and mammals. *Ecology* 95: 2027,
#'   \doi{10.1890/13-1917.1} (CC0).
#' * Myhrvold, N. P. et al. (2015) An amniote life-history database to
#'   perform comparative analyses with birds, mammals, and reptiles.
#'   *Ecology* 96: 3109, \doi{10.1890/15-0846R.1} (CC0).
#' * IUCN (2026) The IUCN Red List of Threatened Species,
#'   <https://www.iucnredlist.org>. Categories are provided for
#'   demonstration; use the current Red List for research.
#'
#' @seealso [birds_traits_all] (all 8,413 species of the paper),
#'   [birds_trait_info], [birds_tree]
#' @examples
#' data(birds_traits)
#' table(birds_traits$iucn)
#' round(colMeans(is.na(birds_traits[, birds_trait_info$trait])), 2)
"birds_traits"

#' Full example dataset: the 8,413 bird species of Toussaint et al. (2026)
#'
#' Same columns as [birds_traits] (see that page for the description of
#' variables and sources) for all 8,413 species analysed by Toussaint et al.
#' (2026, *Ecography*, \doi{10.1002/ecog.08578}), with observed
#' (non-imputed) trait values: 22 % of trait values are missing (all in
#' reproduction traits; 632 species have none of the eight reproduction
#' traits), as reported in the paper (23 %). Use it with [birds_trait_info]
#' to reproduce the published analysis: missing values imputed with
#' [impute_traits()] using `Eigen.1`-`Eigen.10` as predictors, then
#' [build_spaces()], [compare_spaces()] and [species_metrics_all()]. The
#' script `system.file("scripts", "ecography_8413_species.R", package =
#' "DisentFunSpaces")` runs the full analysis (about 4 GB of RAM and several
#' minutes), so this dataset is not used in examples or in the vignette.
#'
#' @format A data frame with 8413 rows and 40 columns (see [birds_traits]).
#'   BirdLife taxonomy (AVONET): 36 orders, 238 families. IUCN categories:
#'   6785 LC, 686 NT, 503 VU, 284 EN, 155 CR.
#' @source See [birds_traits].
#' @seealso [birds_traits] (2000 species),
#'   [birds_trait_info].
#' @examples
#' data(birds_traits_all)
#' dim(birds_traits_all)
#' round(colMeans(is.na(birds_traits_all[, birds_trait_info$trait])), 2)
"birds_traits_all"

#' Phylogeny of the `birds_traits` species
#'
#' Time-calibrated phylogeny (class `phylo`, branch lengths in million
#' years) of the 2000 species of [birds_traits], pruned from the tree used
#' in Toussaint et al. (2026): global bird phylogeny of Jetz et al. (2012)
#' accessed through the megatrees R package (Li 2024; `bird_20221117`),
#' made ultrametric with `phytools::force.ultrametric(method = "extend")`,
#' tip names matched to the trait data through the GBIF backbone, and
#' species still absent added at the root of their genus
#' (`phytools::add.species.to.genus()`).
#'
#' @format An object of class `phylo` with 2000 tips.
#' @source Jetz, W. et al. (2012) The global diversity of birds in space and
#'   time. *Nature* 491: 444-448, \doi{10.1038/nature11631}. Li, D. (2024)
#'   megatrees: mega-trees of plants and animals (R package),
#'   <https://github.com/megatrees>.
"birds_tree"

#' Trait-information table for `birds_traits`
#'
#' Describes the 26 trait columns of [birds_traits] as in Toussaint et al.
#' (2026): their `type` (`Q` or `F`), functional `dimension` (Locomotion,
#' Reproduction, Diet; Beak for the validation space), `block` (`diet` for
#' the fuzzy-coded diet), `transform` (`log10` for all continuous traits),
#' `weight`, `levels`, `circ_min`/`circ_max` and `use` (`FALSE` for the beak
#' traits: set it to `TRUE` to build the beak space). See
#' [check_trait_info()].
#' @format Data frame with 26 rows and 10 columns.
#' @examples
#' data(birds_trait_info)
#' table(birds_trait_info$dimension, birds_trait_info$use)
"birds_trait_info"
