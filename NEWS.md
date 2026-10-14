# DisentFunSpaces 0.5.1

* `guess_trait_info()`: columns sharing a prefix are proposed as a
  fuzzy-coded block only if their rows sum to a near-constant (as
  documented); e.g. `body_mass` and `body_length` are no longer grouped.

* The simulated example data (`toy_traits`, `toy_trait_info`, `toy_tree`)
  are replaced by real data: `birds_traits`, 1000 bird species of
  Toussaint et al. (2026, Ecography) with the fewest missing values (all
  281 complete species plus species missing one or two traits),
  described by the traits of that paper (locomotion: 5 continuous traits,
  AVONET; reproduction: 8 continuous traits, Amniote database; diet: 1
  fuzzy-coded trait with 10 categories, EltonTraits; beak: 3 traits of the
  validation space, `use = FALSE`), with real missing values, IUCN
  categories and 10 phylogenetic eigenvectors; `birds_trait_info` and
  `birds_tree` (script `data-raw/birds_traits.R`).
* Large datasets: `trait_distance()` accumulates per-variable distances one
  at a time (memory no longer grows with the number of traits), and
  `PCoAGraph()` uses a memory-lean PCoA (same result as `ape::pcoa()`) for
  3000 species or more without correction (`options(DisentFunSpaces.lean_pcoa_n)`).
  The 8,413 bird species run in about 6 min with < 4 GB of RAM.
* `birds_traits_all`: the same traits for all 8,413 species of the paper,
  with the script `inst/scripts/ecography_8413_species.R` reproducing the
  published analysis.
* Examples, vignette, README, Shiny application ("Example (1000 bird
  species, real data)") and the worked example of `inst/paper/paper_analyses.R`
  use the bird data. Simulated data remain available with `simulate_traits()`
  and are used by the tests that need circular and multichoice traits.

# DisentFunSpaces 0.4.0

* New module for non-native species, the counterpart of FUn/FSp/FUSE for
  threatened species:
  * `nonnative_species()`: for each non-native species, distance to the
    nearest native (functional novelty, `FNov`), distance to the native
    centroid (`FSp_native`), position inside/outside the native convex hull
    and, with EICAT categories or other weights, `FUSI`
    (`log(1 + W x FNov + W x FSp_native)`, built as FUSE); for each native
    species, relative loss of functional uniqueness once non-native species
    are added (`dFUn_rel`).
  * `nonnative_assemblage()`: historical (natives) vs current (natives -
    extirpated + non-native) assemblages: S, FRic, FDis, FOri, MPD, with
    `"pool"` and `"permute"` null models (SES, p-values); per-assemblage
    status and extirpations supported.
  * `nonnative_all()`, `summary_nonnative()`, `compare_nonnative()`: the
    same across all spaces, to identify the functional dimensions most
    affected by introductions.
  * `plot_nonnative()`, `plot_nonnative_change()`, `simulate_introductions()`,
    `eicat_weights()`.
* Shiny app: status column in tab 1, new tab "7. Non-native species"
  (status from the data or simulated introductions); export renumbered 8.
* `analysis_script()` adds the non-native analysis when a status column is set.
* `geometry` added to Imports (convex hulls).

# DisentFunSpaces 0.3.0

* `simulate_traits()` and `default_sim_design()`: simulate traits of all seven
  types from latent niche axes with known congruence between dimensions
  (exact target correlations, phylogenetic signal, missing values).
* New weighting `weighting = "dimension_scaled"` in `trait_distance()`,
  `build_space()` and the Shiny app: each dimension is rescaled to a mean
  squared distance of 1 before combining, so that dimensions of categorical
  traits no longer dominate combined spaces.
* `inst/paper/paper_analyses.R` reproduces all results of the software paper.

# DisentFunSpaces 0.2.0

* `build_space()` is a wrapper around `PCoAGraph()` (published core function:
  `ktab.list.df()` -> `dist.ktab()` -> `ape::pcoa()` -> `envfit()` + `CorrelType()`);
  `fs_space` objects keep `PCoA`, `varPlot` and `PCoACor`.
* `build_space(engine = "fast")` and `trait_distance(engine = "fast")`: same
  distances as `dist.ktab()`, much lower memory use for thousands of species.
* `impute_traits(covariates = )`: columns of the data table (e.g. `Eigen.1`...)
  as imputation predictors, also selectable in the Shiny app.
* Multichoice blocks with no modality recorded are set to NA; imputed
  multichoice blocks keep at least one modality.
* Undefined distances (species pairs sharing no trait) stop with a clear
  message before the PCoA.
* `species_metrics()`: exact zeros for identical profiles, new `n_identical`;
  `compare_metrics()` handles ties in top sets.
* `plot_space()`: bubbles for discrete spaces.

# DisentFunSpaces 0.1.0

* First version.
