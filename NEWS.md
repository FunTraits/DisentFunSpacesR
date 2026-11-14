# DisentFunSpaces 0.8.1

## Changes requested by CRAN

* DESCRIPTION: references written as `authors (year) <doi:...>`.
* `compare_nonnative()`: the returned value (structure and meaning) is
  documented in full.
* `write_analysis_project()` no longer changes the working directory: the
  generated project is run in a separate R process with
  `source(chdir = TRUE)`, and the archive is written with `zip::zip(root = )`
  ('zip' moved from Suggests to Imports).
* Generated projects install nothing: `utils/library.R` stops with the list of
  missing packages, `scripts/setup_environment.R` only lists them, and
  `run_all.R` neither changes the working directory nor calls
  `renv::restore()`.
* At most 2 cores are used when `_R_CHECK_LIMIT_CORES_` is set (Procrustes
  permutations, default of the application, `data.table::fread()`).

# DisentFunSpaces 0.8.0

## Traits shared by several functional dimensions

* A trait can now be assigned to several functional dimensions: its
  `dimension` in the trait-information table lists them separated by `;`
  (e.g. `"Diet;Locomotion"`). The trait enters the space of each of its
  dimensions and only once a space combining them (`weighting =
  "variable"`); with `weighting = "dimension"` / `"dimension_scaled"` it
  contributes to the within-dimension mean of each of its dimensions. The
  columns of a fuzzy-coded or multichoice block must share the same set of
  dimensions.
* New `trait_dimensions()` (individual dimensions of a set of traits) and
  `dimension_traits()` (one row per trait x dimension, with a `shared`
  flag).
* In the application (tab 2), dimensions are separated by `;` in the
  'dimension' field; the bulk assignment can replace, add or remove
  dimensions for the selected traits, and the validation reports the number
  of traits shared by several dimensions.

## Analysis project

* The folder of the exported analysis project holding `Methods_Results.Rmd`
  (and its rendered `.docx`) is renamed `report/` (previously `manuscript/`).

## Repeated species and group centroids

* A species may now appear in several rows of the trait table (individuals,
  populations, sexes, sites...): each row is a record placed in the spaces
  with its own identifier. New `record_ids()` numbers the records of repeated
  species (`"Salmo trutta_1"`, `"Salmo trutta_2"`...); `prepare_traits()`
  gains `taxon` (species of each record, used to match the phylogeny in
  `impute_traits()`); `community_from_long()` gains `records` (a species of
  an assemblage is given all its records). The application detects repeated
  species, numbers the records (groups, IUCN categories and imputation
  predictors are read per row) and exports scripts and projects that do the
  same.
* New `group_centroids()`: coordinates of the centroid of each group in one
  or several spaces. `plot_groups()` draws the centroids (`centroids = TRUE`,
  diamonds; with a reference group, the reference centroid too). In the
  application (tab 7), a 'Centroids' check box and a 'Centroids' table
  (also in the .zip of results and as `TableS14_group_centroids` in the
  exported project).

## Interactive figures in the application

* All the figures of the Shiny application are interactive 'plotly' widgets
  (when 'plotly', now suggested, is installed; static plots otherwise): hover
  shows the species name and its group, IUCN category, residual or value;
  zoom and pan; a click on a legend entry hides that group; the camera button
  exports an SVG. The PDF downloads are unchanged.
* New `as_plotly()` converts any `plot_*()` figure of the package (keeps the
  arrow heads of trait correlations, the subtitle and the equal scaling of
  the axes).
* The matrix of Procrustes correlations (tab 5, `plot_comparison()`) is drawn
  directly with plotly instead of `ggplotly()`, whose conversion of
  `geom_tile()` failed in the application ('x is null'); hover gives r, p and
  the number of shared species.
* No more warnings when the application starts: `page_navbar()` uses
  `navbar_options()` (bslib >= 0.9.0 now required) and the layout of the
  'Assemblages' tab declares one column width per element.

## Species groups (generalises the non-native analysis)

* New functions for groups of species defined by any categorical attribute
  that is not a trait (origin, human uses, migratory behaviour, activity
  period, threat status, taxonomy...): `groups_all()`, `group_overlap()`,
  `group_species()`, `group_contribution()`, `summary_groups()`,
  `plot_groups()`, `plot_group_overlap()`. Groups may be exclusive or not
  (several labels per species, `sep`), labels can be merged (`recode`) and
  selected (`levels`).
* Group level: hull volume (FRic), FDis, FOri, MPD, share of the hull not
  shared with the other species, segregation of nearest neighbours; pairwise
  hull overlap (intersection / union, with `geometry::intersectn()`),
  functional Jaccard dissimilarity and its turnover / nestedness components
  (Villeger et al. 2013), proportion of species of one group inside the hull
  of the other, distance between centroids; SES against permutations of the
  labels; PERMANOVA and PERMDISP computed from the Euclidean coordinates
  (same values as `vegan::adonis2()` / `betadisper()`, checked in the tests,
  permutations in O(n p)).
* Species level, relative to a reference group: the metrics of
  `nonnative_species()` / `nonnative_assemblage()` for each focal group
  (`FNov`, `FSp_ref`, `inside_ref_hull`, `dFUn_rel`, contribution of the
  focal group), with optional weights (`FUSW`). `nonnative_*()` are
  unchanged and give the same values with `reference = "native"`.

* `groups_from_columns()`: groups from indicator (0 / 1) columns, one group
  per column (e.g. one column per human use: fisheries, aquaculture,
  aquarium, game fish), species in several groups, and a group (`none`) for
  species set to 0 in all columns, the natural reference. App (tab 7): new
  source "Indicator columns (0 / 1, one group per column)", with detection
  of a column that is the union of the others ("All uses"); the exported
  script, analysis project and Methods text follow. Tab 1 proposes common
  use columns as non-trait variables, and warns that undeclared 0 / 1
  columns are used as binary traits.

## Densities, extinction scenarios, continuous variables, assemblages

* `group_overlap()` / `groups_all()`: `density = TRUE` adds the overlap of
  kernel densities of the groups (`overlap_kde`, sum of the minimum of the
  two densities, as the TPD overlap of Carmona et al. 2016; robust to outlying
  species, unlike hulls) and `FRic_kde`, on a grid of 2-4 axes with one
  bandwidth for all groups of a space, optionally weighted by abundance
  (`abundance`), with SES under the null model. `plot_groups(contours = TRUE)`
  draws density contours.
* `extinction_scenarios()`, `plot_extinction()`: species removed by
  decreasing extinction probability (IUCN categories -> EDGE2 probabilities,
  or any probability / priority such as FUSE) vs at random, up to the number
  of threatened species, with FRic, FDis, FOri and MPD relative to all
  species; stochastic extinctions drawn with the probabilities. Which
  dimension loses most, and more than expected from species richness?
* `covariate_space()`, `plot_covariate()`: a continuous variable that is not
  a trait (range size, abundance, trend, economic value...) in every space:
  linear fit on the first axes (vector fitting) with permutation test,
  functional autocorrelation (correlation with the mean of the k nearest
  neighbours, also for non-linear structure), Spearman correlations with FUn
  and FSp.
* `group_assemblage()`, `community_from_long()`, `plot_group_sites()`,
  `groups_all(comm =)`: contribution of focal groups within each assemblage,
  with groups that may differ between assemblages (native in one basin,
  introduced in another); same values as `nonnative_assemblage()` for
  non-native species (checked in the tests).

## Shiny application

* Tab 6: "Extinction scenarios" card (needs the IUCN column).
* Tab 7: three modes (categories, continuous variable, simulated status);
  kernel-density overlap with optional abundance weights; density contours
  on the plot; assemblages from a long table (assemblage, species, optional
  group per assemblage) with an "Assemblages" result tab.
* The exported project archives the assemblages (`data/raw/assemblages.csv`)
  and runs the new analyses (Table8-9, TableS10-S13, Fig6-8), described in
  the Methods paragraph of the report.
* Tab 7 "Non-native species" becomes "Species groups": grouping variables
  are declared in tab 1 (several columns, not used as traits; the IUCN column
  is also offered), one is analysed at a time, with an optional separator of
  multiple labels, a box to merge labels (`threatened = VU, EN, CR`), the
  groups analysed and an optional reference group. Results: separation of
  the groups in each space, focal groups relative to the reference, plot of
  the groups (or of the focal species against the reference hull), heat map
  of pairwise overlaps across spaces, and tables (groups, pairs, species
  mixing, focal and reference species, contribution). The simulated
  native / non-native status remains available for exploration.
* Exported project: `06_ANALYSIS_groups.R` replaces `06_ANALYSIS_nonnative.R`
  (step `"groups"`; `"nonnative"` and the 0.7.0 settings `status_col` /
  `nonnative` are still accepted and translated); new tables Table5-7 and
  TableS6-S9, figures Fig4_groups_* and Fig5_group_overlap, and the Methods
  paragraph of the report.
* Data tab: the separator of text files is now detected automatically
  ("Auto-detect", default; comma, semicolon, tab or pipe) together with the
  decimal mark (decimal comma of European CSV files). A file read as a
  single column with the separator chosen by hand is read again with the
  detected one; the preview reports the separator used, which is also
  written into the exported script.

# DisentFunSpaces 0.7.0

## Shiny application

* New look, shared with the other FunTraits applications (intraitR,
  Rfishmorph): white header and cards on a pale page, IBM Plex Sans shipped
  with the package (works offline), numbered steps in the header that turn
  green once done, and a "Help" menu.
* New home page: what the framework does, the eight steps of the workflow
  (click to open one; done / next / optional), a summary of the session and
  the citation.
* Export tab rebuilt: it now writes a **reproducible project** (see
  `write_analysis_project()`), runs its scripts in a background R process
  (package 'callr', when installed: the app stays usable) and offers one
  `.zip`. The tab previews the project before it is built (architecture,
  README, `parameters.R`) and shows the run log and the Methods & Results
  report once it has run. The single files remain available.

## Functions

* `write_analysis_project()`: writes an analysis as a research repository
  (`scripts/` numbered by step, `utils/parameters.R` with every setting,
  `data/raw`, `data/processed`, `results/figures`, `results/tables`,
  `manuscript/Methods_Results.Rmd`, `run_all.R`, `.Rproj`, `README.md`,
  `data/metadata.csv`, `session_info.txt`), runs it and zips it.

# DisentFunSpaces 0.6.0

Faster computations and a Shiny application for large datasets. All the
changes below give the same results as before (checked against the previous
implementations in `tests/testthat/test-performance.R`), except where stated.

## Shiny application

* New "Performance" panel (tab 4). By default ("auto"), from 1000 species,
  distances are computed with the fast engine (`dist.ktab()` exceeded 6 GB
  of RAM for the 8,413 bird species) and the PCoA with a partial
  eigendecomposition of all positive axes (same results as the published
  settings). From 3000 species (`options(DisentFunSpaces.app_large_n)`) the
  app shows the memory needed and advice.
* Procrustes permutations can use several CPU cores (macOS / Linux) with
  identical results.
* Text files are read with `data.table::fread()` when installed; upload limit
  raised to 1 GB (`run_app(max_upload_mb =)`).
* Plots are cached (`shiny::bindCache()`): changing tab or going back to a
  space does not redraw it; the arrow threshold slider is debounced.
* The exported R script reproduces the engine, partial PCoA and core settings.
* Hint when R uses its reference BLAS (an optimised BLAS speeds up PCoA and
  Procrustes tests several times).

## Functions

* `build_space()` / `build_spaces()`: new argument `eigen = "partial"`:
  only the axes kept are computed (Lanczos algorithm, package 'RSpectra', in
  Suggests): all positive-eigenvalue axes by default (falls back to the full
  decomposition beyond 500 positive axes, `options(DisentFunSpaces.partial_k_max)`),
  or the first `max_axes`. Same coordinates (up to sign) and relative
  inertia; broken-stick values and negative-eigenvalue diagnostics are not
  available. Default unchanged (`"full"`). Combined space of the 8,413 bird
  species: 5 s instead of 11 min.
* `PCoAGraph()`: new arguments `partial` and `n_axes` (partial eigendecomposition).
* `CorrelType()`: vectorised over axes (one matrix operation per trait instead
  of one `cor()` / `lm()` call per trait x axis).
* `trait_distance(engine = "fast")`: about 1.5x faster and less memory
  (fewer temporary vectors of length n(n-1)/2); values unchanged (identical).
* `compare_spaces(method = "ade4")`: exact re-implementation of
  `ade4::procuste.rtest()` computing singular values only (same permutations,
  same p-values); new argument `cores`.
* `nonnative_species()` / `nonnative_assemblage()`: the inside/outside test
  uses the facets of a single convex hull instead of one hull per non-native
  species; FOri and MPD no longer build an n x n matrix; with the
  `"permute"` null model the distances of the current assemblage are
  computed once and MPD / FOri of each replicate are updated from them; with
  the `"pool"` null model the historical assemblage is computed once. One
  space of 8,413 species, 9 replicates: 4 s instead of 225 s.

# DisentFunSpaces 0.5.1

* `guess_trait_info()`: columns sharing a prefix are proposed as a
  fuzzy-coded block only if their rows sum to a near-constant (as
  documented); e.g. `body_mass` and `body_length` are no longer grouped.

* The simulated example data (`toy_traits`, `toy_trait_info`, `toy_tree`)
  are replaced by real data: `birds_traits`, 2000 bird species of
  Toussaint et al. (2026, Ecography) (all 281 complete species plus a random
  sample of species with at least one trait informed in each dimension,
  adult body mass required so that all pairs of species share traits),
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
  published analysis; also available in the Shiny application (tab 1).
* Examples, vignette, README, Shiny application ("Example: 2000 bird
  species") and the worked example of `inst/paper/paper_analyses.R`
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
