# DisentFunSpaces

[![R-CMD-check](https://github.com/FunTraits/DisentFunSpacesR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/FunTraits/DisentFunSpacesR/actions/workflows/R-CMD-check.yaml)

**Build and compare functional spaces from any mixed-type trait dataset.**

R package and Shiny application generalising the framework of

> Toussaint A., Tedesco P., Grenouillet G., Kasari-Toussaint L., Brosse S. (2026).
> Disentangling functional spaces: a multidimensional framework for ecology and
> conservation. *Ecography*, e08578. https://doi.org/10.1002/ecog.08578

The published scripts ([FunTraits/DisentFunSpaces](https://github.com/FunTraits/DisentFunSpaces))
were written for >10,000 bird species and three trait dimensions (locomotion,
diet, reproduction). This package makes the same workflow available for any
taxon, any number of functional dimensions and any trait type.

## Installation

```r
# install.packages("remotes")
remotes::install_github("FunTraits/DisentFunSpacesR", build_vignettes = TRUE)
install.packages(c("missForest", "readxl"))               # optional: imputation, Excel input
```

## Workflow

| Step | Function | Published analysis (`utils/functions.R`, scripts) |
|---|---|---|
| Describe traits (type, dimension, block) | `guess_trait_info()`, `check_trait_info()` | hard-coded trait vectors |
| Prepare traits (transformations, fuzzy rescaling) | `prepare_traits()` | `log10`, `prep.fuzzy()` |
| Missing values | `impute_traits()` (missForest; predictors: columns of the data table such as `Eigen.1`… via `covariates`, and/or eigenvectors computed from a phylogeny via `phylo`), `drop_incomplete()` | `missForest()` with 10 phylogenetic eigenvectors |
| Distance + PCoA (core) | `PCoAGraph()` (same arguments and output as in the paper) | `PCoAGraph()`: `ktab.list.df()` → `dist.ktab(scann = FALSE, "scaledBYrange")` → `ape::pcoa()` → `envfit()` + `CorrelType()` |
| Functional spaces | `build_space()`, `build_spaces()`: build the trait tables (one per dimension × type) from the trait-information table and call `PCoAGraph()` | `PCoAGraph(list(LHT, Morpho, Diet), groups = c("Q","Q","F"))` |
| Distance only | `trait_distance()` | `ade4::dist.ktab()` |
| Space comparison | `compare_spaces()`, `procrustes_residuals()` | `run_procrustes()` (`ade4::procuste.rtest`) |
| Species metrics | `species_metrics_all()`, `compare_metrics()` | `calc_fuse_noINS()` (FUn, FSp, FUSE) |
| Non-native species | `nonnative_species()`, `nonnative_assemblage()`, `nonnative_all()`, `summary_nonnative()`, `compare_nonnative()`, `simulate_introductions()` | Toussaint et al. 2018 (*Ecol. Lett.*) approach, per functional space |
| Plots | `plot_space()`, `plot_comparison()`, `plot_procrustes()`, `plot_quality()`, `plot_eigen()` | |

### Trait types

| Code | Type | Columns | Per-trait distance (`ade4::dist.ktab` defaults) |
|---|---|---|---|
| `Q` | continuous | 1 | \|xi − xj\| / range |
| `O` | ordinal | 1 (numeric ranks, or text + `levels`) | rank-based |
| `N` | categorical | 1 | simple matching |
| `D` | binary | 1 (0/1) | Jaccard |
| `B` | multichoice | ≥ 2 (0/1, same `block`) | Jaccard |
| `F` | fuzzy-coded | ≥ 2 (≥ 0, same `block`), rows rescaled to proportions | Manly overlap |
| `C` | circular | 1 (`circ_min`/`circ_max`) | angular |

Per-trait distances are combined as d = √(Σ w·d² / Σ w) over the traits
available for both species (Pavoine et al. 2009). Two weightings:
`"variable"` (each trait, or each fuzzy/multichoice block, weighs 1 — identical
to the published analysis) and `"dimension"` (each functional dimension weighs
the same in combined spaces). Per-trait weights can also be set in the
`weight` column.

Each `fs_space` keeps the `PCoAGraph()` elements `PCoA`, `varPlot` and
`PCoACor`, so objects can be passed to the published scripts
(e.g. `run_procrustes()` uses `PCA_list[[j]]$PCoA$vectors`).

## Quick start

With the example data:

```r
library(DisentFunSpaces)
traits <- prepare_traits(birds_traits, birds_trait_info, species_col = "species")
preds <- birds_traits[, paste0("Eigen.", 1:10)]; rownames(preds) <- birds_traits$species
traits <- impute_traits(traits, covariates = preds)     # life-history gaps (missForest)
spaces <- build_spaces(traits, engine = "fast")         # Locomotion, Reproduction, Diet, Combined
compare_spaces(spaces, nrepet = 199)
m <- species_metrics_all(spaces, iucn = setNames(birds_traits$iucn, birds_traits$species))
```

With your own data:

```r
library(DisentFunSpaces)
traits <- prepare_traits(my_traits, my_trait_info, species_col = "species")
preds  <- my_traits[, paste0("Eigen.", 1:10)]           # predictors (not traits)
rownames(preds) <- my_traits$species
traits <- impute_traits(traits, covariates = preds)     # and/or phylo = my_tree
spaces <- build_spaces(traits)                          # each dimension + Combined
cmp    <- compare_spaces(spaces, nrepet = 999)          # Procrustes tests
plot_comparison(cmp)
metrics <- species_metrics_all(spaces, iucn = my_iucn)  # FUn, FSp, FUSE
nn <- nonnative_all(spaces, my_status, null = "permute") # non-native species
summary_nonnative(nn)
```

A trait-information template can be generated with
`write_trait_template(my_traits, "trait_info.csv", species_col = "species")`.

## Shiny application

```r
DisentFunSpaces::run_app()
```

Tabs: 1. Data (CSV/TSV/Excel or the bird example) — 2. Traits (type,
dimension, block, transformation, weight, ordinal levels, circular limits) —
3. Missing values — 4. Functional spaces (per dimension, combined, custom
combinations) — 5. Procrustes comparisons and species residuals — 6. Species
metrics — 7. Non-native species — 8. Export (tables, R objects, trait-information CSV and an R script
reproducing the analysis).

### Large datasets

From 1000 species the app uses large-data settings by default ("Performance"
panel of tab 4), which give the same results as the published analysis:

* distances computed with the fast engine (same values as
  `ade4::dist.ktab()`, which needs more than 6 GB of RAM for 8,413 species);
* partial PCoA (package `RSpectra`): only the positive-eigenvalue axes (or
  the "Max. PCoA axes" set) are computed, without the full
  eigendecomposition; same coordinates (up to sign) and inertia; negative
  eigenvalues are not reported;
* Procrustes permutations on several CPU cores (macOS / Linux, when many
  axes are compared), with identical p-values;
* text files read with `data.table::fread()`; uploads up to 1 GB.

Timings for the 8,413 bird species (2-core virtual machine, reference BLAS):

| Step | Before (0.5.1) | Now (0.6.0) |
|---|---|---|
| Distances, `dist.ktab()` (app default before) | out of memory (> 6 GB) | — |
| Distances, fast engine, Combined space | 41 s | 32 s |
| PCoA, Combined space (23 positive axes) | 11 min | 5 s |
| Non-native analysis, 1 space, 9 null replicates | 225 s | 4 s |

```r
install.packages(c("RSpectra", "data.table"))   # recommended for large data
DisentFunSpaces::run_app()
```

In scripts: `build_spaces(traits, engine = "fast", eigen = "partial")`.

PCoA and Procrustes tests rely on matrix algebra: an optimised BLAS
(Accelerate on macOS, see the *R for macOS FAQ*; OpenBLAS on Linux) makes
them several times faster. The app shows a hint when R uses its reference
BLAS.

## Notes

* Default settings reproduce the published analysis (all positive PCoA axes,
  no correction for negative eigenvalues, `ade4::procuste.rtest`). Options:
  first *k* axes, Cailliez/Lingoes corrections, symmetric `vegan::protest`,
  Mantel tests.
* `weighting = "dimension_scaled"` rescales each dimension before combining:
  in simulations it was the only weighting that balanced dimensions of
  different trait types (categorical distances otherwise dominate).
* `simulate_traits()` generates mixed-type traits with known congruence
  between dimensions, to test the workflow or plan a study.
* `trait_distance(engine = "fast")` gives the same distances as
  `dist.ktab()` (checked in the tests) with vectorised matrix algebra, for
  very large datasets.
* `space_quality()` reports the mean absolute deviation between functional
  distances and distances in the first *k* axes (Maire et al. 2015) to choose
  the dimensionality of spaces.
* Non-native species: `nonnative_species()` gives, per non-native species,
  the distance to the nearest native (`FNov`), to the native centroid
  (`FSp_native`), the position inside/outside the native convex hull and,
  with EICAT categories, `FUSI = log(1 + W·FNov + W·FSp_native)` (built as
  FUSE; default EICAT weights are ranks, not empirical estimates); per native
  species, the loss of uniqueness once non-natives are added (`dFUn_rel`).
  `nonnative_assemblage()` compares historical and current assemblages (S,
  FRic, FDis, FOri, MPD) with `"permute"` or `"pool"` null models.
* Example data: `birds_traits`, 2000 bird species of Toussaint et al.
  (2026) (all 281 complete species plus a random sample of species with at
  least one trait informed in each dimension, including adult body mass),
  described by the traits of that
  paper — locomotion (tarsus length, wing length, Kipp's distance, tail
  length, body mass; AVONET, Tobias et al. 2022, CC BY 4.0), reproduction
  (clutch size, clutches per year, adult body mass, incubation time,
  longevity, fledging age, egg mass, body length; Amniote database,
  Myhrvold et al. 2015, CC0), diet (10 categories, fuzzy-coded; EltonTraits
  1.0, Wilman et al. 2014, CC0) and the beak traits of the validation space
  (`use = FALSE`) — with real missing values, IUCN categories and 10
  phylogenetic eigenvectors; `birds_trait_info` and `birds_tree` (see
  `?birds_traits` and `data-raw/birds_traits.R`). `birds_traits_all`: the
  same traits for all 8,413 species of the paper (22 % missing values); the
  full analysis is run by
  `system.file("scripts", "ecography_8413_species.R", package = "DisentFunSpaces")`. Simulated data with
  known congruence: `simulate_traits()`.
