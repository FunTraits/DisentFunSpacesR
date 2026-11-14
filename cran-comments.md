## Resubmission

This is a resubmission. In this version I have:

* Written the references in the Description field as
  `authors (year) <doi:...>`, with years in parentheses.
* Added a complete \value section to compare_nonnative.Rd (class, structure
  and meaning of the output).
* Removed every change of the working directory: R/project.R runs the
  generated project in a separate R process with `source(chdir = TRUE)` and
  writes the archive with `zip::zip(root = )`; the template
  inst/templates/analysis_project/run_all.R no longer calls setwd().
* Removed all package installation from the package code and from the
  project templates (no install.packages(), remotes::install_github() or
  renv::restore() calls; missing packages are only reported).
* Removed installed.packages() (replaced by requireNamespace()).
* Limited the number of cores to 2 when _R_CHECK_LIMIT_CORES_ is set
  (parallel Procrustes permutations, default number of cores of the shiny
  application and data.table::fread() threads).

DisentFunSpaces generalises the framework of Toussaint et al. (2026,
Ecography, doi:10.1002/ecog.08578) to any mixed-type trait dataset: one
functional space per functional dimension (generalised Gower distance +
PCoA), Procrustes comparisons of spaces, species-level metrics, non-native
species, and a 'shiny' application.

## Test environments

* local: macOS (aarch64), R 4.5.2, devtools::check(args = "--as-cran"):
  0 errors | 0 warnings | 0 notes
* R-hub v2 (GitHub Actions): windows (R-devel), macos-arm64 (R-devel),
  linux (R-devel): 0 errors | 0 warnings
* GitHub Actions R-CMD-check: macOS, Windows and Ubuntu (R-release): OK

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Aurèle Toussaint <aurele.toussaint@cnrs.fr>'
  New submission.

  Possibly misspelled words in DESCRIPTION are author names (Gower, Pavoine,
  Tobias, Wilman, Myhrvold, Toussaint; "et al."), names of databases and of a
  journal (AVONET, EltonTraits, Amniote, Ecography), metric abbreviations
  (FUn: functional uniqueness, FSp: functional specialisation) and the trait
  type "multichoice" (several categories recorded per species).

## Notes for CRAN

* The package writes files only where the user asks: `write_analysis_project()`
  has no default `dir` and refuses a non-empty folder; its example writes in
  `tempdir()` and does not run the generated project (`run = FALSE`).
* Tests that launch a separate R process (running a generated project) are
  skipped on CRAN (`skip_on_cran()`).
* `inst/app/www/fonts` holds subsets of the IBM Plex fonts (SIL Open Font
  License 1.1, licence text included) so that the shiny application displays
  correctly offline; see `inst/COPYRIGHTS` and the `Copyright` field.
* Example data are derived from published, openly available datasets (AVONET,
  EltonTraits, Amniote database, bird phylogeny of Jetz et al. 2012), with
  their sources documented in the help pages and `inst/COPYRIGHTS`.
