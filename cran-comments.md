## Submission

This is a new submission.

DisentFunSpaces generalises the framework of Toussaint et al. (2026,
Ecography, doi:10.1002/ecog.08578) to any mixed-type trait dataset: one
functional space per functional dimension (generalised Gower distance +
PCoA), Procrustes comparisons of spaces, species-level metrics, non-native
species, and a 'shiny' application.

## Test environments

* local: Ubuntu 24.04, R 4.3.3 (R CMD check --as-cran): 0 errors, 0 warnings
  (once LaTeX and missForest are available)
<!-- TODO before submitting: run, then fill in the results -->
* local: macOS (aarch64), R-release: devtools::check(args = "--as-cran")
* win-builder, R-release: R 4.6.1 (2026-06-24 ucrt), Windows Server 2022 x64:
  Status: 1 NOTE (new submission; possibly misspelled words are proper names)
* win-builder, R-devel: R Under development (unstable) (2026-10-28 r90783 ucrt),
  Windows Server 2022 x64: Status: 1 NOTE (new submission; possibly misspelled words)
* macOS builder: devtools::check_mac_release()

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
