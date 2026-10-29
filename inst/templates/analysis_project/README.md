# {{TITLE}}

This repository contains the code, data, figures and result tables of a
functional-space analysis of {{N_SPECIES}} species described by {{N_TRAITS}}
traits in {{N_DIMS}} functional dimension(s) ({{DIMENSIONS}}), produced with the
R package [DisentFunSpaces](https://github.com/FunTraits/DisentFunSpacesR)
({{DFS_VERSION}}) on {{DATE}}.

**Author(s):** {{AUTHOR}}

------------------------------------------------------------------------

## 🗂️ Repository structure

``` text
{{PROJECT}}/
├── scripts/                  # Analysis pipeline (sourced in numeric order by run_all.R)
│   ├── 00_START_GeneralScript.R          # Packages, parameters and helpers (sourced by every script)
{{SCRIPT_TREE}}
│   └── setup_environment.R               # Installs the packages (not run by run_all.R)
│
├── utils/                    # Shared code
│   ├── library.R                         # Packages
│   ├── parameters.R                      # Every analysis parameter (one place)
│   └── functions.R                       # Helpers (tables, figures, data reading)
│
├── data/
│   ├── raw/                  # Trait table and trait description as analysed{{PHYLO_TREE}}
│   ├── processed/            # Objects passed between scripts (.rds), written by the pipeline
│   └── metadata.csv          # Definition of every variable of the data and result tables
│
├── results/
│   ├── figures/              # Main (Fig*) and supplementary (FigS*) figures, .png and .pdf
│   └── tables/               # Main (Table*) and supplementary (TableS*) tables, .csv
│
├── manuscript/
│   ├── Methods_Results.Rmd   # Methods and Results, numbers read from results/tables
│   └── Methods_Results.docx  # The same, rendered by run_all.R
│
├── run_all.R                 # Runs the full pipeline, then renders the report
├── {{PROJECT}}.Rproj         # RStudio project file
├── session_info.txt          # R and package versions that produced the results
├── README.md                 # This file
└── .gitignore
```

> **Note:** `data/processed/` and the rendered report are listed in
> `.gitignore`: they are rebuilt by `run_all.R`. Archive `data/raw/`
> with the code (e.g. Zenodo) if the trait data can be shared.

------------------------------------------------------------------------

## 📦 Data

`data/raw/traits.csv` is the trait table analysed (source: *{{DATA_SOURCE}}*),
one row per species (column `{{SPECIES_COL}}`) and one column per trait, or per
modality for fuzzy-coded and multichoice traits.{{EXTRA_COLS}}

`data/raw/trait_info.csv` describes each trait: `type` (Q continuous, O
ordinal, N nominal, D binary, B multichoice, F fuzzy-coded, C circular),
functional `dimension`, `block` (columns of one fuzzy-coded or multichoice
trait), `transform`, `weight`, `levels` (ordinal traits stored as text),
`circ_min` / `circ_max` (circular traits) and `use`.{{PHYLO_SENTENCE}}

`data/metadata.csv` defines every variable of the data and result tables.

------------------------------------------------------------------------

## 🧪 Reproducing the analysis

1.  Open `{{PROJECT}}.Rproj` in RStudio.
2.  Install the packages once: `source("scripts/setup_environment.R")`
    (or let `utils/library.R` install what is missing).
3.  Run the full pipeline:

    ``` r
    source("run_all.R")
    ```

Every script reads its inputs from `data/` and writes its outputs to
`data/processed/` and `results/`, so each can also be run on its own once the
previous ones have been run. All parameters (missing values, spaces,
weighting, PCoA, permutations, number of neighbours, random seed...) are set
in **`utils/parameters.R`**: change them there, not in the scripts.

| Step | Script | Main outputs |
|------|--------|--------------|
{{STEP_TABLE}}

------------------------------------------------------------------------

## ⚙️ Settings used

{{SETTINGS}}

------------------------------------------------------------------------

## 📚 Citation

If you use these results or the code, please cite:

> Toussaint A., Tedesco P., Grenouillet G., Kasari-Toussaint L., Brosse S.
> (2026) Disentangling functional spaces: a multidimensional framework for
> ecology and conservation. *Ecography* e08578.
> <https://doi.org/10.1002/ecog.08578>

and the DisentFunSpaces package ({{DFS_VERSION}}),
<https://github.com/FunTraits/DisentFunSpacesR>.
