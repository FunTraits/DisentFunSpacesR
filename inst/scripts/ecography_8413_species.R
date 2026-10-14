## Full example: the 8,413 bird species of Toussaint et al. (2026, Ecography) --
## Reproduces the workflow of the published analysis with the package:
## traits of `birds_traits_all` (observed values), missForest imputation with
## the first ten phylogenetic eigenvectors, locomotion / reproduction / diet /
## combined spaces, Procrustes tests (all axes and first two axes), FUn, FSp
## (k = 5, first two axes as in the paper) and FUSE, overlap of the 10 % most
## distinct species across spaces.
## Run: source(system.file("scripts", "ecography_8413_species.R", package = "DisentFunSpaces"))
## Requirements: missForest; about 4 GB of RAM; about 6 min on a 2-core Linux
## machine with OpenBLAS (longer with the reference BLAS).
## Outputs are written to ./ecography_outputs
suppressPackageStartupMessages(library(DisentFunSpaces))
out <- "ecography_outputs"; dir.create(out, showWarnings = FALSE)
SEED <- 2026
tic <- function() proc.time()[["elapsed"]]
t0 <- tic()

## 1. traits and imputation ------------------------------------------------------
data(birds_traits_all); data(birds_trait_info)
tr <- prepare_traits(birds_traits_all, birds_trait_info, species_col = "species")
cat("Missing trait values:", round(100 * mean(is.na(tr$data)), 1), "%\n")
preds <- birds_traits_all[, paste0("Eigen.", 1:10)]
rownames(preds) <- birds_traits_all$species
tr <- impute_traits(tr, covariates = preds, ntree = 100, seed = SEED)
print(tr$imputation$oob)
cat("imputation:", round(tic() - t0), "s\n")

## 2. spaces ---------------------------------------------------------------------
spaces <- build_spaces(tr, envfit_perm = 0, engine = "fast", keep_dist = FALSE)
print(spaces)
cat("spaces:", round(tic() - t0), "s\n")
cc <- round(spaces$Combined$PCoACor[, 1:2]^2, 2)          # R2 trait ~ axis
write.csv(cc, file.path(out, "combined_trait_axis_R2.csv"))
print(cc)

## 3. Procrustes tests -----------------------------------------------------------
pr_all <- compare_spaces(spaces, nrepet = 999)
pr_2 <- compare_spaces(spaces, n_axes = 2, nrepet = 999)
tab <- data.frame(pr_all$table[, c("space1", "space2", "n_axes", "r", "p_value")],
                  r_2axes = pr_2$table$r, p_2axes = pr_2$table$p_value)
write.csv(tab, file.path(out, "procrustes.csv"), row.names = FALSE)
print(tab, digits = 3)
cat("Procrustes:", round(tic() - t0), "s\n")

## 4. species-level metrics (first two axes, as in the paper) -------------------------
iucn <- stats::setNames(birds_traits_all$iucn, birds_traits_all$species)
m <- species_metrics_all(spaces, n_axes = 2, k = 5, iucn = iucn)
write.csv(m, file.path(out, "species_metrics_2axes.csv"), row.names = FALSE)
overlap <- function(metric, top = 0.1) {
  sets <- lapply(split(m, m$space), function(x) x$species[x[[metric]] >= stats::quantile(x[[metric]], 1 - top)])
  n_in <- table(unlist(sets))
  data.frame(metric = metric, space = names(sets), n_top = lengths(sets),
             exclusive = vapply(sets, function(v) sum(n_in[v] == 1), integer(1)),
             in_all_spaces = sum(n_in == length(sets)))
}
ov <- rbind(overlap("FUn"), overlap("FSp"))
write.csv(ov, file.path(out, "top10_overlap.csv"), row.names = FALSE)
print(ov)
print(compare_metrics(m, "FUn", top = 0.1)$pairs, digits = 3)
print(compare_metrics(m, "FSp", top = 0.1)$pairs, digits = 3)
topF <- do.call(rbind, lapply(split(m, m$space), function(x) head(x[order(-x$FUSE), c("space", "species", "iucn", "FUSE")], 10)))
write.csv(topF, file.path(out, "top10_FUSE.csv"), row.names = FALSE)
print(topF)
saveRDS(list(procrustes = tab, overlap = ov, metrics = m, oob = tr$imputation$oob,
             eigen = lapply(spaces, `[[`, "eigen")), file.path(out, "ecography_results.rds"))
cat("Total:", round(tic() - t0), "s\n")
