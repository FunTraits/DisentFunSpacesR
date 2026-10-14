## Analyses of the DisentFunSpaces software paper ------------------------------
## Reproduces every number and figure of the paper: worked example on the
## bird data shipped with the package (A), simulation experiments (B, C).
## Run: source(system.file("paper", "paper_analyses.R", package = "DisentFunSpaces"))
## Outputs are written to ./paper_outputs
suppressPackageStartupMessages({library(DisentFunSpaces); library(ggplot2)})
out <- "paper_outputs"; dir.create(out, showWarnings = FALSE)
R_SEED <- 2026

## ---- A. Worked example: 1000 bird species (real data, ?birds_traits) ---------
## Traits of Toussaint et al. (2026): locomotion (5), reproduction (8), diet
## (10 fuzzy categories). Missing values are imputed with missForest and the
## first ten phylogenetic eigenvectors, as in the published analysis.
data(birds_traits); data(birds_trait_info)
tr <- prepare_traits(birds_traits, birds_trait_info, species_col = "species")
preds <- birds_traits[, paste0("Eigen.", 1:10)]
rownames(preds) <- birds_traits$species
tr <- impute_traits(tr, covariates = preds, ntree = 100, seed = R_SEED)
spaces <- build_spaces(tr, envfit_perm = 0, engine = "fast")
spaces$Combined_eqdim <- build_space(tr, name = "Combined_eqdim", weighting = "dimension",
                                     envfit_perm = 0, engine = "fast")
spaces$Combined_scaled <- build_space(tr, name = "Combined_scaled", weighting = "dimension_scaled",
                                      envfit_perm = 0, engine = "fast")
dims <- c("Locomotion", "Reproduction", "Diet")
## validation space of the paper: beak traits (use = FALSE by default)
ti_beak <- birds_trait_info; ti_beak$use <- TRUE
beak <- build_space(prepare_traits(birds_traits, ti_beak, "species"), "Beak", name = "Beak",
                    envfit_perm = 0, engine = "fast")
pr_beak <- compare_spaces(structure(list(Locomotion = spaces$Locomotion, Diet = spaces$Diet,
                                         Beak = beak), class = "fs_spaces"), nrepet = 999)
write.csv(pr_beak$table, file.path(out, "TableA2_beak_validation.csv"), row.names = FALSE)
print(pr_beak$table[, c("space1", "space2", "r", "p_value")], digits = 3)
pr <- compare_spaces(structure(spaces[c(dims, "Combined")], class = "fs_spaces"), nrepet = 999)
pr_eq <- compare_spaces(structure(spaces[c(dims, "Combined_eqdim")], class = "fs_spaces"), nrepet = 999)
pr_sc <- compare_spaces(structure(spaces[c(dims, "Combined_scaled")], class = "fs_spaces"), nrepet = 999)
tabA <- pr$table[, c("space1", "space2", "n_axes", "r", "p_value")]
eqt <- pr_eq$table; eqt$space2[eqt$space2 == "Combined_eqdim"] <- "Combined"
tabA$r_equal_dims <- eqt$r[match(paste(tabA$space1, tabA$space2), paste(eqt$space1, eqt$space2))]
sct <- pr_sc$table; sct$space2[sct$space2 == "Combined_scaled"] <- "Combined"
tabA$r_scaled_dims <- sct$r[match(paste(tabA$space1, tabA$space2), paste(sct$space1, sct$space2))]
write.csv(tabA, file.path(out, "TableA_procrustes_example.csv"), row.names = FALSE)
print(tabA, digits = 3)
iucn <- stats::setNames(birds_traits$iucn, birds_traits$species)
met <- species_metrics_all(structure(spaces[c(dims, "Combined")], class = "fs_spaces"), k = 5,
                           iucn = iucn)
red <- do.call(rbind, lapply(split(met, met$space), function(x)
  data.frame(space = x$space[1], pct_unique = 100 * mean(x$n_identical == 0),
             pct_FUn_zero = 100 * mean(x$FUn_raw == 0), median_identical = median(x$n_identical))))
wide <- reshape(met[, c("space", "species", "n_identical", "FUn")], idvar = "species",
                timevar = "space", direction = "wide")
uq <- wide[wide$n_identical.Combined == 0, ]
only_comb <- sum(rowSums(uq[, paste0("n_identical.", dims)] == 0) == 0)
cm <- suppressWarnings(compare_metrics(met, "FUn", top = 0.1))
write.csv(red, file.path(out, "TableB_redundancy_example.csv"), row.names = FALSE)
write.csv(cm$pairs, file.path(out, "TableC_FUn_congruence_example.csv"), row.names = FALSE)
print(red, digits = 3); print(cm$pairs, digits = 3)
cat("Unique in Combined:", nrow(uq), "| unique in no single dimension:", only_comb, "\n")

pal <- c(Locomotion = "#2E7D32", Reproduction = "#E07B00", Diet = "#1565C0")
ps <- lapply(c("Combined", dims), function(nm)
  suppressMessages(plot_space(spaces[[nm]], min_cor = 0.5, point_size = 0.6, point_alpha = 0.4,
                              title = paste(nm, "space")) +
                     scale_colour_manual(values = pal, limits = names(pal), name = "Dimension")))
if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  figA <- (plot_spacer() + ps[[1]] + plot_spacer() + plot_layout(widths = c(0.5, 1, 0.5))) /
    (ps[[2]] + ps[[3]] + ps[[4]] + plot_layout(nrow = 1)) + plot_layout(guides = "collect")
  ggsave(file.path(out, "Fig2_example_spaces.png"), figA, width = 14, height = 9.5, dpi = 200, bg = "white")
  figB <- plot_comparison(pr) + labs(title = "a. Each trait = 1 (published)", subtitle = NULL) |
    plot_comparison(pr_eq) + labs(title = "b. Each dimension = 1", subtitle = NULL) |
    plot_comparison(pr_sc) + labs(title = "c. Each dimension = 1, rescaled", subtitle = NULL)
  figB <- figB + plot_layout(guides = "collect")
  ggsave(file.path(out, "Fig3_example_procrustes.png"), figB, width = 17, height = 5.5, dpi = 200, bg = "white")
}

## ---- B. Experiment 1: does Procrustes r track the imposed congruence? ---------
## Two dimensions of mixed types; correlation between their first latent axes
## varied from 0 to 0.9; 20 replicates of 300 species.
design1 <- list(A = list(list(type = "Q", n = 4)),
                B = list(list(type = "F", n = 1, modalities = 6), list(type = "B", n = 1, modalities = 5)))
grid1 <- expand.grid(rho = seq(0, 0.9, 0.1), rep = 1:20)
res1 <- do.call(rbind, lapply(seq_len(nrow(grid1)), function(i) {
  s <- simulate_traits(300, design1, cor_dims = grid1$rho[i], seed = 1000 + i)
  t <- prepare_traits(s$traits, s$trait_info, "species")
  sp <- build_spaces(t, spaces = list(A = "A", B = "B"), envfit_perm = 0, keep_dist = FALSE)
  data.frame(grid1[i, ], r = compare_spaces(sp, nrepet = 9)$table$r)
}))
write.csv(res1, file.path(out, "Experiment1_congruence.csv"), row.names = FALSE)
s1 <- aggregate(r ~ rho, res1, function(v) c(mean = mean(v), lo = quantile(v, 0.05), hi = quantile(v, 0.95)))
s1 <- data.frame(rho = s1$rho, s1$r); names(s1) <- c("rho", "mean", "lo", "hi")
print(s1, digits = 3)
cat("Spearman rho(imposed, r):", cor(res1$rho, res1$r, method = "spearman"), "\n")

## ---- C. Experiment 2: which dimension dominates the combined space? -----------
## Two independent dimensions (rho = 0). A: k continuous traits (k = 1..8);
## B: 3 continuous traits. Combined space under both weightings.
design2 <- function(k) list(A = list(list(type = "Q", n = k)), B = list(list(type = "Q", n = 3)))
grid2 <- expand.grid(k = 1:8, rep = 1:20)
res2 <- do.call(rbind, lapply(seq_len(nrow(grid2)), function(i) {
  s <- simulate_traits(300, design2(grid2$k[i]), cor_dims = 0, seed = 5000 + i)
  t <- prepare_traits(s$traits, s$trait_info, "species")
  do.call(rbind, lapply(c("variable", "dimension", "dimension_scaled"), function(w) {
    sp <- build_spaces(t, spaces = list(A = "A", B = "B", Combined = c("A", "B")),
                       weighting = w, envfit_perm = 0, keep_dist = FALSE)
    tb <- compare_spaces(sp, nrepet = 9)$table
    data.frame(grid2[i, ], weighting = w,
               r_A = tb$r[tb$space1 == "A" & tb$space2 == "Combined"],
               r_B = tb$r[tb$space1 == "B" & tb$space2 == "Combined"])
  }))
}))
write.csv(res2, file.path(out, "Experiment2_dominance.csv"), row.names = FALSE)
s2 <- aggregate(cbind(r_A, r_B) ~ k + weighting, res2, mean)
print(s2, digits = 3)

## Experiment 2b: same number of traits (3 vs 3), different trait types for A
types2b <- list(continuous = list(list(type = "Q", n = 3)),
                binary = list(list(type = "N", n = 3, levels = 2)),
                multichoice = list(list(type = "B", n = 3, modalities = 4)),
                nominal = list(list(type = "N", n = 3, levels = 3)))
grid2b <- expand.grid(typeA = names(types2b), rep = 1:20, stringsAsFactors = FALSE)
res2b <- do.call(rbind, lapply(seq_len(nrow(grid2b)), function(i) {
  s <- simulate_traits(300, list(A = types2b[[grid2b$typeA[i]]], B = list(list(type = "Q", n = 3))),
                       cor_dims = 0, seed = 9000 + i)
  t <- prepare_traits(s$traits, s$trait_info, "species")
  do.call(rbind, lapply(c("variable", "dimension", "dimension_scaled"), function(w) {
    sp <- build_spaces(t, spaces = list(A = "A", B = "B", Combined = c("A", "B")),
                       weighting = w, envfit_perm = 0, keep_dist = FALSE)
    tb <- compare_spaces(sp, nrepet = 9)$table
    data.frame(grid2b[i, ], weighting = w,
               r_A = tb$r[tb$space1 == "A" & tb$space2 == "Combined"],
               r_B = tb$r[tb$space1 == "B" & tb$space2 == "Combined"])
  }))
}))
write.csv(res2b, file.path(out, "Experiment2b_trait_types.csv"), row.names = FALSE)
s2b <- aggregate(cbind(r_A, r_B) ~ typeA + weighting, res2b, mean)
print(s2b, digits = 3)

## ---- Figure 4: simulation experiments -------------------------------------------
f1 <- ggplot(s1, aes(rho, mean)) + geom_abline(linetype = 2, colour = "grey60") +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.02, colour = "grey40") +
  geom_point(size = 2, colour = "#1F4E79") + geom_line(colour = "#1F4E79") +
  coord_cartesian(xlim = c(0, 0.95), ylim = c(0, 1)) +
  labs(x = "Imposed correlation between latent niche axes", y = "Procrustes r (A vs B)",
       title = "a. Procrustes r tracks the imposed congruence",
       subtitle = "mean and 90% interval of 20 simulations, 300 species") + theme_classic(base_size = 11)
s2l <- rbind(data.frame(s2[, c("k", "weighting")], dim = "A (k traits)", r = s2$r_A),
             data.frame(s2[, c("k", "weighting")], dim = "B (3 traits)", r = s2$r_B))
s2l$weighting <- factor(s2l$weighting, c("variable", "dimension", "dimension_scaled"),
                        c("each trait = 1 (published)", "each dimension = 1", "each dimension = 1, rescaled"))
f2 <- ggplot(s2l, aes(k, r, colour = dim, linetype = weighting, shape = weighting)) + geom_line() + geom_point(size = 1.6) +
  scale_colour_manual(values = c("#C62828", "#1565C0"), name = NULL) +
  scale_x_continuous(breaks = 1:8) + coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Number of traits in dimension A", y = "Procrustes r with the combined space",
       linetype = NULL, shape = NULL, title = "b. The dimension with more traits dominates the combined space",
       subtitle = "two independent dimensions, continuous traits") +
  theme_classic(base_size = 11) + theme(legend.position = "bottom", legend.box = "vertical")
s2bl <- rbind(data.frame(s2b[, c("typeA", "weighting")], dim = "A (3 traits, type varies)", r = s2b$r_A),
              data.frame(s2b[, c("typeA", "weighting")], dim = "B (3 continuous traits)", r = s2b$r_B))
s2bl$weighting <- factor(s2bl$weighting, c("variable", "dimension", "dimension_scaled"),
                         c("each trait = 1 (published)", "each dimension = 1", "each dimension = 1, rescaled"))
s2bl$typeA <- factor(s2bl$typeA, names(types2b))
f3 <- ggplot(s2bl, aes(typeA, r, fill = dim)) + geom_col(position = position_dodge(0.8), width = 0.7) +
  facet_wrap(~weighting) + scale_fill_manual(values = c("#C62828", "#1565C0"), name = NULL) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Type of the traits of dimension A", y = "Procrustes r with the combined space",
       title = "c. Trait type also changes which dimension dominates",
       subtitle = "same number of traits (3 vs 3), two independent dimensions; with 3 vs 3 traits the first two weightings coincide") +
  theme_classic(base_size = 11) + theme(legend.position = "bottom")
if (requireNamespace("patchwork", quietly = TRUE)) {
  fig4 <- (f1 | f2) / f3
  ggsave(file.path(out, "Fig4_simulation_experiments.png"), fig4, width = 13, height = 10, dpi = 200, bg = "white")
}
saveRDS(list(example = tr, procrustes = tabA, beak_validation = pr_beak$table, redundancy = red, fun_congruence = cm$pairs,
             unique_combined = nrow(uq), unique_only_combination = only_comb,
             exp1 = res1, exp1_summary = s1, exp2 = res2, exp2_summary = s2,
             exp2b = res2b, exp2b_summary = s2b),
        file.path(out, "paper_results.rds"))
message("Done: outputs in ", normalizePath(out))
