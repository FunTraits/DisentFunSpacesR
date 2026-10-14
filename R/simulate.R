#' Default design for simulated trait data
#'
#' Four functional dimensions with different numbers and types of traits,
#' used by [simulate_traits()]. Each dimension is a list of trait
#' specifications: `type` (see [trait_types()]), `n` (number of such traits)
#' and, depending on the type, `levels` (ordinal, nominal), `modalities`
#' (fuzzy-coded, multichoice) or `period` (circular).
#'
#' @return A named list of dimensions.
#' @export
#' @examples
#' str(default_sim_design(), max.level = 2)
default_sim_design <- function() {
  list(
    Locomotion = list(list(type = "Q", n = 5), list(type = "N", n = 1, levels = 3)),
    Diet = list(list(type = "F", n = 1, modalities = 6)),
    Reproduction = list(list(type = "Q", n = 2), list(type = "O", n = 1, levels = 4),
                        list(type = "C", n = 1, period = 12)),
    Habitat = list(list(type = "B", n = 2, modalities = 5), list(type = "D", n = 1))
  )
}

#' Simulate mixed-type trait data with known congruence between dimensions
#'
#' Generates species traits of any type (see [trait_types()]) from latent
#' niche axes, so that the congruence between functional dimensions is known.
#' Each dimension has two latent axes. The first axes of the dimensions are
#' correlated according to `cor_dims`; the second axes are independent. Each
#' latent axis combines a phylogenetic component (Brownian motion on a
#' coalescent tree) and an independent component, weighted by
#' `phylo_signal`. Each trait is a noisy function of the two latent axes of
#' its dimension, with random loadings:
#' * `Q`: log-normal continuous values (transform `log10` in `trait_info`);
#' * `O`: the linear predictor cut into `levels` ordered classes;
#' * `N`: one of `levels` categories drawn from softmax probabilities;
#' * `D`: presence drawn from a logistic model (prevalence ~ 0.3);
#' * `B`: `modalities` binary columns, each with its own logistic model,
#'   at least one modality per species;
#' * `F`: `modalities` affinity scores (0-10) from softmax probabilities;
#' * `C`: the angle of the two predictors mapped on `0:(period - 1)`.
#'
#' @param n_species Number of species.
#' @param design Named list of dimensions (see [default_sim_design()]).
#' @param cor_dims Correlation between the first latent axes of the
#'   dimensions: a single value (same for all pairs) or a correlation matrix
#'   with dimension names. By default, Locomotion and Reproduction are
#'   correlated at 0.6 and all other pairs at 0.2.
#' @param phylo_signal Share (0-1) of latent variance due to phylogeny.
#' @param noise Standard deviation of the noise added to each trait's
#'   linear predictor (higher = traits less tied to their dimension).
#' @param missing Proportion of missing values per trait (whole blocks for
#'   fuzzy and multichoice traits).
#' @param tree Optional phylogeny (class `phylo`) with `n_species` tips;
#'   simulated with [ape::rcoal()] if `NULL`.
#' @param seed Optional random seed.
#'
#' @return An object of class `fs_simulation`: a list with `traits` (data
#'   frame, first column `species`), `trait_info` (ready for
#'   [prepare_traits()]), `tree`, `latent` (species x dimensions matrix of
#'   first latent axes), `cor_dims` (the target correlation matrix) and
#'   `settings`.
#' @export
#' @examples
#' sim <- simulate_traits(200, seed = 1)
#' sim
#' tr <- prepare_traits(sim$traits, sim$trait_info, "species")
#' spaces <- build_spaces(tr, envfit_perm = 0)
#' compare_spaces(spaces, nrepet = 99)
simulate_traits <- function(n_species = 500, design = default_sim_design(), cor_dims = NULL,
                            phylo_signal = 0.5, noise = 0.5, missing = 0, tree = NULL,
                            seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  dims <- names(design)
  if (is.null(dims) || any(dims == "")) stop("`design` must be a named list of dimensions.")
  nd <- length(dims)
  R <- .sim_cor(cor_dims, dims)
  if (phylo_signal < 0 || phylo_signal > 1) stop("`phylo_signal` must be between 0 and 1.")
  sp <- sprintf("sp%04d", seq_len(n_species))
  if (is.null(tree)) {
    tree <- ape::rcoal(n_species, tip.label = sp)
  } else {
    if (length(tree$tip.label) != n_species) stop("`tree` must have `n_species` tips.")
    sp <- tree$tip.label
  }
  L <- chol(R)
  std <- function(m) scale(m)[, , drop = FALSE]
  bm <- function(k) std(sapply(seq_len(k), function(i) ape::rTraitCont(tree, model = "BM")[sp]))
  iid <- function(k) matrix(stats::rnorm(n_species * k), n_species, k)
  mix <- function(k) sqrt(phylo_signal) * bm(k) + sqrt(1 - phylo_signal) * iid(k)
  # whiten the sample, then impose exactly the target correlations, so that
  # the congruence between dimensions is known in the simulated sample
  exact <- function(M, L) {
    M <- scale(M)
    W <- M %*% solve(chol(stats::cor(M)))
    std(W %*% L)
  }
  z1 <- exact(mix(nd), L)                       # correlated first axes
  z2 <- exact(mix(nd), diag(nd))                # independent second axes
  colnames(z1) <- colnames(z2) <- dims
  rownames(z1) <- sp

  cols <- list(); info <- list()
  for (d in seq_along(dims)) {
    zz <- cbind(z1[, d], z2[, d])
    k <- 0
    for (spec in design[[d]]) {
      for (j in seq_len(if (is.null(spec$n)) 1 else spec$n)) {
        k <- k + 1
        base <- sprintf("%s_%s%d", substr(dims[d], 1, 4), spec$type, k)
        res <- .sim_trait(spec, zz, noise, n_species)
        for (nm in names(res$values)) cols[[paste0(base, if (nm == "x") "" else nm)]] <- res$values[[nm]]
        info[[base]] <- data.frame(
          trait = paste0(base, ifelse(names(res$values) == "x", "", names(res$values))), type = spec$type, dimension = dims[d],
          block = base, transform = if (spec$type == "Q") "log10" else "none", weight = 1,
          levels = NA_character_,
          circ_min = if (spec$type == "C") 0 else NA_real_,
          circ_max = if (spec$type == "C") spec$period - 1 else NA_real_,
          use = TRUE, stringsAsFactors = FALSE)
      }
    }
  }
  traits <- data.frame(species = sp, cols, stringsAsFactors = FALSE, check.names = FALSE)
  trait_info <- do.call(rbind, info); rownames(trait_info) <- NULL

  if (missing > 0) {
    for (b in unique(trait_info$block)) {
      cc <- trait_info$trait[trait_info$block == b]
      idx <- sample(n_species, round(missing * n_species))
      traits[idx, cc] <- NA
    }
  }
  structure(list(traits = traits, trait_info = trait_info, tree = tree, latent = z1,
                 cor_dims = R,
                 settings = list(n_species = n_species, phylo_signal = phylo_signal,
                                 noise = noise, missing = missing, seed = seed,
                                 design = design)),
            class = "fs_simulation")
}

.sim_cor <- function(cor_dims, dims) {
  nd <- length(dims)
  if (is.null(cor_dims)) {
    R <- matrix(0.2, nd, nd, dimnames = list(dims, dims)); diag(R) <- 1
    if (all(c("Locomotion", "Reproduction") %in% dims))
      R["Locomotion", "Reproduction"] <- R["Reproduction", "Locomotion"] <- 0.6
  } else if (length(cor_dims) == 1) {
    R <- matrix(cor_dims, nd, nd, dimnames = list(dims, dims)); diag(R) <- 1
  } else {
    R <- as.matrix(cor_dims)
    if (!all(dim(R) == nd)) stop("`cor_dims` must be ", nd, " x ", nd, ".")
    if (!is.null(rownames(R))) R <- R[dims, dims] else dimnames(R) <- list(dims, dims)
  }
  ev <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  if (any(abs(R) > 1) || min(ev) <= 0) stop("`cor_dims` is not a valid correlation matrix.")
  R
}

.sim_trait <- function(spec, zz, noise, n) {
  lin <- function() {
    a <- c(sample(c(-1, 1), 1) * stats::runif(1, 0.6, 1), stats::runif(1, -0.6, 0.6))
    as.numeric(zz %*% a) + stats::rnorm(n, sd = noise)
  }
  softmax <- function(m) { e <- exp(m - apply(m, 1, max)); e / rowSums(e) }
  switch(spec$type,
    Q = list(values = list(x = round(exp(1 + 0.5 * lin()), 3))),
    O = {
      e <- lin()
      list(values = list(x = as.integer(cut(e, stats::quantile(e, seq(0, 1, length.out = spec$levels + 1)),
                                            include.lowest = TRUE))))
    },
    N = {
      m <- sapply(seq_len(spec$levels), function(l) 1.5 * lin())
      p <- softmax(m)
      x <- apply(p, 1, function(pr) sample(seq_along(pr), 1, prob = pr))
      list(values = list(x = paste0("cat", x)))
    },
    D = list(values = list(x = stats::rbinom(n, 1, stats::plogis(1.5 * lin() - 1)))),
    B = {
      m <- sapply(seq_len(spec$modalities), function(j) 1.5 * lin() - 0.5)
      x <- matrix(stats::rbinom(length(m), 1, stats::plogis(m)), n)
      none <- rowSums(x) == 0
      x[cbind(which(none), max.col(m[none, , drop = FALSE]))] <- 1
      v <- lapply(seq_len(ncol(x)), function(j) x[, j])
      list(values = stats::setNames(v, paste0(".m", seq_len(ncol(x)))))
    },
    F = {
      m <- sapply(seq_len(spec$modalities), function(j) 1.5 * lin())
      p <- softmax(m)
      x <- round(10 * p)
      x[rowSums(x) == 0, 1] <- 1
      v <- lapply(seq_len(ncol(x)), function(j) x[, j])
      list(values = stats::setNames(v, paste0(".m", seq_len(ncol(x)))))
    },
    C = {
      ang <- atan2(lin(), lin())
      list(values = list(x = floor(((ang + pi) / (2 * pi)) * spec$period) %% spec$period))
    },
    stop("Unknown trait type in design: ", spec$type)
  )
}

#' @export
print.fs_simulation <- function(x, ...) {
  ti <- x$trait_info
  cat("<fs_simulation>", x$settings$n_species, "species,", length(unique(ti$block)), "traits (",
      nrow(ti), "columns ) in", ncol(x$latent), "dimensions\n")
  print(table(dimension = ti$dimension[!duplicated(ti$block)], type = ti$type[!duplicated(ti$block)]))
  cat("Target correlation between first latent axes:\n")
  print(round(x$cor_dims, 2))
  cat("Phylogenetic signal:", x$settings$phylo_signal, "| noise:", x$settings$noise,
      "| missing:", x$settings$missing, "\n")
  invisible(x)
}
