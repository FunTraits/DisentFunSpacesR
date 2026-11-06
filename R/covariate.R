## ---------------------------------------------------------------------------
## Continuous attributes of species that are not traits (range size,
## population trend, economic value, abundance, body condition index...):
## how are they structured in each functional space?
## ---------------------------------------------------------------------------

## Internal: transformed values of x for the species of a space
.covariate_values <- function(x, species, transform = "none") {
  if (is.null(names(x))) stop("'x' must be named by species.", call. = FALSE)
  v <- suppressWarnings(as.numeric(x[match(species, names(x))]))
  v <- switch(transform,
    none = v,
    log10 = { if (any(v <= 0, na.rm = TRUE)) stop("log10 needs positive values.", call. = FALSE); log10(v) },
    log1p = { if (any(v < 0, na.rm = TRUE)) stop("log1p needs values >= 0.", call. = FALSE); log10(1 + v) },
    rank = rank(v, na.last = "keep") / sum(is.finite(v)),
    stop("Unknown transform.", call. = FALSE))
  stats::setNames(v, species)
}

#' A continuous attribute of species (not a trait) in functional spaces
#'
#' Relates a continuous variable that is not a trait (e.g. range size,
#' population trend, abundance, economic value, a vulnerability score) to the
#' position of species in every functional space, so that the dimensions in
#' which the variable is structured can be identified, without cutting it
#' into classes.
#'
#' Per space (`summary`):
#' * `R2`, `R2_adj`, `p_fit`: linear fit of the variable on the first
#'   `fit_axes` axes (vector fitting, as `vegan::envfit()`), `p_fit` from
#'   `nperm` permutations of the values;
#' * `knn_cor`, `p_knn`: correlation between the value of each species and
#'   the mean value of its `k` nearest neighbours in the space (functional
#'   autocorrelation: > 0 = functionally similar species have similar
#'   values, also when the relation is not linear), permutation p-value;
#' * `rho_FUn`, `p_FUn`, `rho_FSp`, `p_FSp`: Spearman correlations with
#'   functional uniqueness and specialisation ([species_metrics()]): are
#'   functionally distinct species also, e.g., rarer or more exploited?
#'
#' Species without a value are left out.
#'
#' @param spaces An `fs_spaces` object, a named list of `fs_space`, or one
#'   `fs_space`.
#' @param x Named numeric vector (names = species).
#' @param transform `"none"`, `"log10"`, `"log1p"` (log10(1 + x)) or
#'   `"rank"`.
#' @inheritParams species_metrics
#' @param fit_axes Axes used for the linear fit (default 4).
#' @param nperm Permutations.
#' @param seed Optional random seed.
#' @return An object of class `fs_covariate`: list with `summary`,
#'   `species` (space, species, value, knn_mean, FUn_raw, FSp_raw) and the
#'   fitted directions `directions` (space, axis, coefficient).
#' @export
#' @examples
#' \donttest{
#' data(birds_traits); data(birds_trait_info)
#' tr <- drop_incomplete(prepare_traits(birds_traits, birds_trait_info, "species"))
#' sp <- build_spaces(tr, spaces = list(Diet = "Diet", Locomotion = "Locomotion"))
#' # here: the first phylogenetic eigenvector, as an example of a non-trait variable
#' cv <- covariate_space(sp, setNames(birds_traits$Eigen.1, birds_traits$species), nperm = 99)
#' cv$summary
#' }
covariate_space <- function(spaces, x, transform = c("none", "log10", "log1p", "rank"),
                            n_axes = NULL, k = 5, fit_axes = 4, nperm = 999, seed = NULL) {
  transform <- match.arg(transform)
  if (inherits(spaces, "fs_space")) spaces <- list(space = spaces)
  if (!is.null(seed)) set.seed(seed)
  res <- lapply(names(spaces), function(nm) {
    X <- .space_coords(spaces[[nm]], n_axes)
    v <- .covariate_values(x, rownames(X), transform)
    ok <- is.finite(v)
    if (sum(ok) <= max(k, fit_axes) + 2) stop("Too few species with a value in space '", nm, "'.", call. = FALSE)
    m <- species_metrics(spaces[[nm]], n_axes = n_axes, k = k)
    X <- X[ok, , drop = FALSE]; v <- v[ok]; n <- nrow(X)
    fa <- min(fit_axes, ncol(X))
    Q <- qr.Q(qr(cbind(1, X[, seq_len(fa), drop = FALSE])))
    r2 <- function(y) { yc <- y - mean(y); f <- Q %*% crossprod(Q, y) - mean(y); sum(f^2) / sum(yc^2) }
    R2 <- r2(v)
    nnk <- FNN::get.knn(X, k = min(k, n - 1))$nn.index
    kc <- function(y) stats::cor(y, rowMeans(matrix(y[nnk], nrow = n)))
    KC <- kc(v)
    pr <- if (nperm > 0) replicate(nperm, { y <- sample(v); c(r2(y), kc(y)) }) else NULL
    pv <- function(obs, sim) if (is.null(sim)) NA else (sum(sim >= obs - 1e-12) + 1) / (length(sim) + 1)
    mm <- m[match(rownames(X), m$species), ]
    ct <- function(a) suppressWarnings(stats::cor.test(v, a, method = "spearman", exact = FALSE))
    c1 <- ct(mm$FUn_raw); c2 <- ct(mm$FSp_raw)
    beta <- stats::coef(stats::lm(v ~ X[, seq_len(fa), drop = FALSE]))[-1]
    list(summary = data.frame(space = nm, n = n, R2 = R2,
                              R2_adj = 1 - (1 - R2) * (n - 1) / (n - fa - 1),
                              p_fit = pv(R2, pr[1, ]), knn_cor = KC, p_knn = pv(KC, pr[2, ]),
                              rho_FUn = unname(c1$estimate), p_FUn = c1$p.value,
                              rho_FSp = unname(c2$estimate), p_FSp = c2$p.value,
                              stringsAsFactors = FALSE),
         species = data.frame(space = nm, species = rownames(X), value = v,
                              knn_mean = rowMeans(matrix(v[nnk], nrow = n)),
                              FUn_raw = mm$FUn_raw, FSp_raw = mm$FSp_raw,
                              row.names = NULL, stringsAsFactors = FALSE),
         directions = data.frame(space = nm, axis = seq_len(fa), coefficient = unname(beta)))
  })
  bind <- function(w) { d <- do.call(rbind, lapply(res, `[[`, w)); rownames(d) <- NULL; d }
  structure(list(summary = bind("summary"), species = bind("species"),
                 directions = bind("directions"),
                 settings = list(transform = transform, k = k, fit_axes = fit_axes, nperm = nperm)),
            class = "fs_covariate")
}

#' @export
print.fs_covariate <- function(x, ...) {
  cat(sprintf("Continuous variable in functional spaces (transform: %s)\n", x$settings$transform))
  print(x$summary, row.names = FALSE, digits = 3)
  invisible(x)
}

#' Plot a continuous attribute of species in a functional space
#'
#' Species coloured by the value of the variable, with the direction of the
#' linear fit on the plotted axes (from [covariate_space()]).
#'
#' @param space An `fs_space`.
#' @inheritParams covariate_space
#' @param axes Two axes to plot.
#' @param name Name of the variable (legend).
#' @return A ggplot.
#' @export
plot_covariate <- function(space, x, axes = c(1, 2), transform = c("none", "log10", "log1p", "rank"),
                           name = "value") {
  transform <- match.arg(transform)
  X <- space$coords
  v <- .covariate_values(x, rownames(X), transform)
  d <- data.frame(a1 = X[, axes[1]], a2 = X[, axes[2]], x_value = v, species = rownames(X))
  d$hover <- sprintf("%s<br>%s = %s", d$species, name, signif(d$x_value, 3))
  d <- d[order(is.finite(d$x_value), d$x_value), ]
  ok <- is.finite(d$x_value)
  b <- stats::coef(stats::lm(x_value ~ a1 + a2, data = d[ok, ]))[2:3]
  r <- 0.4 * max(diff(range(d$a1)), diff(range(d$a2)))
  b <- if (all(is.finite(b)) && sum(b^2) > 0) b / sqrt(sum(b^2)) * r else c(0, 0)
  c0 <- c(mean(d$a1[ok]), mean(d$a2[ok]))
  ggplot2::ggplot(d, ggplot2::aes(a1, a2)) +
    .geom_point_hover(data = d[!ok, ], ggplot2::aes(text = hover), colour = "grey80", size = 0.9) +
    .geom_point_hover(data = d[ok, ], ggplot2::aes(colour = x_value, text = hover), size = 1.4,
                      alpha = 0.85) +
    ggplot2::annotate("segment", x = c0[1], y = c0[2], xend = c0[1] + b[1], yend = c0[2] + b[2],
                      arrow = ggplot2::arrow(length = ggplot2::unit(0.25, "cm")), linewidth = 0.8) +
    ggplot2::scale_colour_viridis_c(name = if (transform == "none") name else sprintf("%s (%s)", name, transform)) +
    ggplot2::labs(x = paste0("PCoA ", axes[1]), y = paste0("PCoA ", axes[2]),
                  subtitle = "Arrow: direction of increase (linear fit on the plotted axes)") +
    ggplot2::coord_equal() + ggplot2::theme_classic(base_size = 11)
}
