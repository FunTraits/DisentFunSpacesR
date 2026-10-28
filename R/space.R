#' Build a functional space (PCoA) from one or several functional dimensions
#'
#' Wrapper around [PCoAGraph()], the core function of Toussaint et al.
#' (2026): from the trait-information table, the traits of the selected
#' dimensions are assembled into one table per dimension x type (fuzzy and
#' multichoice traits prepared with [ade4::prep.fuzzy()] /
#' [ade4::prep.binary()]), combined with [ade4::ktab.list.df()], and the
#' generalised Gower distance computed by [ade4::dist.ktab()]
#' (`option = "scaledBYrange"`) is ordinated with [ape::pcoa()]. Trait-axis
#' relationships are given by [vegan::envfit()] and [CorrelType()].
#'
#' With `weighting = "variable"` and equal trait weights (the published
#' analysis), the distance is that of a single `dist.ktab()` call on all
#' tables. With `weighting = "dimension"` or unequal weights in
#' `trait_info$weight`, the distance is computed with [trait_distance()]
#' (also based on `dist.ktab()`) and passed to [PCoAGraph()].
#'
#' @inheritParams trait_distance
#' @param name Name of the space (default: dimensions pasted with "+").
#' @param correction Correction for negative eigenvalues passed to
#'   [ape::pcoa()]: `"none"` (as in the original analysis), `"cailliez"` or
#'   `"lingoes"`.
#' @param max_axes Number of PCoA axes kept in `coords` (`NULL`, default as
#'   in the original analysis: all positive-eigenvalue axes). The full
#'   `ape::pcoa()` output is kept in `PCoA`, unless `eigen = "partial"`.
#' @param eigen `"full"` (default): complete eigendecomposition, as in the
#'   original analysis. `"partial"` (without correction, package 'RSpectra'
#'   installed): only the axes kept are computed, i.e. the first `max_axes`
#'   axes or, with `max_axes = NULL`, all positive-eigenvalue axes (see
#'   [PCoAGraph()]): same coordinates (signs of axes are arbitrary) and
#'   relative inertia, seconds instead of about ten minutes per space for
#'   8,413 bird species; `PCoA` then holds these axes only and the
#'   negative-eigenvalue diagnostics are `NA`.
#' @param envfit_perm,envfit_axes Passed to [PCoAGraph()].
#' @param n_cor_axes Number of axes kept in the tidy correlation table
#'   `cor` (`PCoACor` keeps all axes).
#' @param keep_dist Keep the distance in the output (needed for
#'   [space_quality()] and Mantel tests). Default: `TRUE` if fewer than 5000
#'   species.
#' @param saveFile Optional `.rds` path where the [PCoAGraph()] output is
#'   saved.
#' @param engine Distance computation: `"ade4"` (default, `dist.ktab()` as in
#'   the published analysis) or `"fast"` (same values, see
#'   [trait_distance()]; much lower memory use for several thousand species
#'   with missing values in multichoice traits, where `dist.ktab()` builds
#'   many n x n matrices). The PCoA always uses [PCoAGraph()].
#'
#' @return An object of class `fs_space`, a list holding the
#'   [PCoAGraph()] elements `PCoA` (`ape::pcoa()` output), `varPlot`
#'   (envfit scores) and `PCoACor` (trait x axis correlations), so that it
#'   can be used by the published scripts (e.g. `space$PCoA$vectors`),
#'   plus: `name`, `dimensions`, `coords` (species x axes matrix used by the
#'   other functions of the package), `eigen` (eigenvalues, relative and
#'   cumulative inertia), `cor` (tidy correlation table with trait metadata),
#'   `blocks` (tables and types passed to `dist.ktab()`), `variables`,
#'   `weighting`, `correction`, `excluded` (species with no value for any
#'   trait of the space, left out of it) and optionally `dist`.
#' @export
#' @examples
#' data(birds_traits); data(birds_trait_info)
#' tr <- prepare_traits(birds_traits, birds_trait_info, species_col = "species")
#' tr <- drop_incomplete(tr)
#' sp <- build_space(tr, dimensions = c("Diet", "Locomotion"), envfit_perm = 0)
#' sp
#' sp$blocks$groups
#' head(sp$PCoA$vectors[, 1:3])
build_space <- function(x, dimensions = NULL, name = NULL,
                        weighting = c("variable", "dimension", "dimension_scaled"),
                        correction = c("none", "cailliez", "lingoes"),
                        max_axes = NULL, envfit_perm = 999, envfit_axes = c(1, 2),
                        n_cor_axes = 10, keep_dist = NULL, saveFile = NA,
                        engine = c("ade4", "fast"), eigen = c("full", "partial")) {
  stopifnot(inherits(x, "fs_traits"))
  weighting <- match.arg(weighting)
  correction <- match.arg(correction)
  engine <- match.arg(engine)
  eigen <- match.arg(eigen)
  if (is.null(dimensions)) dimensions <- unique(x$trait_info$dimension)
  unk <- setdiff(dimensions, x$trait_info$dimension)
  if (length(unk)) stop("Unknown dimension(s): ", paste(unk, collapse = ", "),
                        ". Available: ", paste(unique(x$trait_info$dimension), collapse = ", "))
  if (is.null(name)) name <- paste(dimensions, collapse = "+")
  if (is.null(keep_dist)) keep_dist <- length(x$species) < 5000

  # species without any value for the traits of this space are excluded
  ti <- x$trait_info[x$trait_info$dimension %in% dimensions, , drop = FALSE]
  empty <- rowSums(!is.na(x$data[, ti$trait, drop = FALSE])) == 0
  excluded <- x$species[empty]
  if (length(excluded)) {
    message(sprintf("Space '%s': %d species without any trait value excluded (e.g. %s).",
                    name, length(excluded), excluded[1]))
    x <- subset_species(x, x$species[!empty])
  }

  bl <- .make_blocks(x$data, ti)
  equal_w <- length(unique(bl$variables$weight)) == 1
  d <- NULL
  if (weighting != "variable" || !equal_w || engine == "fast") {
    d <- trait_distance(x, dimensions = dimensions, weighting = weighting, engine = engine)
  }
  res <- PCoAGraph(x = bl$x, datax = bl$datax, namesx = x$species, groups = bl$groups,
                   saveFile = saveFile, correction = correction,
                   envfit_perm = envfit_perm, envfit_axes = envfit_axes, dist = d,
                   partial = eigen == "partial", n_axes = max_axes)

  pc <- res$PCoA
  vals <- pc$values
  use_cor <- !is.null(pc$vectors.cor)
  vec <- if (use_cor) pc$vectors.cor else pc$vectors
  eig <- if (use_cor) vals$Corr_eig else vals$Eigenvalues
  rel <- if (use_cor) vals$Rel_corr_eig else vals$Relative_eig
  if (!is.null(max_axes)) vec <- vec[, seq_len(min(max_axes, ncol(vec))), drop = FALSE]
  keep <- seq_len(ncol(vec))
  colnames(vec) <- paste0("PCoA", keep)
  eigen_tab <- data.frame(axis = colnames(vec), eigenvalue = eig[keep],
                          relative = rel[keep], cumulative = cumsum(rel[keep]))
  negv <- vals$Eigenvalues[vals$Eigenvalues < -sqrt(.Machine$double.eps)]
  neg <- if (isTRUE(pc$partial)) list(n = NA_integer_, ratio = NA_real_) else
    list(n = length(negv),
         ratio = if (length(negv)) sum(abs(negv)) / sum(vals$Eigenvalues[vals$Eigenvalues > 0]) else 0)

  # tidy correlation table (first n_cor_axes axes) with trait metadata
  k <- min(n_cor_axes, ncol(res$PCoACor))
  cm <- res$PCoACor[, seq_len(k), drop = FALSE]
  colnames(cm) <- paste0("PCoA", seq_len(k))
  m <- match(rownames(cm), bl$trait_info$trait)
  cor_tab <- data.frame(trait = rownames(cm), block = bl$trait_info$block[m],
                        dimension = bl$trait_info$dimension[m], type = bl$trait_info$type[m],
                        cm, row.names = NULL, check.names = FALSE, stringsAsFactors = FALSE)

  structure(list(
    PCoA = pc,
    varPlot = res$varPlot,
    PCoACor = res$PCoACor,
    name = name,
    dimensions = dimensions,
    coords = vec,
    eigen = eigen_tab,
    negative_eigen = neg,
    cor = cor_tab,
    blocks = list(x = bl$x, groups = bl$groups),
    variables = bl$variables,
    weighting = weighting,
    correction = correction,
    dist = if (keep_dist) res$dist else NULL,
    excluded = excluded
  ), class = "fs_space")
}

#' @export
print.fs_space <- function(x, ...) {
  cat("<fs_space> '", x$name, "' - ", nrow(x$coords), " species, ",
      ncol(x$coords), " PCoA axes (", x$weighting, " weighting, correction: ",
      x$correction, ")\n", sep = "")
  cat("Dimensions:", paste(x$dimensions, collapse = ", "), "\n")
  e <- utils::head(x$eigen, 4)
  cat("Inertia of first axes (%):", paste(sprintf("%s %.1f", e$axis, 100 * e$relative), collapse = " | "), "\n")
  if (isTRUE(x$negative_eigen$n > 0))
    cat(sprintf("Negative eigenvalues: %d (|sum| = %.1f%% of positive inertia)\n",
                x$negative_eigen$n, 100 * x$negative_eigen$ratio))
  invisible(x)
}

#' Build several functional spaces at once
#'
#' By default builds one space per functional dimension plus a combined
#' space with all dimensions (the "disaggregated" and "aggregated" spaces of
#' Toussaint et al. 2026). Any other combination can be requested.
#'
#' @inheritParams build_space
#' @param spaces Named list of character vectors; each element gives the
#'   dimensions of one space, e.g. `list(Diet = "Diet", Combined =
#'   c("Diet", "Locomotion"))`. Default: every single dimension plus
#'   `Combined` (all dimensions) when there are several.
#' @param ... Further arguments to [build_space()].
#' @param progress Optional function called as `progress(i, n, name)` (used
#'   by the Shiny app).
#'
#' @return An object of class `fs_spaces` (named list of `fs_space`).
#' @export
#' @examples
#' data(birds_traits); data(birds_trait_info)
#' tr <- drop_incomplete(prepare_traits(birds_traits, birds_trait_info, "species"))
#' spaces <- build_spaces(tr)
#' spaces
build_spaces <- function(x, spaces = NULL, ..., progress = NULL) {
  stopifnot(inherits(x, "fs_traits"))
  if (is.null(spaces)) spaces <- default_spaces(x)
  if (is.null(names(spaces)) || any(names(spaces) == ""))
    stop("`spaces` must be a named list.")
  out <- list()
  for (i in seq_along(spaces)) {
    if (!is.null(progress)) progress(i, length(spaces), names(spaces)[i])
    out[[names(spaces)[i]]] <- build_space(x, dimensions = spaces[[i]],
                                           name = names(spaces)[i], ...)
  }
  structure(out, class = "fs_spaces")
}

#' Default list of spaces: each dimension alone plus all combined
#' @param x An `fs_traits` object.
#' @return Named list of dimension vectors.
#' @export
default_spaces <- function(x) {
  dims <- unique(x$trait_info$dimension)
  sp <- stats::setNames(as.list(dims), dims)
  if (length(dims) > 1) sp$Combined <- dims
  sp
}

#' @export
print.fs_spaces <- function(x, ...) {
  cat("<fs_spaces>", length(x), "functional spaces\n")
  tab <- space_summary(x)
  print(tab, row.names = FALSE)
  invisible(x)
}

#' Summary table of several functional spaces
#' @param x An `fs_spaces` object.
#' @param n_axes Number of axes for which inertia is reported.
#' @return Data frame: one row per space.
#' @export
space_summary <- function(x, n_axes = 4) {
  do.call(rbind, lapply(x, function(s) {
    e <- s$eigen
    k <- min(n_axes, nrow(e))
    res <- data.frame(space = s$name, dimensions = paste(s$dimensions, collapse = " + "),
                      n_species = nrow(s$coords), n_variables = nrow(s$variables),
                      stringsAsFactors = FALSE)
    for (j in seq_len(n_axes))
      res[[paste0("PCoA", j, "_pct")]] <- if (j <= k) round(100 * e$relative[j], 2) else NA
    res[[paste0("cum", n_axes, "_pct")]] <- round(100 * e$cumulative[k], 2)
    res
  }))
}

#' Species coordinates in all spaces, as one long table
#' @param spaces An `fs_spaces` object.
#' @param n_axes Maximum number of axes exported (spaces with fewer axes are
#'   padded with `NA`).
#' @return Data frame with `space`, `species` and `PCoA1`... columns.
#' @export
coords_table <- function(spaces, n_axes = 10) {
  do.call(rbind, lapply(spaces, function(s) {
    m <- matrix(NA_real_, nrow(s$coords), n_axes, dimnames = list(NULL, paste0("PCoA", seq_len(n_axes))))
    k <- min(n_axes, ncol(s$coords))
    m[, seq_len(k)] <- s$coords[, seq_len(k)]
    data.frame(space = s$name, species = rownames(s$coords), m, row.names = NULL,
               stringsAsFactors = FALSE)
  }))
}

#' Quality of a functional space as a function of the number of axes
#'
#' Mean absolute deviation (mAD) between the initial functional distance and
#' the Euclidean distance in the first k PCoA axes, after scaling the latter
#' to the same maximum (Maire et al. 2015, <doi:10.1111/geb.12299>). Lower is
#' better. The Pearson correlation between both distances is also returned.
#'
#' @param space An `fs_space` built with `keep_dist = TRUE`.
#' @param max_axes Maximum number of axes evaluated.
#' @param n_sample If the space has more species, a random subsample of
#'   this size is used.
#' @param seed Random seed for subsampling.
#' @return Data frame: `n_axes`, `mAD`, `cor`.
#' @export
space_quality <- function(space, max_axes = 10, n_sample = 2000, seed = 1) {
  stopifnot(inherits(space, "fs_space"))
  if (is.null(space$dist)) stop("Distance not stored: rebuild the space with keep_dist = TRUE.")
  dm <- space$dist
  X <- space$coords
  n <- nrow(X)
  if (n > n_sample) {
    set.seed(seed)
    idx <- sort(sample(n, n_sample))
    dm <- stats::as.dist(as.matrix(dm)[idx, idx])
    X <- X[idx, , drop = FALSE]
  }
  d0 <- as.vector(dm)
  ks <- seq_len(min(max_axes, ncol(X)))
  do.call(rbind, lapply(ks, function(k) {
    de <- as.vector(stats::dist(X[, seq_len(k), drop = FALSE]))
    de <- de / max(de) * max(d0)
    data.frame(n_axes = k, mAD = mean(abs(d0 - de)), cor = stats::cor(d0, de))
  }))
}

