#' Compare functional spaces with Procrustes tests
#'
#' Pairwise Procrustes analyses between functional spaces, as in Toussaint
#' et al. (2026, `run_procrustes()` of the original scripts). For every pair
#' of spaces, species shared by both are aligned and their PCoA coordinates
#' compared with a permutation test of the Procrustes correlation:
#' * `method = "ade4"`: [ade4::procuste.rtest()] (original analysis);
#' * `method = "vegan"`: symmetric [vegan::protest()], which also returns
#'   the Procrustes sum of squares m^2.
#'
#' A Mantel test between the two functional distance matrices can be added
#' (`mantel = TRUE`) when the spaces were built with `keep_dist = TRUE`.
#'
#' @param spaces An `fs_spaces` object (or a named list of `fs_space`).
#' @param n_axes Number of PCoA axes retained in each space. `NULL` (default,
#'   as in the original analysis) uses all positive-eigenvalue axes. A
#'   vector of length 2 or more is allowed to explore sensitivity, e.g.
#'   `n_axes = 4`.
#' @param method `"ade4"` or `"vegan"`.
#' @param nrepet Number of permutations.
#' @param mantel Logical; also run Mantel tests (Pearson) between distances.
#' @param seed Random seed.
#' @param progress Optional function `progress(i, n, label)`.
#' @param cores Number of CPU cores used for the permutations of
#'   `method = "ade4"` (forked processes, Unix/macOS only; ignored on
#'   Windows), when at least 50 axes are compared (below that, forking costs
#'   more than it saves). Results do not depend
#'   on `cores`: the permutations are drawn beforehand, in the same order as
#'   `ade4::procuste.rtest()`.
#'
#' @details With `method = "ade4"` the test is an exact re-implementation of
#'   [ade4::procuste.rtest()] (same statistic, same random permutations for a
#'   given seed, hence the same p-values) that computes singular values only,
#'   which is several times faster when many axes are used.
#'
#' @return An object of class `fs_comparison`: `table` (one row per pair:
#'   `space1`, `space2`, `n_species`, `n_axes`, `r` Procrustes correlation,
#'   `m2` (vegan only), `p_value`, `signif`, and Mantel statistics if
#'   requested) and `matrix` (symmetric matrix of Procrustes correlations).
#' @export
#' @examples
#' data(birds_traits); data(birds_trait_info)
#' tr <- drop_incomplete(prepare_traits(birds_traits, birds_trait_info, "species"))
#' spaces <- build_spaces(tr)
#' cmp <- compare_spaces(spaces, nrepet = 99)
#' cmp
compare_spaces <- function(spaces, n_axes = NULL, method = c("ade4", "vegan"),
                           nrepet = 999, mantel = FALSE, seed = 1, progress = NULL,
                           cores = 1) {
  method <- match.arg(method)
  if (length(spaces) < 2) stop("At least two spaces are needed.")
  nm <- names(spaces)
  pairs <- utils::combn(nm, 2, simplify = FALSE)
  set.seed(seed)
  rows <- list()
  for (i in seq_along(pairs)) {
    a <- spaces[[pairs[[i]][1]]]; b <- spaces[[pairs[[i]][2]]]
    if (!is.null(progress)) progress(i, length(pairs), paste(pairs[[i]], collapse = " vs "))
    al <- .align(a, b, n_axes)
    if (method == "ade4") {
      tst <- .procuste_rtest(al$X, al$Y, nrepet = nrepet, cores = cores)
      r <- tst$obs; p <- tst$pvalue; m2 <- NA_real_
    } else {
      tst <- vegan::protest(al$X, al$Y, permutations = nrepet, symmetric = TRUE)
      r <- tst$t0; p <- tst$signif; m2 <- tst$ss
    }
    row <- data.frame(space1 = pairs[[i]][1], space2 = pairs[[i]][2],
                      n_species = nrow(al$X),
                      n_axes = paste0(ncol(al$X), "/", ncol(al$Y)),
                      r = r, m2 = m2, p_value = p, signif = .stars(p),
                      stringsAsFactors = FALSE)
    if (mantel) {
      if (is.null(a$dist) || is.null(b$dist)) {
        row$mantel_r <- NA; row$mantel_p <- NA
      } else {
        da <- as.matrix(a$dist)[al$species, al$species]
        db <- as.matrix(b$dist)[al$species, al$species]
        mt <- vegan::mantel(stats::as.dist(da), stats::as.dist(db), permutations = nrepet)
        row$mantel_r <- mt$statistic; row$mantel_p <- mt$signif
      }
    }
    rows[[i]] <- row
  }
  tab <- do.call(rbind, rows)
  m <- matrix(1, length(nm), length(nm), dimnames = list(nm, nm))
  for (i in seq_len(nrow(tab))) {
    m[tab$space1[i], tab$space2[i]] <- m[tab$space2[i], tab$space1[i]] <- tab$r[i]
  }
  structure(list(table = tab, matrix = m,
                 settings = list(method = method, nrepet = nrepet,
                                 n_axes = if (is.null(n_axes)) "all" else n_axes)),
            class = "fs_comparison")
}

# Exact re-implementation of ade4::procuste.rtest() (ade4 1.7): same
# centring/scaling, statistic sum(singular values of X'Y), and permutations
# drawn with sample(n) one after the other, as in ade4, so that results are
# identical for a given seed. Faster: singular values only (svd(nu = 0,
# nv = 0)), crossprod(), optional forked parallelism over permutations.
.procuste_rtest <- function(X, Y, nrepet = 999, cores = 1) {
  X <- scale(as.matrix(X), scale = FALSE)
  Y <- scale(as.matrix(Y), scale = FALSE)
  X <- X / sqrt(sum(apply(X, 2, function(x) sum(x^2))))
  Y <- Y / sqrt(sum(apply(Y, 2, function(x) sum(x^2))))
  attr(X, "scaled:center") <- attr(Y, "scaled:center") <- NULL
  stat <- function(Yp) sum(svd(crossprod(X, Yp), nu = 0, nv = 0)$d)
  obs <- stat(Y)
  if (nrepet == 0) return(list(obs = obs, pvalue = NA_real_))
  n <- nrow(X)
  cores <- max(1L, as.integer(cores))
  # forking only pays off when each SVD is costly (many axes)
  big <- min(ncol(X), ncol(Y)) >= getOption("DisentFunSpaces.parallel_min_axes", 50)
  if (cores > 1L && big && .Platform$OS.type == "unix") {
    perms <- lapply(seq_len(nrepet), function(i) sample(n))
    chunks <- split(seq_len(nrepet), cut(seq_len(nrepet), min(cores * 4L, nrepet), labels = FALSE))
    sim <- unlist(parallel::mclapply(chunks, function(ii)
      vapply(ii, function(i) stat(Y[perms[[i]], , drop = FALSE]), numeric(1)),
      mc.cores = cores), use.names = FALSE)
  } else {
    sim <- vapply(seq_len(nrepet), function(i) stat(Y[sample(n), , drop = FALSE]), numeric(1))
  }
  ade4::as.randtest(obs = obs, sim = sim)
}

.align <- function(a, b, n_axes = NULL) {
  sp <- intersect(rownames(a$coords), rownames(b$coords))
  if (length(sp) < 3) stop("Spaces '", a$name, "' and '", b$name, "' share fewer than 3 species.")
  X <- a$coords[sp, , drop = FALSE]; Y <- b$coords[sp, , drop = FALSE]
  if (!is.null(n_axes)) {
    X <- X[, seq_len(min(n_axes, ncol(X))), drop = FALSE]
    Y <- Y[, seq_len(min(n_axes, ncol(Y))), drop = FALSE]
  }
  list(X = X, Y = Y, species = sp)
}

.stars <- function(p) {
  ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))
}

#' @export
print.fs_comparison <- function(x, digits = 3, ...) {
  cat("<fs_comparison> Procrustes tests (", x$settings$method, ", ", x$settings$nrepet,
      " permutations, axes: ", paste(x$settings$n_axes, collapse = ","), ")\n", sep = "")
  tab <- x$table
  num <- vapply(tab, is.numeric, logical(1))
  tab[num] <- lapply(tab[num], round, digits)
  print(tab, row.names = FALSE)
  invisible(x)
}

#' Species-level Procrustes residuals between two spaces
#'
#' After a symmetric Procrustes superimposition ([vegan::procrustes()]),
#' the residual of a species measures how differently it is positioned in
#' the two spaces: species with large residuals are those whose functional
#' role is perceived differently depending on the traits considered.
#'
#' @param space1,space2 Two `fs_space` objects.
#' @inheritParams compare_spaces
#' @return Data frame with `species`, `residual` and `rank` (1 = largest),
#'   sorted by decreasing residual; the `procrustes` fit is attached as
#'   attribute.
#' @export
procrustes_residuals <- function(space1, space2, n_axes = NULL) {
  al <- .align(space1, space2, n_axes)
  fit <- suppressWarnings(vegan::procrustes(al$X, al$Y, symmetric = TRUE))
  res <- stats::residuals(fit)
  out <- data.frame(species = al$species, residual = as.numeric(res), stringsAsFactors = FALSE)
  out <- out[order(-out$residual), ]
  out$rank <- seq_len(nrow(out))
  rownames(out) <- NULL
  attr(out, "procrustes") <- fit
  attr(out, "spaces") <- c(space1$name, space2$name)
  out
}
