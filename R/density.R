## ---------------------------------------------------------------------------
## Kernel densities of groups of species on a grid of the functional space
##
## Convex hulls depend on the most extreme species only. Densities describe
## how the species of a group are distributed in the space: their overlap
## (sum over the grid of the minimum of the two probability densities; the
## TPD overlap of Carmona et al. 2016, <doi:10.1016/j.tree.2016.02.003>) is
## robust to outliers and can be weighted by abundance.
##
## Densities are binned on a regular grid of the first `axes` axes (2-4) and
## smoothed with a separable Gaussian kernel. The grid and the bandwidth are
## shared by all groups of a space (computed from all species), so densities
## of different groups are comparable.
## ---------------------------------------------------------------------------

## Internal: grid settings for the first d axes of X (all species of a space)
.kde_setup <- function(X, axes = 2, n_grid = NULL, bw = 1) {
  d <- min(axes, ncol(X), 4)
  if (d < 1) stop("No axis for the density.", call. = FALSE)
  X <- X[, seq_len(d), drop = FALSE]
  if (is.null(n_grid)) n_grid <- c(200, 100, 40, 20)[d]
  n <- nrow(X)
  sdv <- apply(X, 2, stats::sd)
  iqr <- apply(X, 2, stats::IQR) / 1.349
  s <- pmin(sdv, ifelse(iqr > 0, iqr, sdv))
  h <- bw * 1.06 * s * n^(-1 / (d + 4))          # Scott / Silverman rule, per axis
  h[!is.finite(h) | h <= 0] <- 1e-6
  lo <- apply(X, 2, min) - 3 * h; hi <- apply(X, 2, max) + 3 * h
  centres <- lapply(seq_len(d), function(j) seq(lo[j], hi[j], length.out = n_grid))
  K <- lapply(seq_len(d), function(j) {
    c <- centres[[j]]
    stats::dnorm(outer(c, c, "-") / h[j])
  })
  list(d = d, n_grid = n_grid, lo = lo, hi = hi, h = h, centres = centres, K = K)
}

## Internal: probability density (vector over the grid cells, sums to 1) of
## the points P (rows) with optional weights w
.kde_density <- function(P, setup, w = NULL) {
  d <- setup$d; g <- setup$n_grid
  P <- P[, seq_len(d), drop = FALSE]
  if (is.null(w)) w <- rep(1, nrow(P))
  w[!is.finite(w) | w < 0] <- 0
  if (!nrow(P) || sum(w) <= 0) return(rep(NA_real_, g^d))
  idx <- vapply(seq_len(d), function(j) {
    i <- floor((P[, j] - setup$lo[j]) / (setup$hi[j] - setup$lo[j]) * (g - 1) + 0.5) + 1
    pmin(pmax(i, 1), g)
  }, numeric(nrow(P)))
  if (!is.matrix(idx)) idx <- matrix(idx, nrow = nrow(P))
  lin <- 1 + as.vector((idx - 1) %*% (g^(seq_len(d) - 1)))
  v <- numeric(g^d)
  agg <- rowsum(w, lin)
  v[as.integer(rownames(agg))] <- agg[, 1]
  a <- array(v, rep(g, d))
  for (j in seq_len(d)) {                 # separable Gaussian smoothing, axis by axis
    perm <- c(j, setdiff(seq_len(d), j))
    b <- aperm(a, perm)
    b <- array(setup$K[[j]] %*% matrix(b, nrow = g), dim(b))
    a <- aperm(b, order(perm))
  }
  v <- as.vector(a)
  v / sum(v)
}

## Internal: share of the grid needed to hold `mass` of the density
.kde_cells <- function(p, mass = 0.99) {
  if (anyNA(p)) return(NA_real_)
  s <- sort(p, decreasing = TRUE)
  sum(cumsum(s) < mass) + 1
}

## Internal: overlap of two densities (1 = identical, 0 = disjoint)
.kde_overlap <- function(p1, p2) if (anyNA(p1) || anyNA(p2)) NA_real_ else sum(pmin(p1, p2))
