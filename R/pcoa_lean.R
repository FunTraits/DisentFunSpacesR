# Memory-lean principal coordinates analysis ---------------------------------
# Same algorithm and output as ape::pcoa(D, correction = "none") (Paradis &
# Schliep 2019; code adapted from ape, GPL >= 2), but the Gower double
# centring is done with row means instead of two n x n matrix products, and
# only the eigenvectors of positive eigenvalues are kept. Peak memory is
# about three n x n matrices instead of about ten, which makes PCoA of
# several thousand species possible on a laptop. Used by PCoAGraph() when
# correction = "none" and the number of species is at least
# getOption("DisentFunSpaces.lean_pcoa_n", 3000).
.pcoa_lean <- function(D) {
  bstick.def <- function(n, tot.var = 1) {
    res <- rev(cumsum(tot.var / n:1) / n)
    names(res) <- paste("Stick", seq(len = n), sep = "")
    res
  }
  names <- labels(D)
  invisible(gc(verbose = FALSE))
  D <- as.matrix(D)
  n <- nrow(D)
  epsilon <- sqrt(.Machine$double.eps)
  D <- D * D
  D <- D * -0.5
  rm <- rowMeans(D)
  gm <- mean(rm)
  D <- D - rm                                   # subtracts rm[i] from row i
  for (j in seq_len(n)) D[, j] <- D[, j] - (rm[j] - gm)
  trace <- sum(diag(D))
  D.eig <- eigen(D, symmetric = TRUE)
  rm(D); invisible(gc(verbose = FALSE))
  min.eig <- min(D.eig$values)
  zero.eig <- which(abs(D.eig$values) < epsilon)
  D.eig$values[zero.eig] <- 0
  eig <- D.eig$values
  if (min.eig > -epsilon) {
    k <- length(which(eig > epsilon))
    rel.eig <- eig[1:k] / trace
    cum.eig <- cumsum(rel.eig)
    vectors <- sweep(D.eig$vectors[, 1:k, drop = FALSE], 2, sqrt(eig[1:k]), FUN = "*")
    bs <- bstick.def(k)
    cum.bs <- cumsum(bs)
    res <- data.frame(eig[1:k], rel.eig, bs, cum.eig, cum.bs)
    colnames(res) <- c("Eigenvalues", "Relative_eig", "Broken_stick", "Cumul_eig", "Cumul_br_stick")
    note <- "There were no negative eigenvalues. No correction was applied"
  } else {
    k <- n
    rel.eig <- eig / trace
    rel.eig.cor <- (eig - min.eig) / (trace - (n - 1) * min.eig)
    if (length(zero.eig)) rel.eig.cor <- c(rel.eig.cor[-zero.eig[1]], 0)
    cum.eig.cor <- cumsum(rel.eig.cor)
    k2 <- length(which(eig > epsilon))
    k3 <- length(which(rel.eig.cor > epsilon))
    vectors <- sweep(D.eig$vectors[, 1:k2, drop = FALSE], 2, sqrt(eig[1:k2]), FUN = "*")
    note <- "No correction was applied to the negative eigenvalues"
    bs <- c(bstick.def(k3), rep(0, (k - k3)))
    cum.bs <- cumsum(bs)
    res <- data.frame(eig[1:k], rel.eig, rel.eig.cor, bs, cum.eig.cor, cum.bs)
    colnames(res) <- c("Eigenvalues", "Relative_eig", "Rel_corr_eig", "Broken_stick",
                       "Cum_corr_eig", "Cumul_br_stick")
  }
  rownames(res) <- 1:nrow(res)
  rownames(vectors) <- names
  colnames(vectors) <- paste0("Axis.", seq_len(ncol(vectors)))
  out <- list(correction = c("none", 1), note = note, values = res, vectors = vectors,
              trace = trace)
  class(out) <- "pcoa"
  out
}

# Partial principal coordinates analysis --------------------------------------
# Leading principal coordinates only, from the same Gower-centred matrix as
# ape::pcoa(), computed with the Lanczos solver of RSpectra (exact leading
# eigenpairs, to numerical precision) instead of a full eigendecomposition
# (O(n^3)): seconds instead of about ten minutes per space for ~8,000
# species. The trace (sum of all eigenvalues) is computed exactly, so the
# relative inertia of the axes is the same as with ape::pcoa(); the
# broken-stick model and the diagnostics of negative eigenvalues, which need
# the whole spectrum, are not available (NA). Signs of axes are arbitrary, as
# with any eigensolver.
#  * k = number of axes: the first k positive-eigenvalue axes;
#  * k = NULL: ALL positive-eigenvalue axes (as ape::pcoa), found by
#    computing 50 axes, then 100, 200... until a non-positive eigenvalue is
#    reached. Fails (error) beyond `k_max` positive axes, in which case the
#    caller falls back to the full decomposition. Gower distances of
#    continuous and fuzzy-coded traits without missing values are Euclidean
#    with about one positive axis per trait column, so that this is fast for
#    many trait datasets.
# Used by PCoAGraph(partial = TRUE) when 'RSpectra' is installed.
.pcoa_partial <- function(D, k = NULL,
                          k_max = getOption("DisentFunSpaces.partial_k_max", 500)) {
  names <- labels(D)
  D <- as.matrix(D)
  n <- nrow(D)
  eps <- sqrt(.Machine$double.eps)
  D <- D * D
  D <- D * -0.5
  rm <- rowMeans(D)
  gm <- mean(rm)
  D <- D - rm
  for (j in seq_len(n)) D[, j] <- D[, j] - (rm[j] - gm)
  trace <- sum(diag(D))
  solve_k <- function(kk) {
    e <- RSpectra::eigs_sym(D, k = kk, which = "LA", opts = list(maxitr = 5000, retvec = TRUE))
    if (is.null(e$nconv) || e$nconv < kk) stop("RSpectra did not converge.", call. = FALSE)
    o <- order(e$values, decreasing = TRUE)
    list(values = e$values[o], vectors = e$vectors[, o, drop = FALSE])
  }
  cap <- min(as.integer(k_max), floor(n / 4))
  if (is.null(k)) {
    kk <- min(50L, cap)
    repeat {
      e <- solve_k(kk)
      if (min(e$values) <= eps) break                 # all positive axes found
      if (kk >= cap) stop("More than ", cap, " positive eigenvalues.", call. = FALSE)
      kk <- min(2L * kk, cap)
    }
  } else {
    e <- solve_k(min(as.integer(k), n - 2L))
  }
  rm(D); invisible(gc(verbose = FALSE))
  pos <- e$values > eps
  eig <- e$values[pos]
  vectors <- sweep(e$vectors[, pos, drop = FALSE], 2, sqrt(eig), FUN = "*")
  rel.eig <- eig / trace
  res <- data.frame(Eigenvalues = eig, Relative_eig = rel.eig, Broken_stick = NA_real_,
                    Cumul_eig = cumsum(rel.eig), Cumul_br_stick = NA_real_)
  rownames(res) <- seq_len(nrow(res))
  rownames(vectors) <- names
  colnames(vectors) <- paste0("Axis.", seq_len(ncol(vectors)))
  out <- list(correction = c("none", 1),
              note = if (is.null(k))
                sprintf("Partial eigendecomposition (RSpectra): all %d positive axes", length(eig)) else
                sprintf("Partial eigendecomposition (RSpectra): first %d axes", length(eig)),
              values = res, vectors = vectors, trace = trace, partial = TRUE)
  class(out) <- "pcoa"
  out
}
