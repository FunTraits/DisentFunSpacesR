#' Functional uniqueness and specialisation of species in a space
#'
#' Species-level metrics computed as in Toussaint et al. (2026):
#' * **FUn** (functional uniqueness): sum of the Euclidean distances to the
#'   `k` nearest neighbours in the space;
#' * **FSp** (functional specialisation): Euclidean distance to the centroid
#'   of the space.
#' Both are min-max rescaled to \[0, 1\] within the space.
#'
#' With coarse categorical traits many species can share exactly the same
#' trait profile: distances below `tol` are set to 0 so that such species
#' get exactly FUn = 0 (and identical FSp) instead of values driven by
#' numerical noise. `n_identical` gives the number of other species with the
#' same position in the space; when more than `k` species share a profile,
#' FUn is 0 for all of them.
#'
#' @param space An `fs_space`.
#' @param n_axes Number of PCoA axes used (`NULL` = all, as in the original
#'   analysis).
#' @param k Number of nearest neighbours for FUn (5 in the original
#'   analysis).
#' @param tol Distances below `tol` are treated as 0 (identical profiles).
#' @return Data frame with `species`, `FUn`, `FSp`, `FUn_raw`, `FSp_raw`,
#'   `n_identical`.
#' @export
species_metrics <- function(space, n_axes = NULL, k = 5, tol = 1e-8) {
  stopifnot(inherits(space, "fs_space"))
  X <- space$coords
  if (!is.null(n_axes)) X <- X[, seq_len(min(n_axes, ncol(X))), drop = FALSE]
  if (nrow(X) <= k) stop("Fewer species than k + 1.")
  cen <- colMeans(X)
  fsp <- sqrt(rowSums(sweep(X, 2, cen)^2))
  nn <- FNN::get.knnx(data = X, query = X, k = k + 1)
  nd <- nn$nn.dist[, -1, drop = FALSE]
  nd[nd < tol] <- 0                       # identical profiles: exact zeros, not numerical noise
  fun <- rowSums(nd)
  fsp <- round(fsp, -floor(log10(tol)))  # identical profiles share exactly the same FSp
  # number of other species with exactly the same position (functional redundancy)
  key <- do.call(paste, as.data.frame(round(X, -floor(log10(tol)) - 2)))
  n_identical <- as.integer(stats::ave(seq_along(key), key, FUN = length)) - 1L
  data.frame(species = rownames(X), FUn = .mm(fun), FSp = .mm(fsp),
             FUn_raw = fun, FSp_raw = fsp, n_identical = n_identical,
             row.names = NULL, stringsAsFactors = FALSE)
}

.mm <- function(x) {
  r <- range(x, na.rm = TRUE)
  if (!all(is.finite(r)) || diff(r) == 0) return(rep(0, length(x)))
  (x - r[1]) / diff(r)
}

#' Extinction probabilities associated with IUCN categories
#'
#' Fixed probabilities used to compute FUSE in Toussaint et al. (2026):
#' median extinction probabilities per IUCN category of the EDGE2 protocol
#' (Gumbs et al. 2023, <doi:10.1371/journal.pbio.3001991>); DD species are
#' given the VU value. FUSE itself follows Pimiento et al. (2020,
#' <doi:10.1126/sciadv.aay7650>).
#' @return Named numeric vector.
#' @export
iucn_probabilities <- function() {
  c(LC = 0.06, NT = 0.12, VU = 0.24, EN = 0.49, CR = 0.97, DD = 0.24)
}

#' Species metrics in every space, with optional FUSE
#'
#' @param spaces An `fs_spaces` object.
#' @inheritParams species_metrics
#' @param iucn Optional named character vector of IUCN categories (names =
#'   species) or a data frame with columns `species` and `iucn`. When
#'   supplied, FUSE = log(1 + P x FSp + P x FUn) is added, with P from
#'   `probabilities`.
#' @param probabilities Named vector of extinction probabilities per
#'   category.
#' @return Long data frame: `space`, `species`, `FUn`, `FSp` (+ `iucn`,
#'   `P`, `FUSE`).
#' @export
species_metrics_all <- function(spaces, n_axes = NULL, k = 5, iucn = NULL,
                                probabilities = iucn_probabilities()) {
  if (is.data.frame(iucn)) iucn <- stats::setNames(as.character(iucn$iucn), iucn$species)
  out <- do.call(rbind, lapply(names(spaces), function(nm) {
    m <- species_metrics(spaces[[nm]], n_axes = n_axes, k = k)
    cbind(space = nm, m, stringsAsFactors = FALSE)
  }))
  if (!is.null(iucn)) {
    out$iucn <- toupper(trimws(unname(iucn[out$species])))
    out$P <- unname(probabilities[out$iucn])
    out$FUSE <- log(1 + out$P * out$FSp + out$P * out$FUn)
  }
  out
}

#' Do different spaces identify the same distinct species?
#'
#' For a species metric computed in every space, returns pairwise Spearman
#' correlations and the overlap (Jaccard index) between the sets of the
#' `top` most extreme species of each space - the comparison displayed as
#' Venn diagrams in Toussaint et al. (2026).
#'
#' @param metrics Output of [species_metrics_all()].
#' @param metric Column to compare (`"FUn"`, `"FSp"` or `"FUSE"`).
#' @param top Proportion of species defining the top set.
#' @return A list with `pairs` (data frame: `space1`, `space2`, `spearman`,
#'   `jaccard_top`, `n_shared_top`), `top_sets` (list of species per
#'   space; species tied with the cutoff value are all included) and
#'   `top_info` (size of each top set; `informative = FALSE` when ties make
#'   it larger than twice the requested proportion).
#' @export
compare_metrics <- function(metrics, metric = "FUn", top = 0.1) {
  if (!metric %in% names(metrics)) stop("Metric '", metric, "' not found.")
  sps <- unique(metrics$space)
  wide <- lapply(sps, function(s) {
    m <- metrics[metrics$space == s, ]
    stats::setNames(m[[metric]], m$species)
  })
  names(wide) <- sps
  # top set: the ceiling(top * n) highest values, extended to all species tied
  # with the cutoff value (ties are frequent with categorical traits)
  top_sets <- lapply(wide, function(v) {
    v <- v[!is.na(v)]
    cut <- sort(v, decreasing = TRUE)[ceiling(length(v) * top)]
    names(v)[v >= cut - 1e-12]
  })
  top_info <- data.frame(space = names(top_sets),
                         n_species = vapply(wide, function(v) sum(!is.na(v)), integer(1)),
                         n_top = lengths(top_sets),
                         stringsAsFactors = FALSE)
  top_info$prop_top <- top_info$n_top / top_info$n_species
  top_info$informative <- top_info$prop_top <= 2 * top
  if (any(!top_info$informative))
    warning("Ties make the top-", 100 * top, "% set uninformative in: ",
            paste(top_info$space[!top_info$informative], collapse = ", "),
            " (", metric, " is identical for most species).", call. = FALSE)
  pairs <- utils::combn(sps, 2, simplify = FALSE)
  tab <- do.call(rbind, lapply(pairs, function(p) {
    common <- intersect(names(wide[[p[1]]]), names(wide[[p[2]]]))
    a <- top_sets[[p[1]]]; b <- top_sets[[p[2]]]
    data.frame(space1 = p[1], space2 = p[2],
               spearman = stats::cor(wide[[p[1]]][common], wide[[p[2]]][common],
                                     method = "spearman", use = "complete.obs"),
               jaccard_top = length(intersect(a, b)) / length(union(a, b)),
               n_shared_top = length(intersect(a, b)),
               stringsAsFactors = FALSE)
  }))
  list(pairs = tab, top_sets = top_sets, top_info = top_info, metric = metric, top = top)
}
