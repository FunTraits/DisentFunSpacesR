utils::globalVariables(c("a1", "a2", "position", "FUSI", "species", "space", "value"))

## ---------------------------------------------------------------------------
## Non-native species in functional spaces
##
## Counterpart of the FUSE analysis for threatened species: where do
## non-native species fall in each functional space relative to the native
## fauna, which native species lose functional uniqueness after introductions,
## and how much each functional space (dimension) is expanded or reshaped.
## ---------------------------------------------------------------------------

## Internal: status vector -> logical (TRUE = non-native), names = species
.as_nonnative <- function(status, species, native = "native") {
  if (is.data.frame(status)) {
    if (!all(c("species", "status") %in% names(status)))
      stop("A data frame 'status' needs columns 'species' and 'status'.")
    status <- stats::setNames(status$status, status$species)
  }
  if (is.null(names(status))) stop("'status' must be named by species.")
  miss <- setdiff(species, names(status))
  if (length(miss))
    stop(length(miss), " species of the space have no status (e.g. ",
         paste(utils::head(miss, 3), collapse = ", "), ").")
  s <- status[species]
  if (is.logical(s)) return(stats::setNames(s, species))
  s <- trimws(as.character(s))
  out <- !is.na(s) & tolower(s) != tolower(native)
  out[is.na(s)] <- NA
  stats::setNames(out, species)
}

## Internal: coordinates of an fs_space, first n_axes axes
.space_coords <- function(space, n_axes = NULL) {
  stopifnot(inherits(space, "fs_space"))
  X <- space$coords
  if (!is.null(n_axes)) X <- X[, seq_len(min(n_axes, ncol(X))), drop = FALSE]
  X
}

## Internal: convex hull volume (NA when not computable)
.hull_vol <- function(X) {
  if (!requireNamespace("geometry", quietly = TRUE))
    stop("Package 'geometry' is needed for convex hulls: install.packages('geometry').")
  X <- unique(round(X, 10))
  if (nrow(X) <= ncol(X)) return(NA_real_)
  tryCatch(geometry::convhulln(X, options = "FA")$vol, error = function(e) NA_real_)
}

## Internal: are points P inside the convex hull of X? (NA when not computable)
## A point is outside when adding it increases the hull volume: exact and
## robust for points lying on the hull boundary (frequent with discrete traits).
.in_hull <- function(X, P, rel_tol = 1e-7) {
  v0 <- .hull_vol(X)
  if (is.na(v0)) return(rep(NA, nrow(P)))
  # Fast path: the hull of X is computed once and each point is tested
  # against its facets (half-space test, signed distance to every facet
  # <= tol). Points on the boundary (within a relative tolerance of the
  # extent of the data) are inside, as with the volume criterion below.
  Xu <- unique(round(X, 10))
  h <- tryCatch(geometry::convhulln(Xu, options = "Tv", output.options = "n"),
                error = function(e) NULL)
  if (!is.null(h) && !is.null(h$normals) && ncol(h$normals) == ncol(X) + 1) {
    tol <- rel_tol * max(apply(Xu, 2, function(v) diff(range(v))))
    sd <- cbind(round(P, 10), 1) %*% t(h$normals)
    return(apply(sd, 1, max) <= tol)
  }
  vapply(seq_len(nrow(P)), function(i) {
    v1 <- .hull_vol(rbind(X, P[i, , drop = FALSE]))
    if (is.na(v1)) NA else v1 <= v0 * (1 + rel_tol)
  }, logical(1))
}

#' Weights associated with EICAT impact categories
#'
#' Default weights used to compute FUSI in [nonnative_species()]. The
#' categories are those of the IUCN Environmental Impact Classification for
#' Alien Taxa (EICAT; Blackburn et al. 2014, <doi:10.1371/journal.pbio.1001850>;
#' IUCN 2020): MC (minimal concern), MN (minor), MO (moderate), MR (major),
#' MV (massive). Unlike the extinction probabilities used for FUSE
#' ([iucn_probabilities()]), these values are **not** empirical estimates:
#' they are equally spaced ranks rescaled to \[0.2, 1\]. Supply your own
#' weights (e.g. invasion stage, establishment probability, abundance or
#' spread) when a justified scale is available.
#' @return Named numeric vector.
#' @export
eicat_weights <- function() {
  c(MC = 0.2, MN = 0.4, MO = 0.6, MR = 0.8, MV = 1)
}

#' Position of non-native species relative to the native fauna
#'
#' For one functional space, describes (i) each non-native species relative
#' to the native species and (ii) the loss of functional uniqueness of each
#' native species once non-native species are added. This is the
#' non-native counterpart of [species_metrics()] / FUSE.
#'
#' Non-native species:
#' * `FNov_raw`: distance to the nearest native species (functional novelty;
#'   the trait-based analogue of the phylogenetic distance used to test
#'   Darwin's naturalisation hypothesis). `FNov` is the same, min-max rescaled
#'   over non-native species.
#' * `mean_dist_native`: mean distance to all native species.
#' * `FSp_native_raw`: distance to the centroid of the native species
#'   (`FSp_native` rescaled), i.e. how peripheral the species is with respect
#'   to the native fauna.
#' * `inside_native_hull`: is the species inside the convex hull of the
#'   native species (first `hull_axes` axes)? `FALSE` = the species expands
#'   the native functional space (adding it increases the hull volume).
#' * `nearest_native`: name of the nearest native species.
#' * `FUSI` (when `impact` is supplied): Functionally Unique, Specialised and
#'   Impactful, `log(1 + W x FNov + W x FSp_native)`, built as FUSE
#'   (Pimiento et al. 2020) with the impact weight `W` replacing the
#'   extinction probability.
#'
#' Native species:
#' * `FUn_before_raw`, `FUn_after_raw`: FUn (sum of distances to the
#'   `k` nearest neighbours, as in [species_metrics()]) computed among native
#'   species only, then among native + non-native species.
#' * `dFUn_rel`: relative loss of uniqueness, `(after - before) / before`
#'   (0 = unaffected, -1 = a non-native species now occupies exactly the
#'   same position); natives with high losses are those for which a
#'   non-native species is functionally redundant (potential competition or
#'   functional replacement).
#' * `dist_nearest_native`, `dist_nearest_nonnative`: distance to the
#'   nearest other native and to the nearest non-native species.
#' * `n_nonnative_knn`: number of non-native species among the `k` nearest
#'   neighbours after introduction; `nearest_nonnative`: nearest non-native
#'   species.
#'
#' @param space An `fs_space`.
#' @param status Named vector (names = species) giving the status of every
#'   species of the space: character (`native` label vs. anything else, e.g.
#'   `"non-native"`, `"introduced"`, `"established"`) or logical
#'   (`TRUE` = non-native). A data frame with columns `species` and `status`
#'   is also accepted.
#' @param native Label of native species in `status`.
#' @inheritParams species_metrics
#' @param hull_axes Number of axes used for convex hulls (hull computations
#'   become unstable above ~6 axes or with few species).
#' @param impact Optional named vector (names = non-native species): EICAT
#'   categories (character, matched to `impact_weights`) or numeric weights
#'   in \[0, 1\]. Adds `impact`, `W` and `FUSI`.
#' @param impact_weights Named vector of weights per category, used when
#'   `impact` is character.
#' @return An object of class `fs_nonnative`: a list with `nonnative` and
#'   `native` data frames and the settings.
#' @export
nonnative_species <- function(space, status, native = "native", n_axes = NULL,
                              k = 5, hull_axes = 4, tol = 1e-8, impact = NULL,
                              impact_weights = eicat_weights()) {
  X <- .space_coords(space, n_axes)
  sp <- rownames(X)
  nn_status <- .as_nonnative(status, sp, native)
  keep <- !is.na(nn_status)
  X <- X[keep, , drop = FALSE]; nn_status <- nn_status[keep]
  Xn <- X[!nn_status, , drop = FALSE]
  Xi <- X[nn_status, , drop = FALSE]
  if (nrow(Xi) == 0) stop("No non-native species in 'status'.")
  if (nrow(Xn) <= k) stop("Fewer native species than k + 1.")

  ## --- non-native species vs natives
  nn1 <- FNN::get.knnx(data = Xn, query = Xi, k = 1)
  fnov <- nn1$nn.dist[, 1]; fnov[fnov < tol] <- 0
  cen <- colMeans(Xn)
  fspn <- sqrt(rowSums(sweep(Xi, 2, cen)^2))
  mdn <- vapply(seq_len(nrow(Xi)), function(i)
    mean(sqrt(colSums((t(Xn) - Xi[i, ])^2))), numeric(1))
  ha <- min(hull_axes, ncol(X))
  inside <- .in_hull(Xn[, seq_len(ha), drop = FALSE], Xi[, seq_len(ha), drop = FALSE])
  nonnat <- data.frame(species = rownames(Xi),
                       FNov = .mm(fnov), FSp_native = .mm(fspn),
                       FNov_raw = fnov, FSp_native_raw = fspn,
                       mean_dist_native = mdn, inside_native_hull = inside,
                       nearest_native = rownames(Xn)[nn1$nn.index[, 1]],
                       row.names = NULL, stringsAsFactors = FALSE)
  if (!is.null(impact)) {
    imp <- impact[nonnat$species]
    W <- if (is.numeric(impact)) unname(imp) else
      unname(impact_weights[toupper(trimws(as.character(imp)))])
    nonnat$impact <- unname(imp)
    nonnat$W <- W
    nonnat$FUSI <- log(1 + W * nonnat$FNov + W * nonnat$FSp_native)
  }

  ## --- native species: uniqueness before vs after introductions
  kb <- FNN::get.knnx(data = Xn, query = Xn, k = k + 1)
  db <- kb$nn.dist[, -1, drop = FALSE]; db[db < tol] <- 0
  ka <- FNN::get.knnx(data = X, query = Xn, k = k + 1)
  ## drop the species itself (first neighbour at distance 0 = itself or a duplicate)
  self <- ka$nn.index == match(rownames(Xn), rownames(X))
  da <- t(vapply(seq_len(nrow(Xn)), function(i) {
    d <- ka$nn.dist[i, ]; s <- self[i, ]
    if (any(s)) d[-which(s)[1]] else d[-1]
  }, numeric(k)))
  if (k == 1) da <- matrix(da, ncol = 1)
  da[da < tol] <- 0
  idx_a <- t(vapply(seq_len(nrow(Xn)), function(i) {
    id <- ka$nn.index[i, ]; s <- self[i, ]
    if (any(s)) id[-which(s)[1]] else id[-1]
  }, numeric(k)))
  if (k == 1) idx_a <- matrix(idx_a, ncol = 1)
  fb <- rowSums(db); fa <- rowSums(da)
  nn_i <- FNN::get.knnx(data = Xi, query = Xn, k = 1)
  nat <- data.frame(species = rownames(Xn),
                    FUn_before_raw = fb, FUn_after_raw = fa,
                    dFUn_rel = ifelse(fb > 0, (fa - fb) / fb, 0),
                    n_nonnative_knn = rowSums(matrix(nn_status[idx_a], ncol = k)),
                    nearest_nonnative = rownames(Xi)[nn_i$nn.index[, 1]],
                    dist_nearest_nonnative = nn_i$nn.dist[, 1],
                    dist_nearest_native = db[, 1],
                    row.names = NULL, stringsAsFactors = FALSE)
  structure(list(nonnative = nonnat, native = nat,
                 settings = list(n_axes = n_axes, k = k, hull_axes = ha,
                                 n_native = nrow(Xn), n_nonnative = nrow(Xi))),
            class = "fs_nonnative")
}

#' @export
print.fs_nonnative <- function(x, ...) {
  s <- x$settings
  cat("Non-native species in a functional space\n")
  cat(sprintf("  %d native, %d non-native species; k = %d; hull on %d axes\n",
              s$n_native, s$n_nonnative, s$k, s$hull_axes))
  cat(sprintf("  non-native species outside the native hull: %d (%.0f%%)\n",
              sum(!x$nonnative$inside_native_hull, na.rm = TRUE),
              100 * mean(!x$nonnative$inside_native_hull, na.rm = TRUE)))
  cat(sprintf("  median relative loss of FUn of natives: %.1f%%\n",
              100 * stats::median(x$native$dFUn_rel)))
  invisible(x)
}

## Internal: assemblage-level metrics for a set of coordinates
.assemblage_metrics <- function(X, ha, vol_pool) {
  n <- nrow(X)
  if (n < 2) return(c(S = n, FRic = NA, FDis = NA, FOri = NA, MPD = NA))
  cen <- colMeans(X)
  # FOri: nearest-neighbour distances (exact k-NN search, no n x n matrix);
  # MPD: mean of the n(n-1)/2 pairwise distances (half the memory of the
  # full matrix; computed by blocks for very large assemblages)
  nn <- FNN::get.knn(X, k = 1, algorithm = if (ncol(X) > 10) "brute" else "kd_tree")
  c(S = n,
    FRic = .hull_vol(X[, seq_len(ha), drop = FALSE]) / vol_pool,
    FDis = mean(sqrt(rowSums(sweep(X, 2, cen)^2))),
    FOri = mean(nn$nn.dist[, 1]),
    MPD = .mean_pairwise_dist(X))
}

## Internal: "permute" null model. The current assemblage is the same in
## every replicate and the historical one is the current minus a random set R
## of species: its distance matrix is computed once, and MPD and FOri of each
## replicate are updated from it (MPD: sums of distances minus the rows of R;
## FOri: nearest neighbours recomputed only for species whose nearest
## neighbour is in R). Same values as .assemblage_metrics(), without
## recomputing O(n^2) distances in every replicate.
.permute_cache <- function(Xc) {
  D <- as.matrix(stats::dist(Xc))
  rs <- rowSums(D)
  diag(D) <- Inf
  list(D = D, rs = rs, S = sum(rs), nn = max.col(-D, ties.method = "first"))
}
.permute_metrics <- function(pc, R, Xh, ha, vol_pool) {
  n <- nrow(pc$D); h <- n - length(R)
  if (h < 2) return(.assemblage_metrics(Xh, ha, vol_pool))
  keep <- rep(TRUE, n); keep[R] <- FALSE
  dRR <- pc$D[R, R, drop = FALSE]; diag(dRR) <- 0
  mpd <- (pc$S - 2 * sum(pc$rs[R]) + sum(dRR)) / (h * (h - 1))
  idx <- which(keep)
  nnd <- pc$D[cbind(idx, pc$nn[idx])]
  lost <- !keep[pc$nn[idx]]                         # nearest neighbour removed
  if (any(lost)) {
    sub <- pc$D[idx[lost], idx, drop = FALSE]
    nnd[lost] <- sub[cbind(seq_len(nrow(sub)), max.col(-sub, ties.method = "first"))]
  }
  cen <- colMeans(Xh)
  c(S = h,
    FRic = .hull_vol(Xh[, seq_len(ha), drop = FALSE]) / vol_pool,
    FDis = mean(sqrt(rowSums(sweep(Xh, 2, cen)^2))),
    FOri = mean(nnd),
    MPD = mpd)
}

## Internal: mean Euclidean distance between all pairs of rows of X
.mean_pairwise_dist <- function(X, max_n = 10000, block = 1000) {
  n <- nrow(X)
  if (n <= max_n) return(mean(stats::dist(X)))
  tot <- 0
  for (s in seq(1, n - 1, by = block)) {
    i <- s:min(s + block - 1, n)
    j <- (s + 1):n
    D2 <- matrix(0, length(i), length(j))
    for (c in seq_len(ncol(X))) D2 <- D2 + outer(X[i, c], X[j, c], "-")^2
    D2[outer(i, j, ">=")] <- 0                     # keep pairs i < j only
    tot <- tot + sum(sqrt(D2))
  }
  tot / (n * (n - 1) / 2)
}

#' Change in functional structure of assemblages after introductions
#'
#' For each assemblage (or for the whole species set when `comm = NULL`),
#' compares the native ("historical") assemblage with the current assemblage
#' (natives - extirpated natives + non-native species) in one functional
#' space, following the approach of Toussaint et al. (2018,
#' <doi:10.1111/ele.13141>) for world freshwater fish faunas.
#'
#' Metrics (computed in the first `n_axes` axes; hulls in the first
#' `hull_axes` axes): `S` richness, `FRic` convex hull volume relative to the
#' hull of all species of the space (Villéger et al. 2008), `FDis` mean
#' distance to the assemblage centroid (unweighted functional dispersion),
#' `FOri` mean nearest-neighbour distance (functional originality; Mouillot
#' et al. 2013), `MPD` mean pairwise functional distance. For each metric
#' the output gives the historical (`_hist`) and current (`_cur`) values and
#' the change `d_` (current - historical); `p_nn_outside` is the proportion
#' of non-native species outside the native hull.
#'
#' Null models (`null`):
#' * `"pool"`: the observed non-native species are replaced by the same
#'   number of species drawn at random from `pool` (default: all species of
#'   the space that are not native to the assemblage). Tests whether
#'   introduced species change functional structure more than a random set
#'   of colonists from the pool. Requires a pool wider than the observed
#'   non-native species (e.g. the regional or global species pool with
#'   `comm`); with a single assemblage and no `pool`, use `"permute"`.
#' * `"permute"`: statuses are permuted within the current assemblage
#'   (same numbers of natives and non-natives). Tests whether non-native
#'   species are more functionally distinct than the same number of
#'   randomly chosen species of the assemblage.
#' For each `d_` metric, `ses_` = (obs - mean null) / sd null and `p_` is the
#' two-sided rank-based p-value.
#'
#' @param space An `fs_space`.
#' @param status Status of the species: named vector (as in
#'   [nonnative_species()]) or, when a species can be native in one
#'   assemblage and non-native in another, a matrix/data frame with the same
#'   dimensions and dimnames as `comm`.
#' @param comm Optional assemblage x species presence (0/1) matrix. `NULL` =
#'   one assemblage with all species having a status.
#' @param extirpated Optional extirpated native species: character vector of
#'   species (all assemblages), or logical matrix with the dimensions of
#'   `comm`. Extirpated natives belong to the historical but not to the
#'   current assemblage.
#' @inheritParams nonnative_species
#' @param null `"none"`, `"pool"` or `"permute"`.
#' @param nrep Number of null-model replicates.
#' @param pool Character vector of candidate colonists for `null = "pool"`.
#' @param seed Optional random seed.
#' @return Data frame, one row per assemblage.
#' @export
nonnative_assemblage <- function(space, status, comm = NULL, extirpated = NULL,
                                 native = "native", n_axes = NULL,
                                 hull_axes = 4,
                                 null = c("none", "pool", "permute"),
                                 nrep = 199, pool = NULL, seed = NULL) {
  null <- match.arg(null)
  if (!is.null(seed)) set.seed(seed)
  X <- .space_coords(space, n_axes)
  sp_all <- rownames(X)
  ha <- min(hull_axes, ncol(X))
  vol_pool <- .hull_vol(X[, seq_len(ha), drop = FALSE])

  if (is.null(comm)) {
    st <- .as_nonnative(status, sp_all, native)
    comm <- matrix(as.integer(!is.na(st)), nrow = 1,
                   dimnames = list("all", sp_all))
    status <- stats::setNames(st[!is.na(st)], sp_all[!is.na(st)])
  }
  comm <- as.matrix(comm)
  miss <- setdiff(colnames(comm)[colSums(comm > 0) > 0], sp_all)
  if (length(miss))
    stop(length(miss), " species of 'comm' are not in the space (e.g. ",
         paste(utils::head(miss, 3), collapse = ", "), ").")
  per_site <- is.matrix(status) || is.data.frame(status) && !all(c("species", "status") %in% names(status))
  if (per_site) status <- as.matrix(status)
  if (is.character(extirpated)) {
    ext_vec <- extirpated; extirpated <- NULL
  } else ext_vec <- character(0)

  dnames <- c("S", "FRic", "FDis", "FOri", "MPD")
  rows <- lapply(seq_len(nrow(comm)), function(j) {
    present <- colnames(comm)[comm[j, ] > 0]
    s_j <- if (per_site) {
      v <- status[j, present]; names(v) <- present
      .as_nonnative(v, present, native)
    } else .as_nonnative(status, present, native)
    nat <- present[!s_j %in% TRUE & !is.na(s_j)]
    nn  <- present[s_j %in% TRUE]
    ext <- if (!is.null(extirpated)) present[as.logical(extirpated[j, present])] else
      intersect(ext_vec, nat)
    ext <- intersect(ext, nat)
    hist <- nat
    cur <- c(setdiff(nat, ext), nn)
    mh <- .assemblage_metrics(X[hist, , drop = FALSE], ha, vol_pool)
    mc <- .assemblage_metrics(X[cur, , drop = FALSE], ha, vol_pool)
    out <- c(n_native = length(nat), n_nonnative = length(nn),
             n_extirpated = length(ext),
             stats::setNames(mh, paste0(dnames, "_hist")),
             stats::setNames(mc, paste0(dnames, "_cur")),
             stats::setNames(mc - mh, paste0("d_", dnames)))
    out["p_nn_outside"] <- if (length(nn) && length(nat) > ha)
      mean(!.in_hull(X[nat, seq_len(ha), drop = FALSE], X[nn, seq_len(ha), drop = FALSE]))
      else NA
    if (null != "none" && length(nn) > 0) {
      keep_nat <- setdiff(nat, ext)
      cand <- if (is.null(pool)) setdiff(sp_all, nat) else setdiff(intersect(pool, sp_all), nat)
      if (null == "pool" && length(cand) <= length(nn))
        warning("Assemblage '", rownames(comm)[j], "': the pool contains no species other than ",
                "the observed non-native species; the 'pool' null model is uninformative ",
                "(supply a wider 'pool' or use null = \"permute\").", call. = FALSE)
      pcache <- if (null == "permute" && !length(ext) && length(cur) >= 3 &&
                    length(cur) <= getOption("DisentFunSpaces.max_n_dist_matrix", 10000))
        .permute_cache(X[cur, , drop = FALSE]) else NULL
      sims <- replicate(nrep, {
        if (null == "pool") {
          if (length(cand) < length(nn)) return(rep(NA_real_, length(dnames)))
          cur_r <- c(keep_nat, sample(cand, length(nn)))
          hist_r <- hist
        } else {
          ## same current assemblage; a random subset of the same size plays
          ## the role of the non-native species
          fake_nn <- sample(cur, length(nn))
          cur_r <- cur
          hist_r <- c(setdiff(cur, fake_nn), ext)
        }
        m_cur <- if (null == "permute") mc else .assemblage_metrics(X[cur_r, , drop = FALSE], ha, vol_pool)
        m_hist <- if (null == "pool") mh else if (!is.null(pcache))
          .permute_metrics(pcache, match(fake_nn, cur), X[hist_r, , drop = FALSE], ha, vol_pool) else
          .assemblage_metrics(X[hist_r, , drop = FALSE], ha, vol_pool)
        m_cur - m_hist
      })
      sims <- matrix(sims, nrow = length(dnames), dimnames = list(dnames, NULL))
      for (m in dnames[-1]) {
        obs <- out[paste0("d_", m)]; v <- sims[m, ]; v <- v[is.finite(v)]
        out[paste0("ses_", m)] <- if (length(v) > 1 && stats::sd(v) > 0)
          (obs - mean(v)) / stats::sd(v) else NA
        out[paste0("p_", m)] <- if (length(v) > 1 && is.finite(obs))
          min(1, 2 * min((sum(v <= obs) + 1) / (length(v) + 1),
                         (sum(v >= obs) + 1) / (length(v) + 1))) else NA
      }
    }
    out
  })
  nm <- unique(unlist(lapply(rows, names)))
  res <- as.data.frame(do.call(rbind, lapply(rows, function(r) r[nm])))
  names(res) <- nm
  res$FRic_change_rel <- ifelse(res$FRic_hist > 0, res$d_FRic / res$FRic_hist, NA)
  cbind(assemblage = rownames(comm), res, row.names = NULL, stringsAsFactors = FALSE)
}

#' Non-native species in every functional space
#'
#' Runs [nonnative_species()] and, optionally, [nonnative_assemblage()] in
#' each space of an `fs_spaces` object, so that the functional dimensions
#' most affected by introductions can be identified - the non-native
#' counterpart of [species_metrics_all()].
#'
#' @param spaces An `fs_spaces` object (or a named list of `fs_space`).
#' @inheritParams nonnative_species
#' @inheritParams nonnative_assemblage
#' @param assemblage Also compute assemblage-level changes?
#' @param ... Passed to [nonnative_assemblage()] (`comm`, `extirpated`,
#'   `null`, `nrep`, `pool`, `seed`).
#' @return An object of class `fs_nonnative_all`: list with long data frames
#'   `nonnative`, `native` (column `space` added) and, when requested,
#'   `assemblage`.
#' @export
nonnative_all <- function(spaces, status, native = "native", n_axes = NULL,
                          k = 5, hull_axes = 4, impact = NULL,
                          impact_weights = eicat_weights(),
                          assemblage = TRUE, ...) {
  res <- lapply(names(spaces), function(nm) {
    r <- nonnative_species(spaces[[nm]], status, native = native, n_axes = n_axes,
                           k = k, hull_axes = hull_axes, impact = impact,
                           impact_weights = impact_weights)
    a <- if (assemblage) nonnative_assemblage(spaces[[nm]], status, native = native,
                                              n_axes = n_axes, hull_axes = hull_axes, ...)
    list(nonnative = cbind(space = nm, r$nonnative, stringsAsFactors = FALSE),
         native = cbind(space = nm, r$native, stringsAsFactors = FALSE),
         assemblage = if (!is.null(a)) cbind(space = nm, a, stringsAsFactors = FALSE))
  })
  out <- list(nonnative = do.call(rbind, lapply(res, `[[`, "nonnative")),
              native = do.call(rbind, lapply(res, `[[`, "native")),
              assemblage = do.call(rbind, lapply(res, `[[`, "assemblage")))
  rownames(out$nonnative) <- rownames(out$native) <- NULL
  if (!is.null(out$assemblage)) rownames(out$assemblage) <- NULL
  structure(out, class = "fs_nonnative_all")
}

#' @export
print.fs_nonnative_all <- function(x, ...) {
  cat("Non-native species across functional spaces\n")
  print(summary_nonnative(x), row.names = FALSE, digits = 3)
  invisible(x)
}

#' Which functional spaces are most affected by non-native species?
#'
#' One row per space: proportion of non-native species outside the native
#' hull (`p_outside_hull`), median novelty (`FNov_median`), mean novelty
#' relative to the mean nearest-neighbour distance among natives of the same
#' space (`FNov_rel`; > 1 = non-native species are on average further from
#' natives than natives are from each other; comparable across spaces),
#' median relative loss of native uniqueness (`dFUn_native_median`), the
#' proportion of natives losing at least 10% of their uniqueness
#' (`p_native_lose_10pc`) and, when
#' assemblage results are present, median relative change in FRic and
#' median SES across assemblages.
#'
#' @param x Output of [nonnative_all()].
#' @return Data frame.
#' @export
summary_nonnative <- function(x) {
  sps <- unique(x$nonnative$space)
  do.call(rbind, lapply(sps, function(s) {
    a <- x$nonnative[x$nonnative$space == s, ]
    b <- x$native[x$native$space == s, ]
    ## reference scale: mean nearest-neighbour distance among natives
    ref <- mean(b$dist_nearest_native, na.rm = TRUE)
    row <- data.frame(space = s,
                      p_outside_hull = mean(!a$inside_native_hull, na.rm = TRUE),
                      FNov_median = stats::median(a$FNov_raw),
                      FNov_rel = if (ref > 0) mean(a$FNov_raw) / ref else NA,
                      dFUn_native_median = stats::median(b$dFUn_rel),
                      p_native_lose_10pc = mean(b$dFUn_rel <= -0.1),
                      stringsAsFactors = FALSE)
    if (!is.null(x$assemblage)) {
      z <- x$assemblage[x$assemblage$space == s, ]
      row$FRic_change_rel_median <- stats::median(z$FRic_change_rel, na.rm = TRUE)
      if ("ses_FRic" %in% names(z)) {
        row$ses_FRic_median <- stats::median(z$ses_FRic, na.rm = TRUE)
        row$ses_FDis_median <- stats::median(z$ses_FDis, na.rm = TRUE)
        row$ses_FOri_median <- stats::median(z$ses_FOri, na.rm = TRUE)
      }
    }
    row
  }))
}

#' Do different spaces identify the same functionally novel non-natives?
#'
#' Wrapper around [compare_metrics()] for non-native species metrics
#' (`"FNov"`, `"FSp_native"`, `"FUSI"`) or native losses (`"dFUn_rel"`,
#' sign reversed so that the top set = natives losing most uniqueness).
#'
#' @param x Output of [nonnative_all()].
#' @param metric Metric to compare.
#' @param top Proportion of species defining the top set.
#' @return A list, as [compare_metrics()]: `pairs` (data frame: `space1`,
#'   `space2`, `spearman`, `jaccard_top`, `n_shared_top`), `top_sets`
#'   (species of the top set of each space), `top_info`, `metric` and `top`.
#' @export
compare_nonnative <- function(x, metric = "FNov", top = 0.1) {
  if (metric == "dFUn_rel") {
    m <- x$native; m$dFUn_loss <- -m$dFUn_rel
    return(compare_metrics(m, "dFUn_loss", top = top))
  }
  compare_metrics(x$nonnative, metric, top = top)
}

#' Plot native and non-native species in a functional space
#'
#' Native species in grey with their convex hull (on the plotted axes),
#' non-native species coloured by whether they are inside the native hull
#' in the first `hull_axes` axes; point size can show FUSI.
#'
#' @param space An `fs_space`.
#' @param status As in [nonnative_species()].
#' @param axes Two axes to plot.
#' @param x Optional output of [nonnative_species()] for the same space
#'   (otherwise computed).
#' @inheritParams nonnative_species
#' @param label Number of non-native species with highest `FNov` to label.
#' @return A ggplot.
#' @export
plot_nonnative <- function(space, status, axes = c(1, 2), x = NULL,
                           native = "native", hull_axes = 4, label = 0) {
  if (is.null(x)) x <- nonnative_species(space, status, native = native,
                                         hull_axes = hull_axes)
  X <- space$coords
  st <- .as_nonnative(status, rownames(X), native)
  d <- data.frame(species = rownames(X), a1 = X[, axes[1]], a2 = X[, axes[2]],
                  nonnative = st, stringsAsFactors = FALSE)
  d <- d[!is.na(d$nonnative), ]
  nat <- d[!d$nonnative, ]
  hull <- nat[grDevices::chull(nat$a1, nat$a2), ]
  nn <- merge(d[d$nonnative, ], x$nonnative, by = "species")
  nn$position <- ifelse(is.na(nn$inside_native_hull), "unknown",
                        ifelse(nn$inside_native_hull, "inside native space",
                               "outside native space"))
  ax <- paste0("PCoA ", axes)
  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(data = hull, ggplot2::aes(a1, a2),
                          fill = "grey85", colour = "grey60", alpha = 0.6) +
    ggplot2::geom_point(data = nat, ggplot2::aes(a1, a2),
                        colour = "grey45", size = 1, alpha = 0.6)
  if ("FUSI" %in% names(nn)) {
    p <- p + ggplot2::geom_point(data = nn, ggplot2::aes(a1, a2,
                                 colour = position, size = FUSI), alpha = 0.85) +
      ggplot2::scale_size_area(max_size = 5, name = "FUSI")
  } else {
    p <- p + ggplot2::geom_point(data = nn, ggplot2::aes(a1, a2,
                                 colour = position), size = 2, alpha = 0.85)
  }
  if (label > 0) {
    lab <- utils::head(nn[order(-nn$FNov), ], label)
    p <- p + ggplot2::geom_text(data = lab, ggplot2::aes(a1, a2,
                                label = species), size = 2.8, vjust = -0.8)
  }
  p + ggplot2::scale_colour_manual(values = c("inside native space" = "#2C7FB8",
                                              "outside native space" = "#D7301F",
                                              "unknown" = "grey20"), name = NULL) +
    ggplot2::labs(x = ax[1], y = ax[2],
                  subtitle = sprintf("Hull position assessed on the first %d axes",
                                     x$settings$hull_axes)) +
    ggplot2::coord_equal() + ggplot2::theme_classic(base_size = 11)
}

#' Plot the effect of non-native species across functional spaces
#'
#' Distribution, across assemblages, of the relative change in FRic or of a
#' standardised effect size (`ses_FRic`, `ses_FDis`, `ses_FOri`, `ses_MPD`)
#' in each space.
#'
#' @param x Output of [nonnative_all()] with assemblage results.
#' @param metric Column of `x$assemblage`.
#' @return A ggplot.
#' @export
plot_nonnative_change <- function(x, metric = "FRic_change_rel") {
  a <- x$assemblage
  if (is.null(a) || !metric %in% names(a)) stop("Metric '", metric, "' not available.")
  a$value <- a[[metric]]
  a$space <- factor(a$space, levels = unique(a$space))
  ref <- if (grepl("^ses_", metric)) c(-1.96, 1.96) else numeric(0)
  p <- ggplot2::ggplot(a, ggplot2::aes(space, value)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey50")
  if (length(ref)) p <- p + ggplot2::geom_hline(yintercept = ref, linetype = "dashed",
                                                colour = "grey60")
  p <- if (nrow(a) > length(unique(a$space))) p + ggplot2::geom_boxplot(outlier.size = 0.8, fill = "grey90") else
    p + ggplot2::geom_point(size = 3)
  p + ggplot2::labs(x = NULL, y = metric) + ggplot2::theme_classic(base_size = 11)
}

#' Simulate introductions of non-native species
#'
#' Assigns a non-native status to a proportion of the species of a
#' functional space, at random or biased towards the periphery or the centre
#' of the space (distance to centroid in the first `n_axes` axes). Intended
#' for simulation studies and for the built-in examples, as
#' [simulate_traits()].
#'
#' @param space An `fs_space` (the status is assigned from its coordinates).
#' @param prop Proportion of species set as non-native.
#' @param bias `"random"`, `"peripheral"` or `"central"`.
#' @param strength Strength of the bias (exponent on the rank of the
#'   distance to the centroid).
#' @param n_axes Number of axes used for the distance to the centroid.
#' @param seed Optional random seed.
#' @return Named character vector `"native"` / `"non-native"`.
#' @export
simulate_introductions <- function(space, prop = 0.1,
                                   bias = c("random", "peripheral", "central"),
                                   strength = 2, n_axes = 4, seed = NULL) {
  bias <- match.arg(bias)
  if (!is.null(seed)) set.seed(seed)
  X <- .space_coords(space, n_axes)
  d <- sqrt(rowSums(sweep(X, 2, colMeans(X))^2))
  r <- rank(d) / length(d)
  w <- switch(bias, random = rep(1, length(d)), peripheral = r^strength,
              central = (1 - r + 1 / length(d))^strength)
  n <- max(1, round(prop * nrow(X)))
  nn <- sample(rownames(X), n, prob = w)
  stats::setNames(ifelse(rownames(X) %in% nn, "non-native", "native"), rownames(X))
}
