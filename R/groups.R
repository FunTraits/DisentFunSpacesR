utils::globalVariables(c("group", "pair", "label", "a1", "a2", "position", "n_removed", "lower", "upper", "scenario", "x_value"))

## ---------------------------------------------------------------------------
## Species groups in functional spaces
##
## Generalisation of the non-native analysis to any categorical attribute of
## species that is not a trait (origin, human uses, migratory behaviour,
## activity period, threat status, taxonomy, ...). Two levels:
##  * groups   : where does each group sit in each space, how much do groups
##               overlap (convex hull intersection, functional beta diversity),
##               are they more segregated than random labels?
##  * species  : species of focal groups relative to a reference group
##               (novelty, position in the reference hull, loss of uniqueness
##               of reference species, contribution to the assemblage).
## ---------------------------------------------------------------------------

## Internal: parse the group membership of species.
## `groups`: named vector (character / factor / logical; names = species), a
## data frame with columns `species` and `group`, or a species x groups
## logical / 0-1 matrix or data frame (row names = species; non-exclusive
## groups). `sep` splits multi-label values ("food;aquarium"). `recode`
## (named list: new label = old labels) merges labels; `levels` keeps a
## subset of the groups. Returns a logical membership matrix (species x
## groups) and whether groups are mutually exclusive.
.as_groups <- function(groups, species, sep = NULL, recode = NULL, levels = NULL) {
  if (is.null(sep) && !is.null(attr(groups, "sep"))) sep <- attr(groups, "sep")  # groups_from_columns()
  if (is.data.frame(groups) && all(c("species", "group") %in% names(groups)))
    groups <- stats::setNames(as.character(groups$group), groups$species)
  if (is.matrix(groups) || is.data.frame(groups)) {
    M <- as.matrix(groups)
    if (is.null(rownames(M))) stop("A membership matrix needs species as row names.", call. = FALSE)
    if (is.null(colnames(M))) colnames(M) <- paste0("group", seq_len(ncol(M)))
    if (!is.logical(M)) {
      M2 <- suppressWarnings(matrix(as.numeric(M), nrow(M), dimnames = dimnames(M)))
      if (anyNA(M2[!is.na(M)])) stop("A membership matrix must be logical or 0/1.", call. = FALSE)
      M <- M2 > 0
    }
    miss <- setdiff(species, rownames(M))
    M <- rbind(M, matrix(NA, length(miss), ncol(M), dimnames = list(miss, colnames(M))))
    M <- M[species, , drop = FALSE]
    M[is.na(M)] <- FALSE
    labs <- lapply(seq_len(nrow(M)), function(i) colnames(M)[M[i, ]])
  } else {
    if (is.null(names(groups))) stop("'groups' must be named by species.", call. = FALSE)
    g <- groups[match(species, names(groups))]
    if (is.logical(g)) g <- ifelse(g, "TRUE", "FALSE")
    g <- trimws(as.character(g))
    g[g %in% c("", "NA")] <- NA
    labs <- if (!is.null(sep) && nzchar(sep)) {
      lapply(strsplit(ifelse(is.na(g), "", g), sep, fixed = TRUE), function(z) {
        z <- trimws(z); unique(z[nzchar(z)])
      })
    } else lapply(g, function(z) if (is.na(z)) character(0) else z)
  }
  if (!is.null(recode) && length(recode)) {
    map <- stats::setNames(rep(names(recode), lengths(recode)), trimws(unlist(recode)))
    labs <- lapply(labs, function(z) unique(ifelse(z %in% names(map), map[z], z)))
  }
  lv <- sort(unique(unlist(labs)))
  if (!is.null(levels)) {
    bad <- setdiff(levels, lv)
    if (length(bad)) stop("Unknown group(s): ", paste(bad, collapse = ", "), call. = FALSE)
    lv <- levels
  }
  M <- matrix(FALSE, length(species), length(lv), dimnames = list(species, lv))
  for (j in seq_along(lv)) M[, j] <- vapply(labs, function(z) lv[j] %in% z, logical(1))
  if (ncol(M) < 1) stop("No group found.", call. = FALSE)
  list(M = M, exclusive = all(rowSums(M) <= 1), levels = lv)
}

#' Species groups from indicator (0 / 1) columns
#'
#' Many non-trait attributes come as one presence / absence column per
#' category: one column per human use (fisheries, aquaculture, aquarium,
#' game fish...), per habitat, per country... `groups_from_columns()` turns
#' such columns into the group labels used by [groups_all()],
#' [group_overlap()], [group_species()] and [plot_groups()]: one label per
#' column set to 1 (several labels per species, separated by `sep`), and
#' `none` for the species set to 0 in all columns (e.g. species without any
#' use: the natural reference group).
#'
#' Values read as 1: `1`, `TRUE`, `"yes"`, `"y"`, `"x"`; as 0: `0`, `FALSE`,
#' `"no"`, `"n"`, empty; any other value is an error. A species with missing
#' values only gets no group (`NA`); a missing value next to a 1 is ignored.
#'
#' @param data Data frame with the species column and the indicator columns.
#' @param columns Names of the indicator columns (one group per column; the
#'   group takes the name of the column).
#' @param species_col Name of the species column.
#' @param none Label of the species set to 0 in all columns; `NULL` or `""`
#'   to leave them without a group.
#' @param sep Separator of the labels of species in several groups (give the
#'   same `sep` to [groups_all()] and the other functions).
#' @return A character vector named by species, with attribute `sep`.
#' @examples
#' d <- data.frame(species = c("a", "b", "c", "d"), Fisheries = c(1, 0, 1, 0),
#'                 Aquarium = c(0, 1, 1, 0))
#' groups_from_columns(d, c("Fisheries", "Aquarium"), "species")
#' @export
groups_from_columns <- function(data, columns, species_col = "species", none = "none", sep = ";") {
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  if (!length(columns)) stop("Give at least one indicator column.", call. = FALSE)
  miss <- setdiff(c(species_col, columns), names(data))
  if (length(miss)) stop("Column(s) not found: ", paste(miss, collapse = ", "), call. = FALSE)
  if (any(grepl(sep, columns, fixed = TRUE)))
    stop("Column names must not contain the separator '", sep, "'.", call. = FALSE)
  M <- vapply(columns, function(cl) .as_indicator(data[[cl]], cl), numeric(nrow(data)))
  M <- matrix(M, nrow(data), dimnames = list(NULL, columns))
  lab <- apply(M, 1, function(z) {
    if (all(is.na(z))) return(NA_character_)
    z <- columns[!is.na(z) & z > 0]
    if (length(z)) paste(z, collapse = sep) else NA_character_
  })
  allzero <- is.na(lab) & rowSums(!is.na(M)) > 0
  if (!is.null(none) && nzchar(none)) {
    if (none %in% columns) stop("`none` ('", none, "') is also the name of a column.", call. = FALSE)
    lab[allzero] <- none
  }
  out <- stats::setNames(lab, as.character(data[[species_col]]))
  attr(out, "sep") <- sep
  out
}

## Internal: read a presence / absence column as 0 / 1 (NA kept); error if not binary
.as_indicator <- function(v, name = "") {
  if (is.logical(v)) return(as.numeric(v))
  if (is.numeric(v)) {
    if (!all(v[!is.na(v)] %in% c(0, 1))) stop("Column '", name, "' is not a 0 / 1 indicator.", call. = FALSE)
    return(as.numeric(v))
  }
  x <- tolower(trimws(as.character(v)))
  out <- rep(NA_real_, length(x))
  out[x %in% c("1", "true", "yes", "y", "x", "oui")] <- 1
  out[x %in% c("0", "false", "no", "n", "", "non")] <- 0
  if (any(is.na(out) & !is.na(x) & x != "na"))
    stop("Column '", name, "' is not a 0 / 1 indicator.", call. = FALSE)
  out
}

## Internal: is a column a 0 / 1 indicator (at least one 1)?
.is_indicator <- function(v) {
  x <- tryCatch(.as_indicator(v), error = function(e) NULL)
  !is.null(x) && any(x == 1, na.rm = TRUE) && length(unique(x[!is.na(x)])) <= 2
}

## Internal: indicator columns that are (almost) the union of the others
## ("All uses"): redundant with the other groups, left out of the default
## selection. Returns the number of species disagreeing with the union, named
## by column (`tol`: maximum proportion of disagreeing species).
.union_columns <- function(data, columns, tol = 0.02) {
  if (length(columns) < 3) return(stats::setNames(integer(0), character(0)))
  M <- vapply(columns, function(cl) .as_indicator(data[[cl]], cl), numeric(nrow(data)))
  M <- matrix(M, nrow(data), dimnames = list(NULL, columns))
  M[is.na(M)] <- 0
  n_bad <- vapply(seq_along(columns), function(j)
    sum(M[, j] != as.numeric(rowSums(M[, -j, drop = FALSE]) > 0)), numeric(1))
  n_one <- colSums(M)
  # the union holds the most species
  keep <- n_bad <= tol * nrow(M) & n_one > 0 & n_one == max(n_one)
  stats::setNames(as.integer(n_bad[keep]), columns[keep])
}

## Internal: vertices of the convex hull of X (fewer points for intersections)
.hull_vertices <- function(X) {
  X <- unique(round(X, 10))
  if (nrow(X) <= ncol(X) + 1) return(X)
  h <- tryCatch(geometry::convhulln(X), error = function(e) NULL)
  if (is.null(h)) X else X[sort(unique(as.vector(h))), , drop = FALSE]
}

## Internal: volume of the intersection of the hulls of A and B (0 = disjoint)
.hull_inter <- function(A, B) {
  if (nrow(unique(round(A, 10))) <= ncol(A) || nrow(unique(round(B, 10))) <= ncol(B)) return(NA_real_)
  A <- .hull_vertices(A); B <- .hull_vertices(B)
  tryCatch({
    r <- geometry::intersectn(A, B, tol = 1e-10)
    v <- r$ch$vol
    if (is.null(v) || !is.finite(v)) 0 else v
  }, error = function(e) NA_real_)
}

## Internal: pairwise overlap of two point sets (hulls on the first `ha` axes)
.pair_overlap <- function(XA, XB, ha, vol_all) {
  A <- XA[, seq_len(ha), drop = FALSE]; B <- XB[, seq_len(ha), drop = FALSE]
  vA <- .hull_vol(A); vB <- .hull_vol(B); vI <- .hull_inter(A, B)
  if (is.finite(vI)) vI <- min(vI, vA, vB, na.rm = TRUE)
  a <- vI; b <- vA - vI; c <- vB - vI
  jac <- if (all(is.finite(c(a, b, c))) && (a + b + c) > 0) (b + c) / (a + b + c) else NA
  jtu <- if (all(is.finite(c(a, b, c))) && (a + 2 * min(b, c)) > 0) 2 * min(b, c) / (a + 2 * min(b, c)) else NA
  c(V1 = vA / vol_all, V2 = vB / vol_all, V_inter = vI / vol_all,
    overlap = if (is.finite(jac)) 1 - jac else NA,
    prop_V1_shared = if (is.finite(vA) && vA > 0) vI / vA else NA,
    prop_V2_shared = if (is.finite(vB) && vB > 0) vI / vB else NA,
    beta_jac = jac, beta_jtu = jtu, beta_jne = jac - jtu,
    centroid_dist = sqrt(sum((colMeans(XA) - colMeans(XB))^2)))
}

## Internal: PERMANOVA and PERMDISP on Euclidean coordinates (exclusive groups).
## With Euclidean distances, the pseudo-F of PERMANOVA (Anderson 2001) is the
## classical between / within sum of squares ratio computed from coordinates,
## so permutations cost O(n p) instead of O(n^2).
.group_tests <- function(X, g, nperm = 999) {
  g <- factor(g); n <- nrow(X); G <- nlevels(g)
  if (G < 2 || n <= G) return(c(permanova_F = NA, permanova_R2 = NA, permanova_p = NA,
                                permdisp_F = NA, permdisp_p = NA))
  sx2 <- sum(X^2); SST <- sx2 - sum(colSums(X)^2) / n
  ssw <- function(gg) sx2 - sum(rowSums(rowsum(X, gg)^2) / as.vector(table(gg)))
  fstat <- function(w) ((SST - w) / (G - 1)) / (w / (n - G))
  W <- ssw(g); F0 <- fstat(W)
  disp <- function(gg) {
    cen <- rowsum(X, gg) / as.vector(table(gg))
    z <- sqrt(rowSums((X - cen[as.character(gg), , drop = FALSE])^2))
    m <- tapply(z, gg, mean); nn <- table(gg)
    sb <- sum(nn * (m - mean(z))^2); sw <- sum((z - m[as.character(gg)])^2)
    if (sw > 0) (sb / (G - 1)) / (sw / (n - G)) else NA
  }
  D0 <- disp(g)
  pf <- pd <- NA
  if (nperm > 0) {
    perms <- replicate(nperm, { gg <- sample(g); c(fstat(ssw(gg)), disp(gg)) })
    pf <- (sum(perms[1, ] >= F0 - 1e-12) + 1) / (nperm + 1)
    pd <- if (is.finite(D0)) (sum(perms[2, ] >= D0 - 1e-12, na.rm = TRUE) + 1) / (nperm + 1) else NA
  }
  c(permanova_F = F0, permanova_R2 = (SST - W) / SST, permanova_p = pf,
    permdisp_F = D0, permdisp_p = pd)
}

#' Functional structure and overlap of species groups in a functional space
#'
#' Describes how groups of species defined by any categorical attribute that
#' is not a trait (native / non-native status, human uses, migratory
#' behaviour, activity period, threat status, taxonomy, ...) occupy one
#' functional space, how much they overlap and whether they are more
#' segregated than expected from random labels. Groups may be mutually
#' exclusive (one label per species) or not (several labels per species,
#' e.g. `"food;aquarium"` with `sep = ";"`).
#'
#' Per group (`groups`): `n` species; `FRic` volume of the convex hull of
#' the group (first `hull_axes` axes) relative to the hull of all species
#' with a group; `FDis` mean distance to the group centroid; `FOri` mean
#' nearest-neighbour distance within the group; `MPD` mean pairwise
#' distance; `p_outside_others` proportion of the species of the group
#' outside the hull of the species not in the group; `prop_V_exclusive`
#' share of the hull of the group not shared with the hull of the other
#' species; `knn_same` mean proportion of the `k` nearest neighbours of its
#' species that belong to the group, `knn_expected` its expectation under
#' random labels and `segregation` their ratio (> 1: species of the group
#' are functionally clustered).
#'
#' Per pair of groups (`pairs`), from the volumes of the hulls of the two
#' groups and of their intersection (Villéger et al. 2013,
#' <doi:10.1111/geb.12021>): `V1`, `V2`, `V_inter` (relative to the hull of
#' all species), `overlap` (intersection / union), `prop_V1_shared`,
#' `prop_V2_shared` (intersection / volume of each group), the functional
#' Jaccard dissimilarity `beta_jac` and its turnover (`beta_jtu`) and
#' nestedness-resultant (`beta_jne`) components (Baselga 2012), `p1_in_hull2`
#' and `p2_in_hull1` (proportion of the species of one group inside the hull
#' of the other) and `centroid_dist` (distance between centroids). With
#' `null = "permute"`, labels are permuted among the species of the two
#' groups (sizes and shared species kept): `ses_` = (obs - mean null) / sd
#' null and `p_` two-sided rank-based p-values for `overlap` and
#' `centroid_dist`. Hull sizes grow with the number of species: compare
#' groups of different sizes with the SES rather than the raw overlap.
#'
#' Global test (`test`, exclusive groups only): PERMANOVA (pseudo-F, R2;
#' Anderson 2001) and PERMDISP (homogeneity of dispersions; Anderson 2006),
#' computed on the Euclidean coordinates of the space (`n_axes` axes), with
#' `nperm` label permutations.
#'
#' Densities (`density = TRUE`): hulls depend on the most extreme species
#' only. The probability density of each group is also estimated on a grid
#' of the first `density_axes` axes (2-4; Gaussian kernel, bandwidth from
#' Scott's rule x `bw`, the same grid and bandwidth for all groups of the
#' space; optionally weighted by `abundance`). Per group, `FRic_kde` is the
#' share of the grid holding 99% of the density of the group relative to
#' that of all species; per pair, `overlap_kde` is the sum over the grid of
#' the minimum of the two densities (the overlap of trait probability
#' densities, Carmona et al. 2016, <doi:10.1016/j.tree.2016.02.003>; 1 =
#' identical distributions, 0 = disjoint), with SES under the null model.
#'
#' Per species and group (`species`): `knn_same` (proportion of the `k`
#' nearest neighbours in the same group), `dist_nearest_same`,
#' `dist_nearest_other` (nearest species not in the group) and
#' `dist_centroid` (distance to the group centroid).
#'
#' @param space An `fs_space`.
#' @param groups Group of every species: named vector (names = species;
#'   character, factor or logical), data frame with columns `species` and
#'   `group`, or a species x groups logical (or 0/1) matrix / data frame with
#'   species as row names. Species without a group (`NA`, `""`) are left
#'   out.
#' @param sep Separator of multiple labels in a character `groups` (e.g.
#'   `";"`); `NULL` = one label per species.
#' @param recode Optional named list merging labels, e.g.
#'   `list(threatened = c("VU", "EN", "CR"), other = c("LC", "NT"))`; labels
#'   not listed are kept.
#' @param levels Optional groups to keep (after `recode`).
#' @inheritParams species_metrics
#' @param hull_axes Number of axes used for convex hulls.
#' @param null `"none"` or `"permute"` (pairwise null model).
#' @param nrep Replicates of the pairwise null model.
#' @param nperm Permutations of the global tests (0 = no p-value).
#' @param seed Optional random seed.
#' @param density Also estimate the overlap of kernel densities?
#' @param density_axes Axes of the density grid (2 to 4).
#' @param bw Multiplier of the bandwidth.
#' @param n_grid Cells per axis (default 200, 100, 40, 20 for 1-4 axes).
#' @param abundance Optional named numeric vector (names = species) weighting
#'   species in the densities (e.g. abundance, biomass, occupancy).
#' @return An object of class `fs_groups`: list with data frames `groups`,
#'   `pairs`, `test`, `species` and `settings`.
#' @seealso [group_species()], [groups_all()], [plot_groups()]
#' @export
#' @examples
#' \donttest{
#' data(birds_traits); data(birds_trait_info)
#' tr <- drop_incomplete(prepare_traits(birds_traits, birds_trait_info, "species"))
#' sp <- build_space(tr, dimensions = c("Diet", "Locomotion"))
#' ord <- stats::setNames(birds_traits$order, birds_traits$species)
#' g <- group_overlap(sp, ord, levels = c("Passeriformes", "Charadriiformes", "Accipitriformes"),
#'                    hull_axes = 3, nperm = 99)
#' g$pairs
#' }
group_overlap <- function(space, groups, sep = NULL, recode = NULL, levels = NULL,
                          n_axes = NULL, k = 5, hull_axes = 4,
                          null = c("none", "permute"), nrep = 99, nperm = 999,
                          seed = NULL, tol = 1e-8, density = FALSE, density_axes = 2,
                          bw = 1, n_grid = NULL, abundance = NULL) {
  null <- match.arg(null)
  if (!is.null(seed)) set.seed(seed)
  X0 <- .space_coords(space, n_axes)
  gr <- .as_groups(groups, rownames(X0), sep = sep, recode = recode, levels = levels)
  keep <- rowSums(gr$M) > 0
  if (sum(keep) < 3) stop("Fewer than three species have a group.", call. = FALSE)
  X <- X0[keep, , drop = FALSE]; M <- gr$M[keep, , drop = FALSE]
  lv <- colnames(M); n <- nrow(X)
  ha <- min(hull_axes, ncol(X))
  vol_all <- .hull_vol(X[, seq_len(ha), drop = FALSE])
  kk <- min(k, n - 1)
  knn <- FNN::get.knn(X, k = kk, algorithm = if (ncol(X) > 10) "brute" else "kd_tree")
  if (density) {
    ks <- .kde_setup(X, density_axes, n_grid, bw)
    wts <- if (is.null(abundance)) rep(1, n) else {
      w <- unname(abundance[rownames(X)]); w[is.na(w)] <- 0; w
    }
    dens <- function(i) .kde_density(X[i, , drop = FALSE], ks, wts[i])
    p_all <- dens(seq_len(n)); c_all <- .kde_cells(p_all)
    p_grp <- lapply(stats::setNames(lv, lv), function(g) dens(which(M[, g])))
  }

  ## --- species x group
  sp_rows <- lapply(lv, function(g) {
    inn <- M[, g]; idx <- which(inn)
    same <- rowMeans(matrix(inn[knn$nn.index[idx, , drop = FALSE]], ncol = kk))
    Xg <- X[idx, , drop = FALSE]; Xo <- X[!inn, , drop = FALSE]
    dsame <- if (length(idx) > 1) {
      d <- FNN::get.knn(Xg, k = 1)$nn.dist[, 1]; d[d < tol] <- 0; d
    } else NA_real_
    doth <- if (nrow(Xo)) FNN::get.knnx(Xo, Xg, k = 1)$nn.dist[, 1] else NA_real_
    data.frame(group = g, species = rownames(Xg), knn_same = same,
               dist_nearest_same = dsame, dist_nearest_other = doth,
               dist_centroid = sqrt(rowSums(sweep(Xg, 2, colMeans(Xg))^2)),
               row.names = NULL, stringsAsFactors = FALSE)
  })
  species <- do.call(rbind, sp_rows)

  ## --- per group
  grp <- do.call(rbind, lapply(lv, function(g) {
    inn <- M[, g]; ng <- sum(inn)
    Xg <- X[inn, , drop = FALSE]; Xo <- X[!inn, , drop = FALSE]
    am <- .assemblage_metrics(Xg, ha, vol_all)
    vg <- .hull_vol(Xg[, seq_len(ha), drop = FALSE])
    vi <- if (nrow(Xo) > ha) .hull_inter(Xg[, seq_len(ha), drop = FALSE],
                                         Xo[, seq_len(ha), drop = FALSE]) else NA
    pout <- if (ng && nrow(Xo) > ha)
      mean(!.in_hull(Xo[, seq_len(ha), drop = FALSE], Xg[, seq_len(ha), drop = FALSE]), na.rm = TRUE) else NA
    ks <- mean(species$knn_same[species$group == g])
    ke <- (ng - 1) / (n - 1)
    data.frame(group = g, n = ng, FRic = am[["FRic"]], FDis = am[["FDis"]],
               FOri = am[["FOri"]], MPD = am[["MPD"]], p_outside_others = pout,
               prop_V_exclusive = if (is.finite(vg) && vg > 0 && is.finite(vi)) 1 - min(vi, vg) / vg else NA,
               knn_same = ks, knn_expected = ke,
               segregation = if (ke > 0) ks / ke else NA,
               FRic_kde = if (density) .kde_cells(p_grp[[g]]) / c_all else NA,
               stringsAsFactors = FALSE)
  }))
  if (!density) grp$FRic_kde <- NULL

  ## --- pairs of groups
  pr <- if (length(lv) > 1) utils::combn(lv, 2, simplify = FALSE) else list()
  pairs <- do.call(rbind, lapply(pr, function(p) {
    a <- M[, p[1]]; b <- M[, p[2]]
    o <- .pair_overlap(X[a, , drop = FALSE], X[b, , drop = FALSE], ha, vol_all)
    p12 <- if (sum(b) > ha) mean(.in_hull(X[b, seq_len(ha), drop = FALSE], X[a, seq_len(ha), drop = FALSE]), na.rm = TRUE) else NA
    p21 <- if (sum(a) > ha) mean(.in_hull(X[a, seq_len(ha), drop = FALSE], X[b, seq_len(ha), drop = FALSE]), na.rm = TRUE) else NA
    row <- data.frame(group1 = p[1], group2 = p[2], n1 = sum(a), n2 = sum(b),
                      n_shared = sum(a & b), t(o), p1_in_hull2 = p12, p2_in_hull1 = p21,
                      stringsAsFactors = FALSE, check.names = FALSE)
    mets <- c("overlap", "centroid_dist")
    if (density) {
      row$overlap_kde <- .kde_overlap(p_grp[[p[1]]], p_grp[[p[2]]])
      o <- c(o, overlap_kde = row$overlap_kde); mets <- c(mets, "overlap_kde")
    }
    if (null == "permute" && nrep > 1) {
      U <- which(a | b); s <- sum(a & b); na <- sum(a); nb <- sum(b)
      sims <- replicate(nrep, {
        u <- sample(U)
        ia <- u[seq_len(na)]
        ib <- c(u[seq_len(s)], u[na + seq_len(nb - s)])
        oo <- .pair_overlap(X[ia, , drop = FALSE], X[ib, , drop = FALSE], ha, vol_all)
        if (density) oo <- c(oo, overlap_kde = .kde_overlap(dens(ia), dens(ib)))
        oo[mets]
      })
      sims <- matrix(sims, nrow = length(mets), dimnames = list(mets, NULL))
      for (m in mets) {
        obs <- o[[m]]; v <- sims[m, ]; v <- v[is.finite(v)]
        row[[paste0("ses_", m)]] <- if (length(v) > 1 && stats::sd(v) > 0 && is.finite(obs))
          (obs - mean(v)) / stats::sd(v) else NA
        row[[paste0("p_", m)]] <- if (length(v) > 1 && is.finite(obs))
          min(1, 2 * min((sum(v <= obs) + 1) / (length(v) + 1),
                         (sum(v >= obs) + 1) / (length(v) + 1))) else NA
      }
    }
    row
  }))
  if (is.null(pairs)) pairs <- data.frame()

  ## --- global tests
  test <- if (gr$exclusive && length(lv) > 1) {
    gl <- lv[max.col(M, ties.method = "first")]
    as.data.frame(t(.group_tests(X, gl, nperm)))
  } else data.frame(permanova_F = NA, permanova_R2 = NA, permanova_p = NA,
                    permdisp_F = NA, permdisp_p = NA)
  test <- cbind(n_groups = length(lv), n_species = n, exclusive = gr$exclusive, test)

  structure(list(groups = grp, pairs = pairs, test = test, species = species,
                 settings = list(n_axes = n_axes, k = kk, hull_axes = ha, null = null,
                                 nrep = nrep, nperm = nperm, exclusive = gr$exclusive,
                                 levels = lv, n_unassigned = sum(!keep),
                                 density = if (density) list(axes = ks$d, n_grid = ks$n_grid,
                                                             bandwidth = ks$h,
                                                             weighted = !is.null(abundance)))),
            class = "fs_groups")
}

#' @export
print.fs_groups <- function(x, ...) {
  s <- x$settings
  cat(sprintf("Species groups in a functional space: %d groups (%s), %d species; hulls on %d axes\n",
              length(s$levels), if (s$exclusive) "exclusive" else "non-exclusive",
              x$test$n_species, s$hull_axes))
  if (is.finite(x$test$permanova_R2))
    cat(sprintf("  PERMANOVA R2 = %.3f (p = %s); PERMDISP p = %s\n", x$test$permanova_R2,
                format(x$test$permanova_p, digits = 3), format(x$test$permdisp_p, digits = 3)))
  if (nrow(x$pairs))
    cat(sprintf("  pairwise hull overlap: median %.2f (range %.2f-%.2f)\n",
                stats::median(x$pairs$overlap, na.rm = TRUE),
                suppressWarnings(min(x$pairs$overlap, na.rm = TRUE)),
                suppressWarnings(max(x$pairs$overlap, na.rm = TRUE))))
  invisible(x)
}

## Internal: logical focal status (TRUE = focal, FALSE = reference, NA = other)
.focal_status <- function(M, reference, focal) {
  ref <- rowSums(M[, reference, drop = FALSE]) > 0
  foc <- rowSums(M[, focal, drop = FALSE]) > 0
  st <- rep(NA, nrow(M)); st[foc] <- TRUE; st[ref] <- FALSE   # shared -> reference
  stats::setNames(st, rownames(M))
}

## Internal: resolve reference / focal groups
.ref_focal <- function(gr, reference, focal, pool_focal) {
  lv <- gr$levels
  if (is.null(reference) || !length(reference)) stop("Give the reference group(s).", call. = FALSE)
  bad <- setdiff(reference, lv)
  if (length(bad)) stop("Unknown reference group(s): ", paste(bad, collapse = ", "), call. = FALSE)
  if (is.null(focal)) focal <- setdiff(lv, reference)
  focal <- setdiff(focal, reference)
  if (!length(focal)) stop("No focal group (all groups are reference).", call. = FALSE)
  foc <- if (pool_focal) stats::setNames(list(focal), paste(focal, collapse = " + ")) else
    stats::setNames(as.list(focal), focal)
  list(reference = reference, focal = foc,
       ref_label = paste(reference, collapse = " + "))
}

#' Species of focal groups relative to a reference group
#'
#' Generalises [nonnative_species()] to any grouping of species: the
#' reference group plays the role of the native fauna and each focal group
#' that of the non-native species (e.g. reference = resident birds, focal =
#' migratory birds; reference = species without human use, focal = species
#' used for food, for the aquarium trade, ...; reference = non-threatened,
#' focal = threatened species). Each focal group is analysed separately
#' against the reference (or all together with `pool_focal = TRUE`).
#'
#' Focal species (`focal`): `FNov` distance to the nearest reference species
#' (novelty; `FNov_raw` before min-max rescaling within the focal group),
#' `FSp_ref` distance to the reference centroid, `mean_dist_ref`,
#' `inside_ref_hull` (inside the convex hull of the reference species on the
#' first `hull_axes` axes), `nearest_ref`; with `weight`, `W` and
#' `FUSW = log(1 + W x FNov + W x FSp_ref)` (built as FUSE / FUSI).
#'
#' Reference species (`reference`, one block per focal group): FUn among
#' reference species only (`FUn_before_raw`) and once the focal species are
#' added (`FUn_after_raw`), `dFUn_rel` relative change (redundancy of focal
#' species with reference species), `n_focal_knn`, `nearest_focal`,
#' `dist_nearest_focal`, `dist_nearest_ref`.
#'
#' Species belonging to both the reference and a focal group (non-exclusive
#' groups) are counted in the reference only.
#'
#' @inheritParams group_overlap
#' @param reference Label(s) of the reference group (several labels are
#'   merged).
#' @param focal Focal groups (default: all groups other than the
#'   reference).
#' @param pool_focal Analyse all focal groups together?
#' @param weight Optional weights of focal species for FUSW: named numeric
#'   vector in \[0, 1\] (names = species), or a named character vector of
#'   categories matched to `weight_values`.
#' @param weight_values Named numeric vector of weights per category (e.g.
#'   [eicat_weights()]).
#' @return An object of class `fs_group_species`: list with data frames
#'   `focal` and `reference` (column `group` = focal group) and `settings`.
#' @seealso [nonnative_species()], [group_overlap()], [groups_all()]
#' @export
group_species <- function(space, groups, reference, focal = NULL, pool_focal = FALSE,
                          sep = NULL, recode = NULL, levels = NULL, n_axes = NULL,
                          k = 5, hull_axes = 4, tol = 1e-8, weight = NULL,
                          weight_values = NULL) {
  X <- .space_coords(space, n_axes)
  gr <- .as_groups(groups, rownames(X), sep = sep, recode = recode, levels = levels)
  rf <- .ref_focal(gr, reference, focal, pool_focal)
  out <- lapply(names(rf$focal), function(fg) {
    st <- .focal_status(gr$M, rf$reference, rf$focal[[fg]])
    if (!any(st %in% TRUE)) return(NULL)
    imp <- NULL
    if (!is.null(weight)) {
      imp <- if (is.numeric(weight)) weight else {
        w <- weight_values[toupper(trimws(as.character(weight)))]
        stats::setNames(unname(w), names(weight))
      }
    }
    r <- nonnative_species(space, st, n_axes = n_axes, k = k, hull_axes = hull_axes,
                           tol = tol, impact = imp, impact_weights = c(x = 1))
    f <- r$nonnative
    names(f) <- c(species = "species", FNov = "FNov", FSp_native = "FSp_ref",
                  FNov_raw = "FNov_raw", FSp_native_raw = "FSp_ref_raw",
                  mean_dist_native = "mean_dist_ref", inside_native_hull = "inside_ref_hull",
                  nearest_native = "nearest_ref", impact = "impact", W = "W",
                  FUSI = "FUSW")[names(f)]
    f$impact <- NULL
    rr <- r$native
    names(rr) <- c(species = "species", FUn_before_raw = "FUn_before_raw",
                   FUn_after_raw = "FUn_after_raw", dFUn_rel = "dFUn_rel",
                   n_nonnative_knn = "n_focal_knn", nearest_nonnative = "nearest_focal",
                   dist_nearest_nonnative = "dist_nearest_focal",
                   dist_nearest_native = "dist_nearest_ref")[names(rr)]
    list(focal = cbind(group = fg, f, stringsAsFactors = FALSE),
         reference = cbind(group = fg, rr, stringsAsFactors = FALSE),
         n = c(n_ref = r$settings$n_native, n_focal = r$settings$n_nonnative))
  })
  out <- out[!vapply(out, is.null, logical(1))]
  if (!length(out)) stop("No focal species.", call. = FALSE)
  structure(list(focal = do.call(rbind, lapply(out, `[[`, "focal")),
                 reference = do.call(rbind, lapply(out, `[[`, "reference")),
                 settings = list(reference = rf$ref_label, focal = names(rf$focal),
                                 n = do.call(rbind, lapply(out, `[[`, "n")),
                                 k = k, hull_axes = min(hull_axes, ncol(X)))),
            class = "fs_group_species")
}

#' @export
print.fs_group_species <- function(x, ...) {
  cat(sprintf("Focal groups relative to the reference group '%s'\n", x$settings$reference))
  for (g in unique(x$focal$group)) {
    a <- x$focal[x$focal$group == g, ]; b <- x$reference[x$reference$group == g, ]
    cat(sprintf("  %s: %d species, %.0f%% outside the reference hull, median loss of FUn of reference species %.1f%%\n",
                g, nrow(a), 100 * mean(!a$inside_ref_hull, na.rm = TRUE),
                100 * stats::median(b$dFUn_rel)))
  }
  invisible(x)
}

#' Contribution of focal groups to the functional structure of the fauna
#'
#' Generalises [nonnative_assemblage()]: compares the reference group alone
#' (`_ref`) with the reference + focal species (`_all`) for `S`, `FRic`,
#' `FDis`, `FOri` and `MPD`; `d_` = all - reference, i.e. the contribution of
#' the focal species (for non-native species: the effect of introductions;
#' for threatened or exploited species: what the fauna would lose without
#' them, read with the opposite sign). `FRic_change_rel` = relative
#' contribution to FRic. With `null = "permute"`, the same number of species
#' drawn at random within the whole set plays the role of the focal species.
#'
#' @inheritParams group_species
#' @param null `"none"` or `"permute"`.
#' @param nrep Null replicates.
#' @param seed Optional random seed.
#' @return Data frame, one row per focal group.
#' @export
group_contribution <- function(space, groups, reference, focal = NULL, pool_focal = FALSE,
                               sep = NULL, recode = NULL, levels = NULL, n_axes = NULL,
                               hull_axes = 4, null = c("none", "permute"), nrep = 99,
                               seed = NULL) {
  null <- match.arg(null)
  X <- .space_coords(space, n_axes)
  gr <- .as_groups(groups, rownames(X), sep = sep, recode = recode, levels = levels)
  rf <- .ref_focal(gr, reference, focal, pool_focal)
  rows <- lapply(names(rf$focal), function(fg) {
    st <- .focal_status(gr$M, rf$reference, rf$focal[[fg]])
    if (!any(st %in% TRUE)) return(NULL)
    a <- nonnative_assemblage(space, st, n_axes = n_axes, hull_axes = hull_axes,
                              null = null, nrep = nrep, seed = seed)
    a$assemblage <- NULL
    nm <- names(a)
    nm <- sub("_hist$", "_ref", sub("_cur$", "_all", nm))
    nm[nm == "n_native"] <- "n_ref"; nm[nm == "n_nonnative"] <- "n_focal"
    nm[nm == "p_nn_outside"] <- "p_focal_outside"
    names(a) <- nm
    a$n_extirpated <- NULL
    cbind(group = fg, a, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
}

#' Species groups in every functional space
#'
#' Runs [group_overlap()] in each space and, when a `reference` group is
#' given, [group_species()] and [group_contribution()], so that the
#' functional dimensions that best separate the groups can be identified.
#' This generalises [nonnative_all()] to any categorical attribute of
#' species (not a trait).
#'
#' @param spaces An `fs_spaces` object (or a named list of `fs_space`).
#' @inheritParams group_overlap
#' @inheritParams group_species
#' @param reference Optional reference group(s) for the species-level
#'   analysis (`NULL` = groups only).
#' @param contribution Also compute [group_contribution()] (needs
#'   `reference`)?
#' @param comm Optional assemblage x species matrix (or output of
#'   [community_from_long()]): adds [group_assemblage()] in each space (needs
#'   `reference`).
#' @param site_groups Optional assemblage x species matrix of labels when
#'   groups differ between assemblages (default: `groups` everywhere). When
#'   `groups` is `NULL`, the species-level analyses use the most frequent
#'   label of each species across assemblages.
#' @return An object of class `fs_groups_all`: list of long data frames
#'   (column `space` added) `groups`, `pairs`, `test`, `species` and, with a
#'   reference, `focal`, `reference`, `contribution` and, with `comm`,
#'   `sites`; plus `settings`.
#' @export
#' @examples
#' \donttest{
#' data(birds_traits); data(birds_trait_info)
#' tr <- drop_incomplete(prepare_traits(birds_traits, birds_trait_info, "species"))
#' sp <- build_spaces(tr, spaces = list(Diet = "Diet", Locomotion = "Locomotion"))
#' iucn <- stats::setNames(birds_traits$iucn, birds_traits$species)
#' x <- groups_all(sp, iucn, recode = list(threatened = c("VU", "EN", "CR"),
#'                                         `not threatened` = c("LC", "NT")),
#'                 levels = c("threatened", "not threatened"),
#'                 reference = "not threatened", hull_axes = 3, nperm = 99)
#' summary_groups(x)
#' summary_groups(x, "focal")
#' }
groups_all <- function(spaces, groups, reference = NULL, focal = NULL, pool_focal = FALSE,
                       sep = NULL, recode = NULL, levels = NULL, n_axes = NULL, k = 5,
                       hull_axes = 4, null = c("none", "permute"), nrep = 99,
                       nperm = 999, seed = NULL, contribution = TRUE,
                       weight = NULL, weight_values = NULL, density = FALSE,
                       density_axes = 2, bw = 1, abundance = NULL, comm = NULL,
                       site_groups = NULL) {
  null <- match.arg(null)
  if (is.list(comm) && !is.data.frame(comm) && !is.null(comm$comm)) {
    if (is.null(site_groups)) site_groups <- comm$groups
    comm <- comm$comm
  }
  if (is.null(groups) && !is.null(site_groups)) {
    # one label per species: its most frequent label across assemblages
    sg <- as.matrix(site_groups)
    groups <- apply(sg, 2, function(v) {
      v <- v[!is.na(v) & nzchar(v)]
      if (!length(v)) NA_character_ else names(sort(table(v), decreasing = TRUE))[1]
    })
  }
  if (is.null(groups)) stop("Give `groups` (or `site_groups` with `comm`).", call. = FALSE)
  add <- function(d, nm) if (!is.null(d) && nrow(d)) cbind(space = nm, d, stringsAsFactors = FALSE)
  res <- lapply(names(spaces), function(nm) {
    s <- spaces[[nm]]
    o <- group_overlap(s, groups, sep = sep, recode = recode, levels = levels, n_axes = n_axes,
                       k = k, hull_axes = hull_axes, null = null, nrep = nrep, nperm = nperm,
                       seed = seed, density = density, density_axes = density_axes, bw = bw,
                       abundance = abundance)
    r <- if (!is.null(reference) && length(reference))
      group_species(s, groups, reference, focal = focal, pool_focal = pool_focal, sep = sep,
                    recode = recode, levels = levels, n_axes = n_axes, k = k,
                    hull_axes = hull_axes, weight = weight, weight_values = weight_values)
    cc <- if (!is.null(r) && contribution)
      group_contribution(s, groups, reference, focal = focal, pool_focal = pool_focal,
                         sep = sep, recode = recode, levels = levels, n_axes = n_axes,
                         hull_axes = hull_axes, null = null, nrep = nrep, seed = seed)
    st <- if (!is.null(comm) && !is.null(r))
      group_assemblage(s, comm, groups = if (is.null(site_groups)) groups else site_groups,
                       reference = reference, focal = focal, pool_focal = pool_focal, sep = sep,
                       recode = recode, levels = levels, n_axes = n_axes, hull_axes = hull_axes,
                       null = null, nrep = nrep, seed = seed)
    list(groups = add(o$groups, nm), pairs = add(o$pairs, nm), test = add(o$test, nm),
         sites = add(st, nm),
         species = add(o$species, nm), focal = if (!is.null(r)) add(r$focal, nm),
         reference = if (!is.null(r)) add(r$reference, nm), contribution = add(cc, nm),
         settings = o$settings)
  })
  bind <- function(what) {
    d <- do.call(rbind, lapply(res, `[[`, what))
    if (!is.null(d)) rownames(d) <- NULL
    d
  }
  st <- res[[1]]$settings
  st$reference <- if (!is.null(reference)) paste(reference, collapse = " + ")
  structure(list(groups = bind("groups"), pairs = bind("pairs"), test = bind("test"),
                 species = bind("species"), focal = bind("focal"),
                 reference = bind("reference"), contribution = bind("contribution"),
                 sites = bind("sites"), settings = st),
            class = "fs_groups_all")
}

#' @export
print.fs_groups_all <- function(x, ...) {
  cat(sprintf("Species groups across functional spaces (%d groups%s)\n",
              length(x$settings$levels),
              if (!is.null(x$settings$reference)) sprintf("; reference: %s", x$settings$reference) else ""))
  print(summary_groups(x), row.names = FALSE, digits = 3)
  invisible(x)
}

#' Which functional spaces best separate the species groups?
#'
#' `what = "groups"`: one row per space with the PERMANOVA R2 and p-value,
#' PERMDISP p-value, median and minimum pairwise hull overlap, median SES of
#' the overlap (null model), and median segregation of the groups.
#' `what = "focal"`: one row per space and focal group with the proportion
#' of focal species outside the reference hull, median novelty, novelty
#' relative to the mean nearest-neighbour distance among reference species
#' (`FNov_rel`, comparable across spaces), median relative loss of FUn of
#' reference species, and the relative contribution of the focal group to
#' FRic.
#'
#' @param x Output of [groups_all()].
#' @param what `"groups"` or `"focal"`.
#' @return Data frame.
#' @export
summary_groups <- function(x, what = c("groups", "focal")) {
  what <- match.arg(what)
  med <- function(v) if (length(v) && any(is.finite(v))) stats::median(v, na.rm = TRUE) else NA_real_
  if (what == "groups") {
    return(do.call(rbind, lapply(unique(x$groups$space), function(s) {
      t <- x$test[x$test$space == s, ]
      p <- if (!is.null(x$pairs)) x$pairs[x$pairs$space == s, ] else data.frame()
      g <- x$groups[x$groups$space == s, ]
      data.frame(space = s, n_groups = t$n_groups,
                 permanova_R2 = t$permanova_R2, permanova_p = t$permanova_p,
                 permdisp_p = t$permdisp_p,
                 overlap_median = if (nrow(p)) med(p$overlap) else NA,
                 overlap_min = if (nrow(p) && any(is.finite(p$overlap))) min(p$overlap, na.rm = TRUE) else NA,
                 ses_overlap_median = if (nrow(p) && "ses_overlap" %in% names(p)) med(p$ses_overlap) else NA,
                 overlap_kde_median = if (nrow(p) && "overlap_kde" %in% names(p)) med(p$overlap_kde) else NA,
                 segregation_median = med(g$segregation), stringsAsFactors = FALSE)
    })))
  }
  if (is.null(x$focal)) stop("No reference group: run groups_all() with `reference`.", call. = FALSE)
  keys <- unique(x$focal[, c("space", "group")])
  do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    s <- keys$space[i]; g <- keys$group[i]
    a <- x$focal[x$focal$space == s & x$focal$group == g, ]
    b <- x$reference[x$reference$space == s & x$reference$group == g, ]
    ref <- mean(b$dist_nearest_ref, na.rm = TRUE)
    row <- data.frame(space = s, group = g, n_focal = nrow(a),
                      p_outside_ref_hull = mean(!a$inside_ref_hull, na.rm = TRUE),
                      FNov_median = stats::median(a$FNov_raw),
                      FNov_rel = if (ref > 0) mean(a$FNov_raw) / ref else NA,
                      dFUn_ref_median = stats::median(b$dFUn_rel),
                      p_ref_lose_10pc = mean(b$dFUn_rel <= -0.1), stringsAsFactors = FALSE)
    cc <- x$contribution
    if (!is.null(cc)) {
      z <- cc[cc$space == s & cc$group == g, ]
      if (nrow(z)) {
        row$FRic_contribution_rel <- z$FRic_change_rel[1]
        if ("ses_FRic" %in% names(z)) row$ses_FRic <- z$ses_FRic[1]
      }
    }
    row
  }))
}

## Internal: colour palette for groups (Okabe-Ito, then HCL)
.group_palette <- function(lv) {
  oi <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#000000")
  cols <- if (length(lv) <= length(oi)) oi[seq_along(lv)] else grDevices::hcl.colors(length(lv), "Dark 3")
  stats::setNames(cols, lv)
}

#' Plot species groups in a functional space
#'
#' Without `reference`: species coloured by group with the convex hull of
#' each group on the plotted axes. With `reference`: reference species in
#' grey with their hull, focal species coloured by group, filled when outside
#' the reference hull on the first `hull_axes` axes (from `x`, output of
#' [group_species()] for the same space, or computed).
#'
#' @inheritParams group_species
#' @param axes Two axes to plot.
#' @param reference Optional reference group(s).
#' @param x Optional output of [group_species()] for this space.
#' @param hulls Draw the hull of each group (without `reference`)?
#' @param contours Draw density contours of each group (2-D kernel density
#'   on the plotted axes) instead of hulls?
#' @param label Number of focal species with highest `FNov` to label (with
#'   a reference).
#' @param centroids Draw the centroid of each group (mean position of its
#'   species, see [group_centroids()]) as a large diamond outlined in black?
#'   With `reference`, the centroid of the reference species is drawn too.
#' @return A ggplot.
#' @export
plot_groups <- function(space, groups, axes = c(1, 2), sep = NULL, recode = NULL,
                        levels = NULL, reference = NULL, focal = NULL,
                        pool_focal = FALSE, x = NULL, hull_axes = 4, hulls = TRUE,
                        label = 0, contours = FALSE, centroids = TRUE) {
  X <- space$coords
  gr <- .as_groups(groups, rownames(X), sep = sep, recode = recode, levels = levels)
  long <- do.call(rbind, lapply(colnames(gr$M), function(g) {
    i <- gr$M[, g]
    if (!any(i)) return(NULL)
    data.frame(species = rownames(X)[i], group = g, a1 = X[i, axes[1]], a2 = X[i, axes[2]],
               stringsAsFactors = FALSE)
  }))
  ax <- paste0("PCoA ", axes)
  hull_of <- function(d) do.call(rbind, lapply(split(d, d$group), function(z)
    if (nrow(z) >= 3) z[grDevices::chull(z$a1, z$a2), ] else NULL))
  theme <- list(ggplot2::labs(x = ax[1], y = ax[2]), ggplot2::coord_equal(),
                ggplot2::theme_classic(base_size = 11))
  # centroid (mean position) of each group on the plotted axes
  cent_of <- function(d) {
    cc <- do.call(rbind, lapply(split(d, d$group), function(z)
      data.frame(group = z$group[1], a1 = mean(z$a1), a2 = mean(z$a2), n = nrow(z),
                 stringsAsFactors = FALSE)))
    cc$hover <- sprintf("Centroid of %s<br>%d species<br>(%.3g, %.3g)", cc$group, cc$n, cc$a1, cc$a2)
    cc
  }
  if (is.null(reference)) {
    lv <- unique(long$group); pal <- .group_palette(colnames(gr$M))[lv]
    long$group <- factor(long$group, levels = lv)
    p <- ggplot2::ggplot()
    if (contours) {
      p <- p + ggplot2::stat_density_2d(data = long, ggplot2::aes(a1, a2, colour = group),
                                        bins = 5, linewidth = 0.45, alpha = 0.9)
    } else if (hulls) {
      h <- hull_of(long)
      if (!is.null(h)) p <- p + ggplot2::geom_polygon(data = h, ggplot2::aes(a1, a2, colour = group,
                                                         fill = group), alpha = 0.12, linewidth = 0.5)
    }
    long$hover <- paste0(long$species, "<br>", long$group)
    p <- p + .geom_point_hover(data = long, ggplot2::aes(a1, a2, colour = group, text = hover), size = 1.1,
                                 alpha = 0.7) +
      ggplot2::scale_colour_manual(values = pal, name = NULL) +
      ggplot2::scale_fill_manual(values = pal, name = NULL, guide = "none")
    if (centroids) {
      cc <- cent_of(long); cc$group <- factor(cc$group, levels = lv)
      p <- p + .geom_point_hover(data = cc, ggplot2::aes(a1, a2, fill = group, text = hover),
                                 shape = 23, size = 3.6, colour = "black", stroke = 0.9,
                                 inherit.aes = FALSE, show.legend = FALSE)
    }
    return(p + theme)
  }
  rf <- .ref_focal(gr, reference, focal, pool_focal)
  if (is.null(x)) x <- group_species(space, groups, reference, focal = focal, pool_focal = pool_focal,
                                     sep = sep, recode = recode, levels = levels,
                                     hull_axes = hull_axes)
  refsp <- rownames(gr$M)[rowSums(gr$M[, rf$reference, drop = FALSE]) > 0]
  rd <- data.frame(a1 = X[refsp, axes[1]], a2 = X[refsp, axes[2]],
                   hover = paste0(refsp, "<br>reference (", rf$ref_label, ")"))
  hr <- rd[grDevices::chull(rd$a1, rd$a2), ]
  f <- x$focal
  f$a1 <- X[f$species, axes[1]]; f$a2 <- X[f$species, axes[2]]
  f$position <- ifelse(is.na(f$inside_ref_hull), "unknown",
                       ifelse(f$inside_ref_hull, "inside reference space", "outside reference space"))
  lv <- unique(f$group); pal <- .group_palette(lv)
  f$hover <- paste0(f$species, "<br>", f$group, ", ", f$position,
                    if (!is.null(f$FNov_raw)) sprintf("<br>FNov = %.3g", f$FNov_raw) else "")
  f$group <- factor(f$group, levels = lv)
  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(data = hr, ggplot2::aes(a1, a2), fill = "grey88", colour = "grey60", alpha = 0.6) +
    .geom_point_hover(data = rd, ggplot2::aes(a1, a2, text = hover), colour = "grey55", size = 0.9,
                      alpha = 0.5) +
    .geom_point_hover(data = f, ggplot2::aes(a1, a2, colour = group, shape = position, text = hover),
                        size = 1.9, alpha = 0.85, stroke = 0.8) +
    ggplot2::scale_colour_manual(values = pal, name = NULL) +
    ggplot2::scale_shape_manual(values = c("inside reference space" = 1,
                                           "outside reference space" = 16, "unknown" = 4), name = NULL)
  if (centroids) {
    cf <- cent_of(f[, c("group", "a1", "a2")]); cf$group <- factor(cf$group, levels = lv)
    cr <- data.frame(a1 = mean(rd$a1), a2 = mean(rd$a2), n = nrow(rd))
    cr$hover <- sprintf("Centroid of the reference (%s)<br>%d species<br>(%.3g, %.3g)",
                        rf$ref_label, cr$n, cr$a1, cr$a2)
    p <- p + .geom_point_hover(data = cr, ggplot2::aes(a1, a2, text = hover), shape = 23, size = 3.6,
                               fill = "grey45", colour = "black", stroke = 0.9, inherit.aes = FALSE) +
      .geom_point_hover(data = cf, ggplot2::aes(a1, a2, fill = group, text = hover), shape = 23,
                        size = 3.6, colour = "black", stroke = 0.9, inherit.aes = FALSE,
                        show.legend = FALSE) +
      ggplot2::scale_fill_manual(values = pal, name = NULL, guide = "none")
  }
  if (label > 0) {
    lab <- utils::head(f[order(-f$FNov_raw), ], label)
    p <- p + ggplot2::geom_text(data = lab, ggplot2::aes(a1, a2, label = species),
                                size = 2.6, vjust = -0.8)
  }
  p + ggplot2::labs(subtitle = sprintf("Grey: reference (%s); hull position on the first %d axes",
                                       rf$ref_label, x$settings$hull_axes)) + theme
}

#' Centroids of species groups in functional spaces
#'
#' Position of the centroid (mean coordinates of its species, as used for
#' `FDis` and `centroid_dist` in [group_overlap()]) of each group on the
#' axes of one or several functional spaces. A species in several groups
#' counts in each of them.
#'
#' @inheritParams group_species
#' @param space An `fs_space`, or a list of them (output of
#'   [build_spaces()]).
#' @param n_axes Number of axes reported (`NULL`: all axes of the space).
#' @return Data frame with `space`, `group`, `n` (number of species) and
#'   `PCoA1`, `PCoA2`... (centroid coordinates; `NA` beyond the axes of a
#'   space).
#' @examples
#' data(birds_traits); data(birds_trait_info)
#' tr <- drop_incomplete(prepare_traits(birds_traits, birds_trait_info, "species"))
#' s <- build_space(tr, "Locomotion")
#' g <- setNames(birds_traits$iucn, birds_traits$species)
#' group_centroids(s, g, n_axes = 3)
#' @export
group_centroids <- function(space, groups, sep = NULL, recode = NULL, levels = NULL, n_axes = NULL) {
  if (!inherits(space, "fs_space")) {
    out <- lapply(space, group_centroids, groups = groups, sep = sep, recode = recode,
                  levels = levels, n_axes = n_axes)
    nm <- unique(unlist(lapply(out, names)))
    out <- lapply(out, function(z) { z[setdiff(nm, names(z))] <- NA; z[nm] })
    return(do.call(rbind, out))
  }
  X <- .space_coords(space, n_axes)
  gr <- .as_groups(groups, rownames(X), sep = sep, recode = recode, levels = levels)
  rows <- lapply(colnames(gr$M), function(g) {
    i <- gr$M[, g]
    if (!any(i)) return(NULL)
    data.frame(space = space$name, group = g, n = sum(i),
               t(colMeans(X[i, , drop = FALSE])), check.names = FALSE, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  colnames(out)[-(1:3)] <- paste0("PCoA", seq_len(ncol(X)))
  rownames(out) <- NULL
  out
}

#' Plot group overlap or separation across functional spaces
#'
#' Heat map of a pairwise metric of [groups_all()] (`x$pairs`: `overlap`,
#' `beta_jtu`, `ses_overlap`, `centroid_dist`, ...): one column per space,
#' one row per pair of groups, to see which functional dimensions separate
#' which groups.
#'
#' @param x Output of [groups_all()].
#' @param metric Column of `x$pairs`.
#' @return A ggplot.
#' @export
plot_group_overlap <- function(x, metric = "overlap") {
  p <- x$pairs
  if (is.null(p) || !nrow(p)) stop("No pair of groups.", call. = FALSE)
  if (!metric %in% names(p)) stop("Metric '", metric, "' not available.", call. = FALSE)
  p$value <- p[[metric]]
  p$pair <- paste(p$group1, "vs", p$group2)
  p$space <- factor(p$space, levels = unique(p$space))
  p$pair <- factor(p$pair, levels = rev(unique(p$pair)))
  p$label <- ifelse(is.finite(p$value), formatC(p$value, digits = 2, format = "f"), "")
  div <- grepl("^ses_", metric)
  g <- ggplot2::ggplot(p, ggplot2::aes(space, pair, fill = value)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = label), size = if (nlevels(p$space) > 3) 2.5 else 3)
  g <- if (div) g + ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                                  midpoint = 0, name = metric) else
    g + ggplot2::scale_fill_gradient(low = "#F7FBFF", high = "#2171B5", name = metric,
                                     limits = if (metric %in% c("overlap", "beta_jac", "beta_jtu", "beta_jne",
                                                                 "prop_V1_shared", "prop_V2_shared",
                                                                 "overlap_kde")) c(0, 1))
  g + ggplot2::labs(x = NULL, y = NULL) + ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   plot.background = ggplot2::element_rect(fill = "white", colour = NA),
                   legend.position = "bottom",
                   axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
}

## Internal: "new = old1, old2" lines -> named list for `recode`
.parse_recode <- function(txt) {
  if (is.null(txt) || !nzchar(trimws(txt))) return(NULL)
  l <- trimws(unlist(strsplit(txt, "\n", fixed = TRUE)))
  l <- l[grepl("=", l, fixed = TRUE)]
  if (!length(l)) return(NULL)
  nm <- trimws(sub("=.*$", "", l))
  val <- lapply(sub("^[^=]*=", "", l), function(v) {
    v <- trimws(unlist(strsplit(v, "[,;]"))); v[nzchar(v)]
  })
  ok <- nzchar(nm) & lengths(val) > 0
  if (!any(ok)) return(NULL)
  stats::setNames(val[ok], nm[ok])
}
