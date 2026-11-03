## ---------------------------------------------------------------------------
## Species groups within assemblages (sites, basins, countries, plots...)
##
## A species can belong to different groups in different assemblages (native
## in one basin and introduced in another, exploited in one country only,
## migrant in one region and resident in another...). Groups are then given
## per assemblage, and the contribution of focal species is measured within
## each assemblage.
## ---------------------------------------------------------------------------

#' Assemblages from a long table
#'
#' Converts a long table (one row per species occurrence in an assemblage,
#' with an optional group label per occurrence) into the assemblage x
#' species matrix and the assemblage x species matrix of labels used by
#' [group_assemblage()].
#'
#' @param data Data frame.
#' @param site_col,species_col Columns of assemblages and species.
#' @param group_col Optional column of the group of the species in the
#'   assemblage (e.g. native / introduced in that basin).
#' @param abundance_col Optional column of abundances (otherwise presence).
#' @return A list: `comm` (numeric matrix, assemblages x species) and
#'   `groups` (character matrix of the same dimensions, `NA` where absent;
#'   `NULL` without `group_col`).
#' @export
community_from_long <- function(data, site_col = "site", species_col = "species",
                                group_col = NULL, abundance_col = NULL) {
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  for (cl in c(site_col, species_col, group_col, abundance_col))
    if (!cl %in% names(data)) stop("Column '", cl, "' not found.", call. = FALSE)
  s <- as.character(data[[site_col]]); sp <- as.character(data[[species_col]])
  ok <- !is.na(s) & !is.na(sp) & nzchar(s) & nzchar(sp)
  s <- s[ok]; sp <- sp[ok]
  si <- unique(s); spp <- unique(sp)
  comm <- matrix(0, length(si), length(spp), dimnames = list(si, spp))
  ab <- if (is.null(abundance_col)) rep(1, length(s)) else as.numeric(data[[abundance_col]][ok])
  ab[!is.finite(ab)] <- 0
  idx <- cbind(match(s, si), match(sp, spp))
  agg <- stats::aggregate(ab, by = list(i = idx[, 1], j = idx[, 2]), FUN = sum)
  comm[cbind(agg$i, agg$j)] <- agg$x
  grp <- NULL
  if (!is.null(group_col)) {
    grp <- matrix(NA_character_, length(si), length(spp), dimnames = list(si, spp))
    grp[idx] <- as.character(data[[group_col]][ok])
  }
  list(comm = comm, groups = grp)
}

#' Contribution of focal groups within each assemblage
#'
#' Site-level counterpart of [group_contribution()] that generalises
#' [nonnative_assemblage()]: in each assemblage, the reference species alone
#' (`_ref`) are compared with reference + focal species (`_all`) for `S`,
#' `FRic` (convex hull volume relative to the hull of all species of the
#' space), `FDis`, `FOri` and `MPD`; `d_` = contribution of the focal
#' species; `overlap` = hull overlap (intersection / union) of focal and
#' reference species in the assemblage; `p_focal_outside` = proportion of
#' focal species outside the reference hull. Groups can differ between
#' assemblages (`groups` as an assemblage x species matrix of labels).
#' With `null = "permute"`, the same number of species of the assemblage,
#' drawn at random, plays the role of the focal species (`ses_`, `p_`).
#'
#' @param space An `fs_space`.
#' @param comm Assemblage x species matrix (presence / abundance; > 0 =
#'   present) or the output of [community_from_long()].
#' @param groups Groups of species: a named vector (the same in every
#'   assemblage) or an assemblage x species matrix of labels. Ignored when
#'   `comm` comes from [community_from_long()] with groups.
#' @inheritParams group_species
#' @param min_ref Minimum number of reference species in an assemblage.
#' @param null `"none"` or `"permute"`.
#' @param nrep Null replicates.
#' @param seed Optional random seed.
#' @return Data frame, one row per assemblage and focal group.
#' @export
group_assemblage <- function(space, comm, groups = NULL, reference, focal = NULL,
                             pool_focal = FALSE, sep = NULL, recode = NULL, levels = NULL,
                             n_axes = NULL, hull_axes = 4, min_ref = 3,
                             null = c("none", "permute"), nrep = 99, seed = NULL) {
  null <- match.arg(null)
  if (!is.null(seed)) set.seed(seed)
  if (is.list(comm) && !is.data.frame(comm) && !is.null(comm$comm)) {
    if (!is.null(comm$groups)) groups <- comm$groups
    comm <- comm$comm
  }
  comm <- as.matrix(comm)
  if (is.null(groups)) stop("Give the groups of species.", call. = FALSE)
  X <- .space_coords(space, n_axes)
  sp_all <- rownames(X)
  inspace <- intersect(colnames(comm), sp_all)
  if (!length(inspace)) stop("No species of 'comm' in the space.", call. = FALSE)
  nmiss <- length(setdiff(colnames(comm)[colSums(comm > 0, na.rm = TRUE) > 0], sp_all))
  ha <- min(hull_axes, ncol(X))
  vol_pool <- .hull_vol(X[, seq_len(ha), drop = FALSE])
  per_site <- is.matrix(groups) || is.data.frame(groups)
  if (per_site) groups <- as.matrix(groups)
  dn <- c("S", "FRic", "FDis", "FOri", "MPD")
  rows <- list()
  for (j in seq_len(nrow(comm))) {
    site <- rownames(comm)[j]
    present <- inspace[comm[j, inspace] > 0 & !is.na(comm[j, inspace])]
    if (length(present) < 2) next
    lab <- if (per_site) {
      if (!site %in% rownames(groups)) next
      stats::setNames(as.character(groups[site, present]), present)
    } else stats::setNames(as.character(groups[present]), present)
    gr <- tryCatch(.as_groups(lab, present, sep = sep, recode = recode), error = function(e) NULL)
    if (is.null(gr)) next
    lv <- if (is.null(levels)) gr$levels else intersect(levels, gr$levels)
    if (!any(reference %in% lv)) next
    M <- gr$M[, lv, drop = FALSE]
    foc_lv <- if (is.null(focal)) setdiff(lv, reference) else intersect(setdiff(focal, reference), lv)
    if (!length(foc_lv)) next
    fl <- if (pool_focal) stats::setNames(list(foc_lv), paste(foc_lv, collapse = " + ")) else
      stats::setNames(as.list(foc_lv), foc_lv)
    for (fg in names(fl)) {
      st <- .focal_status(M, intersect(reference, lv), fl[[fg]])
      ref <- names(st)[st %in% FALSE]; foc <- names(st)[st %in% TRUE]
      out <- c(n_ref = length(ref), n_focal = length(foc))
      if (length(ref) >= min_ref && length(foc)) {
        mr <- .assemblage_metrics(X[ref, , drop = FALSE], ha, vol_pool)
        ma <- .assemblage_metrics(X[c(ref, foc), , drop = FALSE], ha, vol_pool)
        out <- c(out, stats::setNames(mr, paste0(dn, "_ref")), stats::setNames(ma, paste0(dn, "_all")),
                 stats::setNames(ma - mr, paste0("d_", dn)))
        out["FRic_change_rel"] <- if (is.finite(mr[["FRic"]]) && mr[["FRic"]] > 0)
          (ma[["FRic"]] - mr[["FRic"]]) / mr[["FRic"]] else NA
        out["overlap"] <- if (length(ref) > ha && length(foc) > ha)
          .pair_overlap(X[ref, , drop = FALSE], X[foc, , drop = FALSE], ha, vol_pool)[["overlap"]] else NA
        out["p_focal_outside"] <- if (length(ref) > ha)
          mean(!.in_hull(X[ref, seq_len(ha), drop = FALSE], X[foc, seq_len(ha), drop = FALSE]), na.rm = TRUE) else NA
        if (null == "permute" && nrep > 1) {
          allsp <- c(ref, foc)
          sims <- replicate(nrep, {
            f2 <- sample(allsp, length(foc)); r2 <- setdiff(allsp, f2)
            ma - .assemblage_metrics(X[r2, , drop = FALSE], ha, vol_pool)
          })
          for (m in dn[-1]) {
            obs <- out[[paste0("d_", m)]]; v <- sims[m, ]; v <- v[is.finite(v)]
            out[paste0("ses_", m)] <- if (length(v) > 1 && stats::sd(v) > 0 && is.finite(obs))
              (obs - mean(v)) / stats::sd(v) else NA
            out[paste0("p_", m)] <- if (length(v) > 1 && is.finite(obs))
              min(1, 2 * min((sum(v <= obs) + 1) / (length(v) + 1), (sum(v >= obs) + 1) / (length(v) + 1))) else NA
          }
        }
      }
      rows[[length(rows) + 1]] <- list(site = site, group = fg, v = out)
    }
  }
  if (!length(rows)) stop("No assemblage with both reference and focal species.", call. = FALSE)
  nm <- unique(unlist(lapply(rows, function(r) names(r$v))))
  vals <- do.call(rbind, lapply(rows, function(r) r$v[nm]))
  colnames(vals) <- nm
  res <- data.frame(assemblage = vapply(rows, `[[`, "", "site"), group = vapply(rows, `[[`, "", "group"),
                    vals, row.names = NULL, stringsAsFactors = FALSE, check.names = FALSE)
  attr(res, "n_species_not_in_space") <- nmiss
  res
}

#' Plot the contribution of focal groups across assemblages and spaces
#'
#' Distribution across assemblages of a column of the `sites` table of
#' [groups_all()] (e.g. `FRic_change_rel`, `ses_FRic`, `overlap`), per space
#' and focal group.
#'
#' @param x Output of [groups_all()] with `comm`.
#' @param metric Column of `x$sites`.
#' @return A ggplot.
#' @export
plot_group_sites <- function(x, metric = "FRic_change_rel") {
  a <- x$sites
  if (is.null(a) || !metric %in% names(a)) stop("Metric '", metric, "' not available.", call. = FALSE)
  a$x_value <- a[[metric]]
  a$space <- factor(a$space, levels = unique(a$space))
  p <- ggplot2::ggplot(a, ggplot2::aes(space, x_value, fill = group)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey50")
  if (grepl("^ses_", metric)) p <- p + ggplot2::geom_hline(yintercept = c(-1.96, 1.96),
                                                           linetype = "dashed", colour = "grey60")
  p + ggplot2::geom_boxplot(outlier.size = 0.8, alpha = 0.7) +
    ggplot2::scale_fill_manual(values = .group_palette(unique(a$group)), name = NULL) +
    ggplot2::labs(x = NULL, y = metric) + ggplot2::theme_classic(base_size = 11)
}
