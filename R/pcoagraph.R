#' PCoA of a list of trait tables of different types (published core function)
#'
#' Core function of Toussaint et al. (2026, `utils/functions.R`), kept with
#' the same arguments and output so that the published scripts can use it
#' unchanged. Trait tables of different types are combined with
#' [ade4::ktab.list.df()], the generalised Gower distance between species is
#' computed with [ade4::dist.ktab()] (`option = "scaledBYrange"`), and the
#' distance is ordinated with [ape::pcoa()] (for 3000 species or more without
#' correction, an equivalent memory-lean implementation is used, see
#' `options(DisentFunSpaces.lean_pcoa_n)`). Traits are then related to the
#' axes with [vegan::envfit()] (`varPlot`) and with correlations between each
#' trait and each PCoA axis (`PCoACor`, function `CorrelType` of the original
#' scripts).
#'
#' Differences with the original script: the parallel/progress-bar
#' machinery (`future`, `progressr`) is removed; `CorrelType` also accepts
#' nominal traits (correlation ratio, eta) and missing values (pairwise
#' complete observations); a correction for negative eigenvalues can be
#' requested; a precomputed distance can be supplied (used by
#' [build_space()] for dimension or trait weighting).
#'
#' @param x List of data frames, one per block of traits of the same type
#'   (species in rows, in the same order in every table). Fuzzy-coded (`F`)
#'   and multichoice (`B`) tables must be prepared with
#'   [ade4::prep.fuzzy()] / [ade4::prep.binary()], circular (`C`) tables with
#'   [ade4::prep.circular()].
#' @param datax Data frame of all traits, used by [vegan::envfit()].
#' @param namesx Species names (row names of the PCoA coordinates).
#' @param groups Character vector giving the type of each table of `x`
#'   (`"Q"`, `"O"`, `"N"`, `"D"`, `"F"`, `"B"`, `"C"`; see [trait_types()]).
#' @param saveFile Path of an `.rds` file where the result is saved, or `NA`.
#' @param correction Correction for negative eigenvalues passed to
#'   [ape::pcoa()] (`"none"` in the original analysis).
#' @param envfit_perm Number of permutations of [vegan::envfit()] (999 =
#'   envfit default, as in the original analysis; 0 to skip tests on large
#'   datasets).
#' @param envfit_axes Axes used by [vegan::envfit()] (envfit default:
#'   first two).
#' @param dist Optional precomputed `dist` object; if supplied, `x` and
#'   `groups` are only used for `PCoACor`.
#'
#' @return A list with `PCoA` (output of [ape::pcoa()], with species as row
#'   names of `vectors`), `varPlot` (envfit vector scores and a `Species`
#'   column holding trait names, as in the original), `PCoACor` (trait x axis
#'   correlation matrix) and `dist` (the functional distance).
#' @references Pavoine S., Vallet J., Dufour A.-B., Gachet S., Daniel H.
#'   (2009) Oikos 118: 391-402.
#' @export
#' @examples
#' data(birds_traits)
#' tr <- birds_traits[stats::complete.cases(birds_traits), ]
#' morpho <- log10(tr[, c("body_mass", "wing_length", "tarsus_length")])
#' diet <- ade4::prep.fuzzy(tr[, grep("^diet_", names(tr))], 10, labels = "diet")
#' res <- PCoAGraph(x = list(morpho, diet), datax = cbind(morpho, diet),
#'                  namesx = tr$species, groups = c("Q", "F"), envfit_perm = 0)
#' head(res$PCoA$values)
#' res$PCoACor[, 1:3]
PCoAGraph <- function(x, datax, namesx, groups, saveFile = NA,
                      correction = "none", envfit_perm = 999, envfit_axes = c(1, 2),
                      dist = NULL) {
  # --- Checks ---
  if (!is.list(x) || !all(vapply(x, is.data.frame, logical(1L)))) {
    stop("Argument 'x' must be a list of data frames.")
  }
  if (!is.data.frame(datax)) stop("'datax' must be a data frame.")
  if (length(namesx) != nrow(datax)) {
    stop("Length of 'namesx' must match the number of rows in 'datax'.")
  }
  if (length(groups) != length(x)) stop("'groups' must give one type per table of 'x'.")
  if (!is.null(saveFile) && !is.character(saveFile) && !is.na(saveFile)) {
    stop("'saveFile' must be a character string or NA.")
  }

  # --- Create ktab object and compute inter-block trait distances ---
  if (is.null(dist)) {
    ktabList <- ade4::ktab.list.df(x)
    disTraits <- ade4::dist.ktab(ktabList, groups, scann = FALSE, option = "scaledBYrange")
  } else {
    disTraits <- dist
  }
  attr(disTraits, "Labels") <- namesx
  if (anyNA(disTraits)) {
    m <- as.matrix(disTraits)
    bad <- sum(colSums(is.na(m)) > 0)
    stop(sum(is.na(disTraits)), " species pairs (", bad, " species) share no non-missing trait, ",
         "so their distance is undefined and the PCoA cannot be computed. Impute missing values ",
         "(impute_traits) or remove these species (drop_incomplete).", call. = FALSE)
  }

  # --- Run PCoA ---
  lean_n <- getOption("DisentFunSpaces.lean_pcoa_n", 3000)
  PCoATraits <- if (identical(correction, "none") && length(namesx) >= lean_n) {
    .pcoa_lean(disTraits)        # same result as ape::pcoa, less memory
  } else {
    ape::pcoa(disTraits, correction = correction)
  }
  rownames(PCoATraits$vectors) <- namesx
  if (!is.null(PCoATraits$vectors.cor)) rownames(PCoATraits$vectors.cor) <- namesx

  # --- envfit correlation with traits ---
  vec_env <- if (!is.null(PCoATraits$vectors.cor)) PCoATraits$vectors.cor else PCoATraits$vectors
  choices <- envfit_axes[envfit_axes <= ncol(vec_env)]
  num <- vapply(datax, is.numeric, logical(1))
  spp.scrs <- NULL
  if (any(num) && length(choices) >= 2) {
    plotenvfit <- vegan::envfit(vec_env, datax[, num, drop = FALSE], w = rep(1, nrow(vec_env)),
                                permutations = envfit_perm, choices = choices, na.rm = TRUE)
    spp.scrs <- as.data.frame(vegan::scores(plotenvfit, display = "vectors"))
    spp.scrs$Species <- rownames(spp.scrs)
  }

  # --- Correlation between each trait and each PCoA axis (CorrelType) ---
  PcoaCorrel <- do.call(rbind, lapply(x, function(tab) CorrelType(PCoATraits, tab)))

  res <- list(PCoA = PCoATraits, varPlot = spp.scrs, PCoACor = PcoaCorrel, dist = disTraits)
  if (!is.null(saveFile) && !is.na(saveFile)) saveRDS(res, saveFile)
  res
}

#' Correlations between traits and PCoA axes
#'
#' `CorrelType` of Toussaint et al. (2026): Pearson correlation between each
#' column of a trait table and each PCoA axis. Extended to factors
#' (correlation ratio eta, between 0 and 1) and to missing values (pairwise
#' complete observations).
#'
#' @param pcoa Output of [ape::pcoa()] (uses `vectors.cor` when present).
#' @param type Data frame of traits (species in the same order).
#' @return Matrix traits x axes.
#' @export
CorrelType <- function(pcoa, type) {
  vec <- if (!is.null(pcoa$vectors.cor)) pcoa$vectors.cor else pcoa$vectors
  type <- as.data.frame(type)
  Axis <- matrix(NA_real_, ncol = ncol(vec), nrow = ncol(type),
                 dimnames = list(colnames(type), colnames(vec)))
  for (i in seq_len(ncol(type))) {
    v <- type[[i]]
    for (tr in seq_len(ncol(vec))) {
      y <- vec[, tr]
      Axis[i, tr] <- if (is.factor(v) || is.character(v)) {
        v <- factor(v); ok <- !is.na(v)
        if (nlevels(droplevels(v[ok])) < 2) NA_real_ else
          sqrt(summary(stats::lm(y[ok] ~ v[ok]))$r.squared)
      } else {
        if (isTRUE(stats::sd(v, na.rm = TRUE) > 0))
          stats::cor(y, v, use = "pairwise.complete.obs") else NA_real_
      }
    }
  }
  Axis
}

# Build the list of trait tables (one per dimension x type, as in the
# published analysis), their types and the envfit table, from trait_info.
# Fuzzy (F) / multichoice (B) traits of a dimension share one table prepared
# with prep.fuzzy / prep.binary (one col.block per trait); each circular
# trait gets its own table (prep.circular). Grouping does not change
# dist.ktab distances (each column, or each fuzzy/multichoice block, is one
# variable).
.make_blocks <- function(dat, ti) {
  keep <- vapply(seq_len(nrow(ti)), function(i) {
    if (ti$type[i] %in% c("F", "B")) return(TRUE)
    v <- dat[[ti$trait[i]]]
    ok <- length(unique(stats::na.omit(as.character(v)))) >= 2
    if (!ok) warning("Trait '", ti$trait[i], "' is constant and was ignored.", call. = FALSE)
    ok
  }, logical(1))
  ti <- ti[keep, , drop = FALSE]
  if (nrow(ti) == 0) stop("No usable trait in the selected dimension(s).")

  x <- list(); groups <- character(0)
  for (dm in unique(ti$dimension)) {
    tid <- ti[ti$dimension == dm, , drop = FALSE]
    for (ty in unique(tid$type)) {
      r <- tid[tid$type == ty, , drop = FALSE]
      df <- dat[, r$trait, drop = FALSE]
      if (ty %in% c("F", "B")) {
        bl <- unique(r$block)
        df <- dat[, unlist(lapply(bl, function(b) r$trait[r$block == b])), drop = FALSE]
        sizes <- vapply(bl, function(b) sum(r$block == b), integer(1))
        tab <- if (ty == "F") {
          suppressWarnings(ade4::prep.fuzzy(df, col.blocks = sizes, labels = bl))
        } else {
          ade4::prep.binary(df, col.blocks = sizes, labels = bl)
        }
        x[[paste(dm, ty, sep = "_")]] <- tab; groups <- c(groups, ty)
      } else if (ty == "C") {
        for (j in seq_len(nrow(r))) {
          v <- df[, j, drop = FALSE]
          rmin <- if (is.na(r$circ_min[j])) min(v[[1]], na.rm = TRUE) else r$circ_min[j]
          rmax <- if (is.na(r$circ_max[j])) max(v[[1]], na.rm = TRUE) else r$circ_max[j]
          x[[paste(dm, ty, r$trait[j], sep = "_")]] <- ade4::prep.circular(v, rangemin = rmin, rangemax = rmax)
          groups <- c(groups, ty)
        }
      } else {
        x[[paste(dm, ty, sep = "_")]] <- df; groups <- c(groups, ty)
      }
    }
  }
  datax <- do.call(cbind.data.frame, unname(x))
  names(datax) <- unlist(lapply(x, names))
  vars <- do.call(rbind, lapply(unique(ti$block), function(b) {
    rr <- ti[ti$block == b, ]
    data.frame(variable = b, dimension = rr$dimension[1], type = rr$type[1],
               n_columns = nrow(rr), weight = rr$weight[1], stringsAsFactors = FALSE)
  }))
  list(x = x, groups = groups, datax = datax, trait_info = ti, variables = vars)
}
