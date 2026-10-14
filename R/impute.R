#' Impute missing trait values with missForest
#'
#' Follows the imputation strategy of Toussaint et al. (2026): continuous
#' traits are imputed on their transformed and standardised scale, together
#' with the other traits and, optionally, phylogenetic eigenvectors used as
#' additional predictors (Penone et al. 2014, <doi:10.1111/2041-210X.12232>).
#'
#' Type-specific handling: continuous, ordinal and circular traits are
#' imputed as numeric variables (ordinal values are then rounded to the
#' nearest observed rank); nominal and binary traits as factors; fuzzy-coded
#' blocks as numeric columns, re-normalised to proportions afterwards;
#' multichoice blocks as presence probabilities of each modality (regression
#' forests on 0/1 columns): a modality is retained when its probability is
#' >= 0.5 and, when none is, the most probable modality is retained, so that
#' every imputed species has at least one modality (Jaccard distances are
#' undefined otherwise).
#'
#' @param x An `fs_traits` object from [prepare_traits()].
#' @param phylo Optional phylogeny (class `phylo`, tips named as species).
#'   Its first `n_phylo_axes` principal coordinates (classical MDS of the
#'   cophenetic distance) are used as predictors.
#' @param covariates Optional data frame of additional predictors with
#'   species as row names (numeric columns, or factors with < 53 levels),
#'   e.g. phylogenetic eigenvectors or taxonomic/environmental variables
#'   already present in the trait table (`Eigen.1`, ..., as in the published
#'   bird analysis). Predictors are used to impute traits but are not traits:
#'   they do not enter functional distances. They can be combined with
#'   `phylo`.
#' @param n_phylo_axes Number of phylogenetic eigenvectors (10 in the
#'   original analysis).
#' @param seed Random seed.
#' @param ... Further arguments to [missForest::missForest()] (e.g.
#'   `ntree`, `maxiter`).
#'
#' @return The `fs_traits` object with imputed `data`, `imputed = TRUE` and
#'   an `imputation` element holding the out-of-bag error and a logical
#'   matrix flagging imputed cells.
#' @export
impute_traits <- function(x, phylo = NULL, covariates = NULL, n_phylo_axes = 10,
                          seed = 20251029, ...) {
  stopifnot(inherits(x, "fs_traits"))
  if (!requireNamespace("missForest", quietly = TRUE)) {
    stop("Package 'missForest' is required: install.packages('missForest')")
  }
  ti <- x$trait_info
  dat <- x$data
  was_na <- is.na(dat)
  if (!any(was_na)) {
    message("No missing values: nothing to impute.")
    x$imputation <- list(oob = NA, imputed_cells = was_na)
    return(x)
  }

  # species with no information at all cannot be imputed meaningfully
  empty <- rowSums(!was_na) == 0
  if (any(empty)) stop(sum(empty), " species have no trait value at all; remove them first.")

  xm <- dat
  q_scale <- list()
  for (i in seq_len(nrow(ti))) {
    tr <- ti$trait[i]
    if (ti$type[i] == "Q") {
      m <- mean(xm[[tr]], na.rm = TRUE); s <- stats::sd(xm[[tr]], na.rm = TRUE)
      if (is.na(s) || s == 0) s <- 1
      q_scale[[tr]] <- c(m, s)
      xm[[tr]] <- (xm[[tr]] - m) / s
    } else if (ti$type[i] == "D") {
      xm[[tr]] <- factor(xm[[tr]], levels = c(0, 1))
    }
  }
  # missForest requires syntactic, unique names
  orig_names <- names(xm)
  names(xm) <- make.names(orig_names, unique = TRUE)

  extra <- NULL
  if (!is.null(phylo)) {
    if (!inherits(phylo, "phylo")) stop("`phylo` must be an object of class 'phylo'.")
    common <- intersect(x$species, phylo$tip.label)
    if (length(common) < length(x$species)) {
      warning(length(x$species) - length(common),
              " species absent from the phylogeny get NA eigenvectors (imputed too).")
    }
    tree <- ape::keep.tip(phylo, common)
    ev <- stats::cmdscale(ape::cophenetic.phylo(tree), k = min(n_phylo_axes, length(common) - 1))
    colnames(ev) <- paste0("PhyloEigen.", seq_len(ncol(ev)))
    extra <- as.data.frame(ev)[x$species, , drop = FALSE]
  }
  if (!is.null(covariates)) {
    covariates <- as.data.frame(covariates)
    miss_sp <- setdiff(x$species, rownames(covariates))
    if (length(miss_sp)) stop(length(miss_sp), " species absent from the row names of `covariates` (e.g. '",
                              miss_sp[1], "'). Use species names as row names.")
    both <- intersect(names(covariates), x$trait_info$trait)
    if (length(both)) stop("Column(s) used both as trait and as predictor: ", paste(both, collapse = ", "))
    cv <- covariates[x$species, , drop = FALSE]
    cv[] <- lapply(cv, function(v) if (is.character(v)) factor(v) else v)
    big <- vapply(cv, function(v) is.factor(v) && nlevels(v) > 53, logical(1))
    if (any(big)) {
      warning("Covariate(s) with > 53 levels dropped (randomForest limit): ",
              paste(names(cv)[big], collapse = ", "))
      cv <- cv[, !big, drop = FALSE]
    }
    names(cv) <- make.names(paste0("cov.", names(cv)), unique = TRUE)
    extra <- if (is.null(extra)) cv else cbind(extra, cv)
  }
  if (!is.null(extra)) xm <- cbind(xm, extra)

  set.seed(seed)
  mf <- missForest::missForest(xmis = xm, ...)
  imp <- mf$ximp[, seq_along(orig_names), drop = FALSE]
  names(imp) <- orig_names

  for (i in seq_len(nrow(ti))) {
    tr <- ti$trait[i]
    switch(ti$type[i],
      Q = { imp[[tr]] <- imp[[tr]] * q_scale[[tr]][2] + q_scale[[tr]][1] },
      O = {
        obs <- sort(unique(stats::na.omit(dat[[tr]])))
        imp[[tr]] <- obs[vapply(imp[[tr]], function(v) which.min(abs(obs - v)), integer(1))]
      },
      D = { imp[[tr]] <- as.numeric(as.character(imp[[tr]])) },
      F = { imp[[tr]] <- pmax(imp[[tr]], 0) }
    )
  }
  # multichoice blocks: modalities imputed as presence probabilities
  # (regression forests on 0/1); a modality is retained if p >= 0.5 and, if
  # none is, the most probable one, so that every species keeps >= 1 modality
  for (b in unique(ti$block[ti$type == "B"])) {
    cols <- ti$trait[ti$block == b]
    m <- as.matrix(imp[, cols, drop = FALSE])
    m01 <- (m >= 0.5) * 1
    none <- rowSums(m01) == 0
    if (any(none)) m01[cbind(which(none), max.col(m[none, , drop = FALSE], ties.method = "first"))] <- 1
    imp[, cols] <- m01
  }
  for (b in unique(ti$block[ti$type == "F"])) {
    cols <- ti$trait[ti$block == b]
    m <- as.matrix(imp[, cols, drop = FALSE])
    rs <- rowSums(m); rs[rs == 0] <- 1
    imp[, cols] <- m / rs
  }
  # never overwrite observed values
  for (tr in orig_names) {
    obs <- !was_na[, tr]
    imp[[tr]][obs] <- dat[[tr]][obs]
  }
  rownames(imp) <- rownames(dat)

  x$data <- imp
  x$imputed <- TRUE
  x$missing <- colMeans(is.na(imp))
  x$imputation <- list(oob = mf$OOBerror, imputed_cells = was_na,
                       phylo_axes = if (!is.null(phylo)) n_phylo_axes else 0,
                       covariates = if (!is.null(covariates)) names(covariates) else character(0))
  x
}
