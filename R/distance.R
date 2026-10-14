#' Generalised Gower distance between species for mixed trait types
#'
#' Computes the functional distance used in Toussaint et al. (2026): the
#' generalised Gower distance of Pavoine et al. (2009) as implemented in
#' [ade4::dist.ktab()]. Each single-column trait, and each multi-column
#' fuzzy (`F`) or multichoice (`B`) block, is one *variable*; per-variable
#' distances are combined as
#' \deqn{d_{ij} = \sqrt{\sum_v w_v \delta_{ijv} d_{ijv}^2 / \sum_v w_v \delta_{ijv}}}{d_ij = sqrt( sum_v w_v delta_ijv d_ijv^2 / sum_v w_v delta_ijv )}
#' where \eqn{\delta_{ijv}} is 0 when variable v is missing for species i or j
#' (or undefined, e.g. Jaccard with double absence), 1 otherwise.
#'
#' Two weighting schemes are available:
#' * `"variable"` (default, as in the published analysis): every variable
#'   has weight `w_v` (1 unless set in `trait_info$weight`). A dimension
#'   described by many traits then weighs more in a combined space.
#' * `"dimension"`: the squared distance is first averaged within each
#'   functional dimension, then dimensions are averaged with equal weight,
#'   so that each dimension contributes equally whatever its number of traits
#'   (cf. de Bello et al. 2021, <doi:10.1111/2041-210X.13537>).
#'
#' Species pairs sharing no non-missing trait make the distance undefined;
#' [build_space()] excludes species without any value for the selected
#' dimensions before calling this function.
#'
#' * `"dimension_scaled"`: as `"dimension"`, but the squared distances of
#'   each dimension are first divided by their mean, so that every dimension
#'   contributes equally on average. Distances of categorical traits (simple
#'   matching, Jaccard) are on average much larger than range-scaled
#'   distances of continuous traits, so that, with the two other options, a
#'   dimension of categorical traits can dominate a combined space whatever
#'   its number of traits (in the spirit of de Bello et al. 2021,
#'   <doi:10.1111/2041-210X.13537>). Distances are then no longer bounded by
#'   1.
#'
#' With equal weights and `weighting = "variable"` the result is identical to
#' a single call to `ade4::dist.ktab()` on all traits.
#'
#' Engines: `"ade4"` (default, as in the published analysis) calls
#' `ade4::dist.ktab()` on the same trait tables as [PCoAGraph()] (one table
#' per dimension x type); with `weighting = "dimension"`, `dist.ktab()` is
#' called once per dimension and the squared distances are averaged over
#' dimensions; with unequal trait weights, once per trait. `"fast"` gives
#' the same values with vectorised matrix algebra and is useful for very
#' large datasets (the ordinal and circular branches of `dist.ktab()` loop
#' over all pairs of species in R); the package tests check that both
#' engines agree to numerical precision for every trait type, with and
#' without missing values. As in `ade4`, a
#' species with all modalities of a multichoice (`B`) trait at 0 has an
#' undefined (missing) Jaccard distance for that trait.
#'
#' @param x An `fs_traits` object ([prepare_traits()] / [impute_traits()]).
#' @param dimensions Character vector of dimensions to include (default:
#'   all).
#' @param weighting `"variable"`, `"dimension"` or `"dimension_scaled"`
#'   (see Details).
#' @param option Scaling of quantitative variables passed to
#'   `ade4::dist.ktab()` (`"scaledBYrange"`, as in the original analysis;
#'   other options use the `"ade4"` engine).
#' @param engine `"ade4"` (default) or `"fast"` (see Details).
#'
#' @return A `dist` object (species labels) with attributes `variables`
#'   (table of variables used, their dimension, type and weight) and
#'   `weighting`.
#' @references Pavoine S., Vallet J., Dufour A.-B., Gachet S., Daniel H.
#'   (2009) On the challenge of treating various types of variables:
#'   application for improving the measurement of functional diversity.
#'   Oikos 118: 391-402.
#' @export
#' @examples
#' data(birds_traits); data(birds_trait_info)
#' tr <- drop_incomplete(prepare_traits(birds_traits, birds_trait_info, "species"))
#' d <- trait_distance(tr, dimensions = "Diet")
#' summary(as.vector(d))
trait_distance <- function(x, dimensions = NULL,
                           weighting = c("variable", "dimension", "dimension_scaled"),
                           option = c("scaledBYrange", "scaledBYsd", "noscale"),
                           engine = c("ade4", "fast")) {
  stopifnot(inherits(x, "fs_traits"))
  weighting <- match.arg(weighting)
  option <- match.arg(option)
  engine <- match.arg(engine)
  if (option != "scaledBYrange") engine <- "ade4"
  ti <- x$trait_info
  if (!is.null(dimensions)) {
    unk <- setdiff(dimensions, ti$dimension)
    if (length(unk)) stop("Unknown dimension(s): ", paste(unk, collapse = ", "))
    ti <- ti[ti$dimension %in% dimensions, , drop = FALSE]
  }
  vars <- .build_variables(x$data, ti)
  if (length(vars$tables) == 0) stop("No usable trait in the selected dimension(s).")
  vt <- vars$info
  n <- nrow(x$data)

  equal_w <- length(unique(vt$weight)) == 1
  ktab_dist <- function(t) {
    bl <- .make_blocks(x$data, t)
    ade4::dist.ktab(ade4::ktab.list.df(bl$x), bl$groups, scann = FALSE, option = option)
  }
  equal_within <- all(tapply(vt$weight, vt$dimension, function(w) length(unique(w)) == 1))
  if (engine == "ade4" && weighting == "variable" && equal_w) {
    d <- suppressWarnings(ktab_dist(ti))
  } else if (engine == "ade4" && weighting != "variable" && equal_within) {
    dims <- unique(vt$dimension)
    per_dim <- lapply(dims, function(dm) as.vector(suppressWarnings(ktab_dist(ti[ti$dimension == dm, ]))))
    if (weighting == "dimension_scaled") per_dim <- lapply(per_dim, .unit_ms)
    d <- sqrt(.wmean_sq(per_dim, rep(1, length(dims))))
    attr(d, "Size") <- n; attr(d, "Diag") <- FALSE; attr(d, "Upper") <- FALSE
    class(d) <- "dist"
  } else {
    single <- lapply(seq_along(vars$tables), function(i) {
      if (engine == "fast") return(.var_dist(vars$tables[[i]], vt$type[i]))
      as.vector(ade4::dist.ktab(ade4::ktab.list.df(vars$tables[i]), vt$type[i], option = option))
    })
    if (weighting == "variable") {
      d2 <- .wmean_sq(single, vt$weight)
    } else {
      dims <- unique(vt$dimension)
      per_dim <- lapply(dims, function(dm) {
        idx <- which(vt$dimension == dm)
        sqrt(.wmean_sq(single[idx], vt$weight[idx]))
      })
      if (weighting == "dimension_scaled") per_dim <- lapply(per_dim, .unit_ms)
      d2 <- .wmean_sq(per_dim, rep(1, length(dims)))
    }
    d <- sqrt(d2)
    attr(d, "Size") <- n
    attr(d, "Diag") <- FALSE
    attr(d, "Upper") <- FALSE
    class(d) <- "dist"
  }
  attr(d, "Labels") <- rownames(x$data)
  if (anyNA(d)) {
    m <- as.matrix(d)
    worst <- names(sort(colSums(is.na(m)), decreasing = TRUE))[1]
    stop(sum(is.na(d)), " species pairs share no non-missing trait (e.g. '", worst,
         "'). Impute missing values (impute_traits) or drop incomplete species (drop_incomplete).")
  }
  attr(d, "variables") <- vt
  attr(d, "weighting") <- weighting
  d
}

# Per-variable distance (lower triangle, `dist` order), vectorised.
# Reproduces ade4::dist.ktab(option = "scaledBYrange") with its
# non-interactive defaults; NaN marks undefined / missing pairs.
.var_dist <- function(tab, type, tol = 1e-8) {
  out <- .var_dist0(tab, type)
  # ade4 sets (squared, for F/B/D/C blocks) distances below `tol` to 0
  if (type %in% c("Q", "O", "N")) out[!is.na(out) & out < tol] <- 0
  else out[!is.na(out) & out^2 < tol] <- 0
  out
}

.var_dist0 <- function(tab, type) {
  lt <- function(m) m[lower.tri(m)]
  jac <- function(X) {
    X <- as.matrix(X)
    a <- tcrossprod(X); rs <- rowSums(X)
    abc <- outer(rs, rs, "+") - a
    list(d = sqrt(1 - a / abc), abc = abc, rs = rs)
  }
  switch(type,
    Q = {
      v <- tab[[1]]
      as.vector(stats::dist(v, "manhattan")) / diff(range(v, na.rm = TRUE))
    },
    O = {
      r <- rank(tab[[1]], na.last = "keep")
      r <- (r - min(r, na.rm = TRUE)) / (max(r, na.rm = TRUE) - min(r, na.rm = TRUE))
      as.vector(stats::dist(r, "manhattan"))
    },
    N = {
      f <- as.integer(factor(tab[[1]]))
      as.numeric(as.vector(stats::dist(f, "manhattan")) > 0)
    },
    D = {
      j <- jac(tab)
      j$d[j$abc == 0] <- NaN
      lt(j$d)
    },
    B = {
      j <- jac(tab)
      z <- is.na(j$rs) | j$rs == 0
      j$d[z, ] <- NaN; j$d[, z] <- NaN
      lt(j$d)
    },
    F = {
      P <- as.matrix(tab)
      nr <- sqrt(rowSums(P^2))
      s <- tcrossprod(P) / outer(nr, nr)
      lt(sqrt(pmax(1 - s, 0)))
    },
    C = {
      mx <- attr(tab, "max")
      v <- tab[[1]] / mx
      dd <- as.vector(stats::dist(v, "manhattan"))
      d2 <- 2 * pmin(dd, 1 - dd)
      if (abs(mx / 2 - floor(mx / 2)) > 1e-8) d2 <- d2 * mx / (mx - 1)
      sqrt(d2)
    }
  )
}

# rescale a distance vector to a mean squared distance of 1
.unit_ms <- function(v) v / sqrt(mean(v^2, na.rm = TRUE))

# weighted, NA-aware mean of squared distances; `ds` = list of distance vectors
.wmean_sq <- function(ds, w) {
  num <- 0; den <- 0
  for (i in seq_along(ds)) {
    ok <- !is.na(ds[[i]])
    num <- num + ifelse(ok, w[i] * ds[[i]]^2, 0)
    den <- den + w[i] * ok
  }
  out <- num / den
  out[den == 0] <- NA
  out
}

# Build one ade4-ready table per variable (single trait or block)
.build_variables <- function(dat, ti) {
  tables <- list(); info <- list()
  for (b in unique(ti$block)) {
    r <- ti[ti$block == b, , drop = FALSE]
    type <- r$type[1]
    df <- dat[, r$trait, drop = FALSE]
    if (type %in% c("Q", "O", "C")) {
      if (length(unique(stats::na.omit(df[[1]]))) < 2) {
        warning("Trait '", b, "' is constant and was ignored.", call. = FALSE); next
      }
    }
    if (type %in% c("N", "D")) {
      if (length(unique(stats::na.omit(as.character(df[[1]])))) < 2) {
        warning("Trait '", b, "' is constant and was ignored.", call. = FALSE); next
      }
    }
    tab <- switch(type,
      Q = , O = , N = , D = df,
      F = suppressWarnings(ade4::prep.fuzzy(df, ncol(df), labels = b)),
      B = ade4::prep.binary(df, ncol(df), labels = b),
      C = {
        rmin <- if (is.na(r$circ_min[1])) min(df[[1]], na.rm = TRUE) else r$circ_min[1]
        rmax <- if (is.na(r$circ_max[1])) max(df[[1]], na.rm = TRUE) else r$circ_max[1]
        ade4::prep.circular(df, rangemin = rmin, rangemax = rmax)
      }
    )
    tables[[b]] <- tab
    info[[b]] <- data.frame(variable = b, dimension = r$dimension[1], type = type,
                            n_columns = nrow(r), weight = r$weight[1],
                            stringsAsFactors = FALSE)
  }
  list(tables = tables, info = do.call(rbind, info))
}
