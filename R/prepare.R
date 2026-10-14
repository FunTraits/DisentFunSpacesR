#' Prepare a trait table for functional-space analyses
#'
#' Converts each trait to the representation required by its declared type
#' (see [trait_types()]), applies the requested transformations to
#' continuous traits, and rescales fuzzy-coded blocks to proportions.
#' Species whose fuzzy-coded or multichoice block contains only zeros (no
#' modality recorded) get `NA` for that block, since their distance to other
#' species is undefined for it. The
#' result is the input of [impute_traits()], [trait_distance()] and
#' [build_spaces()].
#'
#' @param data Data frame of traits, species in rows.
#' @param trait_info Trait-information table (see [check_trait_info()]).
#' @param species_col Name of the column holding species identifiers. If
#'   `NULL`, row names of `data` are used.
#'
#' @return An object of class `fs_traits`: a list with `data` (prepared
#'   traits, species as row names), `trait_info` (validated table),
#'   `species`, `imputed` (logical) and `missing` (per-trait proportion of
#'   missing values after preparation).
#' @export
#' @examples
#' data(birds_traits); data(birds_trait_info)
#' tr <- prepare_traits(birds_traits, birds_trait_info, species_col = "species")
#' tr
prepare_traits <- function(data, trait_info, species_col = NULL) {
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  if (!is.null(species_col)) {
    if (!species_col %in% names(data)) stop("Column '", species_col, "' not found.")
    sp <- as.character(data[[species_col]])
  } else {
    sp <- rownames(data)
  }
  if (anyNA(sp) || any(sp == "")) stop("Missing species identifiers.")
  if (anyDuplicated(sp)) {
    stop("Duplicated species identifiers (e.g. '", sp[duplicated(sp)][1],
         "'). Aggregate intraspecific records before building species-level spaces.")
  }

  ti <- check_trait_info(trait_info, data)
  out <- data.frame(row.names = sp)
  warn <- character(0)

  for (i in seq_len(nrow(ti))) {
    tr <- ti$trait[i]
    x <- data[[tr]]
    switch(ti$type[i],
      Q = {
        x <- suppressWarnings(as.numeric(as.character(x)))
        res <- .transform(x, ti$transform[i])
        if (res$n_lost > 0) {
          warn <- c(warn, sprintf("%s: %d value(s) incompatible with '%s' set to NA",
                                  tr, res$n_lost, ti$transform[i]))
        }
        x <- res$x
      },
      O = {
        lev <- ti$levels[i]
        if (!is.na(lev) && nzchar(lev)) {
          lev <- trimws(strsplit(lev, ";", fixed = TRUE)[[1]])
          bad <- setdiff(unique(stats::na.omit(as.character(x))), lev)
          if (length(bad)) stop("Ordinal trait '", tr, "' has values absent from `levels`: ",
                                paste(bad, collapse = ", "))
          x <- as.numeric(factor(as.character(x), levels = lev))
        } else if (is.factor(x)) {
          x <- as.numeric(x)
        } else {
          x <- suppressWarnings(as.numeric(as.character(x)))
          if (all(is.na(x)))
            stop("Ordinal trait '", tr, "' is not numeric: give its ordered `levels` (e.g. 'low;medium;high').")
        }
      },
      N = {
        x <- factor(as.character(x))
      },
      D = , B = {
        x <- .as01(x, tr)
      },
      F = {
        x <- suppressWarnings(as.numeric(as.character(x)))
        if (any(x < 0, na.rm = TRUE)) stop("Fuzzy-coded trait '", tr, "' has negative values.")
      },
      C = {
        x <- suppressWarnings(as.numeric(as.character(x)))
      }
    )
    out[[tr]] <- x
  }

  # Fuzzy blocks -> row proportions (as ade4::prep.fuzzy). Rows summing to 0 -> NA.
  for (b in unique(ti$block[ti$type == "F"])) {
    cols <- ti$trait[ti$block == b]
    m <- as.matrix(out[, cols, drop = FALSE])
    rs <- rowSums(m, na.rm = TRUE)
    allna <- rowSums(!is.na(m)) == 0
    zero <- rs == 0 & !allna
    if (any(zero)) warn <- c(warn, sprintf("fuzzy block %s: %d species with all-zero scores set to NA",
                                            b, sum(zero)))
    m[is.na(m) & !allna] <- 0
    m <- m / rs
    m[zero | allna, ] <- NA
    out[, cols] <- m
  }

  # Multichoice blocks: a species with no modality recorded has an undefined
  # Jaccard distance (NaN in ade4::dist.ktab); recode it as missing so that
  # it is handled as such (and excluded from spaces where it has no data).
  for (b in unique(ti$block[ti$type == "B"])) {
    cols <- ti$trait[ti$block == b]
    m <- as.matrix(out[, cols, drop = FALSE])
    zero <- rowSums(m, na.rm = TRUE) == 0 & rowSums(!is.na(m)) > 0
    if (any(zero)) {
      warn <- c(warn, sprintf("multichoice block %s: %d species with no modality recorded set to NA",
                              b, sum(zero)))
      out[zero, cols] <- NA
    }
  }

  if (length(warn)) warning(paste(warn, collapse = "\n"), call. = FALSE)

  structure(list(
    data = out,
    trait_info = ti,
    species = sp,
    imputed = FALSE,
    missing = colMeans(is.na(out))
  ), class = "fs_traits")
}

.transform <- function(x, how) {
  n0 <- sum(!is.na(x))
  y <- switch(how,
    none = x,
    log10 = { x[x <= 0] <- NA; log10(x) },
    log1p = { x[x <= -1] <- NA; log1p(x) },
    sqrt = { x[x < 0] <- NA; sqrt(x) },
    logit = { x[x <= 0 | x >= 1] <- NA; stats::qlogis(x) }
  )
  list(x = y, n_lost = n0 - sum(!is.na(y)))
}

.as01 <- function(x, tr) {
  if (is.logical(x)) return(as.numeric(x))
  if (is.character(x) || is.factor(x)) {
    v <- tolower(trimws(as.character(x)))
    map <- c("1" = 1, "0" = 0, "yes" = 1, "no" = 0, "true" = 1, "false" = 0,
             "oui" = 1, "non" = 0, "y" = 1, "n" = 0, "present" = 1, "absent" = 0)
    y <- unname(map[v])
    if (any(is.na(y) & !is.na(v))) stop("Binary trait '", tr, "' has values other than 0/1.")
    return(y)
  }
  x <- as.numeric(x)
  if (!all(stats::na.omit(x) %in% c(0, 1))) stop("Binary trait '", tr, "' has values other than 0/1.")
  x
}

#' @export
print.fs_traits <- function(x, ...) {
  ti <- x$trait_info
  cat("<fs_traits>", length(x$species), "species,", nrow(ti), "trait columns",
      if (x$imputed) "(imputed)" else "", "\n")
  tab <- trait_summary(x)
  print(tab, row.names = FALSE)
  invisible(x)
}

#' Summarise traits by dimension and type
#'
#' @param x An `fs_traits` object.
#' @return Data frame with, for each trait (multi-column blocks collapsed),
#'   its dimension, type, number of columns, weight and percentage of
#'   species with missing values.
#' @export
trait_summary <- function(x) {
  stopifnot(inherits(x, "fs_traits"))
  ti <- x$trait_info
  blocks <- unique(ti$block)
  do.call(rbind, lapply(blocks, function(b) {
    r <- ti[ti$block == b, ]
    na <- rowSums(!is.na(x$data[, r$trait, drop = FALSE])) == 0
    data.frame(trait = b, dimension = r$dimension[1], type = r$type[1],
               n_columns = nrow(r), transform = r$transform[1], weight = r$weight[1],
               pct_missing = round(100 * mean(na), 1), stringsAsFactors = FALSE)
  }))
}

#' Remove species with missing trait values
#'
#' @param x An `fs_traits` object.
#' @param dimensions Optional character vector; only traits of these
#'   dimensions are considered.
#' @return The `fs_traits` object restricted to complete species.
#' @export
drop_incomplete <- function(x, dimensions = NULL) {
  stopifnot(inherits(x, "fs_traits"))
  ti <- x$trait_info
  if (!is.null(dimensions)) ti <- ti[ti$dimension %in% dimensions, ]
  keep <- stats::complete.cases(x$data[, ti$trait, drop = FALSE])
  subset_species(x, x$species[keep])
}

#' Subset species of an `fs_traits` object
#' @param x An `fs_traits` object.
#' @param species Character vector of species to keep.
#' @return The subset object.
#' @export
subset_species <- function(x, species) {
  stopifnot(inherits(x, "fs_traits"))
  keep <- x$species %in% species
  x$data <- x$data[keep, , drop = FALSE]
  x$species <- x$species[keep]
  x$missing <- colMeans(is.na(x$data))
  x
}
