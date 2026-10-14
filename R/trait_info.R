#' Trait types supported by DisentFunSpaces
#'
#' Returns the trait types that can be declared in a trait-information table.
#' Codes are those of [ade4::dist.ktab()], which computes the generalised
#' Gower distance of Pavoine et al. (2009) used in Toussaint et al. (2026).
#'
#' @return A data frame with one row per type: `code`, `label`,
#'   `multi_column` (does one trait span several columns?), `input`
#'   (expected data) and `distance` (per-trait distance used by
#'   `ade4::dist.ktab()` with its non-interactive defaults).
#' @export
#' @examples
#' trait_types()
trait_types <- function() {
  data.frame(
    code = c("Q", "O", "N", "D", "B", "F", "C"),
    label = c(
      "Continuous (quantitative)",
      "Ordinal",
      "Categorical (nominal)",
      "Binary (single 0/1 trait)",
      "Multichoice (block of 0/1 columns)",
      "Fuzzy-coded (block of affinity / proportion columns)",
      "Circular (hour, month, aspect...)"
    ),
    multi_column = c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE),
    input = c(
      "numeric",
      "numeric ranks, ordered factor, or text with `levels` given",
      "text or factor",
      "0/1 (or TRUE/FALSE)",
      "several 0/1 columns sharing the same `block`",
      "several non-negative columns sharing the same `block`; rows are rescaled to proportions",
      "numeric on a cyclic scale; `circ_min`/`circ_max` = first and last possible values (ade4::prep.circular), e.g. 1/12 for months, 0/23 for whole hours"
    ),
    distance = c(
      "|xi - xj| / range (option 'scaledBYrange')",
      "rank-based distance for ordinal variables (ade4 method 1)",
      "simple matching: 0 if identical, 1 otherwise",
      "Jaccard, d = sqrt(1 - s)",
      "Jaccard on the block, d = sqrt(1 - s)",
      "Manly overlap index on proportions",
      "angular distance on the circle"
    ),
    stringsAsFactors = FALSE
  )
}

.ti_columns <- c("trait", "type", "dimension", "block", "transform",
                 "weight", "levels", "circ_min", "circ_max", "use")

.transforms <- c("none", "log10", "log1p", "sqrt", "logit")

#' Guess a trait-information table from a trait data frame
#'
#' Produces a first draft of the trait-information table that drives the
#' whole workflow. Every column except the species identifier becomes a
#' trait; its type is guessed from its content. **The draft must be
#' checked**: fuzzy-coded and multichoice blocks, ordinal traits and
#' functional dimensions cannot be inferred reliably and are usually set by
#' hand (or in the Shiny application).
#'
#' Heuristics:
#' * logical, or numeric with only 0/1 values -> `D`
#' * character / factor -> `N` (ordered factor -> `O`)
#' * numeric with <= `max_ordinal_levels` distinct integer values -> `O`
#' * other numeric -> `Q`; strictly positive and strongly right-skewed
#'   (skewness > 1) continuous traits get `transform = "log10"`, as in the
#'   original bird analysis.
#' * columns sharing a prefix before a separator (e.g. `Diet.Inv`,
#'   `Diet.Fruit`) whose rows sum to a constant are proposed as one fuzzy
#'   block (`F`); groups of 0/1 columns with a common prefix are proposed as a
#'   multichoice block (`B`).
#'
#' @param data A data frame of traits (species in rows).
#' @param species_col Name of the species identifier column, or `NULL` if
#'   species are in the row names.
#' @param max_ordinal_levels Integer; numeric traits with at most this many
#'   distinct integer values are proposed as ordinal.
#' @param detect_blocks Logical; try to detect fuzzy and multichoice blocks
#'   from column-name prefixes.
#' @param sep Regular expression separating a block prefix from the modality
#'   name in column names.
#' @param predictor_pattern Regular expression (case-insensitive) matching
#'   columns that are imputation predictors rather than traits (by default
#'   phylogenetic eigenvectors such as `Eigen.1`); they are listed with
#'   `use = FALSE`. See the `covariates` argument of [impute_traits()].
#'
#' @return A data frame with columns `trait`, `type`, `dimension`, `block`,
#'   `transform`, `weight`, `levels`, `circ_min`, `circ_max`, `use`.
#' @export
#' @examples
#' data(birds_traits)
#' # taxonomy and IUCN category are not traits
#' x <- birds_traits[, !names(birds_traits) %in% c("order", "family", "iucn")]
#' guess_trait_info(x, species_col = "species")
guess_trait_info <- function(data, species_col = NULL, max_ordinal_levels = 5,
                             detect_blocks = TRUE, sep = "[._]",
                             predictor_pattern = "eigen") {
  data <- as.data.frame(data)
  traits <- setdiff(names(data), species_col)
  if (length(traits) == 0) stop("No trait column found.")

  type <- vapply(traits, function(tr) {
    x <- data[[tr]]
    if (is.logical(x)) return("D")
    if (is.ordered(x)) return("O")
    if (is.factor(x) || is.character(x)) return("N")
    if (is.numeric(x)) {
      ux <- unique(stats::na.omit(x))
      if (length(ux) <= 2 && all(ux %in% c(0, 1))) return("D")
      if (length(ux) <= max_ordinal_levels && all(ux == round(ux))) return("O")
      return("Q")
    }
    "N"
  }, character(1))

  transform <- vapply(traits, function(tr) {
    x <- data[[tr]]
    if (type[tr] != "Q") return("none")
    x <- stats::na.omit(x)
    if (length(x) > 2 && all(x > 0) && .skewness(x) > 1) "log10" else "none"
  }, character(1))

  block <- traits
  if (detect_blocks) {
    prefix <- vapply(strsplit(traits, sep), `[`, character(1), 1)
    for (p in unique(prefix[duplicated(prefix)])) {
      cols <- traits[prefix == p]
      if (length(cols) < 2) next
      sub <- data[, cols, drop = FALSE]
      if (!all(vapply(sub, is.numeric, logical(1)))) next
      vals <- unique(stats::na.omit(unlist(sub)))
      if (all(vals %in% c(0, 1))) {
        type[cols] <- "B"; block[traits %in% cols] <- p; transform[cols] <- "none"
      } else if (all(vals >= 0)) {
        # fuzzy-coded only if rows (complete ones) sum to a near-constant
        rs <- rowSums(sub[stats::complete.cases(sub), , drop = FALSE])
        if (length(rs) > 1 && mean(rs) > 0 && stats::sd(rs) / mean(rs) < 0.05) {
          type[cols] <- "F"; block[traits %in% cols] <- p; transform[cols] <- "none"
        }
      }
    }
  }

  data.frame(
    trait = traits,
    type = unname(type),
    dimension = "All traits",
    block = block,
    transform = unname(transform),
    weight = 1,
    levels = NA_character_,
    circ_min = NA_real_,
    circ_max = NA_real_,
    use = !grepl(predictor_pattern, traits, ignore.case = TRUE),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

.skewness <- function(x) {
  x <- x[is.finite(x)]
  s <- stats::sd(x)
  if (length(x) < 3 || is.na(s) || s == 0) return(0)
  mean((x - mean(x))^3) / s^3
}

#' Check and complete a trait-information table
#'
#' Validates a trait-information table against a trait data frame and fills
#' optional columns with defaults. Called internally by every function that
#' takes a `trait_info` argument; use it directly to check a table read from
#' a file.
#'
#' Required columns: `trait` (column name in `data`), `type` (one of
#' [trait_types()] codes) and `dimension` (functional dimension, e.g.
#' "Locomotion"). Optional columns: `block` (traits of type `F` or `B` that
#' share a `block` form a single multi-column trait; defaults to `trait`),
#' `transform` (`none`, `log10`, `log1p`, `sqrt`, `logit`; continuous traits
#' only), `weight` (relative weight of the trait in the distance, default 1),
#' `levels` (ordered levels of an ordinal trait, separated by `;`),
#' `circ_min`/`circ_max` (cycle limits of circular traits) and `use`
#' (logical; `FALSE` drops the trait).
#'
#' @param trait_info Data frame describing the traits.
#' @param data Optional trait data frame used to check that traits exist.
#'
#' @return The completed table, restricted to traits with `use == TRUE`.
#' @export
check_trait_info <- function(trait_info, data = NULL) {
  ti <- as.data.frame(trait_info, stringsAsFactors = FALSE)
  miss <- setdiff(c("trait", "type"), names(ti))
  if (length(miss)) stop("trait_info lacks column(s): ", paste(miss, collapse = ", "))
  ti$trait <- as.character(ti$trait)
  ti$type <- toupper(trimws(as.character(ti$type)))
  if (is.null(ti$dimension)) ti$dimension <- "All traits"
  ti$dimension <- as.character(ti$dimension)
  ti$dimension[is.na(ti$dimension) | ti$dimension == ""] <- "Unassigned"
  if (is.null(ti$block)) ti$block <- ti$trait
  ti$block <- as.character(ti$block)
  ti$block[is.na(ti$block) | ti$block == ""] <- ti$trait[is.na(ti$block) | ti$block == ""]
  if (is.null(ti$transform)) ti$transform <- "none"
  ti$transform <- as.character(ti$transform)
  ti$transform[is.na(ti$transform) | ti$transform == ""] <- "none"
  if (is.null(ti$weight)) ti$weight <- 1
  ti$weight <- as.numeric(ti$weight)
  ti$weight[is.na(ti$weight)] <- 1
  if (is.null(ti$levels)) ti$levels <- NA_character_
  ti$levels <- as.character(ti$levels)
  if (is.null(ti$circ_min)) ti$circ_min <- NA_real_
  if (is.null(ti$circ_max)) ti$circ_max <- NA_real_
  ti$circ_min <- as.numeric(ti$circ_min)
  ti$circ_max <- as.numeric(ti$circ_max)
  if (is.null(ti$use)) ti$use <- TRUE
  ti$use <- as.logical(ti$use)
  ti$use[is.na(ti$use)] <- TRUE

  ti <- ti[ti$use, .ti_columns, drop = FALSE]
  if (nrow(ti) == 0) stop("No trait selected (all `use` are FALSE).")
  rownames(ti) <- NULL

  bad <- setdiff(unique(ti$type), trait_types()$code)
  if (length(bad)) {
    stop("Unknown trait type(s): ", paste(bad, collapse = ", "),
         ". Use one of: ", paste(trait_types()$code, collapse = ", "), ".")
  }
  badtr <- setdiff(unique(ti$transform), .transforms)
  if (length(badtr)) stop("Unknown transform(s): ", paste(badtr, collapse = ", "))
  if (any(ti$transform != "none" & ti$type != "Q")) {
    stop("Transformations only apply to continuous (Q) traits: ",
         paste(ti$trait[ti$transform != "none" & ti$type != "Q"], collapse = ", "))
  }
  if (any(ti$weight < 0)) stop("Weights must be non-negative.")
  if (anyDuplicated(ti$trait)) stop("Duplicated trait name(s) in trait_info.")

  # single-column types: block must be the trait itself
  single <- !ti$type %in% c("F", "B")
  ti$block[single] <- ti$trait[single]

  # multi-column blocks must be homogeneous
  for (b in unique(ti$block[!single])) {
    rows <- ti$block == b
    if (length(unique(ti$type[rows])) > 1)
      stop("Block '", b, "' mixes several trait types.")
    if (length(unique(ti$dimension[rows])) > 1)
      stop("Block '", b, "' is assigned to several dimensions.")
    if (length(unique(ti$weight[rows])) > 1)
      stop("Block '", b, "' has several weights; give one weight per block.")
    if (sum(rows) < 2)
      warning("Block '", b, "' (type ", ti$type[rows][1],
              ") has a single column; consider type D or Q instead.")
  }
  if (any(ti$block[single] %in% ti$block[!single]))
    stop("A block name is also used as a single-column trait name.")

  if (!is.null(data)) {
    absent <- setdiff(ti$trait, names(data))
    if (length(absent)) stop("Trait(s) absent from data: ", paste(absent, collapse = ", "))
  }
  ti
}

#' Write an empty trait-information template
#'
#' Writes a CSV template (one row per trait column of `data`, pre-filled by
#' [guess_trait_info()]) that can be completed in a spreadsheet and read
#' back with [utils::read.csv()].
#'
#' @inheritParams guess_trait_info
#' @param file Output CSV path.
#' @return Invisibly, the template data frame.
#' @export
write_trait_template <- function(data, file, species_col = NULL) {
  ti <- guess_trait_info(data, species_col = species_col)
  utils::write.csv(ti, file, row.names = FALSE)
  invisible(ti)
}
