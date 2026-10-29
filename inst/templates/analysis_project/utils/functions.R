# ============================================================================
# Helpers shared by the scripts (reading the data, writing tables and figures)
# ============================================================================

# The trait table as archived in data/raw (one row per species).
read_traits_table <- function(file = PARAMS$data_file) {
  utils::read.csv(file, check.names = FALSE, stringsAsFactors = FALSE,
                  na.strings = c("NA", ""), encoding = "UTF-8")
}

# A named vector (names = species) from a column of the trait table, or NULL.
species_vector <- function(data, column) {
  if (is.null(column) || !nzchar(column) || !column %in% names(data)) return(NULL)
  stats::setNames(as.character(data[[column]]), as.character(data[[PARAMS$species_col]]))
}

# Objects passed from one script to the next.
save_object <- function(x, name) saveRDS(x, file.path("data", "processed", paste0(name, ".rds")))
read_object <- function(name) {
  f <- file.path("data", "processed", paste0(name, ".rds"))
  if (!file.exists(f)) stop(f, " not found: run the previous scripts first (run_all.R).", call. = FALSE)
  readRDS(f)
}
has_object <- function(name) file.exists(file.path("data", "processed", paste0(name, ".rds")))

# Result tables: results/tables/<name>.csv
save_table <- function(x, name) {
  f <- file.path("results", "tables", paste0(name, ".csv"))
  utils::write.csv(x, f, row.names = FALSE, na = "")
  invisible(f)
}

# Figures: results/figures/<name>.png (300 dpi) and .pdf
save_figure <- function(p, name, width = 7, height = 6) {
  for (ext in c("png", "pdf")) {
    f <- file.path("results", "figures", paste0(name, ".", ext))
    if (ext == "png") ggplot2::ggsave(f, p, width = width, height = height, dpi = 300, bg = "white")
    else ggplot2::ggsave(f, p, width = width, height = height)
  }
  invisible(name)
}

# File-name-safe version of a space name.
safe_name <- function(x) gsub("^_|_$", "", gsub("[^A-Za-z0-9]+", "_", x))

# rbind of data frames whose columns differ (missing columns filled with NA).
bind_fill <- function(lst) {
  lst <- Filter(Negate(is.null), lst)
  cols <- unique(unlist(lapply(lst, names)))
  do.call(rbind, lapply(lst, function(d) { d[setdiff(cols, names(d))] <- NA; d[cols] }))
}
