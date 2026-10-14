# Example data (real birds): complete cases (281 species) or, with missing
# values, 333 species with clutch size, egg mass and adult body mass known, so
# that every pair of species shares reproduction traits (keeps tests fast).
ex <- function(complete = TRUE) {
  tr <- prepare_traits(birds_traits, birds_trait_info, species_col = "species")
  if (complete) return(drop_incomplete(tr))
  ok <- birds_traits$species[!is.na(birds_traits$clutch_size) & !is.na(birds_traits$egg_mass) &
                               !is.na(birds_traits$adult_body_mass)]
  subset_species(tr, ok[seq_len(333)])
}

# Simulated data with all seven trait types (Q, N, F, O, C, B, D); the bird
# data have no circular or multichoice trait.
sim_tr <- function(n = 120, missing = 0, seed = 7) {
  s <- simulate_traits(n, seed = seed, missing = missing)
  suppressWarnings(prepare_traits(s$traits, s$trait_info, "species"))
}

# Reference generalised Gower distance: one ade4::dist.ktab() call on all
# traits (block order does not matter with equal weights).
ref_distktab <- function(tr) {
  ti <- tr$trait_info; x <- tr$data
  blocks <- list(); groups <- character(0)
  add <- function(df, g) { blocks[[length(blocks) + 1]] <<- df; groups <<- c(groups, g) }
  for (t in c("Q", "N", "O", "D")) {
    tt <- ti$trait[ti$type == t]
    if (length(tt)) add(x[, tt, drop = FALSE], t)
  }
  for (b in unique(ti$block[ti$type %in% c("F", "B")])) {
    cols <- ti$trait[ti$block == b]
    if (ti$type[ti$block == b][1] == "F") add(ade4::prep.fuzzy(x[, cols], length(cols), labels = b), "F")
    else add(ade4::prep.binary(x[, cols], length(cols), labels = b), "B")
  }
  cc <- ti[ti$type == "C", ]
  if (nrow(cc)) add(ade4::prep.circular(x[, cc$trait, drop = FALSE], cc$circ_min, cc$circ_max), "C")
  ade4::dist.ktab(ade4::ktab.list.df(blocks), groups)
}
