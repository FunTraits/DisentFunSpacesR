test_that("variable weighting reproduces a single ade4::dist.ktab call (published analysis)", {
  tr <- ex()                                      # Q, F (2 blocks), N, O, D
  expect_equal(as.vector(trait_distance(tr, engine = "ade4")), as.vector(ref_distktab(tr)),
               tolerance = 1e-12)
  ts <- sim_tr()                                  # all seven types, incl. B and C
  expect_equal(as.vector(trait_distance(ts, engine = "ade4")), as.vector(ref_distktab(ts)),
               tolerance = 1e-12)
})

test_that("explicit equal weights give the same result as a single dist.ktab call", {
  tr <- ex()
  tr2 <- tr; tr2$trait_info$weight[1] <- 1 + 1e-12   # forces the slow path
  expect_equal(as.vector(trait_distance(tr)), as.vector(trait_distance(tr2)), tolerance = 1e-8)
})

test_that("dimension weighting gives equal weight to dimensions", {
  tr <- ex()
  dims <- unique(tr$trait_info$dimension)
  dd <- lapply(dims, function(dm) as.vector(trait_distance(tr, dimensions = dm)))
  manual <- sqrt(Reduce(`+`, lapply(dd, `^`, 2)) / length(dd))
  expect_equal(as.vector(trait_distance(tr, weighting = "dimension")), manual, tolerance = 1e-12)
})

test_that("missing values are handled pairwise", {
  tr <- ex(FALSE)
  expect_true(anyNA(tr$data$longevity_years))
  d <- trait_distance(tr)
  expect_false(anyNA(d))
  expect_equal(attr(d, "Size"), 333)
})

test_that("fast and ade4 engines agree for every type, with missing values", {
  # simulated data (all types): missing values in a continuous, a circular and
  # a multichoice trait, each dimension keeping at least one complete trait
  ts <- sim_tr(); ti <- ts$trait_info
  set.seed(4); idx <- sample(length(ts$species), 15)
  na_cols <- c(ti$trait[ti$type == "Q"][1], ti$trait[ti$type == "C"],
               ti$trait[ti$block == ti$block[ti$type == "B"][1]])
  ts$data[idx, na_cols] <- NA
  for (tr in list(ex(FALSE), ts)) for (dm in unique(tr$trait_info$dimension)) {
    sub <- subset_species(tr, tr$species[rowSums(!is.na(tr$data[, tr$trait_info$trait[tr$trait_info$dimension == dm], drop = FALSE])) > 0])
    expect_equal(as.vector(trait_distance(sub, dm, engine = "fast")),
                 as.vector(trait_distance(sub, dm, engine = "ade4")), tolerance = 1e-10,
                 label = paste("dimension", dm))
  }
  tr <- ex(FALSE)
  for (w in c("variable", "dimension")) {
    expect_equal(as.vector(trait_distance(tr, weighting = w, engine = "fast")),
                 as.vector(trait_distance(tr, weighting = w, engine = "ade4")), tolerance = 1e-10)
  }
  # extra checks on edge cases: multichoice with all-zero rows, odd circular range
  set.seed(5); n <- 40
  df <- data.frame(sp = paste0("s", 1:n), h1 = sample(0:1, n, TRUE), h2 = sample(0:1, n, TRUE),
                   m = sample(1:12, n, TRUE), o = sample(c("a", "b", "c"), n, TRUE),
                   q = c(NA, rnorm(n - 1)))
  df$h1[1:3] <- 0; df$h2[1:3] <- 0
  ti <- data.frame(trait = c("h1", "h2", "m", "o", "q"), type = c("B", "B", "C", "O", "Q"),
                   dimension = "A", block = c("h", "h", "m", "o", "q"),
                   levels = c(NA, NA, NA, "a;b;c", NA), circ_min = c(NA, NA, 1, NA, NA),
                   circ_max = c(NA, NA, 12, NA, NA))
  t2 <- prepare_traits(df, ti, "sp")
  expect_equal(as.vector(trait_distance(t2, engine = "fast")),
               as.vector(trait_distance(t2, engine = "ade4")), tolerance = 1e-10)
})

test_that("dimension_scaled gives each dimension a mean squared distance of 1", {
  tr <- ex()
  dims <- c("Diet", "Reproduction")
  d <- as.vector(trait_distance(tr, dims, weighting = "dimension_scaled"))
  dd <- lapply(dims, function(dm) as.vector(trait_distance(tr, dm)))
  dd <- lapply(dd, function(v) v / sqrt(mean(v^2)))
  expect_equal(d, sqrt((dd[[1]]^2 + dd[[2]]^2) / 2), tolerance = 1e-10)
  expect_equal(d, as.vector(trait_distance(tr, dims, weighting = "dimension_scaled", engine = "fast")),
               tolerance = 1e-10)
})
