tax <- c("order", "family", "iucn")

test_that("guess_trait_info detects types and blocks", {
  g <- guess_trait_info(birds_traits[, setdiff(names(birds_traits), tax)], species_col = "species")
  expect_equal(g$type[g$trait == "body_mass"], "Q")
  expect_equal(g$transform[g$trait == "body_mass"], "log10")
  expect_true(all(g$type[grepl("^diet_", g$trait)] == "F"))
  expect_true(all(g$block[grepl("^diet_", g$trait)] == "diet"))
  expect_true(all(g$type[grepl("^beak_", g$trait)] == "Q"))
  # other types (none in the bird data)
  df <- data.frame(sp = paste0("s", 1:6), hab_forest = c(1, 0, 1, 1, 0, 1),
                   hab_open = c(0, 1, 1, 0, 1, 0), x = c(2.1, 3.5, 1.2, 8.4, 5.5, 3.3),
                   mode = c("a", "b", "a", "c", "b", "a"), mig = c(0, 1, 1, 0, 0, 1),
                   dev = c(1, 2, 3, 2, 1, 3))
  gb <- guess_trait_info(df, species_col = "sp")
  expect_equal(gb$type, c("B", "B", "Q", "N", "D", "O"))
  expect_equal(gb$block[1:2], c("hab", "hab"))
})

test_that("check_trait_info catches inconsistent tables", {
  ti <- birds_trait_info
  ti$type[1] <- "X"
  expect_error(check_trait_info(ti), "Unknown trait type")
  ti <- birds_trait_info; ti$transform[ti$trait == "diet_inv"] <- "log10"
  expect_error(check_trait_info(ti), "only apply to continuous")
  ti <- birds_trait_info; ti$dimension[ti$trait == "diet_inv"] <- "Other"
  expect_error(check_trait_info(ti), "several dimensions")
  expect_error(check_trait_info(birds_trait_info, birds_traits[, 1:3]), "absent from data")
  expect_false(any(grepl("^beak_", check_trait_info(birds_trait_info)$trait)))   # use = FALSE
})

test_that("prepare_traits converts types correctly", {
  tr <- prepare_traits(birds_traits, birds_trait_info, species_col = "species")
  expect_s3_class(tr, "fs_traits")
  expect_equal(nrow(tr$data), 2000)
  d <- tr$data
  expect_false(any(grepl("^beak_", names(d))))
  expect_setequal(unique(tr$trait_info$dimension), c("Locomotion", "Reproduction", "Diet"))
  fz <- d[, grep("^diet_", names(d))]
  expect_equal(unname(rowSums(fz)), rep(1, 2000), tolerance = 1e-12)
  expect_equal(d$body_mass, log10(birds_traits$body_mass))
  expect_equal(mean(is.na(d$longevity_years)), mean(is.na(birds_traits$longevity_years)))
  expect_error(prepare_traits(rbind(birds_traits, birds_traits[1, ]), birds_trait_info, "species"),
               "Duplicated")
  # ordinal text traits are converted to ranks with their levels; nominal to factors
  df <- data.frame(sp = paste0("s", 1:4), dev = c("low", "high", "mid", NA), m = c("a", "b", "a", "b"))
  ti <- data.frame(trait = c("dev", "m"), type = c("O", "N"), dimension = "A",
                   levels = c("low;mid;high", NA))
  t2 <- prepare_traits(df, ti, "sp")
  expect_equal(t2$data$dev, c(1, 3, 2, NA))
  expect_true(is.factor(t2$data$m))
})

test_that("ordinal text without levels fails informatively", {
  df <- data.frame(sp = paste0("s", 1:4), dev = c("low", "high", "mid", "low"))
  ti <- data.frame(trait = "dev", type = "O", dimension = "A")
  expect_error(prepare_traits(df, ti, "sp"), "levels")
})

test_that("eigenvector columns are proposed as predictors, not traits", {
  g <- guess_trait_info(birds_traits[, setdiff(names(birds_traits), tax)], species_col = "species")
  expect_false(any(g$use[grepl("^Eigen", g$trait)]))
  expect_true(all(g$use[g$trait == "body_mass"]))
})

test_that("example data are consistent with their phylogeny", {
  expect_setequal(birds_tree$tip.label, birds_traits$species)
  expect_true(ape::is.ultrametric(birds_tree))
  expect_setequal(unique(birds_traits$iucn), c("LC", "NT", "VU", "EN", "CR"))
})

test_that("simulate_traits returns valid data with the target congruence", {
  sim <- simulate_traits(120, seed = 3, missing = 0.05)
  expect_s3_class(sim, "fs_simulation")
  expect_equal(unname(cor(sim$latent)), unname(sim$cor_dims), tolerance = 1e-8)
  tr <- suppressWarnings(prepare_traits(sim$traits, sim$trait_info, "species"))
  expect_equal(length(tr$species), 120)
  expect_setequal(unique(tr$trait_info$type), c("Q", "N", "F", "O", "C", "B", "D"))
  expect_error(simulate_traits(50, cor_dims = 1.5), "valid correlation")
})

test_that("the full dataset has the structure of the example", {
  expect_equal(nrow(birds_traits_all), 8413)
  expect_identical(names(birds_traits_all), names(birds_traits))
  expect_true(all(birds_traits$species %in% birds_traits_all$species))
  sub <- birds_traits_all[match(birds_traits$species, birds_traits_all$species), ]
  rownames(sub) <- NULL
  expect_equal(sub, birds_traits)
  expect_silent(check_trait_info(birds_trait_info, birds_traits_all))
})
