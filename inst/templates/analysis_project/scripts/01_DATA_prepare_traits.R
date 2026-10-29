# ============================================================================
# 01_DATA_prepare_traits.R -- reads the trait table and the description of the
# traits (type, functional dimension, transformation, weight), checks them and
# builds the trait object used by every analysis.
#   Inputs : data/raw/traits.csv, data/raw/trait_info.csv
#   Outputs: data/processed/traits_prepared.rds
#            results/tables/Table1_traits.csv
# ============================================================================
source("scripts/00_START_GeneralScript.R")

traits_raw <- read_traits_table()
trait_info <- utils::read.csv(PARAMS$trait_info_file, stringsAsFactors = FALSE)

traits <- prepare_traits(traits_raw, trait_info, species_col = PARAMS$species_col)
print(traits)

save_object(traits, "traits_prepared")
save_table(trait_summary(traits), "Table1_traits")
