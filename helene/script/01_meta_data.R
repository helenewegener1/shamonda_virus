library(tidyverse)

setwd("~/Documents/courses/MolecularEvolutionPhylogenomics/")

# ------------------------------------------------------------------------------
# Load data 
# ------------------------------------------------------------------------------

meta1 <- read_csv("data/sequences_shamonda.csv") # Shamonda virus (taxid:15950)
meta2 <- read_csv("data/sequences_Schmallenbergvirus1.csv") # Orthonbunyavirus schmallenvergense (taxid: 3052437)
meta3 <- read_csv("data/sequences_Schmallenbergvirus2.csv") # Schmallenberg virus (taxid: 1133363)

# ------------------------------------------------------------------------------
# Check presence of accissions number 
# ------------------------------------------------------------------------------

grep("OZ50247*", meta1$Accession, value = TRUE)

# ------------------------------------------------------------------------------
# Add taxid and search string 
# ------------------------------------------------------------------------------

meta1$taxid <- "15950"
meta2$taxid <- "3052437"
meta3$taxid <- "1133363"

meta1$search_string <- "Shamonda virus (taxid:15950)"
meta2$search_string <- "Orthonbunyavirus schmallenvergense (taxid: 3052437)"
meta3$search_string <- "Schmallenberg virus (taxid: 1133363)"

# ------------------------------------------------------------------------------
# Combine
# ------------------------------------------------------------------------------

meta_all <- bind_rows(meta1, meta2, meta3)

meta_all <- meta_all %>% select(
  Accession, Organism_Name, Isolate, Segment, Species, Length, Nuc_Completeness, 
  Geo_Location, Host,
  Collection_Date, taxid, search_string
)

# ------------------------------------------------------------------------------
# Overview of features + clean
# ------------------------------------------------------------------------------

# meta_all_clean <- 
meta_all$Organism_Name %>% table(useNA = "always")
meta_all$Species %>% table(useNA = "always")
meta_all$Segment %>% table(useNA = "always")
meta_all$Nuc_Completeness %>% table(useNA = "always")
meta_all$Geo_Location %>% table(useNA = "always") # add country 
meta_all$Host %>% table(useNA = "always")
meta_all$Collection_Date %>% table(useNA = "always")

# Clean 
meta_all_clean <- meta_all %>% 
  mutate(
    Country = sub(":.*", "", Geo_Location)
  )

meta_all_clean$Country %>% table()

saveRDS(meta_all_clean, "helene/out/meta_all_clean.rds")



