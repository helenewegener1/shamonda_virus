library(tidyverse)
library(glue)

setwd("~/Documents/courses/MolecularEvolutionPhylogenomics/shamonda_virus/")

# Amanda is aligning using MAFTT 

# ------------------------------------------------------------------------------
# Load data 
# ------------------------------------------------------------------------------

meta_all_clean <- readRDS("helene/out/meta_all_clean.rds")

fasta <- readDNAStringSet("helene/out/M.fasta")


fasta1 <- readDNAStringSet("data/sequences_shamonda.fasta")             # Shamonda, taxid 15950
fasta2 <- readDNAStringSet("data/sequences_Schmallenbergvirus1.fasta")  # Orthobunyavirus schmallenbergense, taxid 3052437
fasta3 <- readDNAStringSet("data/sequences_Schmallenbergvirus2.fasta")  # Schmallenberg virus, taxid 1133363

seqs_all <- c(fasta1, fasta2, fasta3)


# ------------------------------------------------------------------------------
# Subset to the 55 
# ------------------------------------------------------------------------------

# meta_all_clean %>% nrow()
# 
# meta_all_clean %>% 
#   select(Segment, Isolate, Accession) %>% 
#   distinct() %>%
#   # filter(
#   #   Segment == "M"
#   # ) %>% 
#   count(Isolate, sort = TRUE) %>% 
#   filter(n > 2) 
# 
# meta_all_clean %>% filter(Isolate == "F6") %>% distinct()
# meta_all_clean %>% filter(Isolate == "BE/12/2478") %>% distinct() %>% view()
# # meta_all_clean %>% filter(Isolate == "Ib An 5550") %>% distinct() %>% view()
# 
# 
# seqs_all$`NC_043582.1 |Bos taurus|Netherlands|2011` == seqs_all$`KC355459.1 |Bos taurus|Netherlands|2011` 

# ------------------------------------------------------------------------------
# 
# ------------------------------------------------------------------------------

SML <- read_csv("complete_isolates_SML.csv")

SML$Organism_Name %>% table()

seqs_all
names(seqs_all) <- str_extract(names(seqs_all), "^[^.\\s]+")

# drop duplicates in case the three taxid queries overlap
seqs_all <- seqs_all[!duplicated(names(seqs_all))]

# Filter for accession number in SML 
table(SML$Accession %in% names(seqs_all))

seq_SML <- seqs_all[names(seqs_all) %in% SML$Accession]

# Divide accession numbers in segments 
segment <- list()
segment[["L"]] <- SML %>% filter(Segment == "L") %>% select(Accession) %>% pull()
segment[["S"]] <- SML %>% filter(Segment == "S") %>% select(Accession) %>% pull()
segment[["M"]] <- SML %>% filter(Segment == "M") %>% select(Accession) %>% pull()

# Divide seqs by segments 
seqs_by_seg <- map(segment, ~ seq_SML[names(seq_SML) %in% .x])

# Prep outgroup 
out_seqs <- readDNAStringSet("data/California.fasta")
out_meta <- read_csv("data/California.csv")

names(out_seqs)

# Add sabo S, M and L as out groups
fasta_out_group <- list(
  S = out_seqs["NC_055117.1 |California encephalitis virus|USA: California|Ochlerotatus dorsalis"],
  M = out_seqs["NC_055118.1 |California encephalitis virus|USA: California|Ochlerotatus dorsalis"],
  L = out_seqs["NC_055119.1 |California encephalitis virus|USA: California|Ochlerotatus dorsalis"]
)

# Export 
lapply(names(seqs_by_seg), function(seg){
  
  # seg <- "L"
  out_group <- fasta_out_group[[seg]]
  names(out_group) <- glue("California_{seg}") # Does not have a name 
  comb <- c(seqs_by_seg[[seg]], out_group)
  
  writeXStringSet(comb, glue("helene/out/SML_{seg}.fasta"))
  # writeXStringSet(seqs_by_seg[[seg]], glue("helene/out/{seg}.fasta"))
  
})


