library(tidyverse)
library(glue)

setwd("~/Documents/courses/MolecularEvolutionPhylogenomics/")

# Amanda is aligning using MAFTT 

# ------------------------------------------------------------------------------
# Load data 
# ------------------------------------------------------------------------------

meta_all_clean <- readRDS("helene/out/meta_all_clean.rds")


# Switzerland
# S: 875-900
# M: 4376
# L: 6700-6800
# meta_all_clean$Isolate

meta_all_clean %>% 
  ggplot(aes(x = Length)) + 
  geom_histogram(position = "dodge") +
  facet_wrap(vars(Segment), ncol = 1) + 
  theme_bw()


# meta_all_clean$header <- glue("")
  
# 
segment <- list()
segment[["L"]] <- meta_all_clean %>% filter(Segment == "L") %>% select(Accession) %>% pull()
segment[["S"]] <- meta_all_clean %>% filter(Segment == "S") %>% select(Accession) %>% pull()
segment[["M"]] <- meta_all_clean %>% filter(Segment == "M") %>% select(Accession) %>% pull()
segment[["NA"]] <- meta_all_clean %>% filter(is.na(Segment)) %>% select(Accession) %>% pull()

# ------------------------------------------------------------------------------
# Split fasta by segments 
# ------------------------------------------------------------------------------

library(Biostrings)

fasta1 <- readDNAStringSet("data/sequences_shamonda.fasta")             # Shamonda, taxid 15950
fasta2 <- readDNAStringSet("data/sequences_Schmallenbergvirus1.fasta")  # Orthobunyavirus schmallenbergense, taxid 3052437
fasta3 <- readDNAStringSet("data/sequences_Schmallenbergvirus2.fasta")  # Schmallenberg virus, taxid 1133363

seqs_all <- c(fasta1, fasta2, fasta3)
names(seqs_all) <- str_extract(names(seqs_all), "^[^.\\s]+")

# drop duplicates in case the three taxid queries overlap
seqs_all <- seqs_all[!duplicated(names(seqs_all))]

# Check Accessions
segment[["L"]] %in% names(seqs_all) %>% table()

# 
seqs_by_seg <- map(segment, ~ seqs_all[names(seqs_all) %in% .x])

map_int(seqs_by_seg, length)
map_int(segment, length)

# Export 
names(seqs_by_seg)

lapply(names(seqs_by_seg), function(seg){
  
  writeXStringSet(seqs_by_seg[[seg]], glue("helene/out/{seg}.fasta"))
  
})

# With NA included 
names(seqs_by_seg[["NA"]]) <- glue('{names(seqs_by_seg[["NA"]])}_NA')

lapply(names(seqs_by_seg), function(seg){
  
  # seg <- "L"
  
  names(seqs_by_seg[[seg]]) <- glue("{names(seqs_by_seg[[seg]])}_{seg}")
  
  comb <- c(seqs_by_seg[[seg]], seqs_by_seg[["NA"]])
  
  writeXStringSet(comb, glue("helene/out/{seg}_NA.fasta"))
  
})










seqs_all <- c(fasta1, fasta2, fasta3)

# segment assignment from the FASTA headers, adjust pattern to what yours look like
seg <- case_when(
  str_detect(names(seqs_all), regex("segment S|nucleocapsid|\\bN\\b", ignore_case = TRUE)) ~ "S",
  str_detect(names(seqs_all), regex("segment M|glycoprotein|polyprotein", ignore_case = TRUE)) ~ "M",
  str_detect(names(seqs_all), regex("segment L|polymerase|RdRp",         ignore_case = TRUE)) ~ "L",
  TRUE ~ NA_character_
)
table(seg, useNA = "ifany")

seg <- cut(width(seqs_all),
           breaks = c(0, 2000, 5500, Inf),
           labels = c("S", "M", "L"))


