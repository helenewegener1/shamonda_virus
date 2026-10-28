library(tidyverse)
library(glue)

setwd("~/Documents/courses/MolecularEvolutionPhylogenomics/shamonda_virus/")

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

table(meta_all_clean$Segment, useNA = "always")

meta_all_clean %>%
  mutate(Segment = factor(replace_na(as.character(Segment), "Unassigned"),
                          levels = c("S", "M", "L", "Unassigned"))) %>%
  add_count(Segment, name = "n_seg") %>%
  mutate(Segment_lab = fct_reorder(sprintf("%s (n = %d)", Segment, n_seg),
                                   as.integer(Segment))) %>%
  ggplot(aes(x = Length)) +
  geom_histogram(bins = 60) +
  facet_wrap(vars(Segment_lab), ncol = 1, scales = "free_y") +
  theme_bw() +
  labs(title = "Size of segments", x = "Sequence length (nt)", y = "Count")
  
ggsave("helene/plot/segment_size.png", width = 8, height = 7)

# meta_all_clean$header <- glue("")
  
# Divide accession numbers in segments 
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


# # test
# names(seqs_all)
# 
# meta_all_clean$Accession %>% length()
# meta_all_clean$Accession %>% unique() %>% length()
# 
# meta_all_clean %>% arrange(Accession) %>% view()
# AB183279
# seqs_all$`AB183279.1 ||Japan:Miyazaki|`
# seqs_all$`AB183279.1 ||Japan:Miyazaki|`

grep("AB183279", names(seqs_all), value = TRUE)

# Combine sequnces 
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

# Add sabo S, M and L as out groups
fasta_out_group <- list(
  "S" = readDNAStringSet("data/sequences_sabo_S.fasta"), 
  "M" = readDNAStringSet("data/sequences_sabo_M.fasta"), 
  "L" = readDNAStringSet("data/sequences_sabo_L.fasta")
)

# Export 

lapply(names(fasta_out_group), function(seg){
  
  # seg <- "L"
  out_group <- fasta_out_group[[seg]]
  names(out_group) <- names(out_group) %>% str_split_i("\\|", 1) %>% str_trim()
  comb <- c(seqs_by_seg[[seg]], out_group)

  writeXStringSet(comb, glue("helene/out/{seg}.fasta"))
  # writeXStringSet(seqs_by_seg[[seg]], glue("helene/out/{seg}.fasta"))
  
})

# With NA included 
names(seqs_by_seg[["NA"]]) <- glue('{names(seqs_by_seg[["NA"]])}_NA')

lapply(names(seqs_by_seg), function(seg){
  
  # seg <- "L"
  
  names(seqs_by_seg[[seg]]) <- glue("{names(seqs_by_seg[[seg]])}_{seg}")
  
  comb <- c(seqs_by_seg[[seg]], seqs_by_seg[["NA"]])
  
  writeXStringSet(comb, glue("helene/out/{seg}_NA.fasta"))
  
})

length(seqs_by_seg[["NA"]])



