library(tidyverse)
library(glue)

setwd("~/Documents/courses/MolecularEvolutionPhylogenomics/shamonda_virus/")

# Amanda is aligning using MAFTT 

# ------------------------------------------------------------------------------
# Load data 
# ------------------------------------------------------------------------------

meta_all_clean <- readRDS("helene/out/meta_all_clean.rds")

fasta <- readDNAStringSet("helene/out/M.fasta")


# fasta_S <- readDNAStringSet("helene/out/SML_S_noout.fasta") 
fasta_M <- readDNAStringSet("helene/out/SML_M_noout.fasta") 
fasta_L <- readDNAStringSet("helene/out/SML_L_noout.fasta") 




