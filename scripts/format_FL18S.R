### Setup env
library(tidyverse)

setwd("/mfd_diversity_euk")

## Import OTU table
operon <- data.table::fread("data/MFD_18S_3nfd-15810R_OTUtab_MFD_only.txt", sep = "\t", header = TRUE) %>%
  rename(OTU = 1) %>%
  rename_with(., ~str_remove(., "_pool"),
              starts_with("MFD")) %>%
  # select(-MFD10339) %>%
  filter(rowSums(across(where(is.numeric)))!=0)

samples <- colnames(operon) %>%
  str_subset(., "MFD")

OTUs <- operon$OTU %>%
  unique()

## Import OTU taxonomy and format
taxonomy <- data.table::fread("data/dada2_taxonomy_18s_PR2.csv") %>%
  select(-Supergroup) %>%
  mutate(Phylum = ifelse(!is.na(Subdivision), 
                         paste(Division, Subdivision, sep = "-"), 
                         Division)) %>%
  select(-Division, -Subdivision) %>% # Drop old columns
  rename(OTU = ASV,
         Kingdom = Domain) %>%
  select(OTU, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
  mutate(across(everything(), ~replace_na(., "Unclassified")))

## Combine tables
df.operon <- operon %>%
  left_join(taxonomy, by = "OTU") %>%
  select(OTU, starts_with("MFD"), Kingdom:Species) %>%
  filter(rowSums(across(where(is.integer))) != 0) %>%
  column_to_rownames(var = "OTU") %>%
  filter(Kingdom == "Eukaryota") %>%
  mutate(across(Kingdom:Species, ~na_if(., "")))


## Write file to output directory
data.table::fwrite(df.operon, sep = ",", row.names = TRUE, col.names = TRUE, quote = FALSE,
                   paste0("data/", format(Sys.time(), "%Y-%m-%d"), "_MFD_FL18S_OTU.csv"))
