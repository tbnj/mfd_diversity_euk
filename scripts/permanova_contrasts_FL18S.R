### Setup env
library(tidyverse)
library(vegan)

setwd("/mfd_diversity_euk")

### Import data
meta <- data.table::fread("output/2025-04-16_MFD_FL18S_metadata_filtered.tsv") %>%
  mutate(across(complex, ~factor(.))) %>%
  mutate(across(complex, ~droplevels(.)))

dist <- data.table::fread("output/2025-04-16_MFD_FL18S_BC_distance.csv") %>%
  column_to_rownames(var = "V1") %>%
  mutate(across(where(is.numeric), ~round(., 4)))

# dist <- dist %>%
#   select(all_of(meta$fieldsample_barcode)) %>%
#   filter(rownames(.) %in% meta$fieldsample_barcode)

## Extract levels
levels <- meta %>%
  pull(complex) %>%
  unique() %>%
  as.character()

meta <- meta %>%
  mutate(across(complex, ~factor(., levels = levels)))

## Group summary
group.summary <- meta %>%
  group_by(complex) %>%
  reframe(n = n()) %>%
  arrange(factor(complex, levels = levels))

### Run PERMANOVA
set.seed(123)
permanova.total <- adonis2(dist ~ complex, data = meta, method = "jaccard", by = "terms", permutations = 999, parallel = 10)

df.permanova.total <- data.frame(permanova.total) %>%
  rownames_to_column("Model") %>%
  mutate(across(Model, ~str_replace(., "Model", "complex")),
         across(R2, ~ as.numeric(round(. * 100, digits = 3)))) %>%
  mutate(pval = ifelse(Pr..F. < 0.001, "<0.001", round(Pr..F., 3))) %>%
  select(Model, R2, pval)

gc()


## Write to output
data.table::fwrite(df.permanova.total, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_FL18S_PERMANOVA.tsv"))

## Perform contrasts with PERMANOVA
contrasts.dataframe <- meta %>%
  select(fieldsample_barcode, complex) %>%
  mutate(n = 1) %>%
  pivot_wider(names_from = complex, values_from = n, values_fill = list(n = -1)) %>%
  column_to_rownames(var = "fieldsample_barcode") %>%
  as.matrix()

# set.seed(123)
# permanova.test <- adonis2(dist ~ contrasts.dataframe[, 9], by = "terms", permutations = 999, parallel = 10)
# 
# df.permanova.test <- data.frame(permanova.test) %>%
#   rownames_to_column("Model") %>%
#   mutate(across(Model, ~str_replace(., "Model", "complex")),
#          across(R2, ~ as.numeric(round(. * 100, digits = 3)))) %>%
#   mutate(pval = ifelse(Pr..F. < 0.001, "<0.001", round(Pr..F., 3))) %>%
#   select(Model, R2, pval)
# 
# gc()

perm.list = lst()

levels.test <- levels

for (i in 1:length(levels.test)){
  
  test <- contrasts.dataframe[,i]
  
  set.seed(123)
  permanova<- vegan::adonis2(dist ~ test, method = "jaccard",
                             by = "terms", permutations = 999, parallel = 10)
  
  df <- permanova %>%
    data.frame() %>%
    rownames_to_column("Model") %>%
    mutate(across(Model, ~str_replace(., "test", levels.test[i])),
           across(R2, ~ as.numeric(round(. * 100, digits = 3)))) %>%
    mutate(pval = ifelse(Pr..F. < 0.001, "<0.001", round(Pr..F., 3))) %>%
    select(Model, R2, pval)
  
  perm.list[[levels.test[i]]] <- df
  
  gc()
}


## Combine outputs
df.permanova.contrasts <- bind_rows(perm.list) %>%
  filter(Model %in% levels) %>%
  arrange(factor(Model, levels = levels))


## Write to output
data.table::fwrite(df.permanova.contrasts, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_FL18S_PERMANOVA_contrasts.tsv"))

### Save image
save.image(paste0('output/', format(Sys.time(), "%Y-%m-%d"), "_MFD_FL18S_PERMANOVA.RData"))
