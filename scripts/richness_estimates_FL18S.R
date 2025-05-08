### Setup env
library(tidyverse)
library(vegan)
library(iNEXT)
library(patchwork)
library(rstatix)

### Import color palettes
source('scripts/MFD_colors.R')

setwd("/mfd_diversity_euk")


### Format data
## Import ASV data aggregated to OTUs (98.7% similarity ~ Species representatives)
fl.asv <- data.table::fread("data/2025-04-16_MFD_FL18S_OTU.csv", sep = ",", header = TRUE) %>%
  rename(OTU = V1) %>%
  filter(rowSums(across(where(is.numeric)))!=0)

tax <- fl.asv %>%
  select(OTU, Kingdom:Species)

## Summarise read counts per sample
read.count <- fl.asv %>%
  select(where(is.numeric)) %>%
  colSums() %>%
  data.frame("reads" = .) %>%
  rownames_to_column(var = "fieldsample_barcode") %>%
  arrange(reads)

## Import metdata file and create "complex" correponding to full MFDO1 string
metadata <- readxl::read_excel('data/2025-04-14_mfd_db.xlsx') %>%
  filter(fieldsample_barcode %in% colnames(fl.asv)) %>%
  relocate(coords_reliable, .after = "longitude") %>%
  mutate(across(mfd_hab1, ~str_replace(., "Sclerophyllous scrub", "Temperate heath and scrub")),
         across(mfd_hab2, ~str_replace(., "Scrub", "Sclerophyllous scrub"))) %>%
  mutate(complex = str_c(mfd_sampletype, mfd_areatype, mfd_hab1, sep = ", "),
         across(complex, ~factor(.))) %>%
  left_join(read.count) %>%
  mutate(across(mfd_sampletype, ~factor(., levels = c("Soil", "Sediment", "Water", "Other"))),
         across(mfd_areatype, ~factor(., levels = c("Natural", "Subterranean", "Agriculture", "Urban")))) %>%
  arrange(mfd_sampletype, mfd_areatype)

### Filter data
## Create summary of read count per MFDO1 level
read.count.summary <- metadata %>%
  group_by(complex) %>%
  summarise(n = n(), 
            effort = round(mean(reads), 0))

## Select samples based on minimum read count (<10,000 reads)
min.reads <- read.count %>%
  filter(reads > 6000) %>%
  pull(reads) %>%
  min()

## Create summary on MFDO1 before filtering
groups.summary <- metadata %>%
  group_by(complex) %>%
  summarise(size = n()) %>%
  arrange(desc(size))

## Create summary on MFDO1 after filtering
groups.summary.subset <- metadata %>%
  filter(reads >= min.reads) %>%
  group_by(complex) %>%
  summarise(size = n()) %>%
  arrange(desc(size))

## Remove groups with less than 10 representative samples
groups.summary.filt <- groups.summary.subset %>%
  filter(size >= 9) %>%
  arrange(desc(size))

## Create metadata subset based on read counts
metadata.subset <- metadata %>%
  filter(reads >= min.reads) %>%
  mutate(across(complex, ~factor(., levels = names(mfdo1.palette)))) %>%
  filter(!is.na(complex))

## Filter metadata based on representative samples
metadata.filt <- metadata.subset %>%
  filter(complex %in% c(groups.summary.filt %>% pull(complex)))

### Write metadata files to output
data.table::fwrite(metadata.subset, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_FL18S_metadata_subset.tsv"))
data.table::fwrite(metadata.filt, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_FL18S_metadata_filtered.tsv"))


## Pull fieldsample_barcode from metadata subset
samples.subset <- metadata.subset %>%
  pull(fieldsample_barcode)

## Pull fieldsample_barcode from filtered metadata
samples.filt <- metadata.filt %>%
  pull(fieldsample_barcode)

## Filter observational table based on subsetting
fl.asv.subset <- fl.asv %>%
  select(OTU, any_of(samples.subset), Kingdom:Species) %>%
  filter(rowSums(across(where(is.numeric)))!=0)

## Filter observational table based on filtering
fl.asv.filt <- fl.asv %>%
  select(OTU, any_of(samples.filt), Kingdom:Species) %>%
  filter(rowSums(across(where(is.numeric)))!=0)

## Evaluate differences based on subsetting and filtering
colSums(Filter(is.numeric, fl.asv)) %>% sum()
colSums(Filter(is.numeric, fl.asv.subset)) %>% sum()
colSums(Filter(is.numeric, fl.asv.filt)) %>% sum()

nrow(fl.asv)-nrow(fl.asv.subset) # difference of 169 species
nrow(fl.asv)-nrow(fl.asv.filt) # difference of 5,466 species

### Write metadata files to output
data.table::fwrite(fl.asv.subset, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_FL18S_OTU_subset.csv"))
data.table::fwrite(fl.asv.filt, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_FL18S_OTU_filtered.csv"))


### Perform random subsample without replacement
## Random subsample without replacement for subsetted data
set.seed(123)
fl.asv.subset.ra <- fl.asv.subset %>%
  column_to_rownames(var = "OTU") %>%
  select(where(is.numeric)) %>%
  t() %>%
  rrarefy(., sample = min.reads) %>%
  t() %>%
  data.frame() %>%
  filter(rowSums(across(where(is.numeric)))!=0) %>%
  rownames_to_column(var = "OTU") %>%
  left_join(tax)

## Random subsample without replacement for filtered data
set.seed(123)
fl.asv.filt.ra <- fl.asv.filt %>%
  column_to_rownames(var = "OTU") %>%
  select(where(is.numeric)) %>%
  t() %>%
  rrarefy(., sample = min.reads) %>%
  t() %>%
  data.frame() %>%
  filter(rowSums(across(where(is.numeric)))!=0) %>%
  rownames_to_column(var = "OTU") %>%
  left_join(tax)


### Reduce color palettes
## Pull levels from different metadata files
levels <- metadata %>%
  pull(complex) %>%
  levels()

levels.subset <- metadata.subset %>%
  pull(complex) %>%
  droplevels() %>%
  levels()

levels.filt <- metadata.filt %>%
  pull(complex) %>%
  droplevels() %>%
  levels()

## Filter color palettes
mfdo1.palette.subset <- mfdo1.palette[levels.subset]
mfdo1.palette.filt <- mfdo1.palette[levels.filt]

## Evaluate dropped levels in the different sets
setdiff(levels, levels.subset)
setdiff(levels.subset, levels)

setdiff(levels, levels.filt)
setdiff(levels.filt, levels)


### Write files to output
data.table::fwrite(fl.asv.subset.ra, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_FL18S_OTU_subset_ra.csv"))
data.table::fwrite(fl.asv.filt.ra, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_FL18S_OTU_filtered_ra.csv"))


## Presence/Absence
fl.asv.subset.ra.pa <- fl.asv.subset.ra %>%
  mutate(across(where(is.numeric), ~+as.logical(.x)))

fl.asv.filt.ra.pa <- fl.asv.filt.ra %>%
  mutate(across(where(is.numeric), ~+as.logical(.x)))

### Write metadata files to output
data.table::fwrite(fl.asv.subset.ra.pa, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_FL18S_OTU_subset_ra_pa.csv"))
data.table::fwrite(fl.asv.filt.ra.pa, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_FL18S_OTU_filtered_ra_pa.csv"))


### Estimate DK species richness

## Format data for iNEXT (change to presence/absence format)
data.total <- fl.asv.filt %>%
  column_to_rownames(var = "OTU") %>%
  select(-c(Kingdom:Species)) %>%
  mutate(across(everything(), ~+as.logical(.x)))

list <- lst()

list[["data"]] <- data.total

## Run iNEXT for all data
## Endpoint of extrapolation is by default 2x sample size (n)
set.seed(123)
object.iNEXT.total <- iNEXT(list, q = c(0, 1, 2), datatype = "incidence_raw")

## Inspect iNEXT summary
## Observed and estimates of total richness is given with 95% CI 
## for total as well as Shannon and Simpson diversity
object.iNEXT.total

data.table::fwrite(as.data.frame(object.iNEXT.total$AsyEst), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_FL18S_iNEXT_richness_total.tsv"))

## Create function to format data for iNEXT by MFDO1 categories
format.iNEXT <- function(df, levels, metadata) {
  data.iNEXT <- list()
  
  abund <- df %>%
    column_to_rownames(var = "OTU") %>%
    select(-c(Kingdom:Species)) %>%
    mutate(across(everything(), ~+as.logical(.x)))
  
  for (i in 1:length(levels)) {
    filter <- metadata %>%
      filter(complex == levels[i]) %>%
      pull(fieldsample_barcode)
    
    data.filt <- abund %>%
      select(any_of(filter)) %>%
      filter(rowSums(across(where(is.numeric)))!=0) %>%
      as.matrix()
    
    data.iNEXT[[levels[i]]] <- data.filt
  }
  return(data.iNEXT)
}

## Run formatting function
data.iNEXT <- format.iNEXT(fl.asv.filt, levels.filt, metadata.filt)
data.iNEXT.ra <- format.iNEXT(fl.asv.filt.ra, levels.filt, metadata.filt)

## Run iNEXT for MFDO1 categories
## Endpoint for extrapolation fixed at n=100
set.seed(123)
object.iNEXT <- iNEXT(data.iNEXT, q = c(0, 1, 2), datatype = "incidence_raw", endpoint = 100)

set.seed(123)
object.iNEXT.ra <- iNEXT(data.iNEXT.ra, q = c(0, 1, 2), datatype = "incidence_raw", endpoint = 100)

## Inspect iNEXT summary
## Observed and estimates of total richness is given with 95% CI 
## for total as well as Shannon and Simpson diversity
richness.measure <- object.iNEXT$DataInfo %>% select(1:5) %>%
  left_join(object.iNEXT$AsyEst)

data.table::fwrite(richness.measure, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_FL18S_iNEXT_richness_MFDO1.tsv"))


richness.measure.ra <- object.iNEXT.ra$DataInfo %>% select(1:5) %>%
  left_join(object.iNEXT.ra$AsyEst)

data.table::fwrite(richness.measure.ra, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_FL18S_iNEXT_richness_MFDO1_ra.tsv"))

## Evaluate effect of sampling effort
df.cor <- object.iNEXT.ra$DataInfo %>%
  filter(Assemblage %in% names(mfdo1.palette.filt))

## Visualise data with data and fit GAM model
ggplot(df.cor, aes(T, SC)) + geom_point() + 
  geom_smooth(se = FALSE, method = "lm", formula = y ~ log(x))

## Test for normality and measure correlation
df.cor %>%
  shapiro_test(T)

df.cor %>%
  mutate(across(T, ~log(.))) %>%
  shapiro_test(T)

df.cor %>%
  shapiro_test(SC)

df.cor %>%
  mutate(across(T, ~log(.))) %>%
  cor_test(T, SC, method = "pearson")


## Change grouping value to factor
object.iNEXT$DataInfo$Assemblage <- factor(object.iNEXT.ra$DataInfo$Assemblage, levels = levels.subset)


### Visualise results from iNEXT using internal iNEXT ggplot function
## Species accumulation curve per MFDO1 category
plot.iNEXT <- ggiNEXT(object.iNEXT, type = 1, color.var = "Assemblage", facet.var = "Order.q") +
  scale_color_manual(values = mfdo1.palette.subset, breaks = levels.subset,
                     name = "MFDO1") +
  scale_fill_manual(values = mfdo1.palette.subset,
                    breaks = levels.subset,
                    name = "MFDO1") +
  ggtitle("Species accumulation curves",
          subtitle = "Interpolation and Extrapolation") +
  labs(x = "Number of sampling sites",
       y = "Number of species") +
  theme_bw(base_size = 12) +
  theme(title = element_text(face = "bold"),
        legend.position = "right",
        legend.box = "vetical")

## Remove point layer (curve end-point)
plot.iNEXT$layers[[1]] <- NULL

## Re-plot
plot.iNEXT

### save plots
png(file = 'output/MFD_FL18S_estimated_richness.png',
    width = 1900,
    height = 1000) 
plot.iNEXT
dev.off()

pdf(file = 'output/MFD_FL18S_estimated_richness.pdf',
    width = 19,
    height = 12)
plot.iNEXT
dev.off()

tiff(file = 'output/MFD_FL18S_estimated_richness.tiff',
     width = 1900,
     height = 1200)
plot.iNEXT
dev.off()

ggsave("output/MFD_FL18S_estimated_richness.svg", plot = plot.iNEXT, width = 19, height = 12, units = "in", dpi = "retina")

### Save image
save.image(paste0('output/', format(Sys.time(), "%Y-%m-%d"), "_MFD_FL18S_estimated_richness.RData"))
