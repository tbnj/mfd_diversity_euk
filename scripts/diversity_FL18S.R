### Setup env
library(tidyverse)
library(vegan)
library(rstatix)
library(multcompView)
library(ampvis2)
library(patchwork)

### Import color palettes
setwd("mfd_diversity_euk")

source('scripts/MFD_colors.R')


### Import data
## Load metadata
metadata <- data.table::fread("output/2025-04-16_MFD_FL18S_metadata_subset.tsv", sep = "\t", header = T) %>%
  select(fieldsample_barcode, everything())

## Subset metadata to remove groups with low sample counts
groups.summary <- metadata %>%
  group_by(complex) %>%
  summarise(samples = n()) %>%
  arrange(samples)

metadata.subset <- metadata %>%
  filter(!complex %in% c("Other, Urban, Landfill",
                         "Soil, Urban, Roadside",
                         "Water, Subterranean, Freshwater",
                         "Other, Urban, Drinking water",
                         "Water, Urban, Drinking water",
                         "Soil, Subterranean, Urban",
                         "Sediment, Subterranean, Saltwater"))

## Create new group summary of the subset
groups.summary.subset <- metadata.subset %>%
  group_by(complex) %>%
  summarise(samples = n()) %>%
  arrange(samples)

## Extract MFDO1 levels
levels.subset <- metadata.subset %>%
  mutate(across(complex, ~factor(.))) %>%
  pull(complex) %>%
  droplevels() %>%
  levels()

## Create new group summary of the filtered set
groups.summary.filt <- groups.summary.subset %>%
  filter(samples >= 9)

metadata.filt <- metadata.subset %>%
  filter(complex %in% groups.summary.filt$complex)

## Extract MFDO1 levels
levels.filt <- metadata.filt %>%
  mutate(across(complex, ~factor(.))) %>%
  pull(complex) %>%
  droplevels() %>%
  levels()

## Subset color-palette based on MFDO1 levels
mfdo1.palette.subset <- mfdo1.palette[levels.subset]
mfdo1.palette.filt <- mfdo1.palette[levels.filt]

## Load observational table submitted to random subsampling without replacement
fl.asv.subset.ra <- data.table::fread("output/2025-04-16_MFD_FL18S_OTU_filtered_ra_pa.csv", sep = ",", header = TRUE)

### Alpha diversity
ampvis <- amp_load(metadata = metadata.filt,
                   otutable = fl.asv.subset.ra)

fl.alpha.diversity <- ampvis %>%
  amp_alphadiv(rarefy = NULL,
               measure = c("uniqueotus","shannon","simpson", "invsimpson"),
               richness = T) %>%
  rename(rarefaction_level = Reads,
         observed_species = uniqueOTUs) %>%
  mutate(across(Shannon:Simpson, ~round(., 5)),
         across(Chao1:ACE, ~round(., 0))) %>%
  select(complex, observed_species:ACE) %>%
  rownames_to_column(var = "fieldsample_barcode") %>%
  filter(complex %in% levels.filt) %>%
  mutate(across(complex, ~factor(., levels = rev(levels.filt))))

## Summary of soil MFDO1 categories
fl.alpha.diversity %>%
  group_by(complex) %>%
  summarise(group_size = n(),
            median = median(observed_species),
            mean = mean(observed_species),
            sd = sd(observed_species)) %>%
  arrange(desc(median))

### Statitical test
## Test for normality (Shapiro-Wilk)
fl.alpha.diversity %>%
  group_by(complex) %>%
  shapiro_test(observed_species)

## Test for homoscedasticity
fl.alpha.diversity %>%
  levene_test(observed_species ~ complex)

fl.alpha.diversity %>%
  kruskal_test(observed_species ~ complex)

fl.alpha.diversity %>%
  wilcox_test(observed_species ~ complex, p.adjust.method = "BH")


## Compact letter display
wilcox.pair <- fl.alpha.diversity %>%
  wilcox_test(observed_species ~ complex, alternative = "two.sided", p.adjust.method = "BH", detailed = T) %>%
  mutate(comp = str_c(group1, "-", group2), .keep = "unused", .before = "p.adj") %>%
  select(comp, everything())

wilcox.vec <- setNames(wilcox.pair$p.adj, wilcox.pair$comp)

cld <- multcompLetters(wilcox.vec)

## Attanged table of alpha diversity
dt <- fl.alpha.diversity %>%
  group_by(complex) %>%
  summarise(median = median(observed_species),
            mean=mean(observed_species), 
            sd = sd(observed_species)) %>%
  arrange(desc(median))

## Extract the compact letter display and adding to the table
cld <- as.data.frame.list(cld) %>%
  rownames_to_column(var = "complex")

dt <- dt %>% 
  left_join(cld %>% select(complex, Letters)) %>%
  rename(cld = Letters)

## Investigate
print(dt)

## Write to output
data.table::fwrite(dt, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_FL18S_alpha_diveristy_analysis.csv"))


## Visualise disribution of alpha diversity across MFDO1
plot.order <- fl.alpha.diversity %>%
  group_by(complex) %>%
  reframe(median_alpha = median(observed_species)) %>%
  arrange(median_alpha) %>%
  pull(complex)

plot.alpha <- fl.alpha.diversity %>%
  ggplot(aes(x = complex, y = observed_species)) +
  geom_violin(aes(fill = complex)) +
  geom_boxplot(width = 0.2, alpha = 0.5, fill = 'white', outlier.colour = NA) +
  geom_text(data = dt, aes(x = complex, y = 750, label = cld),
            size = 4, fontface = "bold") +
  scale_x_discrete(limits = plot.order) +
  scale_fill_manual(values = mfdo1.palette.subset,
                    name = "MFDO1") +
  xlab('MFDO1') +
  ylab('Observed species') +
  coord_flip() +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "none")

## Render plot
plot.alpha

## Write alpha table to output
data.table::fwrite(fl.alpha.diversity, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_FL18S_alpha_diversity_MFDO1.tsv"))


### Gamma
## Load gamma table
data.gamma.ra <- data.table::fread("output/2025-04-16_MFD_FL18S_iNEXT_richness_MFDO1_ra.tsv", sep = "\t", header = TRUE)

## Combine into data frame
fl.gamma.diversity <- data.gamma.ra %>%
  mutate(across(Diversity, ~factor(., levels = c("Species richness", "Shannon diversity", "Simpson diversity")))) %>%
  rename(complex = Assemblage,
         Hill_diversity = Diversity) %>% 
  right_join(groups.summary.filt) %>%
  mutate(across(Hill_diversity, ~replace_na(., "Species richness"))) %>%
  complete(complex, Hill_diversity) %>%
  mutate(across(complex, ~factor(., levels = levels.filt))) %>%
  filter(!is.na(complex))

## Plot gamma diversity
plot.gamma.total <- fl.gamma.diversity %>%
  ggplot(aes(x = complex, y = Estimator, group = Hill_diversity, fill = complex)) +
  geom_bar(position=position_dodge(), aes(y=Estimator), stat="identity") +
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0.5, position=position_dodge(width=0.9)) +
  scale_x_discrete(limits = plot.order) +
  scale_fill_manual(values = mfdo1.palette.subset,
                    name = "MFDO1") +
  facet_grid(cols = vars(Hill_diversity), scales = "free") +
  coord_flip() +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

plot.gamma.total

## Plot of Shannon diversity (q=1)
plot.gamma.shannon <- fl.gamma.diversity %>%
  filter(Hill_diversity == "Shannon diversity") %>%
  ggplot(aes(x = complex, y = Estimator, group = Hill_diversity, fill = complex)) +
  geom_col() +
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0.5, position=position_dodge(width=0.9)) +
  # geom_point(shape = 21) +
  geom_text(aes(y = 5000, label = str_c("N = ", samples)), size = 4,
            fontface = "bold") +
  scale_x_discrete(limits = plot.order) +
  scale_fill_manual(values = mfdo1.palette.subset,
                    name = "MFDO1") +
  scale_y_continuous(limits = c(0,6000)) +
  xlab('') +
  ylab('Shannon diversity') +
  coord_flip() +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_blank())

## Render plots
plot.alpha + plot.gamma.shannon

## Write gamma table to output
data.table::fwrite(fl.gamma.diversity, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_FL18S_gamma_diversity_MFDO1.tsv"))


### Combine alpha and gamma table
table.diversity <- fl.alpha.diversity %>%
  group_by(complex) %>%
  summarise(median_alpha = round(median(observed_species), 0),
            mean_alpha = round(mean(observed_species), 0),
            sd_alpha = round(sd(observed_species), 0)) %>%
  right_join(fl.gamma.diversity)

## Write gamma table to output
data.table::fwrite(table.diversity, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_FL18_combined_diversity_MFDO1.tsv"))


### Save image
save.image(paste0('output/', format(Sys.time(), "%Y-%m-%d"), "_MFD_FL18S_diversity.RData"))


