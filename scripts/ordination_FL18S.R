### Setup env
library(tidyverse)
library(vegan)
library(ape)
library(patchwork)

source('scripts/MFD_colors.R')

setwd("/mfd_diversity_euk")

### Load data
data <- data.table::fread('output/2025-04-16_MFD_FL18S_OTU_filtered_ra.csv', na.strings = "")

metadata <- data.table::fread("output/2025-04-16_MFD_FL18S_metadata_filtered.tsv")

### Beta diversity
## Calculate Hellinger-transformed Bray-Curtis dissimilarity in parallel
dist <- data %>%
  select(where(is.numeric)) %>%
  t() %>% 
  parallelDist::parDist(., method = "binary", threads = 20) %>%
  as.matrix() %>%
  data.frame()

## Save Bray-Curtis matrix
data.table::fwrite(dist, sep = ",", row.names = TRUE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_FL18S_BC_distance.csv"))


### Create subset for PCOA
## Filter metadata and create "complex" corresponding to full MFDO1 string
metadata.select <- metadata %>%
  filter(fieldsample_barcode %in% colnames(dist),
    !is.na(mfd_hab1)) %>%
  select(fieldsample_barcode, mfd_sampletype:mfd_hab3) %>%
  mutate(complex = str_c(mfd_sampletype, mfd_areatype, mfd_hab1, sep = ", ")) %>%
  # filter(complex %in% filter) %>%
  mutate(across(complex, ~factor(., levels = names(mfdo1.palette)))) %>%
  group_by(complex) %>%
  mutate(complex_size = n()) %>%
  ungroup()

## Create group summary
groups.summary <- metadata.select %>%
  group_by(complex) %>%
  summarise(size = n()) %>%
  arrange(size)

## Subset metadata to remove groups with low sample counts
metadata.subset <- metadata.select %>%
  filter(!complex %in% c("Other, Urban, Landfill",
                         "Soil, Urban, Roadside",
                         "Water, Subterranean, Freshwater",
                         "Other, Urban, Drinking water",
                         "Water, Urban, Drinking water",
                         "Soil, Subterranean, Urban",
                         "Sediment, Subterranean, Saltwater")) %>%
  mutate(across(complex, ~droplevels(.)))

## Create new group summary of the subset
groups.summary.subset <- metadata.subset %>%
  group_by(complex) %>%
  summarise(size = n()) %>%
  arrange(size)

## Extract sample IDs from samples with MFDO1 information
samples.subset <- metadata.subset %>%
  filter(!is.na(mfd_hab1)) %>%
  pull(fieldsample_barcode)

## Extract MFDO1 levels
levels.subset <- metadata.subset %>%
  pull(complex) %>%
  droplevels() %>%
  levels()

## Subset color-palette based on MFDO1 levels
mfdo1.palette.subset <- mfdo1.palette[levels.subset]

### PCoA - Principal Coordinate Analysis
## Import Bray-Curtis matrix (if uncommented)
# dist <- data.table::fread("output/") %>% column_to_rownames(var = "V1")

## Perform the pcoa
PCOA <- pcoa(dist)

## plot the eigenvalues and interpret 
barplot(PCOA$values$Relative_eig[1:10])

## Get percentage of variance explained by the first 3 principal coordinates
PCO1 <- round(sum(as.vector(PCOA$value$Relative_eig)[1])*100, 1)
PCO2 <- round(sum(as.vector(PCOA$value$Relative_eig)[2])*100, 1)
PCO3 <- round(sum(as.vector(PCOA$value$Relative_eig)[3])*100, 1)

## PCO1 and PCO2 explains 30% of the variation in the communities
sum(as.vector(PCOA$value$Relative_eig)[1:2])*100

## Extract the scores for plotting with ggplot
PCOA.scores <- PCOA$vectors %>%
  as.data.frame() %>%
  select(1:6) %>%
  rename_with(., ~str_replace(., "Axis.", "PCO")) %>%
  cbind(fieldsample_barcode = colnames(dist)) %>%
  relocate(fieldsample_barcode, .before = "PCO1") %>%
  left_join(metadata.subset)

### Stats based on spatial subset
## Betadisper
betadisper <- betadisper(as.dist(dist), metadata.subset$complex)

betadisper.aov <- anova(betadisper)

betadisper.aov

TukeyHSD(betadisper)

## Anosim
set.seed(123)
anosim <- anosim(as.dist(dist), metadata.subset$complex, permutations = 999, parallel = 10)

anosim

## PERMANOVA - import
df.permanova.total <- data.table::fread("output/2025-04-16_MFD_FL18S_PERMANOVA.tsv")

df.permanova.contrasts <- data.table::fread("output/2025-04-16_MFD_FL18S_PERMANOVA_contrasts.tsv")


## Add labels to scores-object
PCOA.scores <- PCOA.scores %>%
  mutate(across(complex, ~factor(., levels = names(mfdo1.palette.subset)))) %>%
  filter(!is.na(complex)) %>%
  arrange(complex) %>%
  left_join(df.permanova.contrasts, by = c("complex" = "Model")) %>%
  mutate(label = str_c(mfd_sampletype, ", ", mfd_areatype, "<br>", mfd_hab1, ",<br>*n* = ", complex_size, sep = "",
                       "<br>", "__PERMANOVA:__ ", "<br>*R^2^* = ", 
                       round(R2, 2), "%, ", "*p* = ", pval)) %>%
  select(-R2, -pval) %>%
  cbind(df.permanova.total %>% filter(!is.na(pval))) %>%
  mutate(label_all = str_c("MFDO1 subset, ", "<br>*n* = ", nrow(.), "<br>", "__ANOSIM:__ ", "*R* = ", 
                           round(anosim$statistic, 2), ", ", "*p* = ", anosim$signif,
                           "<br>", "__PERMANOVA:__ ", "*R^2^* = ", 
                           round(R2, 2), "%, ", "*p* = ", pval)) %>%
  select(-Model, -R2, -pval) %>%
  mutate(across(label, ~str_replace(., "reclaimed lowland", "lowland")),
         across(label, ~factor(.)))

## Write object to output
data.table::fwrite(PCOA.scores, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_FL18S_pcoa_scores.tsv"))


### Visualise results
## Plot of all samples in the subset - colored by MFDO1
p.pcoa.1v2.all <- ggplot() +
  geom_point(data = PCOA.scores, aes(x = PCO1, y = PCO2, fill = complex), 
             size = 4, alpha = 1, color = "black", pch = 21) +
  theme_minimal(base_size = 19) +
  scale_fill_manual(values = mfdo1.palette.subset) +
  guides(fill = guide_legend(override.aes = list(shape = 22, size = 5, alpha = 1))) +
  theme(legend.position = "none",
        aspect.ratio = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_line(colour = "black"),
        strip.text.x = element_text(size = 16),
        plot.margin = margin(0,0,0,0),
        strip.text = ggtext::element_markdown()) +
  xlab(str_c("PCO1 - ", PCO1, "%", sep = "")) +
  ylab(str_c("PCO2 - ", PCO2, "%", sep = "")) +
  facet_wrap(vars(label_all), ncol = 1)

## Render plot
p.pcoa.1v2.all

## Plot of all samples in the subset - add contour to show density
p.pcoa.1v2.all.contour <- ggplot() +
  geom_point(data = PCOA.scores, aes(x = PCO1, y = PCO2, fill = complex),
             size = 5, alpha = 1, color = "black", pch = 21)  +
  # scale_x_continuous(limits = c(-0.625, 0.4)) +
  # scale_y_continuous(limits = c(-0.3, 0.525)) +
  scale_fill_manual(values = mfdo1.palette.subset) +
  ggnewscale::new_scale_fill() +
  stat_density_2d(data = PCOA.scores, aes(x = PCO1, y = PCO2, fill = after_stat(level)), 
                  geom = "polygon", colour = "black", adjust = 2.25, alpha = 0.2) +
  scale_fill_distiller(palette = "Greys", direction = 1) +
  theme_minimal(base_size = 19) +
  guides(fill = guide_legend(override.aes = list(shape = 22, size = 5, alpha = 1))) +
  theme(legend.position = "none",
        aspect.ratio = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_line(colour = "black"),
        strip.text.x = element_text(size = 16),
        plot.margin = margin(0,0,0,0),
        strip.text = ggtext::element_markdown()) +
  xlab(str_c("PCO1 - ", PCO1, "%", sep = "")) +
  ylab(str_c("PCO2 - ", PCO2, "%", sep = "")) +
  facet_wrap(vars(label_all), ncol = 1)

## Render plot
p.pcoa.1v2.all.contour

## Plot of all samples in the subset - faceted panel
p.pcoa.1v2.facet <- PCOA.scores %>%
  ggplot() +
  geom_point(data = PCOA.scores[-15], aes(x = PCO1, y = PCO2), 
             size = 2, alpha = 0.5, color = "black", fill = "white", pch = 21) +
  # scale_x_continuous(limits = c(-0.625, 0.4)) +
  # scale_y_continuous(limits = c(-0.3, 0.525)) +
  geom_point(data = PCOA.scores, aes(x = PCO1, y = PCO2, fill = complex),
             size = 2, alpha = 1, color = "black", pch = 21) +
  theme_minimal(base_size = 19) +
  scale_fill_manual(values = mfdo1.palette.subset) +
  guides(fill = guide_legend(override.aes = list(shape = 22, size = 5, alpha = 1))) +
  theme(legend.position = "none",
        aspect.ratio = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_line(colour = "black"),
        strip.text.x = element_text(size = 16),
        plot.margin = margin(0,0,0,0),
        strip.text = ggtext::element_markdown()) +
  xlab(str_c("PCO1 - ", PCO1, "%", sep = "")) +
  ylab(str_c("PCO2 - ", PCO2, "%", sep = "")) +
  facet_wrap(~factor(label, levels = rev(levels(label))), ncol = 5)

## Render plot
p.pcoa.1v2.facet


## Write arranged plots to output
png(file = 'output/ordination_left_panel.png',
    width = 1900,
    height = 1200) 
p.pcoa.1v2.all
dev.off()

pdf(file = 'output/ordination_left_panel.pdf',
    width = 19,
    height = 12) 
p.pcoa.1v2.all
dev.off()

tiff(file = 'output/ordination_left_panel.tiff',
     width = 1900,
     height = 1200) 
p.pcoa.1v2.all
dev.off()

ggsave("output/ordination_left_panel.svg", 
       plot = p.pcoa.1v2.all, width = 19, height = 12, 
       units = "in", dpi = "retina")

## Write arranged plots to output
png(file = 'output/ordination_righ_panel.png',
    width = 1900,
    height = 1200) 
p.pcoa.1v2.facet
dev.off()

pdf(file = 'output/ordination_right_panel.pdf',
    width = 19,
    height = 12) 
p.pcoa.1v2.facet
dev.off()

tiff(file = 'output/ordination_right_panel.tiff',
     width = 1900,
     height = 1200) 
p.pcoa.1v2.facet
dev.off()

ggsave("output/ordination_right_panel.svg", 
       plot = p.pcoa.1v2.facet, width = 19, height = 12, 
       units = "in", dpi = "retina")

### Save image
save.image(paste0('output/', format(Sys.time(), "%Y-%m-%d"), "_MFD_FL18S_ordinations.RData"))

