### Setup env
library(tidyverse)
library(vegan)
library(parallelDist)
library(ggtree)


set.seed(123)

options(stringsAsFactors = F, gsubfn.engine = "R")

Sys.setenv("LANGUAGE"="En")

source("scripts/MFD_colors.R")

setwd("/mfd_diversity_euk")

wd <- getwd()

data.path <- paste0(wd, '/data')
output.path <- paste0(wd, '/output')

### Load data
## MFD Ontology
mfd_ontology <- readxl::read_excel(paste0(data.path, "/2025-02-11_mfd-habitat-ontology.xlsx"), sheet = 1)

## Bray-Curtis matix
BC.df <- data.table::fread(paste0(output.path, "/2025-04-16_MFD_FL18S_BC_distance.csv")) %>%
  column_to_rownames(var = "V1")

## Set rownames for BC matrix
# BC.df <- BC.df %>%
#   `rownames<-`(colnames(BC.df))

## MFD db
mfd_db <- readxl::read_excel(paste0(data.path, "/2025-04-14_mfd_db.xlsx"), sheet = 1) %>%
  filter(fieldsample_barcode %in% colnames(BC.df))

## 10 km grid representatives
rep.10km.df <- data.table::fread(paste0(output.path, "/2025-04-16_MFD_FL18S_metadata_filtered.tsv")) %>%
  filter(fieldsample_barcode %in% colnames(BC.df))

## Heatmap data frame
heatmap.df <- data.table::fread(paste0(output.path, "/2025-04-16_MFD_FL18S_heatmap_MFDO1.tsv"))

## Alpha diversity data frame
diversity.df <- data.table::fread(paste0(output.path, "/2025-04-16_MFD_FL18S_alpha_diversity_MFDO1.tsv"))

signif <- data.table::fread(paste0(output.path, "/2025-04-23_MFD_FL18S_alpha_diveristy_analysis.csv"))

## Gamma diversity data frame
gamma.df <- data.table::fread(paste0(output.path, "/2025-04-16_MFD_FL18S_gamma_diversity_MFDO1.tsv"))

## Genus aggregated observational table
reduced.otu <-  data.table::fread(paste0(output.path, "/2025-04-16_MFD_FL18S_OTU_filtered_ra_pa.csv")) %>%
  select(OTU, any_of(rep.10km.df$fieldsample_barcode), Kingdom:Genus) %>%
  column_to_rownames(var = "OTU") %>%
  filter(rowSums(across(where(is.numeric)))!=0)

## Perform Hellinger transformation
reduced.df <- reduced.otu %>%
  select(where(is.numeric)) %>%
  t() %>%
  as.data.frame()


### Habitat dissimilarities
## Habitat parsing functions
mfd_db_ontology <- mfd_db %>%
  select(fieldsample_barcode, mfd_sampletype, mfd_areatype, mfd_hab1, mfd_hab2, mfd_hab3) %>%
  mutate(`sample type cumulative` = mfd_sampletype,
         `area type cumulative` = paste0(mfd_sampletype, ", ", mfd_areatype),
         `MFDO1 cumulative` = paste0(mfd_sampletype, ", ", mfd_areatype, ", ", mfd_hab1),
         `MFDO2 cumulative` = paste0(mfd_sampletype, ", ", mfd_areatype, ", ", mfd_hab1, ", ", mfd_hab2),
         `MFDO3 cumulative` = paste0(mfd_sampletype, ", ", mfd_areatype, ", ", mfd_hab1, ", ", mfd_hab2, ", ", mfd_hab3))

rep_10k_ontology <- rep.10km.df %>%
  select(fieldsample_barcode, mfd_sampletype, mfd_areatype, mfd_hab1, mfd_hab2, mfd_hab3) %>%
  mutate(`sample type cumulative` = mfd_sampletype,
         `area type cumulative` = paste0(mfd_sampletype, ", ", mfd_areatype),
         `MFDO1 cumulative` = paste0(mfd_sampletype, ", ", mfd_areatype, ", ", mfd_hab1),
         `MFDO2 cumulative` = paste0(mfd_sampletype, ", ", mfd_areatype, ", ", mfd_hab1, ", ", mfd_hab2),
         `MFDO3 cumulative` = paste0(mfd_sampletype, ", ", mfd_areatype, ", ", mfd_hab1, ", ", mfd_hab2, ", ", mfd_hab3))

get_av_dist <- function(mat, lab.df){
  labs <- unique(lab.df$lab)
  mat.to_return <- matrix(-1, nrow = length(labs), ncol = length(labs)) %>%
    `colnames<-`(labs) %>%
    `rownames<-`(labs)
  for(i in labs){
    samples.i <- lab.df %>% filter(lab == i) %>% pull(fieldsample_barcode)
    for(j in labs){
      samples.j <- lab.df %>% filter(lab == j) %>% pull(fieldsample_barcode)
      #print(mat[samples.i, samples.j])
      mat.to_return[i, j] <- mean(as.matrix(mat[samples.i, samples.j]))
    }
  }
  return(mat.to_return)
}

## Dissimilarity functions
BC.wrap <- function(counts, lab){
  mat <- parallelDist(counts,
                      method = "binary",
                      binary = F,
                      diag = T,
                      upper = T,
                      threads = 50) %>%
    as.matrix()
  mat.to_return <- get_av_dist(mat, lab) %>%
    as.dist()
  return(mat.to_return)
}

### Hab profiles
lab.df <- rep_10k_ontology %>%
  mutate(lab = `MFDO1 cumulative`) %>%
  mutate(across(lab,
                ~str_replace(.,
                             "Soil, Natural, Sclerophyllous scrub",
                             "Soil, Natural, Temperate heath and scrub"))) %>%
  select(fieldsample_barcode, lab)

## MFDO1 summary (mean dissimilarity)
mfdo1.profiles.hell <- reduced.df %>%
  as.data.frame() %>%
  rownames_to_column("fieldsample_barcode") %>%
  left_join(lab.df,
            by = "fieldsample_barcode") %>%
  filter(!is.na(lab)) %>%
  dplyr::select(-fieldsample_barcode) %>%
  pivot_longer(names_to = "Genera", values_to = "value", -lab) %>%
  group_by(lab, Genera) %>%
  summarise(value.mean = mean(value)) %>%
  pivot_wider(names_from = "Genera", values_from = "value.mean") %>%
  column_to_rownames("lab") %>%
  decostand(., method = "hellinger") %>%
  as.data.frame()

## Bootstrap SBHC large
ab.df <- rep_10k_ontology %>%
  mutate(lab = `MFDO1 cumulative`) %>%
  mutate(across(lab,
                ~str_replace(.,
                             "Soil, Natural, Sclerophyllous scrub",
                             "Soil, Natural, Temperate heath and scrub"))) %>%
  select(fieldsample_barcode, lab)

## hclust with BC-wrapper
createHclustObject <- function(x)hclust(BC.wrap(counts = as.matrix(x),
                                                lab = (lab.df %>%
                                                         filter(fieldsample_barcode%in%rownames(reduced.df)))), "ave")

## bootstrap
set.seed(123)
b <- bootstrap::bootstrap(reduced.df,
                          fun = createHclustObject,
                          n = 100L,
                          mc.cores = 60)

## Write to output
data.table::fwrite(as.data.frame(b), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0(output.path, format(Sys.time(), "%Y-%m-%d"), "_MFD_FL18S_bootstrap_weights_MFDO1.csv"))


## Plot results
hc <- createHclustObject(reduced.df)
plot(hc)

## Draw bootstrap values to corresponding node
bootstrap::bootlabels.hclust(hc, b, col="blue")


### Hab summarisation
## Average distance for full dataset
mfdo1.dist.full <- get_av_dist((BC.df %>% as.data.frame() %>% `rownames<-`(colnames(.))), (mfd_db_ontology %>%
                                                                                             mutate(lab = `MFDO1 cumulative`) %>%
                                                                                             mutate(across(lab,
                                                                                                           ~str_replace(.,
                                                                                                                        "Soil, Natural, Sclerophyllous scrub",
                                                                                                                        "Soil, Natural, Temperate heath and scrub"))) %>%
                                                                                             select(fieldsample_barcode, lab)))

## Average distance for 10 km reference subset
mfdo1.dist.rep <- get_av_dist((BC.df %>% as.data.frame() %>% `rownames<-`(colnames(.))), (rep_10k_ontology %>%
                                                                                            mutate(lab = `MFDO1 cumulative`) %>%
                                                                                            mutate(across(lab,
                                                                                                          ~str_replace(.,
                                                                                                                       "Soil, Natural, Sclerophyllous scrub",
                                                                                                                       "Soil, Natural, Temperate heath and scrub"))) %>%
                                                                                            select(fieldsample_barcode, lab)))


### Visualisation - Beta and fingerprint using geom_facet

## Transform phylum level 18S metagenome profiles
heatmap.long <- heatmap.df

# phylum.long <- (phylum.df %>%
#   mutate(hab = complex,
#          rel.ab = abund))[c("hab", "rel.ab", "Phylum")]

## Sample count summary of 18S metagenome profiles
dd.MG.n <- rep_10k_ontology %>%
  mutate(lab = `MFDO1 cumulative`) %>%
  mutate(across(lab,
                ~str_replace(.,
                             "Soil, Natural, Sclerophyllous scrub",
                             "Soil, Natural, Temperate heath and scrub"))) %>%
  group_by(`MFDO1 cumulative`) %>%
  summarise(n=n())

## Dendogram of 18S metagenome profiles
p.rep <- hc %>%
  ggtree(ladderize = T) +
  theme_tree2() +
  scale_x_continuous(labels = abs)

b.df <- data.frame(x = -(hc$height),
                   bootstrap = b) %>%
  left_join((p.rep$data %>% select(x, node) %>% as.data.frame()),
            by = "x")

p.colors.rep <- p.rep %<+% data.frame(label = colnames(mfdo1.dist.rep),
                                      beta.within = diag(mfdo1.dist.rep)) +
  geom_tiplab(aes(label = label), align = T, linesize = .5, as_ylab = T) +
  geom_tippoint(aes(fill = label, size = beta.within), shape = 21) +
  scale_fill_manual(values = mfdo1.palette) +
  scale_size(range = c(3,6)) +
  guides(fill = "none") +
  ggnewscale::new_scale_fill()

p.bootstrap.rep <- p.colors.rep %<+% b.df +
  geom_point(data = (. %>% filter(!isTip)), aes(fill=bootstrap), shape = 23, size = 3, color = "black") +
  scale_fill_distiller(type = "seq", palette = "Greys", direction = 1,  na.value = NA) +
  scale_fill_stepsn(colors = rev(grey.colors(8))[c(1, 8)], na.value = NA,
                    limits = c(0.9, 1),
                    breaks = c(0.95, 1),
                    labels = c("<0.95", ">0.95"),
                    name = "Bootstrap") +
  guides(fill = guide_legend(override.aes = list(shape = 23))) +
  scale_color_continuous(na.value = NA) +
  ggnewscale::new_scale_fill() +
  ggnewscale::new_scale_color()

color_vector <- c("#440154FF",  # 0-24% mapped to dark purple
                  "#2A788EFF",  # 25-49% mapped to blue-green
                  "#22A884FF",  # 50-74% mapped to turquoise
                  "#FDE725FF")  # 75-100% mapped to yellow-green

## Phylum level 18S metagenome profiles
p.fingerprint <- p.bootstrap.rep +
  geom_facet(panel = "Fingerprint",
             data = heatmap.long,
             geom = geom_tile,
             aes(x = Display, fill = category),
             color = "black",
             lwd = .25,
             linetype = 1) +
  scale_x_discrete() +
  scale_fill_manual(
    values = c(
      "0-24%" = color_vector[1],   
      "25-49%" = color_vector[2],  
      "50-74%" = color_vector[3],  
      "75-100%" = color_vector[4]   
    ),
    breaks = c("0-24%", "25-49%", "50-74%", "75-100%"),
    labels = c("0-24%", "25-49%", "50-74%", "75-100%")
  ) +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1)) +
  guides(size = "none") +
  theme(legend.position="bottom")

## Sample count summary of 18S metagenome profiles
p.n.MG <- p.bootstrap.rep

## Alpha diversity based on FL18S amplicons
p.a_div <- p.n.MG +
  ggnewscale::new_scale_fill() +
  geom_facet(panel = "Alpha diversity",
             data = diversity.df[,c("complex", "observed_species")],
             geom = geom_jitter,
             # orientation = "y",
             aes(x = observed_species, color = label),
             size = 1.0) +
  scale_fill_manual(values = mfdo1.palette) +
  scale_color_manual(values = mfdo1.palette) +
  geom_facet(panel = "Alpha diversity",
             data = diversity.df[,c("complex", "observed_species")],
             geom = geom_boxplot,
             orientation = "y",
             aes(x = observed_species, fill = label),
             color = "black",
             width = 0.8,
             alpha = 0.6,
             outlier.size = -1) +
  scale_fill_manual(values = mfdo1.palette) +
  guides(fill = "none",
         color = "none") +
  geom_facet(panel = "Alpha diversity",
             data = signif,
             geom = geom_text,
             aes(x = -100, label = cld)) +
  guides(color = "none")

## Gamma diversity based on FL18S amplicons
p.c_div <- p.a_div +
  geom_facet(pane = "Gamma diversity",
             data = (gamma.df %>% filter(Hill_diversity == "Shannon diversity")),
             geom = geom_errorbar,
             aes(x= Estimator, xmin = LCL, xmax = UCL),
             color = "black",
             width = 0.5,
             orientation = "y") +
  geom_facet(panel = "Gamma diversity",
             data = (gamma.df %>% filter(Hill_diversity == "Shannon diversity")),
             geom = geom_col,
             aes(x = Estimator, fill = label),
             # pch = 21,
             color = "black",
             width = 1,
             alpha = 0.8,
             orientation = "y")

## Sample count summary of FL18S amplicons
p.n.18S <- p.c_div +
  geom_facet(panel = "# FL18S",
             data = (gamma.df %>% filter(Hill_diversity == "Shannon diversity")),
             geom = geom_text,
             aes(x = 1, label = paste0("n = ", samples))) +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1))

## Render and save plots in svg format
p.n.18S

ggsave(p.n.18S, filename = paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_18S_tree_cont_scale.svg"), 
       width = 12, height = 4, dpi = "retina")

p.fingerprint

ggsave(p.fingerprint, filename = paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_18S_tree_disc_scale.svg"), 
       width = 8, height = 6, dpi = "retina")

save.image(paste0('output/', format(Sys.time(), "%Y-%m-%d"), "_MFD_FL18S_combined_diversity.RData"))

