### Setup env
library(tidyverse)
library(ampvis2)

source('scripts/MFD_colors.R')

setwd("/mfd_diversity_euk")


### Import data
filter <- data.table::fread("output/2025-04-16_MFD_FL18_combined_diversity_MFDO1.tsv") %>%
  select(complex) %>%
  distinct() %>%
  pull(complex)

## Filter and create "complex" correponding to full MFDO1 string for 10 km reference grid samples
metadata <- data.table::fread('output/2025-04-16_MFD_FL18S_metadata_filtered.tsv', na.strings = "") %>%
  mutate(across(mfd_hab1, ~str_replace(., "Sclerophyllous scrub", "Temperate heath and scrub")),
         across(mfd_hab2, ~str_replace(., "Scrub", "Sclerophyllous scrub"))) %>%
  filter(complex %in% filter) %>%
  mutate(across(complex, ~factor(., levels = names(mfdo1.palette)))) %>%
  select(fieldsample_barcode, everything())

## Create group summary
groups.summary <- metadata %>%
  group_by(complex) %>%
  summarise(size = n())

## Pull sample IDs
samples.subset <- metadata %>%
  pull(fieldsample_barcode)

## Import count data of the genus-aggregated observational table of 16S dervied from metagenomes
data <- data.table::fread('output/2025-04-16_MFD_FL18S_OTU_filtered_ra_pa.csv', na.strings = "")

## Subset data to 10 km reference samples
data.subset <- data %>%
  select(OTU, all_of(samples.subset), Kingdom:Genus) %>%
  filter(rowSums(across(where(is.numeric)))!=0)

rm(data)

## Create ampvis2 object with data subset
ampvis <- amp_load(otutable = data.subset,
                   metadata = metadata)

genus_all_bin <- amp_heatmap(
  ampvis,
  group_by = "fieldsample_barcode",
  tax_show = 667,
  normalise = FALSE,
  plot_values = TRUE,
  color_vector = c("white", "red"),
  tax_aggregate = "Genus"
)

otutable_binary <- genus_all_bin$data$Abundance
otutable_binary[otutable_binary > 0] <- 1

genus_all_bin$data$Abundance <- otutable_binary

##filter after Presence

#Count the number of MFD01 groups with presence for each taxonomic group (Display)
presence_totals <- genus_all_bin$data %>%
  group_by(Display) %>%  # Group by taxonomic group (Display)
  summarise(total_presence = sum(Abundance)) %>%  # Count presence (sum of 1's) across MFD01 groups
  arrange(desc(total_presence))  # Sort by most present

# Filter out placeholder names
filtered_presence_totals <- presence_totals %>%
  filter(!str_detect(Display, "Unclassified|_X$|_XX$|_XXX$|ASV"))  

# Subset to keep only the top 30 most present genera
top_30 <- filtered_presence_totals %>%
  slice_head(n = 30)  

genus_all_bin$data <- genus_all_bin$data %>%
  filter(Display %in% top_30$Display) 

#Reorder the Display factor based on total presence
genus_all_bin$data$Display <- factor(genus_all_bin$data$Display, levels = top_30$Display)

long_data_gen <- genus_all_bin$data %>%
  select(Display, Sample, Abundance)


# Merge long_data with metadata on Sample column
long_data_with_metadata_gen <- left_join(long_data_gen, metadata, by = c("Sample" = "fieldsample_barcode")) %>%
  select(Display, Sample, Abundance,complex)

# Calculate the average Abundance per taxon in each Habitat
avg_data_gen <- long_data_with_metadata_gen %>%
  group_by(complex, Display) %>%
  summarise(average_presence = mean(Abundance)) %>%
  ungroup()

# Change the x-axis names
avg_data_gen <- long_data_with_metadata_gen %>%
  mutate(Display = recode(Display,
                          "Mortierella" = "Mortierella (Fungi)",
                          "Cryptococcus" = "Cryptococcus (Fungi)",
                          "Pseudogymnoascus" = "Pseudogymnoascus (Fungi)",
                          "Penicillium" = "Penicillium (Fungi)",
                          "Monocystis" = "Monocystis (Apicomplexa)",
                          "Gibberella" = "Gibberella (Fungi)",
                          "Paracercomonas" = "Paracercomonas (Cercozoa)",
                          "Leidyana1" = "Leidyana1 (Apicomplexa)",
                          "Cercomonas" = "Cercomonas (Cercozoa)",
                          "Rhizoscyphus" = "Rhizoscyphus (Fungi)",
                          "Trichoderma" = "Trichoderma (Fungi)",
                          "Plectosphaerella" = "Plectosphaerella (Fungi)",
                          "Arcuospathidium" = "Arcuospathidium (Ciliophora)",
                          "Villosiclava" = "Villosiclava (Fungi)",
                          "Spumella" = "Spumella (Heterokontophyta)",
                          "Eocercomonas" = "Eocercomonas (Cercozoa)",
                          "Spongospora" = "Spongospora (Cercozoa)",
                          "Asterotremella" = "Asterotremella (Fungi)",
                          "Chaetomium" = "Chaetomium (Fungi)",
                          "Paramicrosporidium" = "Paramicrosporidium (Fungi)",
                          "Davidiella" = "Davidiella (Fungi)",
                          "Polymyxa" = "Polymyxa (Cercozoa)",
                          "Exophiala" = "Exophiala (Fungi)",
                          "Syncystis" = "Syncystis (Apicomplexa)",
                          "Neocercomonas" = "Neocercomonas (Cercozoa)",
                          "Flamella" = "Flamella (Evosea)",
                          "Acremonium" = "Acremonium (Fungi)",
                          "Pleosporales" = "Pleosporales (Fungi)",
                          "Group-Te" = "Group-Te (Cercozoa)",
                          "Sordaria" = "Sordaria (Fungi)"  )) %>%
  group_by(complex, Display) %>%
  summarise(average_presence = mean(Abundance), .groups = "drop") %>%
  ungroup()

# Ensure MFD01 is ordered according to mfd_order
avg_data_gen$complex <- factor(avg_data_gen$complex, levels = mfdo1.levels)

# Plot the heatmap
avg_bin_gen <- ggplot(avg_data_gen, aes(x = Display, y = complex, fill = average_presence)) +
  geom_tile(color = "black", width = 1, height = 1) +  # Black borders around tiles
  scale_fill_gradient(low = "white", high = "red", na.value = "white", 
                      labels = scales::percent) +  # Format the legend as percentage
  theme_bw() +  # Use black-and-white theme
  labs(
    title = "Average Taxa Presence/Absence by MFD01 Group",
    x = "Taxonomic Group (Display)",  # x-axis is Display
    y = "MFD01 Group",  # y-axis is MFD01
    fill = "Taxa Presence (%)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    panel.border = element_blank(),  # Remove panel border
    panel.background = element_blank(),  # Remove panel background
    axis.line = element_line(colour = "black"),  # Add black line around the plot
    legend.position = "left",  # Place the legend on the left
    legend.title = element_text(size = 12),  # Customize legend title size
    legend.text = element_text(size = 10)  # Customize legend text size
  )

# Define a custom color vector for four categories
color_vector <- c("#440154FF",  # 0-24% mapped to dark purple
                  "#2A788EFF",  # 25-49% mapped to blue-green
                  "#22A884FF",  # 50-74% mapped to turquoise
                  "#FDE725FF")  # 75-100% mapped to yellow-green

# Modify `mutate()` to ensure proper category assignment
avg_data_gen <- avg_data_gen %>%
  mutate(category = case_when(
    average_presence >= 0.75 ~ "75-100%",   # 75% to 100%
    average_presence >= 0.50 ~ "50-74%",    # 50% to 74.99%
    average_presence >= 0.25 ~ "25-49%",    # 25% to 49.99%
    TRUE ~ "0-24%"                          # 0% to 24.99%
  ))


print(unique(avg_data_gen$category))

# Plotting using the corrected category labels
avg_bin_gen <- ggplot(avg_data_gen, aes(x = Display, y = complex, fill = category)) +
  geom_tile(color = "black") +  # Black borders around tiles
  
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
  
  theme_bw() +  # Use black-and-white theme
  
  labs(
    title = "Average Taxa Presence/Absence by MFD01 Group",
    x = "Genus (Division)",  
    y = "MFD01 Group",         
    fill = "Taxa Presence"
  ) +
  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    panel.border = element_blank(),      
    panel.background = element_blank(),  
    axis.line = element_line(colour = "black"),  
    legend.position = "left",  
    legend.title = element_text(size = 12),  
    legend.text = element_text(size = 10)
  )


# Display the plot
avg_bin_gen


## Write table to output
data.table::fwrite(avg_data_gen, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_FL18S_heatmap_MFDO1.tsv"))


### Save iamge
save.image(paste0('output/', format(Sys.time(), "%Y-%m-%d"), "_MFD_FL18S_heatmap.RData"))

