# mfd_diversity_euk
The scripts in this repository are part of the [Microflora Danica project](https://github.com/cmc-aau/mfd_wiki/wiki). 

Be advised, that for the metagenomic-derived data, the term "OTU" is only used due to format requirements by ampvis2, and they do not represent classical OTUs. 
The generated profiles can be thought of as taxonomic bins. 

All scritps are needed to make the combined figure panel at the end. 

## Scripts
### Operon amplicon 18S data 
`scripts/format_FL18s.R` formats an OTU table from the 18S rRNA gene data extracted from the eukaryotic operons. 


`scripts/richness_estimates_FL18S.R` estimates the total richness, using the near full-length 18S OTUs from the PacBio data as well as Shannon and Simpson diversity - [not to confused with their corresponding indexes](https://johnsonhsieh.github.io/iNEXT/). Estimates are made for the total data and categoty-specific estimates based on MFDO1. 


`scripts/diversity_FL18S.R` calculates alpha diversity using the near full-length 18S OTUs from the PacBio data across the MFDO1 ontology level. It performs the statistical analysis on alpha diveristy and combines the alpha with the gamma diversity calculated for the same data.


`scripts/heatmap_FL18S.R` recreates the eukaryotic fingerprint profile across the MFDO1 ontology level. 


`scripts/permanova_contrasts_FL18S.R` performs the PERMANOVA and contrasts abalysis across the MFDO1 ontology level.


`scripts/ordination_FL18S.R` calculates the Jaccard dissimlarity matrix and recreates the eukaryotic ordination panels with addition of results from ANOSIM and PERMANOVA. 


`scripts/combine_diversity_FL18S.R` combines the different aspects of diveristy and calcualtes the between-group Jaccard dissimilarites across the MFDO1 ontology level. The script also renders the figures used in the manuscript. 


## Data
The scripts rely on data files available from the MFD [github](https://github.com/cmc-aau/mfd_metadata) and the MFD Zenodo [repo](https://zenodo.org/records/12605769). 
