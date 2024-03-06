## Load packages
library(ggplot2)
library(phyloseq)
library(dplyr)
library(vegan)
library(tidyr)
library(reshape2)
library(rbiom)
library(devtools)
library(ranacapa)

getwd()

#Import data 
asv_table <- read.table("raw_feature_table.txt", sep = '\t', row.names = 1, header = T, strip.white = T)
asv.mat <- as.matrix(asv_table)

taxonomy <- read.table("taxonomy_gg.tsv", sep = '\t', row.names = 1, header = T, strip.white = T)
taxonomy$Confidence <- NULL

metadata <- read.table("metadata_daynight.txt", sep = '\t', row.names = 1, header = T, strip.white = T) 

#Separate taxonomy into different columns #you can ignore the error. It just mean that the empty value is treated as NA
tax_sep <- separate(taxonomy, Taxon, c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                    sep = ";", remove = TRUE, convert = FALSE, extra = "warn", fill = "warn")

#There will be warning message but do not worry about it

# Replace all NA's with the taxa indicator only
tax_sep$Phylum[is.na(tax_sep$Phylum)] <- "p__unknown"
tax_sep$Class[is.na(tax_sep$Class)] <- "c__unknown"
tax_sep$Order[is.na(tax_sep$Order)] <- "o__unknown"
tax_sep$Family[is.na(tax_sep$Family)] <- "f__unknown"
tax_sep$Genus[is.na(tax_sep$Genus)] <- "g__unknown"
tax_sep$Species[is.na(tax_sep$Species)] <- "s__unknown"

# 'Phyloseq-ize' the data
asv.tab <- otu_table(as.matrix(asv_table), taxa_are_rows = T)   #this is ASV sequence, but since Phyloseq was initially designed for Otu, hence why they named the data as OTU. But, it does not matter.

tax.tab <- tax_table(as.matrix(tax_sep))

sample.data <- sample_data(metadata)

physeq <- phyloseq(asv.tab, tax.tab, sample.data)
physeq

sample_data(physeq)$day.time[sample_data(physeq)$day.time == "1_0915"] #just a test

#Make some data as factor
sample_data(physeq)$day.time = factor(sample_data(physeq)$day.time,
                                          levels = c("1_0915", "1_1230", "1_1612", "1_1740"))

sample_data(physeq)$sample.name = factor(sample_data(physeq)$sample.name)

#Make Rarefaction curves
rarefaction_plot <- ggrare(physeq, step = 500, color = "sample.name", label = NULL, se = FALSE) # The number of rarefaction depths to include between mindepth and max-depth

rarefaction_plot_sample <- rarefaction_plot + facet_wrap(~sample.name, ncol =5, nrow = 6) + 
  guides(colour = guide_legend(title = "Sample")) + 
  scale_x_continuous(breaks = c(5000, 25000, 50000, 75000, 100000)) +
  xlab("Sampling depth") +  
  ylab("ASV richness") +
  theme_bw()
rarefaction_plot_sample

#ggsave("rarefaction.pdf", width = 15, height = 10, dpi = "print")


