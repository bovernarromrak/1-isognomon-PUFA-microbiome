#Loading library

library(reshape2)
library(vegan)
library(phyloseq)
library(microbiome)
library(tidyr)
library(ggplot2)
library(scales)
library(dplyr)
library(forcats)
library(gridExtra)
library(RColorBrewer)
library(randomcoloR)

##############################################################################################################################
#####Taxonomic profiles of the 10 most abundant ASV ##########################################################################
##############################################################################################################################

###################################################################################################################################################################################################################################################################################
#import data 
asv <- read.table("raw_feature_table.txt", sep = '\t', row.names = 1, header = T, strip.white = T)    
asv <- asv[,colnames(asv) != "BRM3" & colnames(asv) != "DRM3"]

map <- read.table("metadata_daynight.txt", sep = '\t', row.names = 1, header = T, strip.white = T)

tax <- read.table("taxonomy_gg.tsv", sep = '\t', row.names = 1, header = T, strip.white = T)
# Separate taxonomy into different columns #you can ignore the error. It just mean that the empty value is treated as NA
tax <- separate(tax, Taxon, c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                sep = ";" , remove = TRUE, convert = FALSE, extra = "warn", fill = "warn")
tax$Confidence <- NULL

######################################################################################################################################d#######################################################################################################################################d#######################################################################################################################################
#Start with renaming first - Replace all NA's with the taxa indicator only or like in the kingdom case - change the name
#Always check with the unique names first for each column and check again at the end

tax$Kingdom[is.na(tax$Kingdom)] <- " k__unknown" 
tax$Kingdom[tax$Kingdom == " k__"] <- " k__unknown"
tax$Phylum[is.na(tax$Phylum)] <- " p__unknown"
tax$Phylum[tax$Phylum == " p__"] <- " p__unknown"
tax$Class[is.na(tax$Class)] <- " c__unknown"
tax$Class[tax$Class == " c__"] <- " c__unknown"
tax$Order[is.na(tax$Order)] <- " o__unknown"
tax$Order[tax$Order == " o__"] <- " o__unknown"
tax$Family[is.na(tax$Family)] <- " f__unknown"
tax$Family[tax$Family == " f__"] <- " f__unknown"
tax$Genus[is.na(tax$Genus)] <- " g__unknown"
tax$Genus[tax$Genus == " g__"] <- " g__unknown"
tax$Species[is.na(tax$Species)] <- " s__unknown"
tax$Species[tax$Genus == " s__"] <- " s__unknown"
unique(tax$Kingdom)
unique(tax$Phylum)
unique(tax$Class)
unique(tax$Order)
unique(tax$Family)
unique(tax$Genus)
unique(tax$Species)


# 'Phyloseq-ize' the data - following the Savary et al. (2021)
otu.t= otu_table(asv, taxa_are_rows=TRUE)
sam.t= sample_data(data.frame(map))
tax.t= tax_table(as.matrix(tax))

phy.all= phyloseq(otu.t, tax.t,  sam.t)
phy.all

#Group up to Family level
Family_data <- aggregate_taxa(phy.all, "Family") 
Family_level_abundance <- as.data.frame(otu_table(Family_data))
Family_level_abundance_transpose <- t(Family_level_abundance)

total <- colSums(Family_level_abundance_transpose)
total
Family_level_abundance_total <- rbind(Family_level_abundance_transpose, total)

# re-order columns by total
Family_level_abund_sort <- Family_level_abundance_total [,order(-Family_level_abundance_total[which(rownames(Family_level_abundance_total) == 'total'),]) ] #Here essentially it is re-ordering the total by descending order

## Convert counts to percentages on data matrix & remove row total
Family_level_abund_percent <- Family_level_abund_sort[-which(rownames(Family_level_abund_sort) == 'total'),] / rowSums(Family_level_abund_sort[-which(rownames(Family_level_abund_sort) == 'total'),]) * 100

# create new column "other" which is the sum of all the taxa not in top 30, remove other taxa
Other <- rowSums(Family_level_abund_percent[,-c(1:10)])

# Combine data
Family_top10 <- cbind(Family_level_abund_percent[,1:10], Other)
Family_rel_abun_10 <- cbind(Family_top10, map)

#Get the output for the plotting in the PRISM
getwd()

#write.table(Family_rel_abun_10,  "Family_levelTop10percent_results.txt", 
#            sep = "\t", quote = F, row.names = T ) #Create a new file
