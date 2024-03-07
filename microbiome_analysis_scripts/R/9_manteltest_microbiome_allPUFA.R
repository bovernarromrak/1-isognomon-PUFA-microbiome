# Load packages

library(vegan)
library(phyloseq)
library(dplyr)
library(phangorn)
library(textshape)
library(microbiome)
library(ape)
library(geosphere)
library(tidyr)
library(scales)
library(dplyr)
library(ggplot2)
library(gg3D) 

##Sample-grouped rarefied 53527###################################################################################################################################################################################################################################################################################################

asv <- read.delim("table_daytime_samplegroup_rarefied53527.txt", sep = '\t', row.names = 1, header = T, strip.white = T) #grouped by timepoints and rarefied within the QIIME2 platform
asv_t <- t(asv)

# Edit row names
row.names(asv_t) <- gsub("X1_0915", "D1_0915", row.names(asv_t))
row.names(asv_t) <- gsub("X1_1230", "D1_1230", row.names(asv_t))
row.names(asv_t) <- gsub("X1_1612", "D1_1612", row.names(asv_t))
row.names(asv_t) <- gsub("X1_1740", "D1_1740", row.names(asv_t))

#verify if all samples have equal sampling depth
apply(asv_t, 1, sum) #yes

#upload metabolite data
metabolites_all <- read.delim("metabolite_products_rawdata.txt", sep = '\t',row.names = 1 , header = T, strip.white = T) #this version follows the data structure in the file: "asv_table_wide_rarefied5612"

#transform the data by autoscaling following the suggestion by van den Berg et al.(2006): see the details description in "Things_Iamnot_veryconfident_inmyPhDthesisyet.docx"
metabolites_all_autoscaled <- data.frame(scale(metabolites_all[,1:15], center = TRUE, scale = TRUE)) 

#Check variance
apply(metabolites_all_autoscaled[,1:15], 2, var) #the variance is the same now

#Get the average of each timepoint and each metabolite
D1_0915 <- data.frame(colMeans(metabolites_all_autoscaled[1:5,]))
colnames(D1_0915) <- "D1_0915"

D1_1230<- data.frame(colMeans(metabolites_all_autoscaled[6:9,]))
colnames(D1_1230) <- "D1_1230"

D1_1612 <- data.frame(colMeans(metabolites_all_autoscaled[10:14,]))
colnames(D1_1612) <- "D1_1612"

D1_1740 <- data.frame(colMeans(metabolites_all_autoscaled[15:18,]))
colnames(D1_1740) <- "D1_1740"

#combine all the data and transform it 
metabolites_all_ave<- t(cbind(D1_0915, D1_1230, D1_1612, D1_1740))

# Check matching row names
is.element(row.names(asv_t), row.names(metabolites_all_ave)) #Perfect

# align rows between dataframes - if the two data frame is not the same
asv_t <- asv_t[match(row.names(metabolites_all_ave), row.names(asv_t)),]
apply(asv_t, 1, sum)

#1. Getting the distance matrix - ####################################################################################################################################################################################################################################################

#microbiome distance matrix
#Check again the sum
apply(asv_t, 1, sum)

sample_dist <- vegdist(asv_t[,1:4301], method = "bray")

#non-enzymatic products
metabolites_all_ave_dist <- vegdist(metabolites_all_ave[,1:15], method = "euclidean")

#2. Mantel Test ##################################################################################################################################################
#non-enzymatic products
metabolites_all_ave.spearman <- mantel(sample_dist, metabolites_all_ave_dist , method = "spearman", permutations = 9999, na.rm = TRUE)
metabolites_all_ave.spearman

