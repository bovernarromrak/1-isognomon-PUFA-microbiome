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
library(gg3D) #gg3D is a package created to extend ggplot2 to produce 3D plots


#import data (important to import with read.table because the first column can be set as a vector - you can try with read.delim() to import and you will notice that at later section you can't work with this format): 
##Sample-grouped rarefied 53527###################################################################################################################################################################################################################################################################################################

asv <- read.delim("table_daytime_samplegroup_rarefied53527.txt", sep = '\t', row.names = 1, header = T, strip.white = T)
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
metabolites_significantfour <- metabolites_all[, c("X5.F2t.IsoP", "X10.F4.NeuroP", "X13.F4.NeuroP", "X16.F1.PhytoP")]

#transform the data by autoscaling following the suggestion by van den Berg et al.(2006)
metabolites_significantfour_autoscaled <- data.frame(scale(metabolites_significantfour[,1:4], center = TRUE, scale = TRUE)) 

#Check variance
apply(metabolites_significantfour_autoscaled[,1:4], 2, var) #the variance is the same now

#Output data: Print out the values in excel for averaging from Excel
write.csv(metabolites_significantfour_autoscaled, file = "metabolites_significantfour_autoscaled.csv") #as a backup and to compare the result below

#Get the average of each timepoint and each metabolite
D1_0915 <- data.frame(colMeans(metabolites_significantfour_autoscaled[1:5,]))
colnames(D1_0915) <- "D1_0915"

D1_1230<- data.frame(colMeans(metabolites_significantfour_autoscaled[6:9,]))
colnames(D1_1230) <- "D1_1230"

D1_1612 <- data.frame(colMeans(metabolites_significantfour_autoscaled[10:14,]))
colnames(D1_1612) <- "D1_1612"

D1_1740 <- data.frame(colMeans(metabolites_significantfour_autoscaled[15:18,]))
colnames(D1_1740) <- "D1_1740"

#combine significant four the data and transform it 
metabolites_significantfour_ave<- t(cbind(D1_0915, D1_1230, D1_1612, D1_1740))

# Check matching row names
is.element(row.names(asv_t), row.names(metabolites_significantfour_ave)) #Perfect

# align rows between dataframes - if the two data frame is not the same
asv_t <- asv_t[match(row.names(metabolites_significantfour_ave), row.names(asv_t)),]
apply(asv_t, 1, sum)

#1. Getting the distance matrix - ####################################################################################################################################################################################################################################################
#Note: 1. always check if the column is right! it could be not. 2. The result always changes because of the random permutation and so the significance return will always be different. 
#Samples distance matrix - use bray-curtis

#Check again the sum
apply(asv_t, 1, sum)

sample_dist <- vegdist(asv_t[,1:4301], method = "bray")

#Metabolites distance matrix 
metabolites_significantfour_dist <- vegdist(metabolites_significantfour_ave[,1:4], method = "euclidean")


#2. Mantel Test ##################################################################################################################################################
#non-enzymatic products
metabolites_significantfour_spearman <- mantel(sample_dist, metabolites_significantfour_dist , method = "spearman", permutations = 9999, na.rm = TRUE)
metabolites_significantfour_spearman

