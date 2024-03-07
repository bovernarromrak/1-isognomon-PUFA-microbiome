## Load packages

library(permute)
library(dplyr)
library(lattice)
library(vegan)
library(cluster)
library(pairwiseAdonis)
library(ape)
library(ggplot2)
library(gg3D) #gg3D is a package created to extend ggplot2 to produce 3D plots

##rarefied data- 21282###################################################################################################################################################################################################################################################################################################
#Import data 

getwd() 
asv_table_rarefied21282 <- read.delim("rarefied_feature_table-21282.txt", sep = '\t', header = T, row.names = 1, strip.white = T)
asv_table_wide_rarefied21282 <- t(asv_table_rarefied21282) # transpose to wide format
metadata <- read.delim("metadata_daynight.txt", header = T, row.names = 1,comment.char="#")

##########################################################################################################################################################################################################################################################################################################################################################################################################################
## Rarefied21282 data - Statistical analysis on dist matrix ----
#PERMANOVA 

#Remove any taxa that is not present in the any of the samples first
asv_table_wide_rarefied21282_onlytaxapresentineitheronesamples <- asv_table_wide_rarefied21282[,(colSums(asv_table_wide_rarefied21282) > 0 )]

#Verification steps - that only taxa present in either one of the samples
#Calculate column totals
total <- colSums(asv_table_wide_rarefied21282_onlytaxapresentineitheronesamples)
asv_w_total <- rbind(asv_table_wide_rarefied21282_onlytaxapresentineitheronesamples, total)

# re-order columns by total - check the last column
asv_sort <- asv_w_total[,order(-asv_w_total[which(rownames(asv_w_total) == 'total'), ])] #Here essentially it is re-ordering the total by descending order
asv_sort <- asv_sort[-19,] #remove the "total" row

#Know which one to combine, here I didn't use any of asv_sort because I didn't removed any samples
asv_rarefied21282_metadata <- cbind(asv_table_wide_rarefied21282, metadata)

#1. PERMANOVA#######################################################################################################################################################################################################################################################################################
#Note: 1. always check if the column is right! Make sure not to include empty count taxa 2. The result always changes because of the random permutation and so the significance return will always be different. 

adonis.microbiome.noday2 <- adonis(asv_rarefied21282_metadata[,1:4586] ~ asv_rarefied21282_metadata$day.time, method = "bray")
adonis.microbiome.noday2$aov.tab

#Post-hoc test: I plan to use fdr as correction method because bonferroni can be a bit strict.
#But I proceed with fdr which was used in qiime2 
pairwiseAdonis::pairwise.adonis(asv_rarefied21282_metadata[,1:4586], asv_rarefied21282_metadata$day.time, 
                                p.adjust.m ='fdr', sim.method = 'bray') 

#2. PERMDISP######################################################################################################################################################################################################################################################################################
#The code below adapted from (Ezzat et al., 2020) – “parrotfish_paper.R” 
rarefied21282_bray_dist <- vegdist(asv_rarefied21282_metadata[,1:4586], method = "bray") #code adapted from O'Brien et al.(2021): "Calculate_and_plot_beta-diversity.R" file in PERMANOVA section
bdisp <- betadisper(rarefied21282_bray_dist, asv_rarefied21282_metadata$day.time, type=c("centroid"))
bdisp   #This is with PCoA remember! and produce the interesting graph which seem significant and different from nmds above. 
aov.bdisp <-anova(bdisp)
aov.bdisp       #aov and permutest (next code) produced the same result - with minor difference in the p-value -see Pat Schloss video on the interpretation

#pairwise comparison methods is here: https://www.rdocumentation.org/packages/vegan/versions/2.4-2/topics/permutest.betadisper
permutest(bdisp, permutations = 999, pairwise = TRUE) #Significant : How to interpret more: See this http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html


# PCoA plotting #########################################################################################################################################################################################################################################################################################
# calculate principal coordinates analysis (Bray-Curtis)
bray.pcoa.eigen <- cmdscale(rarefied21282_bray_dist, k = 2, eig = T)

# extract axis positions for each sample from cmdscale object and create a dataframe for plotting
bray.pcoa.eigen.plotting <- as.data.frame(bray.pcoa.eigen$points)
colnames(bray.pcoa.eigen.plotting) <- c("axis_1", "axis_2")
bray.pcoa.eigen.plotting$sample <- rownames(bray.pcoa.eigen.plotting)

bray.pcoa.eigen.plotting <- cbind(bray.pcoa.eigen.plotting, metadata)

# calculate the proportion of variance in the data which is explained by the first two PCoA axes
#Understand what is eigenvalue and eigenvector in here: https://www.youtube.com/watch?v=FgakZw6K1QQ
(bray.pcoa.eigen$eig[1]/(sum(bray.pcoa.eigen$eig)))*100 #20.22653%

(bray.pcoa.eigen$eig[2]/(sum(bray.pcoa.eigen$eig)))*100 #15.93464%

# create a PCoA plot - result same as QIIME2!
pcoa.meio.bray.plot <- ggplot(bray.pcoa.eigen.plotting, aes(x = axis_1, y = axis_2, colour = day.time)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95) +
  theme_bw() + 
  xlab("PCoA 1 (20.23%)") +
  ylab("PCoA 2 (15.93%)") +
  labs(color = "Timepoints") +
  scale_color_manual(values = c("1_0915"="blue", "1_1230"="red", "1_1612"="#04CE13", "1_1740"="purple"),labels = c("morning", "noon", 'early-\nafternoon', "late-\nafternoon") ) +
  theme(axis.text.y = element_text(colour = "black", size = 12), 
        axis.text.x = element_text(colour = "black", size = 12), 
        legend.text = element_text(size = 12, colour ="black"), 
        legend.position = "right", axis.title.y = element_text( size = 14), 
        axis.title.x = element_text(size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black"), #changing the title of legend
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank())
pcoa.meio.bray.plot + annotate("text", x = -0.6, y = 0.6, label = "PERMANOVA:\n p-value > 0.05\n PERMDISP:\n p-value < 0.05")

#ggsave(filename = "pcoa_plot-finalised.tiff", width = 20, height = 15, units = "cm", device = "tiff", dpi = "print")

