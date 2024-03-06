getwd()

############Download gg3D#############################################################################################################################################################################################################################################################################################################
devtools::install_github("AckerDWM/gg3D")

######################################################################################################################################################################################################################################################################################################################################
## Load packages
library(vegan)
library(cluster)
library(pairwiseAdonis)
library(dplyr)
library(ape)
library(ggplot2)
library(gg3D) #gg3D is a package created to extend ggplot2 to produce 3D plots

#Import data 
metabolites_all <- read.delim("metabolite_products_rawdata.txt", sep = '\t',row.names = 1 , header = T, strip.white = T) 

metadata <- read.delim("metadata_daynight.txt", header = T, row.names = 1,comment.char="#")


#check the variance of each metabolites (for assessment to transform data or not) : https://www.quantitative-biology.ca/multi.html#multi
apply(metabolites_all[,1:15], 2, var) #the variance is super huge across all the data


#Before transformation (PC1 values)
PCA_beforetransformation <- prcomp(metabolites_all[,1:15])
summary(PCA_beforetransformation) #variance quite variable across metabolites


#PCA with automatically transformed datasets 
PCA_aftertransformation <- prcomp(metabolites_all[,1:15], scale = TRUE)
summary(PCA_aftertransformation) # The variance greatly reduced. 

#scaling the datasets
metabolites_all_autoscaled <- data.frame(scale(metabolites_all[,1:15], center = TRUE, scale = TRUE)) 

##PCA with manually transformed datasets 
PCA_aftertransformation_autoscaled <- prcomp(metabolites_all_autoscaled[,1:15], scale = FALSE)
summary(PCA_aftertransformation_autoscaled) #PC1 reduced
apply(metabolites_all_autoscaled[,1:15], 2, var) #the variance is the same now

#############################################################################################################################################################################################################
## Statistical analysis on dist matrix ----
## Combine the data first ###
metabolites_all_autoscaled_metadata <- cbind(metabolites_all_autoscaled, metadata)
metabolites_all_autoscaled_metadata$date <- factor(metabolites_all_autoscaled_metadata$date)
metabolites_all_autoscaled_metadata$day.time <- factor(metabolites_all_autoscaled_metadata$day.time)

#1. PERMANOVA#######################################################################################################################################################################################################################################################################################
#Note: 1. always check if the column is right! it could be not. 2. The result always changes because of the random permutation and so the significance returned will always be different. 
adonis.nonenzymaticmetabolites_autoscaled <- adonis(metabolites_all_autoscaled_metadata[,1:15] ~ metabolites_all_autoscaled_metadata$day.time, method = "euclidean")
adonis.nonenzymaticmetabolites_autoscaled$aov.tab

#                                                      Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# metabolites_all_autoscaled_metadata$day.time  3    80.865  26.955  2.1671 0.31712  0.011 *
# Residuals                                    14   174.135  12.438         0.68288         
# Total                                        17   255.000                 1.00000   

#Post-hoc test
pairwise.adonis(metabolites_all_autoscaled_metadata[,1:15], metabolites_all_autoscaled_metadata$day.time, 
                p.adjust.m ='fdr', sim.method = 'euclidean') 

#2. PERMDISP######################################################################################################################################################################################################################################################################################
metabolites_all_autoscaled_metadata_euclidean_dist <- vegdist(metabolites_all_autoscaled_metadata[,1:15], method = "euclidean") 

bdisp <- betadisper(metabolites_all_autoscaled_metadata_euclidean_dist, metabolites_all_autoscaled_metadata$day.time, type=c("centroid"))
bdisp   
aov.bdisp <-anova(bdisp)
aov.bdisp       #aov and permutest (next code) produced the same result - with minor difference in the p-value -see Pat Schloss video on the interpretation CC208
permutest(bdisp, pairwise = TRUE) #Significant : How to interpret more: See this http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html

# Response: Distances
#              Df  Sum Sq Mean Sq      F N.Perm Pr(>F)
# Groups     3  1.8921 0.63071 0.6869    999  0.629
# Residuals 14 12.8545 0.91818  

#Simple overview with  PCoA plot
labs <- paste("Dimension", 1:8, "(", 
              round(100*bdisp$eig / sum(bdisp$eig), 2), "%)")

plot(bdisp, cex=1, pch=15:17,
     main="Euclidean Dissimilarity index", cex.lab=1.25,
     xlab=labs[1], ylab=labs[2],
     hull=FALSE, ellipse=TRUE, conf=0.95, lwd=2)  


#PCoA plotting #########################################################################################################################################################################################################################################################################################
euclidean.pcoa.eigen <- cmdscale(metabolites_all_autoscaled_metadata_euclidean_dist, k = 2, eig = T)

# extract axis positions for each sample from cmdscale object and create a dataframe for plotting
euclidean.pcoa.eigen.plotting <- as.data.frame(euclidean.pcoa.eigen$points)
colnames(euclidean.pcoa.eigen.plotting) <- c("axis_1", "axis_2")
euclidean.pcoa.eigen.plotting$sample <- rownames(euclidean.pcoa.eigen.plotting)

euclidean.pcoa.eigen.plotting <- cbind(euclidean.pcoa.eigen.plotting, metadata)

#Calculate the proportion of variance in the data which is explained by the first two PCoA axes
(euclidean.pcoa.eigen$eig[1]/(sum(euclidean.pcoa.eigen$eig)))*100 #38.54946%

(euclidean.pcoa.eigen$eig[2]/(sum(euclidean.pcoa.eigen$eig)))*100 #27.97446%

# create a PCoA plot
pcoa.meio.bray.plot <- ggplot(euclidean.pcoa.eigen.plotting, aes(x = axis_1, y = axis_2, colour = day.time)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95) +
  theme_bw() + 
  xlab("PCoA 1 (38.55%)") +
  ylab("PCoA 2 (27.97%)") +
  labs(color = "Timepoints") +
    scale_color_manual(values = c("T_0915"="blue", "T_1230"="red", 'T_1612'="#04CE13", "T_1740"="purple"),labels = c("morning", "noon", 'early-\nafternoon', "late-\nafternoon") ) +
    theme(axis.text.y = element_text(colour = "black", size = 12), 
        axis.text.x = element_text(colour = "black", size = 12), 
        legend.text = element_text(size = 12, colour ="black"), 
        legend.position = "right", axis.title.y = element_text( size = 14), 
        axis.title.x = element_text(size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black"), #changing the title of legend
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) 
pcoa.meio.bray.plot
pcoa.meio.bray.plot + annotate("text", x =5, y = 9.5, label = "PERMANOVA: p-value < 0.05\n PERMDISP: p-value > 0.05")
#ggsave(filename = "xxx.pdf", width = 20, height = 15, units = "cm", device = "pdf", dpi = "print")

