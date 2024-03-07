#Relative abundance bar plots 

getwd()

#Loading library
library(reshape2)
library(ggplot2)
library(scales)
library(dplyr)
library(forcats)
library(gridExtra)
library(RColorBrewer)
library(randomcoloR)

##############################################################################################################################
#####Taxonomic profiles of the 30 most abundant families######################################################################
##############################################################################################################################


## Plot family level (Top 30) ----
asv_table_f_30 <- read.csv("raw_taxa_fam.csv", sep = ',', header = T, row.names = 1, strip.white = T) # Family table created using qiime2
asv_table_f_30 <- asv_table_f_30[1:20,] #remove Day2 samples, not in analysis
asv_table_f_30 <- asv_table_f_30[-c(8,18),] #remove BRM3 and DRM3

#subset metadata from the .csv file - need to do this because when we merged the table it has to be in the same number of roles and column
meta_data <- subset(asv_table_f_30, select=c("sample.name","date", "day", "time", "day.time", "body.temperature", "body.size.mm", "extraction.batch")) 
meta_data <- subset(meta_data, sample.name != "BRM3" &  sample.name != "DRM3") #rarefied table have no these samples

asv_table_f2_30 <- as.matrix(asv_table_f_30[, -c(400:407)]) # Remove metadata columns from the asv_table_f

# Calculate Top 30 abundant family #################################################################

# Calculate column totals
total_30 <- colSums(asv_table_f2_30)
asv_w_total_30 <- rbind(asv_table_f2_30, total_30)

# re-order columns by total
asv_sort_30 <- asv_w_total_30[,order(-asv_w_total_30[which(rownames(asv_w_total_30) == 'total_30'), ])] #Here essentially it is re-ordering the total by descending order

## Convert counts to percentages on data matrix & remove row total
asv_percent_30 <- asv_sort_30[-which(rownames(asv_sort_30) == 'total_30'),] / rowSums(asv_sort_30[-which(rownames(asv_sort_30) == 'total_30'),]) * 100

# create new column "other" which is the sum of all the taxa not in top 30, remove other taxa
Other <- rowSums(asv_percent_30[,-c(1:30)])

# Combine data
asv_top30 <- cbind(asv_percent_30[,1:30], Other)
asv_rel_abunF_30 <- cbind(asv_top30, meta_data)
asv_table_longF_30 <- melt(asv_rel_abunF_30, variable.name = "Family",id=c("sample.name", "date", "day", "time", 
                                                                           "day.time", "body.temperature", "body.size.mm",
                                                                           "extraction.batch"))
#Just as a reference - the opposite of melt is dcast in reshape2 package or you can youse reshape function to change it back to wide mode####################################################################################################################################################################
#Method 1
asv_table_longF_30_recast <- dcast(asv_table_longF_30, Family   ~ sample.name)
asv_table_longF_30_recast <- t(asv_table_longF_30_recast)
colnames(asv_table_longF_30_recast) <- asv_table_longF_30_recast[1,]
asv_table_longF_30_recast <- asv_table_longF_30_recast[-1,]
asv_table_longF_30_recast_metadata <- cbind(asv_table_longF_30_recast, meta_data) # need this extra step to get the metadata

#method2 - better method
method2 <-reshape(asv_table_longF_30, 
        idvar = c("sample.name", "date", "day", "time", 
                  "day.time", "body.temperature", "body.size.mm",
                  "extraction.batch"),
        timevar = "Family",
        direction = "wide")
colnames(method2) <- gsub("value.", "", colnames(method2))

##########################################################################################################################################################################################################################################

# Remove long names
unique(asv_table_longF_30$Family)

#First example: asv_table_longF$Family <- gsub("d__Bacteria.__.__.__.__", "Bacteria unknown", asv_table_longF$Family) #Always start with renaming first
#Second example: asv_table_longF$Family <- gsub("Unassigned.__.__.__.__", "Unassigned", asv_table_longF$Family)
asv_table_longF_30$Family <- gsub("Rickettsiales.__", "Rickettsiales unknown", asv_table_longF_30$Family)
asv_table_longF_30$Family <- gsub("Bacteroidia.__.__", "Bacteroidia unknown", asv_table_longF_30$Family)
asv_table_longF_30$Family <- gsub("HOC36.f__", "HOC36 unknown", asv_table_longF_30$Family)
asv_table_longF_30$Family <- gsub("Gracilibacteria.o__.f__", "Gracilibacteria unknown", asv_table_longF_30$Family)
asv_table_longF_30$Family <- gsub("Alphaproteobacteria.__.__", "Alphaproteobacteria unknown ", asv_table_longF_30$Family)
asv_table_longF_30$Family <- gsub("Flavobacteriales.__", "Flavobacteriales unknown", asv_table_longF_30$Family)
asv_table_longF_30$Family <- gsub("d__Bac.*f__", "\\1", asv_table_longF_30$Family)  #The order here is super important - do one by one to remove those unwanted parts from the names
asv_table_longF_30$Family <- gsub("d__Bac.*o__", "\\1", asv_table_longF_30$Family)  # The * means "to" : From d__Bac. "to" "p__" | or you can think of it as "times" : "a name"*"a name" as a vector. You can add a fullstop if you dont want to write the full name. e.g. Bac. instead of bacteria before the "*". Only with "*" you can do this
asv_table_longF_30$Family <- gsub("d__Bac.*c__", "\\1", asv_table_longF_30$Family)
asv_table_longF_30$Family <- gsub("d__Bac.*p__", "\\1", asv_table_longF_30$Family)

unique(asv_table_longF_30$Family) #Check again if anything is weird

# Factorizing ,grouping and summarise : this part is super important! Make sure you have factorise all the factor
asv_table_longF_30$day.time <- factor(asv_table_longF_30$day.time)
asv_table_longF_30$date <- factor(asv_table_longF_30$date)
asv_table_longF_30$sample.name <- factor(asv_table_longF_30$sample.name)
asv_table_longF_30$body.temperature <- factor(asv_table_longF_30$body.temperature)
asv_table_longF_30$body.size.mm <- factor(asv_table_longF_30$body.size.mm)
asv_table_longF_30$extraction.batch <- factor(asv_table_longF_30$extraction.batch)
asv_table_longF_30$day <- factor(asv_table_longF_30$day)
asv_table_longF_30$time <- factor(asv_table_longF_30$time)

#Group by day.time
group_taxonF_30 <- asv_table_longF_30 %>% 
  group_by(sample.name, day.time,Family)
group_taxonF_30

taxon_summaryF_30 <- summarise(group_taxonF_30, Percentage = mean(value, na.rm = T))
taxon_summaryF_30 

taxon_summaryF_30 <- taxon_summaryF_30 %>% mutate(Family = fct_relevel(Family, "Other", after = Inf)) #this part for re-leveling; From ?help description: Using 'Inf' allows you to relevel to the end when the number of levels is unknown or variable (e.g. vectorised operations)


#31 colors vector (The colors was chosen from distinctColorPalette() function in the package called randomcoloR) Reference:  https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
#Color choosing (the process is random and I selected the most attractive and distinctive ones) - You can ignore these and choose the one below if you want colors fewer than 26/31
n <- 31
palette <- distinctColorPalette(n)
palette

P31 = c("#7D7187","#5EB7DC","#DBD784","#98E0AC","#E587DA","#D359D8","#6991DF","#DCAE4D","#5E62D9","#D88F86","#819A67","#7BE54F","#DAE9C5",
        "#A372D3","#AAB9E3","#EC854D","#DA3DE8","#7A37D8","#E08AB4","#ADE4E4","#E3D9DA","#DD5772","#E3BEDB","#E4BA98","#AEE679","#5EE195",
        "#DBE247","#C89FDD","#E546A2","#60DBCC","#97AAA4")
new_daytime_label <- c("1_0915" = "morning", "1_12.30" = "noon",  "1_1612" ="early-afternoon", "1_1740" = "late-afternoon")
ggplot() +
  geom_bar(aes(y = Percentage, x = sample.name, fill = Family), data = taxon_summaryF_30, stat="identity", position = "fill") +  
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) + 
  labs(y= "Relative Abundance (Family)\n", x="Sample", fill = "Taxa") +
  scale_fill_manual(values=P31) + 
  theme_bw() +
  theme(axis.text=element_text(size=11),
        axis.text.x=element_text(angle = 90), 
        axis.title=element_text(size = 13),
        strip.background = element_rect(fill = "white", colour = "black"),
        legend.key.size = unit(0.2, "cm"),
        legend.key.width = unit(0.2,"cm"), 
        legend.position = 'right') + 
  guides(fill=guide_legend(ncol=1)) +
  facet_wrap(~day.time, drop = TRUE, scales = "free_x", labeller = as_labeller(new_daytime_label),
           nrow = 1, ncol = 6) 

#ggsave(filename = "raw_top30_abundfamily_finalised.png", device = "png", width = 23, height = 11, units = "cm", dpi = "print")

