#Use an unrarefied table (ANCOM-BC has built-in normalization procedures)
getwd()

#Loading library
library(phyloseq)
library(tidyr)
library(reshape2)
library(dplyr)
library(ggplot2)
library(rbiom)
library(ranacapa)
library(microbiome)
library(DT)
library(ANCOMBC) #constantly updated

#import data 
asv <- read.table("raw_feature_table.txt", sep = '\t', row.names = 1, header = T, strip.white = T)    
asv <- asv[,colnames(asv) != "BRM3" & colnames(asv) != "DRM3"] #removed these samples, because they are not in the rarefied table

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


# 'Phyloseq-ize' the data 
otu.t= otu_table(asv, taxa_are_rows=TRUE)
sam.t= sample_data(data.frame(map))
tax.t= tax_table(as.matrix(tax))

phy.all= phyloseq(otu.t, tax.t,  sam.t)
phy.all

################################################################################################################################################################################################################################################################################################################################################################################################
##ANCOM-BC Analysis ############################################################################################################################################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################################################################################################################################################################
#Extract taxa at family level
family_data <- aggregate_taxa(phy.all, "Family") 
tax_table(family_data)

#Comparison 1: 1_0915 vs 1_1230  (C1)############################################################################################################################################################################################################################################################################################################################################################
C1=subset_samples(family_data, day.time %in% c("1_0915", "1_1230"))
C1

res1 = ancombc(phyloseq = C1,formula = "day.time", 
             p_adj_method = "holm",zero_cut = 0.9, lib_cut=1000,                  
             group = "day.time",struc_zero = TRUE ,neg_lb = FALSE, tol = 1e-05, 
             max_iter = 100,conserve = TRUE,alpha = 0.05, global = TRUE)


res1_df = data.frame(Family= row.names(res1$res$beta), beta=res1[["res"]][["beta"]], se=res1[["res"]][["se"]], W=res1[["res"]][["W"]],pval=res1[["res"]][["p_val"]], 
                     qval=res1[["res"]][["q_val"]],DA=res1[["res"]][["diff_abn"]])  #[[]][[]] means selecting the item(or column) within list - we are creating a new data frame here by extracting each parameters and combine them all rather than doing one by one 
                                                                                    
colnames(res1_df)=c("Family","lfc",	"se", "W",	"pval",	"qval", "Diff_abundant") #give the correct names to each of the column 
res1_df$Comparison="Comparison 1: D1_0915 vs D1_1230" 
res1_sig = subset(res1_df, res1_df$Diff_abundant == "TRUE") # Here selecting all that shows significant difference
res1_sig$Diff_more_abundant=ifelse(res1_sig$W < 0 , "1_0915", "1_1230") #Here it tells that if less than W < 0 means it comes from the first group = 1_0915 , 
                                                                        #otherwise it is the 2nd group = 1_1230. You can refer to the Lin & Peddada (2020) paper on how the statistics was calculated. But based on a query to Lin in Github: https://github.com/FrederickHuangLin/ANCOMBC/issues/25
                                                                        #If <0, it is the reference that is high, if > 0 then the other group is high in abundance.

message("Number of DA ASVs: ", nrow(res1_sig), "\nNumber of DA ASVs enriched in D1_0915: ", nrow(subset(res1_sig, Diff_more_abundant == "1_0915" )), "\nNumber of DA ASVs enriched in  D1_1230: ", nrow(subset(res1_sig, Diff_more_abundant == "1_1230")))

#Comparison 2: 1_0915 vs 1_1612  (C2)############################################################################################################################################################################################################################################################################################################################################################
?subset_samples
C2=subset_samples(family_data, day.time %in% c("1_0915", "1_1612"))
C2
res2 = ancombc(phyloseq = C2,formula = "day.time", 
               p_adj_method = "holm",zero_cut = 0.9,lib_cut=1000,                  
               group = "day.time",struc_zero = TRUE ,neg_lb = FALSE, tol = 1e-05, 
               max_iter = 100,conserve = TRUE,alpha = 0.05, global = TRUE)

res2_df = data.frame(Family= row.names(res2$res$beta),Beta=res2[["res"]][["beta"]], Beta_se=res2[["res"]][["se"]], W=res2[["res"]][["W"]],pval=res2[["res"]][["p_val"]], 
                     qval=res2[["res"]][["q_val"]],DA=res2[["res"]][["diff_abn"]])  

colnames(res2_df)=c("Family","Beta",	"se", "W",	"pval",	"qval", "Diff_abundant") 
res2_df$Comparison="Comparison 1: D1_0915 vs D1_1612" 
res2_sig = subset(res2_df, res2_df$Diff_abundant == "TRUE") 
res2_sig$Diff_more_abundant=ifelse(res2_sig$W < 0 , "1_0915", "1_1612") 
message("Number of DA ASVs: ", nrow(res2_sig), "\nNumber of DA ASVs enriched in D1_0915: ", nrow(subset(res2_sig, Diff_more_abundant == "1_0915" )), "\nNumber of DA ASVs enriched in  D1_1612: ", nrow(subset(res2_sig, Diff_more_abundant == "1_1612")))

#Comparison 3: 1_0915 vs 1_1740  (C3)############################################################################################################################################################################################################################################################################################################################################################
C3=subset_samples(family_data, day.time %in% c("1_0915", "1_1740"))
C3
res3 = ancombc(phyloseq = C3,formula = "day.time", 
               p_adj_method = "holm",zero_cut = 0.9,lib_cut=1000,                  
               group = "day.time",struc_zero = TRUE ,neg_lb = FALSE, tol = 1e-05, 
               max_iter = 100,conserve = TRUE,alpha = 0.05, global = TRUE)

res3_df = data.frame(Family= row.names(res3$res$beta),Beta=res3[["res"]][["beta"]], Beta_se=res3[["res"]][["se"]], W=res3[["res"]][["W"]],pval=res3[["res"]][["p_val"]], 
                     qval=res3[["res"]][["q_val"]],DA=res3[["res"]][["diff_abn"]])  

colnames(res3_df)=c("Family","Beta",	"se", "W",	"pval",	"qval", "Diff_abundant") 
res3_df$Comparison="Comparison 1: D1_0915 vs D1_1740" 
res3_sig = subset(res3_df, res3_df$Diff_abundant == "TRUE") 
res3_sig$Diff_more_abundant=ifelse(res3_sig$W < 0 , "1_0915", "1_1740") 
message("Number of DA ASVs: ", nrow(res3_sig), "\nNumber of DA ASVs enriched in D1_0915: ", nrow(subset(res3_sig, Diff_more_abundant == "1_0915" )), "\nNumber of DA ASVs enriched in  D1_1740: ", nrow(subset(res3_sig, Diff_more_abundant == "1_1740")))

#Comparison 6: 1_1230 vs 1_1612  (C6)############################################################################################################################################################################################################################################################################################################################################################
C6=subset_samples(family_data, day.time %in% c("1_1230", "1_1612"))
C6
res6 = ancombc(phyloseq = C6,formula = "day.time", 
               p_adj_method = "holm",zero_cut = 0.9,lib_cut=1000,                  
               group = "day.time",struc_zero = TRUE ,neg_lb = FALSE, tol = 1e-05, 
               max_iter = 100,conserve = TRUE,alpha = 0.05, global = TRUE)

res6_df = data.frame(Family= row.names(res6$res$beta), Beta=res6[["res"]][["beta"]], Beta_se=res6[["res"]][["se"]], W=res6[["res"]][["W"]],pval=res6[["res"]][["p_val"]], 
                     qval=res6[["res"]][["q_val"]],DA=res6[["res"]][["diff_abn"]])  

colnames(res6_df)=c("Family","Beta",	"se", "W",	"pval",	"qval", "Diff_abundant") 
res6_df$Comparison="Comparison 1: D1_1230 vs D1_1612" 
res6_sig = subset(res6_df, res6_df$Diff_abundant == "TRUE") 
res6_sig$Diff_more_abundant=ifelse(res6_sig$W < 0 , "1_1230", "1_1612") 
message("Number of DA ASVs: ", nrow(res6_sig), "\nNumber of DA ASVs enriched in D1_1230: ", nrow(subset(res6_sig, Diff_more_abundant == "1_1230" )), "\nNumber of DA ASVs enriched in  D1_1612: ", nrow(subset(res6_sig, Diff_more_abundant == "1_1612")))

#Comparison 7: 1_1230 vs 1_1740  (C7)############################################################################################################################################################################################################################################################################################################################################################
C7=subset_samples(family_data, day.time %in% c("1_1230", "1_1740"))
C7
res7 = ancombc(phyloseq = C7,formula = "day.time", 
               p_adj_method = "holm",zero_cut = 0.9,lib_cut=1000,                  
               group = "day.time",struc_zero = TRUE ,neg_lb = FALSE, tol = 1e-05, 
               max_iter = 100,conserve = TRUE,alpha = 0.05, global = TRUE)

res7_df = data.frame(Family= row.names(res7$res$beta),Beta=res7[["res"]][["beta"]], Beta_se=res7[["res"]][["se"]], W=res7[["res"]][["W"]],pval=res7[["res"]][["p_val"]], 
                     qval=res7[["res"]][["q_val"]],DA=res7[["res"]][["diff_abn"]])  

colnames(res7_df)=c("Family","Beta",	"se", "W",	"pval",	"qval", "Diff_abundant") 
res7_df$Comparison="Comparison 1: D1_1230 vs D1_1740" 
res7_sig = subset(res7_df, res7_df$Diff_abundant == "TRUE") 
res7_sig$Diff_more_abundant=ifelse(res7_sig$W < 0 , "1_1230", "1_1740") 
message("Number of DA ASVs: ", nrow(res7_sig), "\nNumber of DA ASVs enriched in D1_1230: ", nrow(subset(res7_sig, Diff_more_abundant == "1_1230" )), "\nNumber of DA ASVs enriched in  D1_1740: ", nrow(subset(res7_sig, Diff_more_abundant == "1_1740")))

#Comparison 10: 1_1612 vs 1_1740  (C10)############################################################################################################################################################################################################################################################################################################################################################
C10=subset_samples(family_data, day.time %in% c("1_1612", "1_1740"))
C10
res10 = ancombc(phyloseq = C10,formula = "day.time", 
               p_adj_method = "holm",zero_cut = 0.9,lib_cut=1000,                  
               group = "day.time",struc_zero = TRUE ,neg_lb = FALSE, tol = 1e-05, 
               max_iter = 100,conserve = TRUE,alpha = 0.05, global = TRUE)

res10_df = data.frame(Family= row.names(res10$res$beta), Beta=res10[["res"]][["beta"]], Beta_se=res10[["res"]][["se"]], W=res10[["res"]][["W"]],pval=res10[["res"]][["p_val"]], 
                     qval=res10[["res"]][["q_val"]],DA=res10[["res"]][["diff_abn"]])  

colnames(res10_df)=c("Family","Beta",	"se", "W",	"pval",	"qval", "Diff_abundant") 
res10_df$Comparison="Comparison 1: D1_1612 vs D1_1740" 
res10_sig = subset(res10_df, res10_df$Diff_abundant == "TRUE") 
res10_sig$Diff_more_abundant=ifelse(res10_sig$W < 0 , "1_1612", "1_1740") 
message("Number of DA ASVs: ", nrow(res10_sig), "\nNumber of DA ASVs enriched in D1_1612: ", nrow(subset(res10_sig, Diff_more_abundant == "1_1612" )), "\nNumber of DA ASVs enriched in  D1_1740: ", nrow(subset(res10_sig, Diff_more_abundant == "1_1740")))

#######################################################################################################################################
#Combining all the result and produce a tabulated file
ANCOMresults_sig <- rbind(res1_sig,res2_sig, res3_sig,res6_sig, res7_sig,res10_sig)
ANCOMresults_sig
getwd()
write.table(ANCOMresults_sig,  "ANCOMBC_Familylevel_results_significant_taxa.txt", 
            sep = "\t", quote = F, row.names = T ) #Create a new file

ANCOMresults_df <- rbind(res1_df,res2_df, res3_df,res6_df, res7_df,res10_df)
ANCOMresults_df

write.table(ANCOMresults_df,  "ANCOMBC_Familylevel_results_alltaxa.txt", 
            sep = "\t", quote = F, row.names = T ) #Create a new file


