## Use SIMPER analysis to identify metabolites that may have more contribution than the others 

library(vegan)

getwd() 

metabolites_all <- read.delim("metabolite_products_rawdata_v2.txt", sep = '\t',row.names = 1 , header = T, strip.white = T) #this version follows the data structure in the file: "asv_table_wide_rarefied5612"

metadata <- read.delim("metadata_daynight.txt", header = T, row.names = 1,comment.char="#")


##########################################################################################################################################################################################################################################################################################################################################################################################################################
#SIMPER analysis

simpler.test <- simper(metabolites_all, metadata$day.time)
sum.simpler <- summary(simp.test)
sum.simpler

# Pull out the comparison you want to look at if you want
morning-noon <- sum.simpler$T_0915_T_1230

View(morning-noon)

write.csv(morning-noon, file = "morning-noon.csv")

