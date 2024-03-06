#Loading library## Load packages
library(ggplot2)
library(FSA)
library(phyloseq)
library(forcats)
library(reshape2)
library(dplyr) 
library(carData)
library(car)
library(multcomp)
library(permute)
library(lattice)
library(vegan)
library(performance)
library(cowplot)
library(vegan)

###Rarefied data-21282###############################################################################################################################################################################################################################################################################################################################################################################################################################################################################
#Import data 
getwd()
asv_table_rarefied21282 <- read.delim("rarefied_feature_table-21282.txt", sep = '\t', header = T, row.names = 1, strip.white = T)
asv_table_rarefied21282 <- asv_table_rarefied21282[,1:18]#remove the balance from samples from D2
asv_table_rarefied21282 <- t(asv_table_rarefied21282)      #t() here means transpose the matrix

metadata <- read.delim("metadata_daynight.txt", header = T, row.names = 1,comment.char="#")

#the %in% operator returns a logical vector that tells us whether or not there is a match between the first object and the second and subset() does only for those that present. This checking is important because sometimes we want to might rarefied to sampling depth that might removed some samples. So we want to know what are those samples. 
metadata  <- subset(metadata, (rownames(metadata) %in%  rownames(asv_table_rarefied21282))) 

### Observed ASV ##############################################################################################################################################################################################################
## Plot and calculate statistics for species richness (in this case, ASV richness) 
ASVrichness <- apply(asv_table_rarefied21282  >0, 1, sum) 

rich_meta_rarefied21282<- cbind(ASVrichness, metadata)  

# Setting some data as categorical factor:
rich_meta_rarefied21282$sample.name <- factor(rich_meta_rarefied21282$sample.name)
rich_meta_rarefied21282$date <- factor(rich_meta_rarefied21282$date)
rich_meta_rarefied21282$extraction.batch <- factor(rich_meta_rarefied21282$extraction.batch)
rich_meta_rarefied21282$day.time <- factor(rich_meta_rarefied21282$day.time)

# Plot as ASVrichness boxplot - By day.time
ASVrichness_rarefied21282_daytime <- ggplot(rich_meta_rarefied21282, aes(x=day.time, y=ASVrichness, fill=day.time)) +
                        geom_boxplot(alpha = 1, outlier.shape = 16, outlier.size = 2) +
                        #scale_fill_discrete(labels = c("D1_0915", "D1_1230", "D1_1612", "D1_1740")) +
                        geom_jitter(color="black", size=2) +
                        scale_x_discrete(name = "Timepoints", labels = c("1_0915" = "morning", "1_1230" = "noon",  "1_1612" ="early-\nafternoon", "1_1740" = "late-\nafternoon")) +
                        scale_y_continuous(name = "ASV richness", limits = c(0, 1000), breaks = c(0, 100,200,300,400,500,600,700,800,900,1000)) +
                        scale_fill_manual(values = c("blue", "red", "#04CE13", "purple")) +
                        #guides(fill = guide_legend(title = "Day_Time")) +
                        theme_bw() +
                        theme(axis.text=element_text(size=11),
                              axis.title=element_text(size = 13),
                              legend.position = "none"
                              #,legend.text=element_text(size=11), 
                              #legend.title=element_text(size=13)
                              )  + annotate("text", x = 3.5, y = 950, label = expression("F"["(3,14)"]*" = 3.67 "*"; p-value < 0.05")  , size= 5)
                        

ASVrichness_rarefied21282_daytime

#Create model - rarefied21282 data
#Calculate statistics By day.time
asv_rarefied21282.lm <- lm(ASVrichness ~ day.time, rich_meta_rarefied21282) #create a model
anova(asv_rarefied21282.lm) #Result: F-value: 3.6655; p-value:  0.03872*

#Normality test - rarefied21282 data 
asv_rarefied21282.lm.residuals <- residuals(object = asv_rarefied21282.lm)
shapiro.test(x = asv_rarefied21282.lm.residuals) #result: W = 0.94686, p-value = 0.3779

#Homogeneity of Variance - raw data
asv_rarefied21282.lm.levene <- leveneTest(asv_rarefied21282.lm, center =mean)
asv_rarefied21282.lm.levene   #Result: Passed: p-value =  0.5579

#Post-hoc test 
glht.rich.sum.xB <- summary(glht(asv_rarefied21282.lm, linfct = mcp(day.time = "Tukey"))) 
glht.rich.sum.xB #The result is the same as the next section post-hoc test for asv_rarefied21282.aov - only to two decimal places.
rich.groups <- cld(glht.rich.sum.xB) #EXTRACT THE INFORMATION FROM THE TEST
rich.groups.df_rarefied21282 <- fortify(rich.groups) #fortify() function convert a complex "data.type" into a data frame. only some information can be made into data frame. 
colnames(rich.groups.df_rarefied21282) <- c("day.time", "letters") #The letters refers to statistical significance. seems like all of it is not significant from time to time

ymax_rarefied21282 <- tapply(rich_meta_rarefied21282$ASVrichness, rich_meta_rarefied21282$day.time, max) #get the MAX ASV richness value by grouping accoridn to day,time
rich.groups.df_rarefied21282$Ymax <- ymax_rarefied21282 # add new column information

ASVrichness_rarefied21282_daytime_finalised <- ASVrichness_rarefied21282_daytime + 
                                      geom_text(data=rich.groups.df_rarefied21282, aes(x = day.time, y = Ymax+25, label = letters), size = 7,vjust=-0.5)
ASVrichness_rarefied21282_daytime_finalised 


### Shannon Diversity ##############################################################################################################################################################################################################
shannon <- diversity(asv_table_rarefied21282, index = "shannon")
shannon_rarefied21282_meta <- cbind(shannon, metadata)

# Setting some data as categorical factor:
shannon_rarefied21282_meta$sample.name <- factor(shannon_rarefied21282_meta$sample.name)
shannon_rarefied21282_meta$date <- factor(shannon_rarefied21282_meta$date)
shannon_rarefied21282_meta$extraction.batch <- factor(shannon_rarefied21282_meta$extraction.batch)
shannon_rarefied21282_meta$day.time <- factor(shannon_rarefied21282_meta$day.time)

# Plot as shannon boxplot - By day.time
shannon_rarefied21282_daytime  <- ggplot(shannon_rarefied21282_meta, aes(x=day.time, y=shannon, fill=day.time)) +
                                geom_boxplot(alpha = 1, outlier.shape = 16, outlier.size = 2) +
                                #scale_fill_discrete(labels = c("D1_0915", "D1_1230", "D1_1612", "D1_1740", "D2_0636", "D2_0945")) +
                                geom_jitter(color="black", size=2) +
                                scale_x_discrete(name = "Timepoints", labels = c("1_0915" = "morning", "1_1230" = "noon",  "1_1612" ="early-\nafternoon", "1_1740" = "late-\nafternoon")) +
                                scale_y_continuous(name = "Shannon-Wiener Index (H')\n", limits = c(2,6.0), breaks = c(2,2.5,3,3.5,4,4.5,5,5.5,6)) +
                                scale_fill_manual(values = c("blue", "red", "#04CE13", "purple")) +
                                #guides(fill = guide_legend(title = "Day_Time")) +
                                theme_bw() +
                                theme(axis.text=element_text(size=11),
                                      axis.title=element_text(size = 13),
                                      legend.position = "none"
                                      #,legend.text=element_text(size=11), 
                                      #legend.title=element_text(size=13)
                                      )  + annotate("text", x = 3.5, y = 5.75, label = expression("F"["(3,14)"]*" = 3.76 "*"; p-value < 0.05")  , size= 5)
shannon_rarefied21282_daytime

#Create shannon model - rarefied21282 data
shannon_rarefied21282_meta.lm <- lm(shannon ~ day.time, shannon_rarefied21282_meta) #create a model
anova(shannon_rarefied21282_meta.lm) #Results: Significant : F-value: 3.7594;  p-value: 0.036*

#Normality test - rarefied21282 data 
par(mfrow=c(2,2))
plot(shannon_rarefied21282_meta.lm) #1. Some data point seem to be an outlier, but there Q-Q plot is quite good; 2. Homogeneity of variance seem to be met

shannon_rarefied21282_meta.lm.residuals <- residuals(object = shannon_rarefied21282_meta.lm)
shapiro.test(x = shannon_rarefied21282_meta.lm.residuals) #result: Passed:W = 0.98585, p-value = 0.9903

#Homogeneity of Variance - rarefied21282 data
shannon_rarefied21282_meta.lm.levene <- leveneTest(shannon_rarefied21282_meta.lm, center =mean)
shannon_rarefied21282_meta.lm.levene  #Result: Passed: p-value =  0.2158

#Post-hoc test - Tukey test + remove Bonferroni adjustment##################################################################################################################################################################################################
glht.shan.sum.xB <- summary(glht(shannon_rarefied21282_meta.lm, linfct = mcp(day.time = "Tukey"))) 
glht.shan.sum.xB 
shan.groups <- cld(glht.shan.sum.xB) #EXTRACT THE INFORMATION FROM THE TEST
shan.groups.df_rarefied21282 <- fortify(shan.groups) #fortify() function convert a complex "data.type" into a data frame. only some information can be made into data frame. 
colnames(shan.groups.df_rarefied21282) <- c("day.time", "letters") #The letters refers to statistical significance. seems like all of it is not significant from time to time
ymax_rarefied21282 <- tapply(shannon_rarefied21282_meta$shannon, shannon_rarefied21282_meta$day.time, max)
shan.groups.df_rarefied21282$Ymax <- ymax_rarefied21282 # add new column information

shannon_rarefied21282_daytime_finalised <- shannon_rarefied21282_daytime + 
              geom_text(data=shan.groups.df_rarefied21282, aes(x = day.time, y = Ymax+0.1, label = letters),size = 7, vjust=-1)
shannon_rarefied21282_daytime_finalised

### Observed ASV -temperature correlation ##############################################################################################################################################################################################################
#Create model (linear regression) - rarefied21282 data 
#Calculate statistics By day.time
asv_rarefied21282_temp.lm <- lm(ASVrichness ~ body.temperature, rich_meta_rarefied21282) #create a model
asv_rarefied21282_temp.lm
summary(asv_rarefied21282_temp.lm) #Adjusted R-squared:  0.2716; p < 0.05
#anova(asv_rarefied21282_temp.lm) #Result: F(1,16): 7.339 p-value: 0.01548*

#Normality test - rarefied21282 data 
asv_rarefied21282_temp.lm.residuals <- residuals(object = asv_rarefied21282_temp.lm)
shapiro.test(x = asv_rarefied21282_temp.lm.residuals) #result: W = 0.9689, p-value = 0.7768

xlabs = expression("Temperature " ( degree*C))

# Plot as ASVrichness boxplot - By day.time
ASVrichness_rarefied21282_daytime_templm <- ggplot(rich_meta_rarefied21282, aes(x=body.temperature, y=ASVrichness, col= day.time)) +
  geom_point(size = 3) +
  geom_smooth(method ="lm", se = TRUE, colour = "black") +
  scale_x_continuous(name = expression("Body Temperature " ( degree*C)), limits = c(32,40)) +
  scale_y_continuous(name = "ASV richness", limits = c(100, 1000), breaks = c(100,200,300,400,500,600,700,800,900,1000)) +
  scale_color_manual(values = c("1_0915"="blue", "1_1230"="red", '1_1612'="#04CE13", "1_1740"="purple"),labels = c("morning", "noon", 'early-\nafternoon', "late-\nafternoon") ) +
  theme_bw() +
  labs(color = "Timepoints") + #change the legend title names according to the labels - works for categorical datasets
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size = 13),
        legend.text = element_text(size = 12, colour ="black"), 
        legend.position = "right", axis.title.y = element_text( size = 14), 
        #legend.position = "none",
        legend.title = element_text(size = 14, colour = "black") #changing the title of legend
  )  + annotate("text", x = 37, y = 950, label = expression(italic(R)^2*"= 0.27; p-value < 0.05")  , size= 4.5)
ASVrichness_rarefied21282_daytime_templm


### Shannon Diversity - temperature correlation ##############################################################################################################################################################################################################

shannon <- diversity(asv_table_rarefied21282, index = "shannon")
shannon_rarefied21282_meta <- cbind(shannon, metadata)

# Setting some data as categorical factor:
shannon_rarefied21282_meta$sample.name <- factor(shannon_rarefied21282_meta$sample.name)
shannon_rarefied21282_meta$date <- factor(shannon_rarefied21282_meta$date)
shannon_rarefied21282_meta$extraction.batch <- factor(shannon_rarefied21282_meta$extraction.batch)
shannon_rarefied21282_meta$day.time <- factor(shannon_rarefied21282_meta$day.time)

#Create shannon model - rarefied21282 data
shannon_rarefied21282_meta_temp.lm <- lm(shannon ~ body.temperature, shannon_rarefied21282_meta) #create a model
summary(shannon_rarefied21282_meta_temp.lm) #Adjusted R-squared:  0.2542 ; p < 0.05
#anova(shannon_rarefied21282_meta_temp.lm) #Results: Significant : F-value: 6.7945;  p-value: 0.01908*

#Normality test - rarefied21282 data 
shannon_rarefied21282_meta_temp.lm.residuals <- residuals(object = shannon_rarefied21282_meta_temp.lm)
shapiro.test(x = shannon_rarefied21282_meta_temp.lm.residuals) #result: Passed:W = W = 0.90551, p-value = 0.07171

# Plot as shannon boxplot - By day.time
shannon_rarefied21282_daytime_templm  <- ggplot(shannon_rarefied21282_meta, aes(x=body.temperature, y=shannon, col=day.time)) +
  geom_point(size = 3) +
  geom_smooth(method ="lm", se = TRUE, colour = "black") +
  scale_x_continuous(name = expression("Body Temperature " ( degree*C)), limits = c(32,40)) +
  scale_y_continuous(name = "Shannon-Wiener Index (H')\n", limits = c(2.5,6.0), breaks = c(2.5,3,3.5,4,4.5,5,5.5,6)) +
  scale_color_manual(values = c("1_0915"="blue", "1_1230"="red", '1_1612'="#04CE13", "1_1740"="purple"),labels = c("morning", "noon", 'early-\nafternoon', "late-\nafternoon") ) +
  labs(color = "Timepoints") + #change the legend title names according to the labels - works for categorical datasets
  theme_bw() +
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size = 13),
        legend.text = element_text(size = 12, colour ="black"), 
        legend.position = "right", axis.title.y = element_text( size = 14), 
        legend.title = element_text(size = 14, colour = "black") #changing the title of legend
  )  + annotate("text", x = 37, y = 5.75, label =  expression(italic(R)^2*"= 0.25; p-value < 0.05")  , size= 4.5)  
shannon_rarefied21282_daytime_templm


## Combine plots into one figure ----

library(cowplot) 
cowplot::plot_grid(ASVrichness_rarefied21282_daytime_finalised,
                   shannon_rarefied21282_daytime_finalised,
                   nrow = 2, ncol = 1, 
                   label_fontface = "plain",
                   label_fontfamily = "Times New Roman",
                   labels = c("A","B"))
ggsave("Alpha_diversity_combined.png", width = 8, height = 8, device = "png", dpi = "print")

cowplot::plot_grid(ASVrichness_rarefied21282_daytime_templm,
                   shannon_rarefied21282_daytime_templm,
                   nrow = 1, ncol = 2,
                   label_fontface = "plain",
                   label_fontfamily = "Times New Roman",
                   labels = c("A","B"))

ggsave("Alpha_diversity_combined_2.png", width = 12, height = 4, device = "png", dpi = "print")







