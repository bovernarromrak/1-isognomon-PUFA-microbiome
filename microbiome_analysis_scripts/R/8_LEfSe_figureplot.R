#Use the alpha0.05_LDAscore2.0_alltoall data produced from the Galaxy Huttenhower website

### library ===============================================================
library(ggplot2)
library(dplyr)


## the original pathway name - MetaCyc
alpha0.05_LDAscore2.0_alltoall <- read.csv("alpha0.05_LDAscore2.0_alltoall_simplified.csv",header=T, sep=",")
version1 <- alpha0.05_LDAscore2.0_alltoall 
version1$pathways<- as.factor(version1$pathways)
version1$pathways <- factor(version1$pathways, levels = c("PWY-7323", "PWY-7332",
  "PWY-5896", "PWY-5850", "PWY-6897", "PWY-5860", "PWY-5529", "PWY-7159", "PWY-5531", "PWY0-845",
  "P441-PWY", "PWY-7255", "PWY-7084", "PWY1G-0", "PWY-6396", "GLCMANNANAUT-PWY", "PWY-7090",
  "RHAMCAT-PWY", "CHLOROPHYLL-SYN", "PWY-6728", "PWY-7187", "ANAEROFRUCAT-PWY", "GLUCARGALACTSUPER-PWY",
  "GALACTARDEG-PWY", "DENITRIFICATION-PWY", "GLUCARDEG-PWY", "PWY-7295"))
version1$perfect_names <- factor(version1$perfect_names, levels =   
c("Superpathway of GDP-mannose-derived O-antigen building blocks biosynthesis", "Superpathway of UDP-N-acetylglucosamine-derived O-antigen building blocks biosynthesis", "Superpathway of menaquinol-10 biosynthesis","Superpathway of menaquinol-6 biosynthesis", "Thiamine diphosphate salvage II",
  "Superpathway of demethylmenaquinol-6 biosynthesis I", "Superpathway of bacteriochlorophyll a biosynthesis", "3,8-divinyl-chlorophyllide a biosynthesis III (aerobic, light independent)", "3,8-divinyl-chlorophyllide a biosynthesis II (anaerobic)",
  "Superpathway of pyridoxal 5'-phosphate biosynthesis and salvage", "Superpathway of N-acetylneuraminate degradation", "Ergothioneine biosynthesis I (bacteria)", "Nitrifier denitrification", "Mycothiol biosynthesis",
"Superpathway of 2,3-butanediol biosynthesis", "Superpathway of N-acetylglucosamine, N-acetylmannosamine and N-acetylneuraminate degradation", "UDP-2,3-diacetamido-2,3-dideoxy-&alpha;-D-mannuronate biosynthesis", "L-rhamnose degradation I",
"3,8-divinyl-chlorophyllide a biosynthesis I (aerobic, light-dependent)", "Methylaspartate cycle", "Pyrimidine deoxyribonucleotides de novo biosynthesis II", "Homolactic fermentation", "Superpathway of D-glucarate and D-galactarate degradation",
"D-galactarate degradation I", "Nitrate reduction I (denitrification)", "D-glucarate degradation I", "L-arabinose degradation IV"))
version1$representation <- factor(version1$representation, levels = c("morning", "noon", "early-afternoon", "late-afternoon"))

version1_plot <-ggplot(version1) +
  geom_bar(aes(x=perfect_names,y=lda_score_log10,fill=representation),stat="identity",position = position_stack(reverse = TRUE)) +
  #ylim(0,3.2) +
  theme_bw() +
  theme(panel.grid.major=element_line(colour="white"),
        panel.grid.minor=element_line(colour="white"),
        strip.background = element_rect(fill="white"),
        axis.text.x = element_text(size = 12), 
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 9), 
        axis.title.y = element_text(size = 14))  +
  guides(fill = guide_legend(title = "Timepoints", keywidth = 1, keyheight = 1)) +
  scale_x_discrete(limits = rev) + #to reverse a categorical variable on axis
  xlab("MetaCyc Pathways") +
  geom_hline(yintercept = 2, alpha = 0.5, linetype = 2) +
  ylab("LDA scores (log 10)") +
  scale_fill_discrete() +
  coord_flip()

version1_plot

#ggsave(filename = "LEfSe_result.png" ,height = 7, width = 12, dpi = "print")

