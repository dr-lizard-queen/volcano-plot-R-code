# anna neustaeter
## volcano plot for papers/presentations
### March 24
```R 
library(MASS)
library(epiDisplay)
library(foreign)
library(ggplot2)
library(tidyr)
library(ggsignif)
library(plyr)
library(dplyr)
library(tableone)
library(geepack)
library(tidyverse)
library(ggpubr)
library(gtsummary)
library(broom)
library(jtools)
library(broom.mixed)
library(effects)
library(rstatix)
library(survminer)
library(jtools)
library(Matrix)
library("pspearman")
library(moonBook)
library(Epi)
library(ztable)
library(ROCR)
library("cowplot")
library(standardize)
library(clogitL1)
library(beeswarm)
library(microbiome)
library(vegan)
library(Maaslin2)
library(magrittr)

#if (!requireNamespace('BiocManager', quietly = TRUE))
#  install.packages('BiocManager')
#BiocManager::install('EnhancedVolcano')

library(EnhancedVolcano)
```

##directory
```R
setwd("C:/Users/Anna neustaeter/OneDrive - SickKids/maaslin_output")
```
##https://forum.biobakery.org/t/trying-to-understand-coef-column-and-how-to-convert-it-to-fold-change/3136
#from a biobakery forum, indicate that 
#'log2foldchange' IS coef from a maaslin output


#load data

res_all <- read.csv("maaslin_2weeks.csv", header = TRUE, 
                as.is = TRUE, na.strings=c(""," ","NA"), strip.white=TRUE)


res <- subset(res_all, value=="LP7")

rownames(res) <- res[,1]
glimpse(res)
str(res)

#turning off scientific notation
options(scipen = 999)

p <-  EnhancedVolcano(res,
                      x = 'coef',
                      y = 'qval',
                      title = 'Maaslin output: LP7 at vs placebo',
                      subtitle = '2 weeks',
                      xlab = bquote(~Log[2]~ 'fold change'),
                      ylab = bquote('-log10(Q-value)'),
                      pCutoff = 0.05,
                      FCcutoff = 0.5,
                      pointSize = 2.0,
                      labSize = 5.0,
                      lab = as.character(res$feature),
                      cutoffLineCol = 'grey34',
                      drawConnectors = TRUE,
                      maxoverlapsConnectors = Inf,
                      widthConnectors =0.5,
                      colConnectors = 'black',
                      lengthConnectors = unit(7, "npc"),
                      #border = "full",
                      col=c('grey44', 'cyan4', 'maroon', 'mediumorchid3'),
                      legendLabels=c('Not sig.','Log (base 2) FC','q-value',
                                     'q-value & Log (base 2) FC'),
                      legendPosition = 'bottom',
                      legendLabSize = 10,
                      legendIconSize = 5.0,
                      arrowheads = FALSE,
                      gridlines.major = FALSE,
                      gridlines.minor = FALSE,
                      colAlpha = 0.8,
                      max.overlaps=Inf)
p 

p + coord_flip()

dev.off()


# define different interventions
res_all_2weeks <- read.csv("maaslin_2weeks_concatname_mar24.csv", header = TRUE, 
                    as.is = TRUE, na.strings=c(""," ","NA"), strip.white=TRUE)


rownames(res_all_2weeks) <- res_all_2weeks[,2]
glimpse(res_all_2weeks)

str(res_all_2weeks)
#LF7 <- c('LP7')
LF7_FOS <- c('FOS_')

# create custom key-value pairs for different interventions

keyvals.shape <- ifelse(grepl(LF7_FOS, rownames(res_all_2weeks)), 16, 17)
#keyvals.shape[is.na(keyvals.shape)] <- 3
names(keyvals.shape)[keyvals.shape == 17] <- 'LP7'
names(keyvals.shape)[keyvals.shape == 16] <- 'LP7 + FOS'


###########################################################################
minimum<-min(res_all_2weeks$coef, na.rm=T)
maximum<-max(res_all_2weeks$coef, na.rm=T)
pdf("volcanoplot_maaslin_2weeks_ .pdf", width=12, height=15)

p2 <- EnhancedVolcano(res_all_2weeks,
                      lab = as.character(res_all_2weeks$feature),
                      x = 'coef',
                      y = 'qval',
                      xlab = bquote('coefficient'),
                      ylab = bquote('-log10(Q-value)'),
                      xlim = c(minimum-1, maximum+1),
                      ylim = c(-1,40),
                      title = 'Maaslin output: LP7 and LP7+FOS vs placebo',
                      subtitle = '2 weeks',
                      #subtitle = '2 weeks',
                      pCutoff = 0.05, #significance threshold
                      FCcutoff = 1.5, # coefficient threshold
                      pointSize = 3.4, #size of points in figure
                      labSize = 4.4, #size of label in figure
                      shapeCustom = keyvals.shape, #calls to custom figure shapes
                      colAlpha = 0.5, #transparency of points
                      col=c('grey44', 'cyan4', 'maroon', 'mediumorchid3'),
                      legendLabels=c('Not sig.','coef>1.5','q-value>0.05',
                                     'q-value & coef>1.5'), #labels for colours
                      legendPosition = 'right', #legend position
                      legendIconSize = 5.0,
                      legendLabSize = 10,
                      arrowheads = FALSE, #if lines should have arrow heads
                      drawConnectors = TRUE,
                      maxoverlapsConnectors = Inf,
                      widthConnectors = 0.5,
                      min.segment.length = 0.05,
                      colConnectors = 'grey5',
                      gridlines.major = FALSE,
                      gridlines.minor = FALSE,
                      borderWidth = 1.0,
                      borderColour = 'black')

p2+ theme(axis.title.x = element_text(size = 15)) +
  theme(axis.text.y = element_text(size = 15)) 



dev.off()

###############################################################################
#from: https://www.nature.com/articles/s41467-023-43688-z#MOESM1 figure 6b
#https://github.com/nathanvannierinrae/Scripts-from-Vannier-et-al.-2023/blob/main/figureS3/script_volcano_plots.R
colnames(res_all_2weeks)
notable_genera <- res_all_2weeks[with(res_all_2weeks, qval <= 0.05 & abs(coef <= 2) ), ]
unique(notable_genera$feature)

strains=c( "Alcaligenes_ammonioxydans" ,    "Alcaligenes_faecalis"        ,  "Arthrobacter_mobilis"         ,
           "Atlantibacter_subterranea" ,    "Bradyrhizobium_frederickii"  ,  "Buttiauxella_agrestis"        ,
           "Citrobacter_amalonaticus"  ,    "Edwardsiella_piscicida"      ,  "Enterococcus_avium"           ,
            "Enterococcus_gallinarum"  ,     "Erwinia_mediterraneensis"   ,   "Limosilactobacillus_fermentum",
            "Morganella_morganii"      ,     "Pasteurella_canis"          ,   "Providencia_huaxiensis"       ,
            "Pseudomonas_alloputida"   ,     "Rahnella_inusitata"         ,   "Serratia_ficaria"             ,
            "Siccibacter_colletis"     ,     "Streptococcus_vicugnae"     ,   "Vibrio_hepatarius"            ,
            "Yersinia_enterocolitica"  ,     "Yersinia_pseudotuberculosis"  )
strains=as.factor(strains)

   
    

min(res_all_2weeks$qval)    
    #create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
    # set the base colour as 'black'
    keyvals.colour <- rep('black', nrow(res_all_2weeks))
    
    # set the base name/label as 'Mid'
    names(keyvals.colour) <- rep('Not enriched', nrow(res_all_2weeks))
    
    # modify keyvals for variables with fold change > 2.5
    keyvals.colour[which(res_all_2weeks$coef > 2 & res_all_2weeks$qval<0.05)] <- 'forestgreen'
    names(keyvals.colour)[which(res_all_2weeks$coef > 2 & res_all_2weeks$qval<0.05)] <- 'Increased abundance'
    
    # modify keyvals for variables with fold change < -2.5
    keyvals.colour[which(res_all_2weeks$coef < -2)] <- 'darkorange3'
    names(keyvals.colour)[which(res_all_2weeks$coef < -2 & res_all_2weeks$qval<0.05)] <- 'Decreased adundance'
    
    
    EnhancedVolcano(res_all_2weeks,
                    lab = as.character(res_all_2weeks$feature),
                    x = 'coef',
                    y = 'qval',
                    xlab = bquote('coefficient'),
                    ylab = bquote('-log10(Q-value)'),
                    xlim = c(minimum-1, maximum+1),
                    ylim = c(0,40),
                    #title = paste(levels(strains)[i]) ,
                    title = 'Maaslin output: LP7 and LP7+FOS vs placebo',
                    subtitle = '2 weeks',
                    pCutoff = 0.05,
                    FCcutoff = 2,
                    pointSize = 2.5,
                    labSize = 5,
                    labCol = "black",
                    cutoffLineType = "blank",
                    vline=c(-2, 2),
                    caption = "coef cutoff, 2 ;q-value cutoff, 0.05",
                    legendPosition = "right",
                    legendLabSize = 12,
                    colAlpha = 0.5, 
                    drawConnectors = TRUE,
                    maxoverlapsConnectors = Inf,
                    arrowheads = FALSE,
                    colCustom = keyvals.colour,
                    shapeCustom = keyvals.shape, #calls to custom figure shapes
                    hline = c(0.05))
    dev.off()
    



################################################################################


res_all_2months <- read.csv("maaslin_2months_concatname_mar24.csv", header = TRUE, 
                    as.is = TRUE, na.strings=c(""," ","NA"), strip.white=TRUE)


rownames(res_all_2months) <- res_all_2months[,2]
glimpse(res_all_2months)

str(res_all_2months)
#LF7 <- c('LP7')
LF7_FOS <- c('FOS_')

# create custom key-value pairs for different interventions

keyvals.shape <- ifelse(grepl(LF7_FOS, rownames(res_all_2months)), 16, 17)
#keyvals.shape[is.na(keyvals.shape)] <- 3
names(keyvals.shape)[keyvals.shape == 17] <- 'LP7'
names(keyvals.shape)[keyvals.shape == 16] <- 'LP7 + FOS'


pdf("volcanoplot_maaslin_2months .pdf", width=10, height=10)
#pdf("volcanoplot_maaslin_2weeks .pdf", width=10, height=10)
p2 <- EnhancedVolcano(res_all_2months,
                      lab = as.character(res_all_2months$feature),
                      x = 'coef',
                      y = 'qval',
                      xlab = bquote('coefficient'),
                      ylab = bquote('-log10(Q-value)'),
                      title = 'Maaslin output: LP7 and LP7+FOS vs placebo',
                      subtitle = '2 months',
                      #subtitle = '2 weeks',
                      pCutoff = 0.05, #significance threshold
                      FCcutoff = 1.0, # coefficient threshold
                      pointSize = 5.5, #size of points in figure
                      labSize = 5.0, #size of label in figure
                      shapeCustom = keyvals.shape, #calls to custom figure shapes
                      colAlpha = 0.9, #transparency of points
                      col=c('grey44', 'cyan4', 'maroon', 'mediumorchid3'), #colour of points, given thresholds above
                      legendLabels=c('Not sig.','coefficient > 1.0','q-value>0.05',
                                     'q-value > 0.05 & coefficient > 1.0'), #labels for colours
                      legendPosition = 'bottom', #legend position
                      legendIconSize = 5.0,
                      legendLabSize = 10,
                      arrowheads = FALSE, #if lines should have arrow heads
                      drawConnectors = TRUE,
                      maxoverlapsConnectors = Inf,
                      widthConnectors = 0.5,
                      colConnectors = 'grey5',
                      gridlines.major = FALSE,
                      gridlines.minor = FALSE,
                      borderWidth = 1.0,
                      borderColour = 'black')

p2 + theme(axis.title.x = element_text(size = 20)) +
  theme(axis.text.y = element_text(size = 20)) 


dev.off()

dev.off()

