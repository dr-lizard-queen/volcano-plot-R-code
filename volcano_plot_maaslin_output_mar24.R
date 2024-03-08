# anna neustaeter
## volcano plot for papers/presentations
### March 24

#installation of libraries
#if (!requireNamespace('BiocManager', quietly = TRUE))
#  install.packages('BiocManager')
#BiocManager::install('EnhancedVolcano')

library(MASS)
library(foreign)
library(ggplot2)
library(tidyr)
library(ggsignif)
library(plyr)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(Maaslin2)
library(EnhancedVolcano)

##directory
setwd("PATH/TO/DIRECTORY/maaslin_output")

##https://forum.biobakery.org/t/trying-to-understand-coef-column-and-how-to-convert-it-to-fold-change/3136
#from a biobakery forum, indicate that 
#'log2foldchange' IS coef from a maaslin output

#load data
res_all_2weeks <- read.csv("maaslin_2weeks_concatname_mar24.csv", header = TRUE, 
                    as.is = TRUE, na.strings=c(""," ","NA"), strip.white=TRUE)

#convert second column to row names
rownames(res_all_2weeks) <- res_all_2weeks[,2]
glimpse(res_all_2weeks)

#turning off scientific notation
options(scipen = 999)


# create custom key-value pairs for different interventions
#LF7 <- c('LP7')
LF7_FOS <- c('FOS_')

keyvals.shape <- ifelse(grepl(LF7_FOS, rownames(res_all_2weeks)), 16, 17)

names(keyvals.shape)[keyvals.shape == 17] <- 'LP7'
names(keyvals.shape)[keyvals.shape == 16] <- 'LP7 + FOS'


###########################################################################
#find minimum and maximum of coefficients
minimum<-min(res_all_2weeks$coef, na.rm=T)
maximum<-max(res_all_2weeks$coef, na.rm=T)

#name output file and type
pdf("volcanoplot_maaslin_2weeks_ .pdf", width=12, height=15)

p2 <- EnhancedVolcano(res_all_2weeks, #file name
                      lab = as.character(res_all_2weeks$feature), #labels
                      x = 'coef', #column for X axis
                      y = 'qval', #column for y-axis
                      xlab = bquote('coefficient'),
                      ylab = bquote('-log10(Q-value)'),
                      xlim = c(minimum-1, maximum+1),
                      ylim = c(-1,40),
                      title = 'Maaslin output: LP7 and LP7+FOS vs placebo',
                      subtitle = '2 weeks',
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

