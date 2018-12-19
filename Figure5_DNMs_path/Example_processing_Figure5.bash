#!/bin/bash
## Example processing pipeline for Figure5 analyses in DDD Splicing paper ("Pathogenicity and selective constraint on variation near splice sites")
## ROC plots with AUC calculations for set of DNMs in dominant DD genes which recieved a score from all 5 of the pathogenicity metrics (AdaBoost, RandomForest, MaxEntScan, Spidex and CADD)

## Full list of de novo mutations (DNMs) in DDD ~8000 trios was matched against the splicing regions of interest to identify DNMs in splicing regions (for figure 2)
## Only included high quality DNMs with pp_dnm score >=0.8
## All DNMs which were scored by all 5 metrics (AdaBoost, RandomForest, Spidex, MaxEntScan and CADD) were included in this analysis
## See DeNovo_variants_for_ROC_plots.xlsx for full list of variants and scores.
## ROC plots and area under the curve calculated as below.

## Plot ROC curves in RStudio
library(ggplot2)
setwd("U:/Splicing/Exons5tag/DDD8K/MAPS/MAPS_v_path/redo_percentiles_mRate_weighted/de_novos_by_path_scores/plots_with_CADD")

data <- read.table("dom_can_ROCplot_data_with_CADD.txt", header = TRUE)
p <- ggplot(data, aes(Bracket, value, colour=score, shape=score)) + geom_point() + geom_line() + scale_x_reverse(lim=c(21,1)) +
     theme(axis.text=element_text(size=18), axis.title=element_text(size=18), legend.text=element_text(size=18), title=element_text(size=16))
  p + theme(
    panel.background = element_rect(fill = "transparent",colour = NA), 
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_blank(),
    plot.background = element_rect(fill = "transparent",colour = NA))+
  theme(axis.line.x = element_line(color="black", size = 0.2),
        axis.line.y = element_line(color="black", size = 0.2))
	 ggsave("dom_can_curve_plot_overlap_CADD_shape.svg", width = 15, height = 18, units = "cm")

data <- read.table("dom_ROCplot_data_NC_with_CADD.txt", header = TRUE)
p <- ggplot(data, aes(Bracket, value, colour=score, shape=score)) + geom_point() + geom_line() + scale_x_reverse(lim=c(21,1)) +
     theme(axis.text=element_text(size=18), axis.title=element_text(size=18), legend.text=element_text(size=18), title=element_text(size=16))
  p + theme(
    panel.background = element_rect(fill = "transparent",colour = NA), 
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_blank(),
    plot.background = element_rect(fill = "transparent",colour = NA))+
  theme(axis.line.x = element_line(color="black", size = 0.2),
        axis.line.y = element_line(color="black", size = 0.2))
	 ggsave("dom_NC_curve_plot_overlap_CADD_shape.svg", width = 15, height = 18, units = "cm")

## Colours for each pathogenicity score adjusted in Illustrator to reflect the colours used in the MAPS v Path plots for consistency with that work

## AUCs calculated in RStudio - 
setwd("U:/Splicing/Exons5tag/DDD8K/MAPS/MAPS_v_path/redo_percentiles_mRate_weighted/de_novos_by_path_scores/plots_with_CADD")
library(zoo)

## For overlapping set in centiles
data <- read.table("overlapping_20ths_incCADD.txt", header=TRUE)
B <- order(data$Bracket)
N <- data$Null
AUC_N <- sum(diff(B[B])*rollmean(N[B], 2))/2000

DAC <- data$D_Ada_C
DRC <- data$D_RF_C
DMC <- data$D_MESPD_C
DSC <- data$D_Spidex_C
DCC <- data$D_CADD_C
DANC <- data$D_Ada_NC
DRNC <- data$D_RF_NC
DMNC <- data$D_MESPD_NC
DSNC <- data$D_Spidex_NC
DCNC <- data$D_CADD_NC

AUC_DAC <- sum(diff(B[B])*rollmean(DAC[B], 2))/2000
AUC_DANC <- sum(diff(B[B])*rollmean(DANC[B], 2))/2000
AUC_DRC <- sum(diff(B[B])*rollmean(DRC[B], 2))/2000
AUC_DRNC <- sum(diff(B[B])*rollmean(DRNC[B], 2))/2000
AUC_DMC <- sum(diff(B[B])*rollmean(DMC[B], 2))/2000
AUC_DMNC <- sum(diff(B[B])*rollmean(DMNC[B], 2))/2000
AUC_DSC <- sum(diff(B[B])*rollmean(DSC[B], 2))/2000
AUC_DSNC <- sum(diff(B[B])*rollmean(DSNC[B], 2))/2000
AUC_DCC <- sum(diff(B[B])*rollmean(DCC[B], 2))/2000
AUC_DCNC <- sum(diff(B[B])*rollmean(DCNC[B], 2))/2000

AUC_N
AUC_DAC
AUC_DRC
AUC_DMC
AUC_DSC
AUC_DCC
AUC_DANC
AUC_DRNC
AUC_DMNC
AUC_DSNC
AUC_DCNC

