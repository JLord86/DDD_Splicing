#!/bin/bash
## Example processing pipeline for Figure4 analyses in DDD Splicing paper ("Pathogenicity and selective constraint on variation near splice sites")
## Example given for AdaBoost score with canonical sites included, but also done for RandomForest, MaxEntScan, Spidex and CADD, with and without canonical sites included 
## (For the non-canonical analyses, acc-2, acc-1, don+1 and don+2 sites excluded)

## Scripts in CodeForFigures/Figure4_MAPS_v_path/Scripts 
## Make files of target sites for this analysis - output is chrom	pos	score bracket (1-20)
bsub -o run_MvP_1.txt python make_sites_files_for_sing_v_path_mRate.py

## Extract allele counts in unaffected DDD parents for sites of interest from multisample VCF using modified version of this code: https://github.com/jeremymcrae/count_singletons
## Example given for chr1, but run for chrs 1-22, and for the other pathogenicity scores
bsub -o run_sing_chr1_Ada.txt -q long -J sing1 -R "select[mem>1000] rusage[mem=1000]" -M 1000 python Farm_count_singletons_MAPS_5tag_DDD8K_Ada.py --chrom 1 --vcf /lustre/scratch115/projects/ddd/8ktrio_multisample_perchr_postQC_vcfs/1_1-249250621_postQC.vcf.gz --singletons 5tag_MAPS_chr1_singletons_Ada.txt --totals 5tag_MAPS_chr1_totals_Ada.txt

## Compile the MAPS input information for Ada score brackets - chr	pos	ref	alt	AC	Ada_Bracket	triplet
## Run separately for each chromosome, then concatenate output files > all_Ada_MAPS_in_FS10_syn.txt
bsub -o compile_Ada.txt python compile_MAPS_infiles_Ada_mRate.py

## Run MAPS analysis
setwd("~/MAPS/dddMAPS")
source("./MAPS.R")

## Ada
parental_gencode_snps = read.table("../data/all_Ada_MAPS_in_FS10_syn.txt", header = TRUE, sep = "\t")
synonymous_parental_vars = subset(parental_gencode_snps, vep_consequence == "synonymous_variant")
maps_lm = maps_fit(synonymous_parental_vars)
out = maps_adjust(variants = parental_gencode_snps, split_factor = parental_gencode_snps$vep_consequence, maps_lm = maps_lm, noncoding = FALSE)
write.table(out, "all_Ada_MAPS_v_path_FS10.txt", sep='\t')
maps_ggplot(names(out$ps_adjusted), out$ps_adjusted, out$standard_error, already_ordered = FALSE)

## Plotting - in RStudio
## Data used are in CodeForFigures/Figure4_MAPS_v_path/plotting_data
setwd("U:/Splicing/Exons5tag/DDD8K/Manuscript/GenomeResearch_revision/FS10_MAPS/MAPS_v_path/plotting")
library(ggplot2)

cbPalette4 <- c("#0072B2","#009E73","#56B4E9",   "#CC79A7",  "#D55E00", "#E69F00",  "#F0E442",  "#999999")

data <- read.table("Ada_FS10_to_plot.txt", header = TRUE)
limits <- aes(ymin = data$ps_adjusted - 1.96*data$standard_error, ymax = data$ps_adjusted + 1.96*data$standard_error )
labelsCQ = c("syn","1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20")

plot <-
   ggplot(data, aes(y=ps_adjusted, x=Order, colour=Set)) + 
   geom_point(aes(size = 0.02)) +
   xlab("Score bracket") + ylim(-0.1,0.2) +
   ylab("MAPS") + theme(legend.position = "none")+
   geom_errorbar(limits, width = 0) +
   scale_colour_manual(values = cbPalette4) +
   scale_x_discrete(limits=labelsCQ)+
       theme(axis.text=element_text(size=20), axis.title=element_text(size=20))+
    theme(legend.text=element_text(size=30)) + theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust = 0.9))+
  theme(axis.line.x = element_line(color="black", size = 0.2),
        axis.line.y = element_line(color="black", size = 0.2))
  plot<- plot + theme(
    panel.background = element_rect(fill = "transparent",colour = NA), 
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_blank(),
    plot.background = element_rect(fill = "transparent",colour = NA)) +
    expand_limits(x=c(0,22))
plot
   ggsave("Ada_FS10_maps_v_path.svg", bg="transparent" ,  width = 20, height = 17, units = "cm")


## Correlation - in RStudio
setwd("U:/Splicing/Exons5tag/DDD8K/Manuscript/GenomeResearch_revision/FS10_MAPS/MAPS_v_path/plotting/correlation")

data <- read.delim("Ada_FS10_to_cor.txt", header=T)  
 cor(data$CQ, data$ps_adjusted, method = "spearman")


