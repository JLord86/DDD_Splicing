#!/bin/bash
## Example processing pipeline for Figure2 analyses in DDD Splicing paper ("Pathogenicity and selective constraint on variation near splice sites")

## Figure 2A - observed and expected de novo mutations across splicing region and Figure 2B - grouped Polypyrimidine tract sites
## Scripts in CodeForFigures/Figure2_DNMs/scripts
## Expected:
## Calculating expected number of de novo mutations based on null model of mutation rate from Samocha et al., adjusting for observed depth of coverage
bsub -o run_exp1.txt python calculate_cov_adj_expected_chr.py --chrom 1
## Run for each chromosome (1-22), then output compiled.
## Calculated separately for DDG2P dominant, recessive, nonDDG2P and all genes

## Observed:
## Find observed de novos that overlap with region of interest
bsub -o run_obs.txt python get_observed_splice_denovos_and_cov_check.py
## Output contained all quality scores - filtered after this to remove anything with ppdmn <0.8
## Added back in a few recurrent sites between individuals where only one sample's DNM was in the output
## Output annotated with whether the gene is in DDG2P as dominant or recessive, or not in the DDG2P gene list (nonDDG2P)

## Poisson test + FDR correction
## Data in CodeForFigures/Figure2_DNMs/stats
## Poisson test used to compare the observed with expected values in R using the poisson.test function.
## Example given below for canonical splice donor site in dominant DD genes, following the format "poisson.test(observed,expected) # set position"
poisson.test(21, 1.04700011912) # dom don+1
poisson.test(5, 0.500733379526) # dom don+2

## To control for multiple testing, an FDR correction was applied at 5% to all p-values
## (57 positions + combined + combined without canonical, in 5 states - dominant, recessive, nonDDG2P, diagnosed and undiagnosed, plus the polypyrimidine split tests in the dominant, recessive and nonDDG2P groups - total 303 tests)
## In RStudio:
data <- read.table("all_DNMs_paper_p_values_to_FDR_correct.txt", header=T)
pvals <- data$poisson_p
valAdj = p.adjust(pvals, method= 'fdr')
data1 <- cbind(data, valAdj)
write.table(data1, "FDR_adj_pvals_splice_all_in_paper.txt", sep='\t')

## Plotting (observed + expected, and FDR corrected Poisson p-values)
## For code and data used in Figure 2A and 2B, see R_code_obs_exp_and_poisson_plots.txt in CodeForFigures/Figure2_DNMs/plotting

## Figure 2C - Positive Predictive Values (PPVs) for different classes of variants
## PPV calculation:
## PPV = (observed-expected)/observed
## Full values for all dominant genes given in CodeForFigures/Figure2_DNMs/plotting/PPVs_data_dominant_genes.txt, derived from the observed and expected values used in Figure2A

## Plotting for Figure 2C in RStudio:
library(ggplot2)
data <- read.table("PPVs_to_plot.txt", header=T) ## in CodeForFigures/Figure2_DNMs/plotting
p <- ggplot(data, aes(x=Order_dom, y=PPV, fill=PPV))
p <- p+ geom_bar(stat="identity") + theme(axis.text=element_text(size=15), axis.title=element_text(size=15)) + theme(legend.text=element_text(size=15)) +
     scale_x_discrete(name = "Annotation", limits=c("Stop\ngained", "Missense", "", "Canonical\nacceptor", "acc-2", "acc-1", "", "Canonical\ndonor", "don+1", "don+2", "", "don+5", "PPT\nPyPu", "Exonic\nnear splice", "PPT\nother", "Other\nnear splice","")) 
   p + theme(
    panel.background = element_rect(fill = "transparent",colour = NA), 
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_blank(),
    plot.background = element_rect(fill = "transparent",colour = NA))
ggsave("PPV_by_splice_annotation_test.svg", width = 50, height = 15, units = "cm")
## Modified in Adobe Illustrator to all be same colour, and to combine with other panels of Figure 2

## Figure 2D - Fold enrichment (observed/expected) for don+5, PolyPy PyPu and synonymous variants split by ascending pLI
## Scripts in CodeForFigures/Figure2_DNMs/scripts
## Example given for don+5 processing.

## Calculating expected rates with coverage adjustment (per chromosome 1-22):
bsub -o run_don5_1.txt python calculate_cov_adj_expected_chr_don+5_pLI_sextiles.py --chrom 1

## Compiling expected rates:
bsub -o run_don5.txt python compile_splice_mRates_don5_pLI_sextiles.py

## Observed de novo mutations by pLI sextile:
## Observed de novo mutations identified for Figure 2A annotated with pLI bracket, and observed/expected calculated for each bracket using the expected values above.

## Data in CodeForFigures/Figure2_DNMs/plotting/obs_exp_d5_PPT_syn_pLI_to_plot.txt
## Rcode for plotting in CodeForFigures/Figure2_DNMs/plotting/Rcode_Fig2D_Enrichment_by_pLI.txt

