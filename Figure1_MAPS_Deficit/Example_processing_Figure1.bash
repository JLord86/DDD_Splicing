#!/bin/bash
## Example processing pipeline for Figure1 analyses in DDD Splicing paper ("Pathogenicity and selective constraint on variation near splice sites")
## Example given for near splice positions, but also applied to grouped and individual PolyPy positions split by PyPu vs other changes, last base of the exon split by nucleotide, and VEP annotated missense, stop_gained and synonymous variants

######## DDD MAPS:
## Example given for chromosome 1, but was run for chr1-22 in parallel
## Also ran this for variants with VEP annotated synonymous consequences only for use as the baseline in the MAPS analysis
## All scripts in CodeForFigures/Figure1_MAPS_Deficit/DDD

## Prepare target sites files
bsub -o create_sites.txt make_sites_files_for_MAPS_input_DDD8K.py --chrom 1

## Extract allele counts in unaffected DDD parents for sites of interest from multisample VCF using modified version of this code: https://github.com/jeremymcrae/count_singletons
bsub -o run_sing_chr1.txt -q long -J sing1 -R "select[mem>1000] rusage[mem=1000]" -M 1000 python Farm_count_singletons_MAPS_5tag_DDD8K_FS10.py --chrom 1 --vcf /lustre/scratch115/projects/ddd/8ktrio_multisample_perchr_postQC_vcfs/1_1-249250621_postQC.vcf.gz --singletons 5tag_MAPS_chr1_singletons.txt --totals 5tag_MAPS_chr1_totals.txt
## The outfile we are interested in is the one specified within the script, in this instance, named /lustre/scratch115/realdata/mdt1/teams/hurles/users/jl18/Splicing/Exons5tag/DDD8K/MAPS/1_DDD8K_ACs_MAPS1.txt"

## Reformat the output of the above file for input into MAPS code
bsub -o run_process_chr1.txt python get_input_files_for_MAPS_5tag_DDD8K_FS10.py --chrom 1

## Combine the results from chr1-22 and add header for use in MAPS
## Format for the file (and the header itself) should be: chr	pos	ref	alt	allele_count	vep_consequence	context
## Now have data for all splicing and synonymous variants on all chromosomes
## Running MAPS analyses in R (via RStudio) using code from https://github.com/pjshort/dddMAPS:

setwd("~/MAPS/dddMAPS")
source("./MAPS.R")

parental_gencode_snps = read.table("../data/all_DDD_MAPS_in_FS10.txt", header = TRUE, sep = "\t")
synonymous_parental_vars = subset(parental_gencode_snps, vep_consequence == "synonymous_variant")
maps_lm = maps_fit(synonymous_parental_vars)
out = maps_adjust(variants = parental_gencode_snps, split_factor = parental_gencode_snps$vep_consequence, maps_lm = maps_lm, noncoding = FALSE)
write.table(out, "all_MAPS_out_DDD8K_FS10.txt", sep='\t')

######## ExAC MAPS:
## All scripts in CodeForFigures/Figure1_MAPS_Deficit/ExAC
## Match up our sites of interest to the ExAC data from ExAC.r0.3.sites.vep.vcf.gz
bsub -o chr1_ExAC_step1.txt python Find_ExAc_vars_5tag_regions_redo.py --chrom 1

## Process the files into desired format for MAPS input
bsub -o ExAC_step2.txt python get_MAPS_vars_ACs_ExAC_syn_redo.py 
bsub -o chr1_ExAC_step3.txt python get_input_files_for_MAPS_5tag_ExAC_redo.py --chrom 1

## Bring consequence filtering in line with DDD data
bsub -o ExAC_step4.txt python filter_ExAC_CQs.py 

## Bring filtering stringency in line with DDD data (apply FS10 filtering)
bsub -o ExAC_step3.txt filter_properly_filtered_ExAC_data_by_qual.py 

## Combine the results from chr1-22 and add header for use in MAPS
## Format for the file (and the header itself) should be: chr	pos	ref	alt	allele_count	vep_consequence	context
## Running MAPS analyses in R (via RStudio) using code from https://github.com/pjshort/dddMAPS:
setwd("~/MAPS/dddMAPS")
source("./MAPS.R")

parental_gencode_snps = read.table("../data/all_ExAC_MAPS_in_FS10.txt", header = TRUE, sep = "\t")
synonymous_parental_vars = subset(parental_gencode_snps, vep_consequence == "synonymous")
maps_lm = maps_fit(synonymous_parental_vars)
out = maps_adjust(variants = parental_gencode_snps, split_factor = parental_gencode_snps$vep_consequence, maps_lm = maps_lm, noncoding = FALSE)
write.table(out, "all_ExAC_MAPS_out_FS10.txt", sep='\t')
maps_ggplot(names(out$ps_adjusted), out$ps_adjusted, out$standard_error, already_ordered = FALSE)

## Output from DDD and ExAC then reformatted to add desired plotting order and combine with ExAC data
## Plotting done using code: Rcode_plot_Figure1A.R and data MAPS_FS10_combined_to_CI.txt



######## DDD deficit of parental variants in high pLI genes:
## Use the output from the first step of the DDD MAPS processing as input for the deficit work
## All scripts in CodeForFigures/Figure1_MAPS_Deficit/DDD

## Annotate the allele counts (ACs) file with geen and strand
bsub -o run_deficit_st1.txt python annotate_ACs_with_gene_and_strand_FS10.py

## Work out counts per site and high pLI (>0.9) gene counts per site, then calculate the proportion in high pLI genes
bsub -o run_deficit_st2.txt python counts_by_splice_annotation_noNA_FS10.py

## Compile data and plot in R using code: Rcode_plot_Figure1B.R and data Rcode_plot_Figure1B.R

