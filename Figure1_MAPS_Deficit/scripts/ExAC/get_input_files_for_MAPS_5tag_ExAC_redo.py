## Script to compile MAPs input for splicing variants
## /lustre/scratch110/sanger/jl18/Exons5tag/AnnotationTables/Full_ann_tables/  for triplets (words[7])
## Run from in /lustre/scratch110/sanger/jl18/Exons5tag/MAPS/DDD/



import gzip
import argparse

def get_input_files():

    parser = argparse.ArgumentParser()
    parser.add_argument("--chrom", required=True, help="chr of interest")
    args = parser.parse_args()
    return args

args = get_input_files()

CQs = ['frameshift', 'stop_lost', 'stop_gained', 'missense_variant', 'inframe_insertion', 'inframe_deletion', 'initiator_codon_variant']    

infile_anns = open("/lustre/scratch115/teams/hurles/users/jl18/Splicing/Exons5tag/AnnotationTables/Full_ann_tables/Path_score_percentiles/Rmethod_050616/all_5tag_annotations_percentiles_mRate_DDG2P_update.txt")

infile_ACs = open("/lustre/scratch115/teams/hurles/users/jl18/Splicing/DDD8K/MAPS/ExAC_redo/MAPS_input_working1.txt")

Outfile = ''.join(("/lustre/scratch115/teams/hurles/users/jl18/Splicing/DDD8K/MAPS/ExAC_redo/infiles_working/", args.chrom, "_MAPS_ExAC_in.txt"))
outfile = open(Outfile, 'w')

sites = {}

for line in infile_ACs:
    line = line.strip()
    words = line.split('\t')
    if words[0] != args.chrom: continue
    variant = '-'.join((words[0], words[1], words[2], words[3]))
    sites[variant] = words[4]

for line in infile_anns:
    line = line.strip()
    words = line.split('\t')
    if words[0] != args.chrom: continue
    variant = '-'.join((words[0], words[1], words[2], words[3]))
    if variant in sites:
        if words[9] in CQs: continue
        ## chr	pos	ref	alt	AC	splice_CQ	triplet
        outline = ''.join((words[0], '\t', words[1], '\t', words[2], '\t', words[3], '\t', sites[variant],'\t', words[4], '\t', words[7], '\n'))
        outfile.write(outline)
outfile.close()











