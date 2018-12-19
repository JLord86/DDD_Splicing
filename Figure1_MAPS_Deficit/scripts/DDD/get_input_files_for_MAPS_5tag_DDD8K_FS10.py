## Script to compile MAPs input for splicing variants
## Run after Farm_count_singletons_MAPS_5tag_DDD8K.py has generated the AC files.
##  

import gzip
import argparse

def get_input_files():

    parser = argparse.ArgumentParser()
    parser.add_argument("--chrom", required=True, help="chr of interest")
    args = parser.parse_args()
    return args

args = get_input_files()

CQs_dict = {"synonymous_variant": 0, "intron_variant": 0, "splice_region_variant": 0, "splice_donor_variant": 0, "splice_acceptor_variant": 0, }

infile_anns = open("/lustre/scratch115/realdata/mdt1/teams/hurles/users/jl18/Splicing/Exons5tag/AnnotationTables/Full_ann_tables/Path_score_percentiles/Rmethod_050616/all_5tag_annotations_percentiles_mRate.txt")

Infile_ACs = ''.join(("/lustre/scratch115/realdata/mdt1/teams/hurles/users/jl18/Splicing/GenRes_redo/MAPS/DDD_qual/FS10/", args.chrom, "_DDD8K_ACs_MAPS1.txt"))
infile_ACs = open(Infile_ACs)

Outfile = ''.join(("/lustre/scratch115/realdata/mdt1/teams/hurles/users/jl18/Splicing/GenRes_redo/MAPS/DDD_qual/MAPS_IN/", args.chrom, "_MAPS_in_FS10.txt"))
outfile = open(Outfile, 'w')

sites = {}

for line in infile_ACs:
    line = line.strip()
    words = line.split('\t')
    variant = '-'.join((words[0], words[1], words[2], words[3]))
    sites[variant] = words[4]

for line in infile_anns:
    line = line.strip()
    words = line.split('\t')
    if words[0] != args.chrom: continue
    variant = '-'.join((words[0], words[1], words[2], words[3]))
    if variant in sites:
        if words[9] not in CQs_dict: continue
        ## chr	pos	ref	alt	AC	splice_CQ	triplet
        outline = ''.join((words[0], '\t', words[1], '\t', words[2], '\t', words[3], '\t', sites[variant],'\t', words[4], '\t', words[7], '\n'))
        outfile.write(outline)
outfile.close()











