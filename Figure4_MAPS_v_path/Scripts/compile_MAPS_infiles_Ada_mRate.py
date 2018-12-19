## Script to compile MAPs input for splicing variants
## Run after Sing_v_path_count_sing_Ada.py has generated the AC files.
## /lustre/scratch110/sanger/jl18/Exons5tag/AnnotationTables/Full_ann_tables/  for triplets (words[7])
##  
## Run from in /lustre/scratch110/sanger/jl18/Exons5tag/Path_scores/MAPS_v_path

import gzip
import argparse

def get_input_files():

    parser = argparse.ArgumentParser()
    parser.add_argument("--chrom", required=True, help="chr of interest")
    args = parser.parse_args()
    return args

args = get_input_files()


Infile_anns = ''.join(("/lustre/scratch110/sanger/jl18/Exons5tag/AnnotationTables/Full_ann_tables/Path_score_percentiles/Rmethod_050616/", "all_5tag_annotations_percentiles_mRate.txt"))
infile_anns = open(Infile_anns)

Infile_ACs = ''.join(("/lustre/scratch109/sanger/jl18/DDD8K/MAPS/MAPS_v_path/", args.chrom, "_DDD8K_ACs_MAPS1_Ada.txt"))
infile_ACs = open(Infile_ACs)                                                  

Outfile = ''.join(("/lustre/scratch109/sanger/jl18/DDD8K/MAPS/MAPS_v_Path_mRate/Compile/MAPS_in_Ada/", args.chrom, "_Ada_MAPS_in.txt"))
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
        ## chr	pos	ref	alt	AC	Ada_Bracket	triplet
        outline = ''.join((words[0], '\t', words[1], '\t', words[2], '\t', words[3], '\t', sites[variant],'\t', words[20], '\t', words[7], '\n'))
        outfile.write(outline)
outfile.close()











