## Script to find ExAc varaints (from ExAC.r0.3.sites.vep.vcf.gz) that coincide with our regions of interest
## Run from in /lustre/scratch110/sanger/jl18/Exons5tag/AnnotationTables/ExAC_5tag_region/
## No VEP CQ filtering at this stage - output all overlapping positions


import gzip
import argparse


def get_options():
    parser = argparse.ArgumentParser(description="Get ExAC variants that overlap with the 5tag filter regions of interest.")
    parser.add_argument("--chrom", required=True, help="chromosome to investigate.")
    args = parser.parse_args()
    return args

args = get_options()

infile_ExAc = gzip.open("/software/hgi/resources/ExAC_release/release0.3/ExAC.r0.3.sites.vep.vcf.gz")   ## ExAC variants file

ann_dict = {}
site_dict = {}
infile_pos = open("/lustre/scratch115/realdata/mdt1/teams/hurles/users/jl18/Splicing/Exons5tag/AnnotationTables/Full_ann_tables/Path_score_percentiles/Rmethod_050616/all_5tag_annotations_percentiles_mRate.txt")

Outfile = ''.join(("/lustre/scratch115/teams/hurles/users/jl18/Splicing/DDD8K/MAPS/ExAC_redo/", args.chrom, "_ExAC_5tag_region_noCQfilt.txt"))
outfile = open(Outfile, 'w')
for line in infile_pos:           ## For each line in the splice sites file, a dictionary entry is made containing the chromosome (i) and position as the key and the site (e.g. acc) as the value
    line = line.strip()
    words = line.split('\t')
    if words[0] != args.chrom: continue
    variant = '-'.join((words[0], words[1]))
    site_dict[variant] = words[4]
    ann_dict[variant] = words[9]

for line in infile_ExAc:           ## Goes through the ExAC sites file and prints out any variants (full line, plus our site annotation) that match our splicing sites of interest.
    line = line.strip()
    if line.startswith('#'): continue
    words = line.split('\t')
    if words[0] != args.chrom: continue
    variant = '-'.join((words[0], words[1]))
    if variant in site_dict:
        outline = ''.join((line, '\t', site_dict[variant], '\t', ann_dict[variant], '\n'))
        outfile.write(outline)
outfile.close()















