## Comparing observed and expected numbers of de novo mutations across ascending pLI sextiles
## Script to calculate expected values using Liu's coverage correction
## 
## for depth >= 50, no adjustment needed (i.e. mRate*(individuals>50 * 2)
## for depth <= 1, exp*0.119 (i.e. mRate*0.119*(individuals<1*2))
## for depth >1 and <50, per individual:  exp*(0.119 + 0.204(log(depth)))*2, summed across all individuals

import math 
import argparse

def get_options():
    """get the command line options"""
    parser = argparse.ArgumentParser(description="###")
    parser.add_argument("--chrom", required=True, help="# Get chromosome to investigate (1-22, X)")
    args = parser.parse_args()
    return args
args = get_options()

## Coverage files
coverage_counts = open("/lustre/scratch115/realdata/mdt1/teams/hurles/users/jl18/Splicing/DDD8K/Denovos/Liu_coverage_correction/coverage_summary_counts.txt")
coverage_values = open("/lustre/scratch115/realdata/mdt1/teams/hurles/users/jl18/Splicing/DDD8K/Denovos/Liu_coverage_correction/coverage_values_by_exon.txt")

## Annotation table
annotation_table = open("/lustre/scratch115/realdata/mdt1/teams/hurles/users/jl18/Splicing/Exons5tag/AnnotationTables/Full_ann_tables/Path_score_percentiles/Rmethod_050616/all_5tag_annotations_percentiles_mRate.txt")

## Output files
Error_file = ''.join(("/lustre/scratch115/realdata/mdt1/teams/hurles/users/jl18/Splicing/DDD8K/Denovos/pLI_5tag/splice_exp/sites_with_cov_not_found_", str(args.chrom), "_sextiles.txt"))
error_file =  open(Error_file, 'w')
Outfile = ''.join(("/lustre/scratch115/realdata/mdt1/teams/hurles/users/jl18/Splicing/DDD8K/Denovos/pLI_5tag/splice_exp/cov_adj_expected_full_", str(args.chrom), "_sextiles.txt"))
outfile = open(Outfile, 'w')

## Set up dictionary of pLI scores to separate mRates by grouped pLI
pLI_file = open("/lustre/scratch115/realdata/mdt1/teams/hurles/users/jl18/Splicing/DDD8K/Denovos/pLI_5tag/cleaned_exac_with_pLI_march16_sextiles.txt")
pLI_dict = {}

for line in pLI_file:
    line = line.strip()
    words = line.split('\t')
    pLI_dict[words[1]] = words[21]

## Make dictionaries with coverage counts and values stored for acc and don positions
low_dict = {}
mid_dict = {}
high_dict = {}
values_dict = {}

for line in coverage_counts:
    line = line.strip()
    words = line.split('\t')
    if line.startswith('X'): continue
    get_chrom = words[0].split(':')
    chrom = get_chrom[0]
    if str(chrom) != str(args.chrom): continue
    if '-' not in words[0]:
        start = words[0]
        end = words[0]
    else:
        get_sites = get_chrom[1].split('-')
        start = ''.join((str(chrom), ":", str(get_sites[0])))
        end = ''.join((str(chrom), ":", str(get_sites[1])))
    low_dict[start] = words[2]
    low_dict[end] = words[2]
    mid_dict[start] = words[3]
    mid_dict[end] = words[3]
    high_dict[start] = words[4]
    high_dict[end] = words[4]

for line in coverage_values:
    line = line.strip()
    words = line.split('\t')
    if line.startswith('X'): continue
    get_chrom = words[0].split(':')
    chrom = get_chrom[0]
    if str(chrom) != str(args.chrom): continue
    if '-' not in words[0]:
        start = words[0]
        end = words[0]
    else:
        get_sites = get_chrom[1].split('-')
        start = ''.join((str(chrom), ":", str(get_sites[0])))
        end = ''.join((str(chrom), ":", str(get_sites[1])))
    values = words[2:]
    values_dict[start] = values
    values_dict[end] = values

## Dictionaries for running totals of mutation rates in 6 groups ascending (G1-G6)
pLI_running = {"G1": 0, "G2": 0, "G3": 0, "G4": 0,"G5": 0, "G6": 0,  "NA": 0,}

## Go through annotation table and use strand and distance from exon to calculate the coordinates of the nearest splice site to use to look up the coverage values
for line in annotation_table:
    line = line.strip()
    words = line.split('\t')
    if line.startswith('X'): continue
    if words[0] != str(args.chrom): continue
    get_distance = words[4].split('-')
    distance = get_distance[0]
    strand = words[6]
    mRate = words[8]
    gene = words[5]

    if words[4] != "don+5": continue ## We're just interested in don+5 sites here
## Work out the position of the nearest splice site to look up
    if words[4] == "acc" or words[4] == "don":
        working_coor = ''.join((words[0], ":", words[1]))
    if strand == "+":
        if '+' in words[4]:
            get_adjustment = words[4].split('+')
            coor = int(words[1]) - int(get_adjustment[1])
            working_coor = ''.join((words[0], ":", str(coor)))
        if "-" in words[4]:
            get_adjustment = words[4].split('-')
            coor = int(words[1]) + int(get_adjustment[1])
            working_coor = ''.join((words[0], ":", str(coor)))
    if strand == "-":
        if '+' in words[4]:
            get_adjustment = words[4].split('+')
            coor = int(words[1]) + int(get_adjustment[1])
            working_coor = ''.join((words[0], ":", str(coor)))
        if "-" in words[4]:
            get_adjustment = words[4].split('-')
            coor = int(words[1]) - int(get_adjustment[1])
            working_coor = ''.join((words[0], ":", str(coor)))
## Work out the pLI class each gene (with a pLI score) falls in to
    pLI = "NA"
    pLI_class = "NA"
    if gene in pLI_dict:
        pLI = pLI_dict[gene]
    if pLI != "NA":
        pLI_class = pLI

## Lookup the coordinate in the tables above and calculate the expected value, adjusting appropriately
    if working_coor not in high_dict: 
        outline = ''.join((working_coor, '\t', line, '\n'))
        error_file.write(outline)
    if working_coor not in high_dict: continue

## High coverage (no adjustment needed, just multiply by individuals and by 2 (for alleles))
    high_individuals = high_dict[working_coor]
    mRate_current = float(mRate) * float(high_individuals) * 2
    pLI_running[pLI_class] = float(pLI_running[pLI_class]) + float(mRate_current) ## pLI group total
## Low coverage (mRate*0.119*(individuals<1*2))
    low_individuals = low_dict[working_coor]
    mRate_current = float(mRate) * 0.119 * (float(low_individuals) * 2)
    pLI_running[pLI_class] = float(pLI_running[pLI_class]) + float(mRate_current) ## pLI group total
## Middle coverage (exp*(0.119 + 0.204(log(depth)))*2, summed across all individuals)
    for i in values_dict[working_coor]:
        mRate_current = float(mRate) * (0.119 + 0.204*(math.log(float(i)))) * 2
        pLI_running[pLI_class] = float(pLI_running[pLI_class]) + float(mRate_current) ## pLI group total

print "starting output"
for i in pLI_running:
    outline = ''.join((str(args.chrom), '\t', i, '\t',  str(pLI_running[i]), '\n'))
    outfile.write(outline)
outfile.close()
error_file.close()
