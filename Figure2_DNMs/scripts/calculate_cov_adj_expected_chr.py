## Script to calculate expected values using Liu's coverage correction
## This is to correct for there being less chance of identifying true de novos when coverage is lower, so our expected value has to be adjusted down
## for depth >= 50, no adjustment needed (i.e. mRate*(individuals>50 * 2)
## for depth <= 1, exp*0.119 (i.e. mRate*0.119*(individuals<1*2))
## for depth >1 and <50, per individual:  exp*(0.119 + 0.204(log(depth)))*2, summed across all individuals
## 
## annotation table: /lustre/scratch115/realdata/mdt1/teams/hurles/users/jl18/Splicing/Exons5tag/AnnotationTables/Full_ann_tables/Path_score_percentiles/Rmethod_050616/all_5tag_annotations_percentiles_mRate.txt 
## Run from in: /lustre/scratch115/realdata/mdt1/teams/hurles/users/jl18/Splicing/DDD8K/Denovos/Liu_coverage_correction/cov_adj_expected

import math 
import argparse


def get_options():
    """get the command line options"""
    parser = argparse.ArgumentParser(description="###")
    parser.add_argument("--chrom", required=True, help="# Get chromosome to investigate (1-22)")
    args = parser.parse_args()
    return args
args = get_options()



## coverage_summary_counts - how many individuals fall into high/medium/low categories (>50/>1 but <50/<= 1)
## coverage_summary_values - the actual coverage values for all the individuals in the medium coverage category, where the mRate calculation takes into account the actual observed depth of coverage
coverage_counts = open("/lustre/scratch115/realdata/mdt1/teams/hurles/users/jl18/Splicing/DDD8K/Denovos/Liu_coverage_correction/coverage_summary_counts.txt")
coverage_values = open("/lustre/scratch115/realdata/mdt1/teams/hurles/users/jl18/Splicing/DDD8K/Denovos/Liu_coverage_correction/coverage_values_by_exon.txt")

annotation_table = open("/lustre/scratch115/realdata/mdt1/teams/hurles/users/jl18/Splicing/Exons5tag/AnnotationTables/Full_ann_tables/Path_score_percentiles/Rmethod_050616/all_5tag_annotations_percentiles_mRate.txt")

Error_file = ''.join(("/lustre/scratch115/realdata/mdt1/teams/hurles/users/jl18/Splicing/DDD8K/DDG2P_update/Update_splice_obs_exp/exp/sites_with_cov_not_found_", str(args.chrom), ".txt"))
error_file =  open(Error_file, 'w')
Outfile = ''.join(("/lustre/scratch115/realdata/mdt1/teams/hurles/users/jl18/Splicing/DDD8K/DDG2P_update/Update_splice_obs_exp/exp/cov_adj_expected_full_", str(args.chrom), ".txt"))
outfile = open(Outfile, 'w')

## Dictionaries for running totals of mutation rates
all_dict = {"don+8": 0, "acc-13": 0, "don+6": 0, "acc-24": 0, "acc-23": 0, "acc+1": 0, "acc-20": 0, "don-2": 0, "don+2": 0, "acc+4": 0, "don+1": 0, "acc+10": 0, "acc+8": 0, "don+7": 0, "acc+5": 0, "acc+2": 0, "acc-18": 0, "acc-22": 0, "acc+7": 0, "acc-2": 0, "acc+6": 0, "don": 0, "acc+9": 0, "don-8": 0, "acc-14": 0, "acc-7": 0, "acc-5": 0, "don-10": 0, "don-1": 0, "acc-1": 0, "acc-3": 0, "acc-12": 0, "acc-17": 0, "don+9": 0, "don-6": 0, "don-5": 0, "acc-15": 0, "acc-21": 0, "don+4": 0, "acc-16": 0, "acc-9": 0, "acc-8": 0, "don-7": 0, "don+5": 0, "don-3": 0, "acc": 0, "don+10": 0, "acc-4": 0, "don-4": 0, "don+3": 0, "acc-11": 0, "don-9": 0, "acc-19": 0, "acc+3": 0, "acc-6": 0, "acc-10": 0, "acc-25": 0,}
dom_dict = {"don+8": 0, "acc-13": 0, "don+6": 0, "acc-24": 0, "acc-23": 0, "acc+1": 0, "acc-20": 0, "don-2": 0, "don+2": 0, "acc+4": 0, "don+1": 0, "acc+10": 0, "acc+8": 0, "don+7": 0, "acc+5": 0, "acc+2": 0, "acc-18": 0, "acc-22": 0, "acc+7": 0, "acc-2": 0, "acc+6": 0, "don": 0, "acc+9": 0, "don-8": 0, "acc-14": 0, "acc-7": 0, "acc-5": 0, "don-10": 0, "don-1": 0, "acc-1": 0, "acc-3": 0, "acc-12": 0, "acc-17": 0, "don+9": 0, "don-6": 0, "don-5": 0, "acc-15": 0, "acc-21": 0, "don+4": 0, "acc-16": 0, "acc-9": 0, "acc-8": 0, "don-7": 0, "don+5": 0, "don-3": 0, "acc": 0, "don+10": 0, "acc-4": 0, "don-4": 0, "don+3": 0, "acc-11": 0, "don-9": 0, "acc-19": 0, "acc+3": 0, "acc-6": 0, "acc-10": 0, "acc-25": 0,}
rec_dict = {"don+8": 0, "acc-13": 0, "don+6": 0, "acc-24": 0, "acc-23": 0, "acc+1": 0, "acc-20": 0, "don-2": 0, "don+2": 0, "acc+4": 0, "don+1": 0, "acc+10": 0, "acc+8": 0, "don+7": 0, "acc+5": 0, "acc+2": 0, "acc-18": 0, "acc-22": 0, "acc+7": 0, "acc-2": 0, "acc+6": 0, "don": 0, "acc+9": 0, "don-8": 0, "acc-14": 0, "acc-7": 0, "acc-5": 0, "don-10": 0, "don-1": 0, "acc-1": 0, "acc-3": 0, "acc-12": 0, "acc-17": 0, "don+9": 0, "don-6": 0, "don-5": 0, "acc-15": 0, "acc-21": 0, "don+4": 0, "acc-16": 0, "acc-9": 0, "acc-8": 0, "don-7": 0, "don+5": 0, "don-3": 0, "acc": 0, "don+10": 0, "acc-4": 0, "don-4": 0, "don+3": 0, "acc-11": 0, "don-9": 0, "acc-19": 0, "acc+3": 0, "acc-6": 0, "acc-10": 0, "acc-25": 0,}
nD_dict = {"don+8": 0, "acc-13": 0, "don+6": 0, "acc-24": 0, "acc-23": 0, "acc+1": 0, "acc-20": 0, "don-2": 0, "don+2": 0, "acc+4": 0, "don+1": 0, "acc+10": 0, "acc+8": 0, "don+7": 0, "acc+5": 0, "acc+2": 0, "acc-18": 0, "acc-22": 0, "acc+7": 0, "acc-2": 0, "acc+6": 0, "don": 0, "acc+9": 0, "don-8": 0, "acc-14": 0, "acc-7": 0, "acc-5": 0, "don-10": 0, "don-1": 0, "acc-1": 0, "acc-3": 0, "acc-12": 0, "acc-17": 0, "don+9": 0, "don-6": 0, "don-5": 0, "acc-15": 0, "acc-21": 0, "don+4": 0, "acc-16": 0, "acc-9": 0, "acc-8": 0, "don-7": 0, "don+5": 0, "don-3": 0, "acc": 0, "don+10": 0, "acc-4": 0, "don-4": 0, "don+3": 0, "acc-11": 0, "don-9": 0, "acc-19": 0, "acc+3": 0, "acc-6": 0, "acc-10": 0, "acc-25": 0,}

exclude_dict = {"missense_variant": 0, "frameshift": 0, "stop_gained": 0, "stop_lost": 0, "initiator_codon_variant": 0,}




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


## Go through annotation table and use strand and distance from exon to calculate the coordinates of the nearest splice site to use to look up the coverage values
for line in annotation_table:
    line = line.strip()
    words = line.split('\t')
    if line.startswith('X'): continue
    if words[9] in exclude_dict: continue
    if words[0] != str(args.chrom): continue
    get_distance = words[4].split('-')
    distance = get_distance[0]
    strand = words[6]
    mRate = words[8]
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

## Lookup the coordinate in the tables above and calculate the expected value, adjusting appropriately
    if working_coor not in high_dict: 
        outline = ''.join((working_coor, '\t', line, '\n'))
        error_file.write(outline)
    if working_coor not in high_dict: continue

## High coverage (no adjustment needed, just multiply by individuals and by 2 (for alleles))
    high_individuals = high_dict[working_coor]
    mRate_current = float(mRate) * float(high_individuals) * 2
    all_dict[words[4]] = float(all_dict[words[4]]) + float(mRate_current)
    if words[10] == "dominant":
        dom_dict[words[4]] = float(dom_dict[words[4]]) + float(mRate_current)
    if words[10] == "recessive":
        rec_dict[words[4]] = float(rec_dict[words[4]]) + float(mRate_current)
    if words[10] == "nonDDG2P":
        nD_dict[words[4]] = float(nD_dict[words[4]]) + float(mRate_current)
## Low coverage (mRate*0.119*(individuals<1*2))
    low_individuals = low_dict[working_coor]
    mRate_current = float(mRate) * 0.119 * (float(low_individuals) * 2)
    all_dict[words[4]] = float(all_dict[words[4]]) + float(mRate_current)
    if words[10] == "dominant":
        dom_dict[words[4]] = float(dom_dict[words[4]]) + float(mRate_current)
    if words[10] == "recessive":
        rec_dict[words[4]] = float(rec_dict[words[4]]) + float(mRate_current)
    if words[10] == "nonDDG2P":
        nD_dict[words[4]] = float(nD_dict[words[4]]) + float(mRate_current)
## Middle coverage (exp*(0.119 + 0.204(log(depth)))*2, summed across all individuals)
    for i in values_dict[working_coor]:
        mRate_current = float(mRate) * (0.119 + 0.204*(math.log(float(i)))) * 2
        all_dict[words[4]] = float(all_dict[words[4]]) + float(mRate_current)
        if words[10] == "dominant":
            dom_dict[words[4]] = float(dom_dict[words[4]]) + float(mRate_current)
        if words[10] == "recessive":
            rec_dict[words[4]] = float(rec_dict[words[4]]) + float(mRate_current)
        if words[10] == "nonDDG2P":
            nD_dict[words[4]] = float(nD_dict[words[4]]) + float(mRate_current)



for i in all_dict:
    outline = ''.join((str(args.chrom), '\t', i, '\t', str(all_dict[i]), '\t', str(dom_dict[i]), '\t', str(rec_dict[i]), '\t', str(nD_dict[i]), '\n'))
    outfile.write(outline)
outfile.close()

error_file.close()

