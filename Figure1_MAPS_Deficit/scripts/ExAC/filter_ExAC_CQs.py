## Script to bring CQ filtering of ExAC data more in line with CQ filtering of DDD data
## run from in: /lustre/scratch115/realdata/mdt1/teams/hurles/users/jl18/Splicing/GenRes_redo/MAPS/finding_origin/


infile = open("/lustre/scratch115/realdata/mdt1/teams/hurles/users/jl18/Splicing/DDD8K/MAPS/ExAC_redo/ExAC_MAPS_infile_syn_and_head.txt")
annfile = open("/lustre/scratch115/realdata/mdt1/teams/hurles/users/jl18/Splicing/Exons5tag/AnnotationTables/Full_ann_tables/Path_score_percentiles/Rmethod_050616/all_5tag_annotations_percentiles_mRate.txt")
vepfile = open("/lustre/scratch115/realdata/mdt1/teams/hurles/users/jl18/Splicing/DDD8K/MAPS/ExAC_redo/check_syn/All_VEP_ExAC_syn.txt")

outfile = open("/lustre/scratch115/realdata/mdt1/teams/hurles/users/jl18/Splicing/GenRes_redo/MAPS/finding_origin/properly_filtered_ExAC_data.txt", 'w')

vep_dict = {}
ann_dict = {}

print("storing_VEP")
for line in vepfile:
    line = line.strip()
    words = line.split('\t')
    if line.startswith('#'): continue
    variant = '-'.join((words[0], words[1], words[3], words[4]))
    vep_dict[variant] = words[7]

print("storing_anns")
for line in annfile:
    line = line.strip()
    words = line.split('\t')
    variant = '-'.join((words[0], words[1], words[2], words[3]))
    ann_dict[variant] = words[9]

print("making_comparisons")
for line in infile:
    line = line.strip()
    words = line.split('\t')
    if line.startswith("chr"): continue
    variant = '-'.join((words[0], words[1], words[2], words[3]))
    
    if words[5] == "synonymous":
        if "synonymous" not in vep_dict[variant]: continue
        outline = ''.join((line, '\n'))
        outfile.write(outline)
    else:
        if ann_dict[variant] == "missense_variant": continue
        outline = ''.join((line, '\n'))
        outfile.write(outline)
outfile.close()
