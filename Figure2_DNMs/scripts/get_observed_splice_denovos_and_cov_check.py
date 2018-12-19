## Identify de novo mutations that overlap with splicing positions of interest
## (output will contain all quality scores, so filtering after this process will be needed)
## 
## Compile sites with coverage not found (output from the calculation of expected values) 
## to remove observed de novos that fall within sites that didn't contribute to the expected values because they didn't have coverage information

denovos = open("/lustre/scratch113/projects/ddd/users/jm33/de_novos.ddd_8k.noncoding_included.2016-06-23.txt") ## 8832 de novos, indels and SNVs
outfile = open("/lustre/scratch115/realdata/mdt1/teams/hurles/users/jl18/Splicing/DDD8K/Denovos/Liu_coverage_correction/cov_adj_expected/compile/observed_splice_denovos_all_ppdnm_noX_CQs_liu_cov_to_cov_check.txt", 'w')

annotation_table = open("/lustre/scratch115/realdata/mdt1/teams/hurles/users/jl18/Splicing/Exons5tag/AnnotationTables/Full_ann_tables/Path_score_percentiles/Rmethod_050616/all_5tag_annotations_percentiles_mRate.txt")

denovo_dict = {}

exclude_dict = {"missense_variant": 0, "frameshift": 0, "stop_gained": 0, "stop_lost": 0, "initiator_codon_variant": 0,}

for line in denovos:
    line = line.strip()
    words = line.split('\t')
    variant = '-'.join((words[1],  words[2], words[3], words[4]))
    denovo_dict[variant] = line

for line in annotation_table:
    line = line.strip()
    words = line.split('\t')
    if words[9] in exclude_dict: continue
    variant = '-'.join((words[0],  words[1], words[2], words[3]))
    if variant in denovo_dict:
        outline = ''.join((denovo_dict[variant], '\t', line, '\n'))
        outfile.write(outline)
outfile.close()


denovos2 = open("/lustre/scratch115/realdata/mdt1/teams/hurles/users/jl18/Splicing/DDD8K/Denovos/Liu_coverage_correction/cov_adj_expected/compile/observed_splice_denovos_all_ppdnm_noX_CQs_liu_cov_to_cov_check.txt")
outfile = open("/lustre/scratch115/realdata/mdt1/teams/hurles/users/jl18/Splicing/DDD8K/Denovos/Liu_coverage_correction/cov_adj_expected/compile/observed_splice_denovos_noX_noNS_cov_adj.txt", 'w')

sites_dict = {}

for i in range(1,23):
    Infile = ''.join(("/lustre/scratch115/realdata/mdt1/teams/hurles/users/jl18/Splicing/DDD8K/Denovos/Liu_coverage_correction/cov_adj_expected/other_CQs/sites_with_cov_not_found_otherCQs_", str(i), ".txt"))
    infile = open(Infile)
    print Infile
    for line in infile:
        line = line.strip()
        words = line.split('\t')
        variant = ''.join((words[1], ":", words[2]))
        sites_dict[variant] = 0

for line in denovos2:
    line = line.strip()
    words = line.split('\t')
    if line.startswith("chr"): continue
    variant = ''.join((words[1], ":", words[2]))
    if variant in sites_dict: continue
    outline = ''.join((line, '\n'))
    outfile.write(outline)
outfile.close()
