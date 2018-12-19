## Create target sites files for DDD MAPS work
## Write directories of A) just synonymous sites for use as baseline, and B) all other sites, excluding exonic non-synonymous consequences

import argparse

def get_options():
    """ get the command line options
    """
    parser = argparse.ArgumentParser(description="removes pathogenic annotated sites")
    parser.add_argument("--chrom", required=True, help="chromosome to investigate.")
    args = parser.parse_args()
    return args

args = get_options()

infile = open("/lustre/scratch115/realdata/mdt1/teams/hurles/users/jl18/Splicing/Exons5tag/AnnotationTables/Full_ann_tables/Path_score_percentiles/Rmethod_050616/all_5tag_annotations_percentiles_mRate.txt")

## Want to have separate sites files for just synonymous (to be used as baseline), and all NOT_nonSyn to be used in splice site MAPS
Outfile_Syn = ''.join(("/lustre/scratch115/realdata/mdt1/teams/hurles/users/jl18/Splicing/GenRes_redo/MAPS/TargetSites_syn/", args.chrom, ".txt"))
outfile_Syn = open(Outfile_Syn, 'w')

Outfile_NoNonSyn = ''.join(("/lustre/scratch115/realdata/mdt1/teams/hurles/users/jl18/Splicing/GenRes_redo/MAPS/TargetSites_noNonSyn/", args.chrom, ".txt"))
outfile_NoNonSyn = open(Outfile_NoNonSyn, 'w')

## These are the only CQs to be included - filters out any exonic missense and nonsense variants which could skew signature of selection without involvement in splicing
CQs_dict = {"synonymous_variant": 0, "intron_variant": 0, "splice_region_variant": 0, "splice_donor_variant": 0, "splice_acceptor_variant": 0, }

## Pull out sites from exons of interest which don't have non-synonymous consequences for the current chromosome
for line in infile:
    line = line.strip()
    words = line.split('\t')
    if words[0] != args.chrom: continue
    if words[9] not in CQs_dict: continue
    outline = ''.join((words[1], '\t', words[2], '\t', words[3], '\t', words[4], '\t', words[5], '\t', words[6], '\t', words[7], '\t', words[9], '\n'))
    outfile_NoNonSyn.write(outline)
    if words[9] == "synonymous_variant":
        outfile_Syn.write(outline)
outfile_Syn.close()
outfile_NoNonSyn.close()