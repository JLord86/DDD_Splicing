## Script to match up the MAPS ACs output to the strand and gene for the work looking for a deficit of variants
## in high pLI genes in parental data
## 
## This input data has already had exonic non-synonymous sites removed
## 



for i in range(1,23):
    Annfile = ''.join(("/lustre/scratch115/realdata/mdt1/teams/hurles/users/jl18/Splicing/Exons5tag/AnnotationTables/DDD8K_VEP_filt_for_singletons/", str(i), ".txt"))
    annfile = open(Annfile)
    Outfile = ''.join(("/lustre/scratch115/teams/hurles/users/jl18/Splicing/GenRes_redo/MAPS/FS10_MAPS/Deficit/Splice_region/Annotated_ACs_DDD/", str(i), ".txt"))
    outfile = open(Outfile, 'w')
    Infile = ''.join(("/lustre/scratch115/realdata/mdt1/teams/hurles/users/jl18/Splicing/GenRes_redo/MAPS/DDD_qual/FS10/", str(i), "_DDD8K_ACs_MAPS1.txt"))
    infile = open(Infile)
    dict_strand = {}
    dict_gene = {}
    for line in annfile:
        line = line.strip()
        words = line.split('\t')
        variant = '-'.join((str(i), words[0]))
        dict_strand[variant] = words[2]
        dict_gene[variant] = words[3]
    for line in infile:
        line = line.strip()
        words = line.split('\t')
        variant = '-'.join((words[0], words[1]))
        if len(words[2]) != 1 or len(words[3]) != 1: continue  ## Take out indels, which are problematic for positional annotations
        if variant in dict_strand:
            outline = ''.join((line, '\t', dict_strand[variant], '\t', dict_gene[variant], '\n'))
            outfile.write(outline)
    outfile.close()









