## Script to start building input file for MAPS - wants chr pos ref alt AC CQ triplet
## Run from in /lustre/scratch110/sanger/jl18/Exons5tag/MAPS/ExAC
## Infile = /lustre/scratch110/sanger/jl18/Exons5tag/AnnotationTables/ExAC_5tag_region/All_5tag_ExAC_variants_noCQfilt.txt which has full ExAC info + CQ



infile_ExAc = open("/lustre/scratch115/teams/hurles/users/jl18/Splicing/DDD8K/MAPS/ExAC_redo/All_5tag_ExAC_variants_noCQfilt.txt")
out = open("/lustre/scratch115/teams/hurles/users/jl18/Splicing/DDD8K/MAPS/ExAC_redo/MAPS_input_working1.txt", 'w')

exonic_dict = {"acc": 0,  "acc+1": 0,  "acc+2": 0,  "acc+3": 0,  "acc+4": 0,  "acc+5": 0,  "acc+6": 0,  "acc+7": 0,  "acc+8": 0,  "acc+9": 0,  "acc+10": 0, "don-10": 0, "don-9": 0, "don-8": 0, "don-7": 0, "don-6": 0, "don-5": 0, "don-4": 0, "don-3": 0, "don-2": 0, "don-1": 0, "don": 0,}

totals_dict = {"acc": 0,  "acc+1": 0,  "acc+2": 0,  "acc+3": 0,  "acc+4": 0,  "acc+5": 0, "bp": 0,  "acc+6": 0,  "acc+7": 0,  "acc+8": 0,  "acc+9": 0,  "acc+10": 0, "don-10": 0, "don-9": 0, "don-8": 0, "don-7": 0, "don-6": 0, "don-5": 0, "don-4": 0, "don-3": 0, "don-2": 0, "don-1": 0, "don": 0, "acc-25": 0,  "acc-24": 0,  "acc-23": 0,  "acc-22": 0,  "acc-21": 0,  "acc-20": 0,  "acc-19": 0,  "acc-18": 0,  "acc-17": 0,  "acc-16": 0,  "acc-15": 0,  "acc-14": 0,  "acc-13": 0,  "acc-12": 0,  "acc-11": 0,  "acc-10": 0,  "acc-9": 0,  "acc-8": 0,  "acc-7": 0,  "acc-6": 0,  "acc-5": 0,  "acc-4": 0,  "acc-3": 0,  "acc-2": 0,  "acc-1": 0,"don+1": 0, "don+2": 0, "don+3": 0, "don+4": 0, "don+5": 0, "don+6": 0, "don+7": 0, "don+8": 0, "don+9": 0, "don+10": 0,}


for line in infile_ExAc:
    line = line.strip()
    words = line.split('\t')
    if words[6] != "PASS": continue ## Filter sites that don't have "PASS" annotation from ExAC
    counts = words[7].split(';')
    counts1 = counts[0].split('=')
    vars = words[4].split(',')
    no_vars = len(vars)
    AC_counts = counts1[1]
    AC = AC_counts.split(',')
    for i in range(0,len(vars)):
        if len(words[3]) != 1 or len(vars[i]) != 1: continue
        count = AC[i]
        outline = ''.join((words[0], '\t', words[1], '\t', words[3], '\t', vars[i], '\t', str(count), '\t', words[8], '\n'))
        out.write(outline)
out.close()
