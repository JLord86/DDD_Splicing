## Script to count variants per splice annotation in pLI >0.9 genes, as well as the total, so I can work out proportions
## 
## Second stage of the process - run after annotate_ACs_with_gene_and_strand_FS10.py

outfile = open("/lustre/scratch115/teams/hurles/users/jl18/Splicing/GenRes_redo/MAPS/FS10_MAPS/Deficit/Splice_region/splice_annotations_proportions_by_site_noNA.txt", 'w')

pLI_file = open("/lustre/scratch115/realdata/mdt1/teams/hurles/users/jl18/Splicing/DDD8K/Denovos/pLI_5tag/Deficit_in_parental_data/MAPS_ACs_approach/cleaned_exac_with_pLI_march16_from_JM_pLIgenes.txt")
high_dict = {}
all_dict = {}

for line in pLI_file:
    line = line.strip()
    if line.startswith("transcript"): continue
    words = line.split('\t')
    gene = words[1]
    pLI = words[19]
    if float(pLI) > 0.9:
        high_dict[gene] = pLI
    all_dict[gene] = pLI

counts_all = {"acc-25": 0, "acc-24": 0, "acc-23": 0, "acc-22": 0, "acc-21": 0, "acc-20": 0, "acc-19": 0, "acc-18": 0, "acc-17": 0, "acc-16": 0, "acc-15": 0, "acc-14": 0, "acc-13": 0, "acc-12": 0, "acc-11": 0, "acc-10": 0, "acc-9": 0,
    "acc-8": 0, "acc-7": 0, "acc-6": 0,
    "acc-5": 0, "acc-4": 0, "acc-3": 0,
    "acc-2": 0, "acc-1": 0, "acc": 0,
    "acc+1": 0, "acc+2": 0,
    "acc+3": 0, "acc+4": 0,
    "acc+5": 0, "acc+6": 0,
    "acc+7": 0, "acc+8": 0,
    "acc+9": 0, "acc+10": 0,
    "don-10": 0, "don-9": 0,
    "don-8": 0, "don-7": 0, "don-6": 0,
    "don-5": 0, "don-4": 0, "don-3": 0,
    "don-2": 0, "don-1": 0, "don": 0,
    "don+1": 0, "don+2": 0,
    "don+3": 0, "don+4": 0,
    "don+5": 0, "don+6": 0,
    "don+7": 0, "don+8": 0,
    "don+9": 0, "don+10": 0}

counts_high = {"acc-25": 0, "acc-24": 0, "acc-23": 0, "acc-22": 0, "acc-21": 0, "acc-20": 0, "acc-19": 0, "acc-18": 0, "acc-17": 0, "acc-16": 0, "acc-15": 0, "acc-14": 0, "acc-13": 0, "acc-12": 0, "acc-11": 0, "acc-10": 0, "acc-9": 0,
    "acc-8": 0, "acc-7": 0, "acc-6": 0,
    "acc-5": 0, "acc-4": 0, "acc-3": 0,
    "acc-2": 0, "acc-1": 0, "acc": 0,
    "acc+1": 0, "acc+2": 0,
    "acc+3": 0, "acc+4": 0,
    "acc+5": 0, "acc+6": 0,
    "acc+7": 0, "acc+8": 0,
    "acc+9": 0, "acc+10": 0,
    "don-10": 0, "don-9": 0,
    "don-8": 0, "don-7": 0, "don-6": 0,
    "don-5": 0, "don-4": 0, "don-3": 0,
    "don-2": 0, "don-1": 0, "don": 0,
    "don+1": 0, "don+2": 0,
    "don+3": 0, "don+4": 0,
    "don+5": 0, "don+6": 0,
    "don+7": 0, "don+8": 0,
    "don+9": 0, "don+10": 0}


for i in range(1,23):
    Infile = ''.join(("/lustre/scratch115/teams/hurles/users/jl18/Splicing/GenRes_redo/MAPS/FS10_MAPS/Deficit/Splice_region/Annotated_ACs_DDD/", str(i), ".txt"))
    infile = open(Infile)
    for line in infile:
        line = line.strip()
        words = line.split('\t')
        if words[7] not in all_dict: continue
        counts_all[words[5]] = int(counts_all[words[5]]) + 1
        if words[7] in high_dict:
            counts_high[words[5]] = int(counts_high[words[5]]) + 1

for i in counts_all:
    print(i, counts_high[i], counts_all[i])
    prop = float(counts_high[i]) / float(counts_all[i])
    outline = ''.join((i, '\t', str(counts_all[i]), '\t', str(counts_high[i]), '\t', str(prop), '\n'))
    outfile.write(outline)
outfile.close()



