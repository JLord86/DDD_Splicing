## Compile mRates from different chromosomes for splicing annotations

all_dict = {"G1": 0, "G2": 0, "G3": 0, "G4": 0,"G5": 0, "G6": 0,  "NA": 0,}

outfile = open("/lustre/scratch115/realdata/mdt1/teams/hurles/users/jl18/Splicing/DDD8K/Denovos/pLI_5tag/don5_expected_coverage_adjusted_compiled_full_int_pLI_sextiles.txt", 'w')

for i in range(1,23):
    Infile = ''.join(("/lustre/scratch115/realdata/mdt1/teams/hurles/users/jl18/Splicing/DDD8K/Denovos/pLI_5tag/splice_exp/cov_adj_expected_full_", str(i), "_sextiles.txt"))
    infile = open(Infile)
    print Infile
    for line in infile:
        line = line.strip()
        words = line.split('\t')
        all_dict[words[1]] = float(all_dict[words[1]]) + float(words[2])

for i in all_dict:
    outline = ''.join((i, '\t', str(all_dict[i]), '\n'))
    outfile.write(outline)
outfile.close()

