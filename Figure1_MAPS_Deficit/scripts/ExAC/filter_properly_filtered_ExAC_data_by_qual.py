## quality filter the properly CQ filtered ExAC data to bring in line with FS10 quality filtering on DDD data
## run from in: /lustre/scratch115/realdata/mdt1/teams/hurles/users/jl18/Splicing/GenRes_redo/MAPS/finding_origin/ExAC_qual
import gzip

infile_filt = open("/lustre/scratch115/realdata/mdt1/teams/hurles/users/jl18/Splicing/GenRes_redo/MAPS/finding_origin/properly_filtered_ExAC_data.txt")
infile_ExAC = gzip.open("/lustre/scratch115/realdata/mdt1/teams/hurles/users/jl18/General/ExAC.r0.3.1.sites.vep.vcf.gz", "rt")

outfile_10 = open("/lustre/scratch115/realdata/mdt1/teams/hurles/users/jl18/Splicing/GenRes_redo/MAPS/finding_origin/ExAC_qual/FS10_ExAC_filt.txt", 'w')


site_dict = {}

for line in infile_filt:
    line = line.strip()
    words = line.split('\t')
    variant = '-'.join((words[0], words[1], words[2], words[3]))
    site_dict[variant] = line
print("filt_variant", variant, line)

for line in infile_ExAC:
    line = line.strip()
    if line.startswith('#'): continue
    words = line.split('\t')
    alts = words[4].split(',')
    for i in alts:
        variant = '-'.join((words[0], words[1], words[3], i))
       # print("ExAC_pre_match", variant)
        if variant in site_dict:
            info = words[7].split(';')
            for j in info:
                if j.startswith('FS='):
                    get_FS = j.split('=')
                    FS = get_FS[1]
            outline = ''.join((site_dict[variant], '\n'))
            if float(FS) < 10:
                outfile_10.write(outline)
outfile_10.close()
