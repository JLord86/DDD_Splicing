## Script to make the sites files for MAPS v pathogenicity scores (path) analysis
## Have score brackets in /lustre/scratch110/sanger/jl18/Exons5tag/AnnotationTables/Full_ann_tables/Path_score_percentiles/Rmethod_050616/all_5tag_annotations_percentiles.txt
## Want for each score to have 1 file per chromosome with just the position and the bracket
## Do in two stages - get chr	pos	bracket in a single file for each score, then split these files by chromosome
## Run from in /lustre/scratch110/sanger/jl18/Exons5tag/Path_scores
## Outputs to go in /lustre/scratch110/sanger/jl18/Exons5tag/Path_scores/Ada/ or /RF/ or /MESPD/ or /Spidex/
## 
## 


infile = open("/lustre/scratch110/sanger/jl18/Exons5tag/AnnotationTables/Full_ann_tables/Path_score_percentiles/Rmethod_050616/all_5tag_annotations_percentiles_mRate.txt")
outfile_A = open("/lustre/scratch110/sanger/jl18/Exons5tag/Path_scores/Ada/Ada_brackets_all_chr.txt", 'w')

## Pull the bracket info from the annotation file
for line in infile:
    line = line.strip()
    words = line.split('\t')
    outline_A = ''.join((words[0], '\t', words[1], '\t', words[20], '\n'))
    outfile_A.write(outline_A)
outfile_A.close()


infile_A = open("/lustre/scratch110/sanger/jl18/Exons5tag/Path_scores/Ada/Ada_brackets_all_chr.txt")
out1 = open("/lustre/scratch110/sanger/jl18/Exons5tag/Path_scores/Ada/1.txt", 'w')
out2 = open("/lustre/scratch110/sanger/jl18/Exons5tag/Path_scores/Ada/2.txt", 'w')
out3 = open("/lustre/scratch110/sanger/jl18/Exons5tag/Path_scores/Ada/3.txt", 'w')
out4 = open("/lustre/scratch110/sanger/jl18/Exons5tag/Path_scores/Ada/4.txt", 'w')
out5 = open("/lustre/scratch110/sanger/jl18/Exons5tag/Path_scores/Ada/5.txt", 'w')
out6 = open("/lustre/scratch110/sanger/jl18/Exons5tag/Path_scores/Ada/6.txt", 'w')
out7 = open("/lustre/scratch110/sanger/jl18/Exons5tag/Path_scores/Ada/7.txt", 'w')
out8 = open("/lustre/scratch110/sanger/jl18/Exons5tag/Path_scores/Ada/8.txt", 'w')
out9 = open("/lustre/scratch110/sanger/jl18/Exons5tag/Path_scores/Ada/9.txt", 'w')
out10 = open("/lustre/scratch110/sanger/jl18/Exons5tag/Path_scores/Ada/10.txt", 'w')
out11 = open("/lustre/scratch110/sanger/jl18/Exons5tag/Path_scores/Ada/11.txt", 'w')
out12 = open("/lustre/scratch110/sanger/jl18/Exons5tag/Path_scores/Ada/12.txt", 'w')
out13 = open("/lustre/scratch110/sanger/jl18/Exons5tag/Path_scores/Ada/13.txt", 'w')
out14 = open("/lustre/scratch110/sanger/jl18/Exons5tag/Path_scores/Ada/14.txt", 'w')
out15 = open("/lustre/scratch110/sanger/jl18/Exons5tag/Path_scores/Ada/15.txt", 'w')
out16 = open("/lustre/scratch110/sanger/jl18/Exons5tag/Path_scores/Ada/16.txt", 'w')
out17 = open("/lustre/scratch110/sanger/jl18/Exons5tag/Path_scores/Ada/17.txt", 'w')
out18 = open("/lustre/scratch110/sanger/jl18/Exons5tag/Path_scores/Ada/18.txt", 'w')
out19 = open("/lustre/scratch110/sanger/jl18/Exons5tag/Path_scores/Ada/19.txt", 'w')
out20 = open("/lustre/scratch110/sanger/jl18/Exons5tag/Path_scores/Ada/20.txt", 'w')
out21 = open("/lustre/scratch110/sanger/jl18/Exons5tag/Path_scores/Ada/21.txt", 'w')
out22 = open("/lustre/scratch110/sanger/jl18/Exons5tag/Path_scores/Ada/22.txt", 'w')
outX = open("/lustre/scratch110/sanger/jl18/Exons5tag/Path_scores/Ada/X.txt", 'w')


## Split into sites file for use with singleton counting script
for line in infile_A:
    line = line.strip()
    words = line.split('\t')
    if words[2] == "NA": continue

    if words[0] == '1':
        outline = ''.join((words[1], '\t', words[2], '\n'))
        out1.write(outline)
    if words[0] == '2':
        outline = ''.join((words[1], '\t', words[2], '\n'))
        out2.write(outline)
    if words[0] == '3':
        outline = ''.join((words[1], '\t', words[2], '\n'))
        out3.write(outline)
    if words[0] == '4':
        outline = ''.join((words[1], '\t', words[2], '\n'))
        out4.write(outline)
    if words[0] == '5':
        outline = ''.join((words[1], '\t', words[2], '\n'))
        out5.write(outline)
    if words[0] == '6':
        outline = ''.join((words[1], '\t', words[2], '\n'))
        out6.write(outline)
    if words[0] == '7':
        outline = ''.join((words[1], '\t', words[2], '\n'))
        out7.write(outline)
    if words[0] == '8':
        outline = ''.join((words[1], '\t', words[2], '\n'))
        out8.write(outline)
    if words[0] == '9':
        outline = ''.join((words[1], '\t', words[2], '\n'))
        out9.write(outline)
    if words[0] == '10':
        outline = ''.join((words[1], '\t', words[2], '\n'))
        out10.write(outline)
    if words[0] == '11':
        outline = ''.join((words[1], '\t', words[2], '\n'))
        out11.write(outline)
    if words[0] == '12':
        outline = ''.join((words[1], '\t', words[2], '\n'))
        out12.write(outline)
    if words[0] == '13':
        outline = ''.join((words[1], '\t', words[2], '\n'))
        out13.write(outline)
    if words[0] == '14':
        outline = ''.join((words[1], '\t', words[2], '\n'))
        out14.write(outline)
    if words[0] == '15':
        outline = ''.join((words[1], '\t', words[2], '\n'))
        out15.write(outline)
    if words[0] == '16':
        outline = ''.join((words[1], '\t', words[2], '\n'))
        out16.write(outline)
    if words[0] == '17':
        outline = ''.join((words[1], '\t', words[2], '\n'))
        out17.write(outline)
    if words[0] == '18':
        outline = ''.join((words[1], '\t', words[2], '\n'))
        out18.write(outline)
    if words[0] == '19':
        outline = ''.join((words[1], '\t', words[2], '\n'))
        out19.write(outline)
    if words[0] == '20':
        outline = ''.join((words[1], '\t', words[2], '\n'))
        out20.write(outline)
    if words[0] == '21':
        outline = ''.join((words[1], '\t', words[2], '\n'))
        out21.write(outline)
    if words[0] == '22':
        outline = ''.join((words[1], '\t', words[2], '\n'))
        out22.write(outline)
    if words[0] == 'X':
        outline = ''.join((words[1], '\t', words[2], '\n'))
        outX.write(outline)
out1.close()
out2.close()
out3.close()
out4.close()
out5.close()
out6.close()
out7.close()
out8.close()
out9.close()
out10.close()
out11.close()
out12.close()
out13.close()
out14.close()
out15.close()
out16.close()
out17.close()
out18.close()
out19.close()
out20.close()
out21.close()
out22.close()
outX.close()

