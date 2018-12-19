"""
"""

from __future__ import print_function
from __future__ import division

import glob
import os
import json
import gzip
import sys
import argparse

IS_PYTHON3 = sys.version[0] == "3"

VCF_DIR = "/lustre/scratch115/realdata/mdt0/projects/ddd/variant_calling/final"
CHROM = "22"

consequence_counts = {"acc-25": 0, "acc-24": 0, "acc-23": 0, "acc-22": 0, "acc-21": 0, "acc-20": 0, "acc-19": 0, "acc-18": 0, "acc-17": 0, "acc-16": 0, "acc-15": 0, "acc-14": 0, "acc-13": 0, "acc-12": 0, "acc-11": 0, "acc-10": 0, "acc-9": 0, "acc-8": 0, "acc-7": 0, "acc-6": 0, "acc-5": 0, "acc-4": 0, "acc-3": 0, "acc-2": 0, "acc-1": 0, "acc": 0, "acc+1": 0, "acc+2": 0, "acc+3": 0, "acc+4": 0, "acc+5": 0, "acc+6": 0, "acc+7": 0, "acc+8": 0, "acc+9": 0, "acc+10": 0, "don-10": 0, "don-9": 0, "don-8": 0, "don-7": 0, "don-6": 0, "don-5": 0, "don-4": 0, "don-3": 0, "don-2": 0, "don-1": 0, "don": 0, "don+1": 0, "don+2": 0, "don+3": 0, "don+4": 0, "don+5": 0, "don+6": 0, "don+7": 0, "don+8": 0, "don+9": 0, "don+10": 0, "bp": 0,}

def get_options():
    """ get the command line options
    """
    parser = argparse.ArgumentParser(description="Counts singletons within \
        multiple sample VCFs.")
    parser.add_argument("--chrom", required=True, help="chromosome to investigate.")
    parser.add_argument("--vcf", default=sys.stdout,
        help="location of vcf file to investigate")
    parser.add_argument("--singletons", default=sys.stdout,
        help="path to send the list of singletons to.")
    parser.add_argument("--totals", default=sys.stdout,
        help="path to send the totals for each consequence class.")
    
    args = parser.parse_args()
    
    return args

def get_ddd_parents():
    """ get a dictionary of unaffected DDD parents, to their sex
    """
    
    DIR = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2016-10-03/"
    family_path = os.path.join(DIR, "family_relationships.ped")
    sanger_id_path = os.path.join(DIR, "decipher_stable_sanger_ega.txt")
    
    # find all the unaffected parental DDD IDs.
    parental_ids = {}
    with open(family_path, "r") as handle:
        for line in handle:
            line = line.strip().split("\t")
            person_id = line[1]
            paternal_id = line[2]
            sex = line[4]
            affected_status = line[5]
            
            if paternal_id != "0" or affected_status != "1":
                continue
            
            parental_ids[person_id] = sex
    
    # get the unaffected parental sanger IDs (since the sanger IDs are the IDs
    # used in the multisample VCFs)
    sanger_ids = {}
    with open(sanger_id_path, "r") as handle:
        for line in handle:
            line = line.strip().split("\t")
            person_id = line[1]
            decipher_id = line[0]
            sanger_id = line[2]
            
            if person_id in parental_ids:
                sanger_ids[sanger_id] = parental_ids[person_id]
    
    return sanger_ids

def get_splice_annotations(chrom):
    
    splice_dir = "/lustre/scratch115/realdata/mdt1/teams/hurles/users/jl18/Splicing/GenRes_redo/MAPS/TargetSites_noNonSyn"  ## VEP filtered so no non-syn - based on AnnotationTable
    splice_path = os.path.join(splice_dir, "{}.txt".format(chrom))
    
    splice_consequences = {}
    with open(splice_path, "r") as handle:
        exclude_header(handle)
        for line in handle:
            line = line.strip().split("\t")
            
            pos = int(line[0])
            consequence = line[3]
            

            splice_consequences[pos] = consequence
    
    return splice_consequences

def get_header(vcf):
    """ get the full VCF header as a list of lines
    
    Args:
        vcf: handle to vcf file
    """
    
    current_pos = vcf.tell()
    vcf.seek(0)
    
    header = []
    for line in vcf:
        if not line.startswith("#"):
            break
        
        header.append(line)
    vcf.seek(current_pos)
    
    return(header)

def exclude_header(vcf):
    """ move the file handle to just beyond the VCF header lines
    """
    
    current_pos = vcf.tell()
    
    while vcf.readline().startswith("#"):
        current_pos = vcf.tell()
    
    vcf.seek(current_pos)

def get_sample_positions(header):
    """ make a dictionary of sample IDs (from the VCF header) vs their position
    
    Args:
        header: list of VCF header lines, the final line defines the VCF columns
    """
    
    samples = {}
    
    sample_ids = header[-1].strip().split("\t")[9:]
    sample_map = dict(zip(sample_ids, range(0, len(sample_ids))))
    
    return sample_map

def get_variant_key(variant):
    """ get the chrom, pos, and alleles for a VCF variant
    
    Args:
        variant: list of data from VCF line for the first 8 columns
    """
    
    chrom = variant[0]
    pos = int(variant[1])
    ref = variant[3]
    alts = variant[4]
    
    
    return (chrom, pos, ref, alts) 

def get_format(format_string):
    """ figure out the format of the sample columns, map data to position
    
    Args:
        format_string: text from format column of VCF, colon separated string
    """
    
    format = format_string.split(":")
    format = dict(zip(format, range(0, len(format))))
    
    return(format)

def get_sample_genotypes(samples, format, sample_pos):
    """ get a dictionary of genotypes for specific sample IDs
    """
    
    samples = samples.strip().split("\t")
    
    genotypes = {}
    for sample_id in sample_pos:
        sample = samples[sample_pos[sample_id]]
        sample = sample.split(":")
        genotype = sample[format["GT"]]
        
        # drop missing genotypes, since they cannot contribute to the total
        # allele count
        if genotype == "./.":
            continue
        
        genotypes[sample_id] = genotype
    
    return genotypes

def reformat_chrX_genotypes(key, genotypes, ddd_parents):
    """ swap male genotypes on chrX to a hemizgous type
    """
    
    # don't alter genotypes not on the X
    if key[0] != "X":
        return genotypes
    
    # define the pseudoautosomal regions on the X chromosome
    x_par = [(60001, 2699520), (154930290, 155260560), (88456802, 92375509)]
    
    # don't alter genotypes within the pseudoautosomal regions
    if any([key[1] >= x[0] and key[1] < x[1] for x in x_par ]):
        return genotypes
    
    exclude_ids = []
    for sample_id in genotypes:
        # don't alter female parents, since they are diploid for chrX
        if ddd_parents[sample_id] != "1":
            continue
        
        geno = genotypes[sample_id].split("/")
        
        # drop genotypes for heterozygous males on chrX
        if geno[0] != geno[1]:
            exclude_ids.append(sample_id)
        
        genotypes[sample_id] = geno[0]
    
    # remove the abberrant heterozygous chrX male genotypes
    for sample_id in exclude_ids:
        del genotypes[sample_id]
    
    return genotypes

def tally_alleles(genotypes, alts):
    """ count the alleles used in the genotypes
    """
    
    # make sure we have entries for all the alt alleles, even if the count
    # ends up as zero, that way when we later match to the consequences for the alt
    # alleles, the allele counts should be consistent with the consequences
    alts = alts.split(",")
    allele_numbers = ["{}".format(x + 1) for x in range(len(alts))]
    allele_numbers.append("0")
    
    # tally each allele found in the genotypes
    counts = dict(zip(allele_numbers, [0]* len(allele_numbers)))
    for genotype in genotypes.values():
        genotype = genotype.split("/")
        
        for allele in genotype:
            counts[allele] += 1
    
    return counts

def check_singletons(key, counts, splice, output, outfile):
    """ checks to see if any of the alleles at a site are singletons
    
    Args:
        key: tuple of (chrom, pos, ref, alt)
        counts: dictionary of number of alleles seen in the unaffecetd DDD parents
        splice: dictionary of consequence strings for each allele, indexed by chrom position

    """    
    consequences = splice[key[1]]
    
    # remove the reference allele, since we don't have a consequence for that
    del counts["0"]
    
    allele_numbers = sorted(counts)
  
    for number in allele_numbers:
        allele_count = counts[number]
        
        # don't include alleles with zero alleles in the unaffected parents,
        # otherwise they will skew the proportion of variants as singletons.
        if allele_count == 0:
            continue
        
        consequence = consequences
        consequence_counts[consequence] += 1
        altse = key[3].split(',')
        num = int(number) -1
        outline = ''.join((str(key[0]), '\t', str(key[1]), '\t', str(key[2]), '\t', str(altse[int(num)]), '\t', str(allele_count), '\t', consequence, '\n'))
        outfile.write(outline)
        
        # we only want to look at singletons
        if allele_count != 1:
            continue
        
        line = "{}\t{}\t{}\n".format(key[0], key[1], consequence)
        
        output.write(line)

def parse_vcf(chrom, ddd_parents, splice,  output_path, outfile):
    """ run through a VCF, counting the alleles, looking for singletons
    """
    
    try:
        output = open(output_path, "w")
    except TypeError:
        output = output_path
    
    args = get_options()
    vcf = gzip.open(args.vcf, "rt") ## JLtest
    # remove the header lines
    header = get_header(vcf)
    exclude_header(vcf)
    
    # figure out where the DDD parents are in the sample list
    index = get_sample_positions(header)
    sample_pos = dict([(x, index[x]) for x in ddd_parents if x in index])
    
    for line in vcf:
        line = line.split("\t", 9)
        variant = line[:9]
        key = get_variant_key(variant)
        
        if key[1] not in splice: continue
        
        format = get_format(variant[8])
        genotypes = get_sample_genotypes(line[9], format, sample_pos)
        genotypes = reformat_chrX_genotypes(key, genotypes, ddd_parents)
        counts = tally_alleles(genotypes, key[3])
        pos = (key[0], key[1])

        ## FS<10 filter added in to combat error rate issue
        info = line[7].split(';')
        for i in info:
            if "FS=" in i:
                FSscore1 = i.split('=')
                FSscore2 = FSscore1[1]
                print(FSscore2)
                if float(FSscore2) < 10: 
                    check_singletons(key, counts, splice, output, outfile)
                if float(FSscore2) > 10: continue

def main():
    args = get_options()
    vcf = gzip.open(args.vcf)
    ddd_parents = get_ddd_parents()
    splice = get_splice_annotations(args.chrom)
    Outfile = ''.join(("/lustre/scratch115/realdata/mdt1/teams/hurles/users/jl18/Splicing/GenRes_redo/MAPS/DDD_qual/FS10/", args.chrom, "_DDD8K_ACs_MAPS1.txt"))
    outfile = open(Outfile, 'w')
    try:
        parse_vcf(args.chrom, ddd_parents, splice, args.singletons, outfile)
    finally:
        # and finally show how many times each consequence was seen in the DDD
        # unaffected parents, so we can determine the proportion of that each
        # consequence is seen as a singleton
        try:
            output = open(args.totals, "w")
        except TypeError:
            output = args.totals
        
        output.write("consequence\toccurrences\n")
        for consequence in sorted(consequence_counts):
            line = "{}\t{}\n".format(consequence, consequence_counts[consequence])
            output.write(line)
        outfile.close()
if __name__ == "__main__":
    main()
