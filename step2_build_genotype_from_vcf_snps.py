####################################################################################################
# This script
# 1) transforms all the snps found in the WGS VCF file into the simplified genotype format,
#    e.g. snp_anme Allel_X Allel_Y (hu1 A C)
# 2) adds in those Biomarks SNPs which are not in the WGS VCF and assign the reference genotype
#####################################################################################################


import sys
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("biomark_ref", help="Biomarks reference SNP file")
parser.add_argument("input_wgs_vcf_text", help="Input WGS VCF text file parsed by the step1 script")
parser.add_argument("output_geno_file", help="Output genotype file")
args = parser.parse_args()


fout = open(args.output_geno_file, "w")
snp_in_vcf = []  # track the snps that show in WGS VCF

reverse_snps = ["rs445251","rs447818","rs525869","rs530501","rs722290","rs870347","rs952718","rs985492","rs1008730","rs1019029","rs1040045","rs1336071","rs1344870","rs1513181","rs1554472","rs1760921","rs1821380","rs1823718","rs1872575","rs1876482","rs2125345","rs3737576","rs3780962","rs10488710"] # list of the reverse strand SNPs in the refSNP database.

def complement(base):
    complement = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    return complement[base]

# write genotype info for all the SNPs that appear in WGS VCF
with open(args.input_wgs_vcf_text) as vcf_snp:
    for line in vcf_snp:
        line = line.strip()
        arr = line.split("\t")
        rsID = arr[7]
        snp_in_vcf.append(rsID)
        name = arr[6]
        info = arr[5]
        geno = ""
        if rsID not in reverse_snps:
            if info[0:3] == "0/1":
                geno = arr[2]+"\t"+arr[3]
            elif info[0:3] == "1/0":
                geno = arr[3]+"\t"+arr[2]
            elif info[0:3] == "1/1":
                geno = arr[3]+"\t"+arr[3]
        else:
            if info[0:3] == "0/1":
                geno = complement(arr[2])+"\t"+complement(arr[3])
            elif info[0:3] == "1/0":
                geno = complement(arr[3])+"\t"+complement(arr[2])
            elif info[0:3] == "1/1":
                geno = complement(arr[3])+"\t"+complement(arr[3])
            
        print >>fout, rsID+"\t"+name+"\t"+geno


# write genotype info for the SNPs that are not present in WGS VCF
with open(args.biomark_ref) as biomarks_in:
    for line in biomarks_in:
        line = line.strip()
        if line[0].isdigit(): # skip the title line, non digit###
            arr = line.split("\t")
            if arr[2] not in snp_in_vcf:
                rsID = arr[2]
                name = arr[1]
                #check if it is a reverse strand snp, if yes, use the ref_allel_2 in the Biomarks ref file
                #otherwise use the ref_allel
                if rsID not in reverse_snps:
                    geno = arr[-2]+"\t"+arr[-2]
                else:
                    geno = arr[-1]+"\t"+arr[-1]
                    
                print >>fout, rsID+"\t"+name+"\t"+geno

            
fout.close()
