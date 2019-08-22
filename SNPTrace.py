import sys
import itertools

file = sys.argv[1]
sample_SNPs = {}

line_no = 0
with open(file) as fn:
    for line in fn:
        line_no = line_no + 1
        #Skip comment lines (assume top 16 lines)
        if line_no > 16:
            line = line.strip()
            arr = line.split(",")
            #print arr[4],arr[1],arr[9]
            # determine 'Sample' lines
            if arr[4].startswith("Samp"):
                if arr[4] not in sample_SNPs.keys():
                    #initiate the snpname/type list
                    sample_SNPs[arr[4]] = [(arr[1], arr[9])]
                else:
                    #populate the snpname/type list
                    sample_SNPs[arr[4]].append((arr[1], arr[9]))

def calculate_Identity_Score(sample1, sample2):
    sample1_SNP_No = 0
    sample2_SNP_No = 0

    # get unique number of snp calls for sample1
    for snp in sample_SNPs[sample1]:
        if snp[1] not in ["No Call", "Invalid"]:
            sample1_SNP_No = sample1_SNP_No + 1

    # get unique number of snp calls for sample2
    for snp in sample_SNPs[sample2]:
        if snp[1] not in ["No Call", "Invalid"]:
            sample2_SNP_No = sample2_SNP_No + 1
            
    # get number of snp calls shared in common between sample1 and sample2
    # assume snp list order is same between sample1 and sample2
    shared_snps = [s1 for s1, s2 in zip(sample_SNPs[sample1], sample_SNPs[sample2]) \
                   if (s1[1] not in ["No Call","Invalid"]) and (s2[1] not in ["No Call","Invalid"]) and (s1 == s2)]

    # calculate identity score as '2*common_SNP_No/(sample1_SNP_No+sample2_SNP_No)'
    ident_score = 2.0*len(shared_snps)/(sample1_SNP_No+sample2_SNP_No)
    
    return ident_score, shared_snps, sample1_SNP_No, sample2_SNP_No


# Calculate the identity scores of all sample pairs
identity_scores = {}

# generate all possible pairs from a list
sample_pairs = list(itertools.permutations(sample_SNPs.keys(),2))

# calculate pairwise sample identity scores and output to a file
ofile = "SNPTrace_result.txt"
with open(ofile, "w") as fout:
    for p in sample_pairs:
        ident_score, shared_snps, sample1_SNP_No, sample2_SNP_No = calculate_Identity_Score(p[0],p[1])
        print>>fout, p[0]+","+p[1]+","+str(ident_score)+","+str(sample1_SNP_No)+","+\
        str(sample2_SNP_No)+","+str(len(shared_snps))
        #print>>fout, shared_snps
    
