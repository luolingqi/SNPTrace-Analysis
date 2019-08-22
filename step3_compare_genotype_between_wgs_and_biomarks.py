############################################################################################
# This script compares the Genotype information between Biomarks and WGS, and generate a
# report which records:
# 1) # of matching SNPs with their rsID and name;
# 2) # of unique SNPs from either WGS or Biomarks result with their rsID and name
# 3) percentage of matching SNPs out of the whole set of 93 SNP (Not including the 3 Y-SNPs
##############################################################################################

import sys
##################################################################
# input two files:
# 1) WGS genotype info from step 2 with the format "rsID\tname\tallelX\tallelY"
# 2) Biomarks genotype info with the format "ID(run id),Assay(name),allelX,allelY". Make sure the file
#    a) is related to only one sample;
#    b) contains only the SNP lines, not other comment lines
##################################################################
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("WGS_Genotype_file", help="The WGS genotype file generated from step2 script ")
parser.add_argument("Biomark_Genotype_file", help="The Biomark genotype file, containing snp calling lines for single sample")
args = parser.parse_args()

Gender_SNPs = ['hu98Y','hu111Y','hu209Y']

# function to get intersection among elements of a list
def getIntersection(s):
    temp = [ set(i) for i in s[0] if i[0] not in Gender_SNPs ] # exclude Y SNPs    
    return [ set(j) for j in s[1] if set(j) in temp ]


def getDifference(s):
    temp0 = [ set(i) for i in s[0] ]
    unique1 = [ set(j) for j in s[1] if set(j) not in temp0 and j[0] not in Gender_SNPs ]
    temp1 = [ set(i) for i in s[1] ]
    unique2 = [ set(j) for j in s[0] if set(j) not in temp1 and j[0] not in Gender_SNPs ]
    return unique1, unique2



##############################################################
# Prepare the input lists for both WGS and Biomarks genotypes
##############################################################
# declare a dictionary to hold genotype from both Sequencer and Biomarks in one place. e.g. {hu1:"A T A T", hu2:"G C G C",....}
# The format should be snpName:"Biomark_allel1 Biomark_allel2 WGS_allel1 WGS_allel2"

dict_geno = {}

# make a list of SNPs from WGS VCF, padded with the unmuted SNPs with Ref allel done from step2: e.g.: ["hu1 G G","hu2 C C", ....]
list_WGS = []
with open(args.WGS_Genotype_file) as WGS_in:
    for line in WGS_in:
        line = line.strip()
        arr = line.split("\t")
        list_WGS.append([arr[1], arr[2], arr[3]])
        # Populate the dictionary with WGS geno calls
        if arr[1] in dict_geno: #Biomarks geno already populated!
            dict_geno[arr[1]] = dict_geno[arr[1]] + " " + arr[2] + " " + arr[3]
        else:
            dict_geno[arr[1]] = arr[2] + " " + arr[3]
        
list_Biomarks = []
nocall = 0
with open(args.Biomark_Genotype_file) as Biomarks_in:
    for line in Biomarks_in:
        line = line.strip()
        arr = line.split(",")
        if arr[10] != "No Call":
            list_Biomarks.append([arr[1], arr[10][0], arr[10][2]])
        else:
            list_Biomarks.append([arr[1], "No Call", "No Call"])
            #print out the no call snps
            if arr[1][-1] != 'Y': # if not Y SNP
                print arr[1],arr[10]
                nocall = nocall + 1

        # Populate the dictionary with Biomarks geno calls
        if arr[1] in dict_geno: # WGS geno already populated!
            dict_geno[arr[1]] = arr[10][0] + " " + arr[10][2] + " " + dict_geno[arr[1]]
        else:
            dict_geno[arr[1]] = arr[10][0] + " " + arr[10][2]
            

################################################
# Perform set operation to gather statistics
################################################
data = [list_WGS, list_Biomarks]
intersec = getIntersection(data)
uniq_WGS, uniq_Biomarks = getDifference(data)

# Calculate the percentage of genotype identity (excluding the Y SNPs!!!)
# because padding Y ref allele is not equal to "No call"
print len(intersec), intersec
#print "The concordance level of genotype between WGS and Biomarks is: ", len(intersec)/93.0, "!"
print("The concordance level is: {0:.2f}%!".format(len(intersec)/93.0*100))

#print "The No-call rate is: ", No_nocall/93.0
print("The call rate is: {0:.2f}%!".format((93-nocall)/93.0*100))

#print len(uniq_WGS), uniq_WGS
#print len(uniq_Biomarks), uniq_Biomarks

#################################################
# Print out the summary file for Cynthe's excel
#################################################

for key, value in dict_geno.iteritems():
    match = "N"
    if value[0:3] == value[-3:] or value[0:3] == value[::-1][0:3]: # if geno between WGS and Biomark match in unphased mode
        match = "Y"
    # Print the line of the snp
    if key[-1] != "Y":
        print key, value[0:3].replace(" ",":"), value[-3:].replace(" ",":"), match
