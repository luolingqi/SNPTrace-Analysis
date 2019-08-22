import argparse
parser = argparse.ArgumentParser()
#parser.add_argument("WGS_Genotype_file", help="The WGS genotype file generated from step2 script ")
parser.add_argument("Biomark_Genotype_file", help="The Biomark genotype file, containing snp calling lines for single sample")
args = parser.parse_args()

Gender_SNPs = ['hu98Y','hu111Y','hu209Y','hu103X','hu107X','hu109X']

def determine_gender(d_sex):
    sex = "NA"
    # determine female
    x = set()
    y = set()
    for key in d_sex:
        if key.endswith("X"):
            x.add(d_sex[key])
    for key in d_sex:
        if key.endswith("Y"):
            y.add(d_sex[key])
    #print "x, ", x
    #print "y, ", y
    
    # determine female
    if len(y) == 1 and list(y)[0] == "No Call:No Call":
        sex = "Female"
    # determine male
    if "No Call:No Call" not in y and all([len(set(i.split(":")))==1 for i in list(x)]): # check all 3 X-SNP are homozygous
        sex = "Male"
    if "No Call:No Call" not in y and any([len(set(i.split(":")))!=1 for i in list(x)]): # check any one of the 3 S-SNPs is heterozygous
        sex = "Klinefelter"
    return sex


list_Biomarks = []
dict_Sex = {}
sample_ID = ""
nocall_snps = []

nocall = 0
with open(args.Biomark_Genotype_file) as Biomarks_in:
    for line in Biomarks_in:
        line = line.strip()
        arr = line.split(",")
        # get sample ID info
        sample_ID = arr[0][0:3]+"_"+arr[4]
        
        if arr[9] != "No Call":
            list_Biomarks.append([arr[1], arr[9][0], arr[9][2]])
            if arr[1][-1] == 'X' or arr[1][-1] == 'Y':
                #print arr[1], arr[10]
                dict_Sex[arr[1]] = arr[9] # e.g.: hu98Y: G:T
        else:
            nocall_snps.append(arr[1])
            list_Biomarks.append([arr[1], "No Call", "No Call"])
            if arr[1][-1] == 'X' or arr[1][-1] == 'Y':
                dict_Sex[arr[1]] = "No Call"+":"+ "No Call"                         
            #print out the no call snps
            if arr[1][-1] != 'Y': # if not Y SNP
                #print arr[1],arr[9]
                nocall = nocall + 1

#print dict_Sex
gender = determine_gender(dict_Sex)
#print "The No-call rate is: ", No_nocall/93.0
#print("The non-Y call rate is: {0:.2f}%!".format((93-nocall)/93.0*100))
#print("The gender is: ", gender)

print sample_ID+ "\t" + "{0:.2f}%".format((93-nocall)/93.0*100) + "\t" + gender + "\t" + str(len(nocall_snps)) + "(" +  "_".join(nocall_snps) + ")"
