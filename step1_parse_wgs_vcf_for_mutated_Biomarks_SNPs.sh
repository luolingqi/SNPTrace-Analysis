#!/bin/bash

if [[ $# -ne 2 ]];then
    echo "Both the human genome version No. and gzipped WGS VCF are required only!!! "
    exit 1
fi

hg_ver=$1
wgs_vcf_gz=$2
snptrace="../Biomarks_96SNPs_w_refallel.txt"
bn=$(basename $wgs_vcf_gz)

zcat $wgs_vcf_gz > ./${bn:0:-3}

if [ $hg_ver = "hg19" ]; then
   awk '{print $4,$6}' $snptrace  # chromosome, position in hg19
elif [ $hg_ver = "hg38" ]; then
   awk '{print $7,$9}' $snptrace  # chromosome, position in hg38
fi   | while read line;

do
    chr=`echo $line | cut -d" " -f1`;
    pos=`echo $line | cut -d" " -f2`;
    echo -e "$chr\t$pos";
done | perl -pe 's/chr//g' | while read pair;

do
    grep "$pair" ./${bn:0:-3};  # grep the 96 SNPs from the ungzipped WGS VCF 
    
done | awk '{OFS="\t"}{print $1,$2,$4,$5,$9,$10}' | while read snp;

do
    echo -ne "$snp\t";grep `echo $snp | awk '{print $2}'` $snptrace | awk '{OFS="\t"}{print $2,$3}';
done > SNPs.from.${bn:0:-3}.txt

