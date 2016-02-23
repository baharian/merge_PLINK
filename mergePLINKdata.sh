#!/bin/bash

# mandatory sanity check!
if [ $# -ne 3 ]
then
        printf "USAGE: mergePLINKdata.sh {dataSet1} {dataSet2} {outputName}\n"
        printf "\n"
        printf "{dataSet1} : prefix for the name of 1st set of .bed/.bim/.fam files\n"
        printf "{dataSet2} : prefix for the name of 2nd set of .bed/.bim/.fam files\n"
        printf "{outputName} : prefix for the output files' name\n"
        printf "\n"
        printf "example: ./mergePLINKdata.sh 1000Genomes HRS mergedData\n"
        exit 1
fi

file1=$1
file2=$2
output=$3

execPath="$(cd "$(dirname "$0")" && pwd)"

if [[ ! (-f ${file1}.bed && -f ${file1}.bim && -f ${file1}.fam && -f ${file2}.bed && -f ${file2}.bim && -f ${file2}.fam) ]]
then
	printf "ERROR: missing input file(s).\n"
	exit 1
fi



printf "\n\n-1) checking for problematic SNPs with different names at the same position...\n\n"

# some SNPs have two (or multiple) names at the same position
# e.g., rs35731280 and SNP1-246974301 on chromosome 1 at bp position 248907678 (data from the 1000 Genomes Project)
# e.g., rs3780548 and SNP9-98565308 on chromosome 9 at bp position 99525487 (data from the 1000 Genomes Project)
# note that they should be on the same chromosome to be considered problematic
# e.g., rs28415373 on chromosome 1 and SNP20-841981 on chromosome 20 appear at the same bp position 893981 (data from the 1000 Genomes Project)
# so the data should be looked at by chromosome
BIM1unnamedSNPs=0
for ((i = 0; i <= 26; i++))
do
	awk -v chr=${i} '{if ($1 == chr) print;}' ${file1}.bim > temp1_-2.bim
	unnamedSNPsInChr=$(cut -f4 temp1_-2.bim | sort -n | uniq -d | wc -l)
	
	if [ ${unnamedSNPsInChr} -ne 0 ]
	then
		(( BIM1unnamedSNPs += unnamedSNPsInChr ))
		cut -f4 temp1_-2.bim | sort -n | uniq -d > unnamedSNP_position_1.list
		python ${execPath}/removeUnnamedSNPs.py temp1_-2.bim unnamedSNP_position_1.list PLINK_exclude_unnamedSNPs_1.list
		rm unnamedSNP_position_1.list
	fi
	rm temp1_-2.bim
done

printf "...in dataset 1: found "${BIM1unnamedSNPs}" SNPs\n"
if [ ${BIM1unnamedSNPs} -ne 0 ]
then
	rename1=0
	plink --bfile ${file1} --exclude PLINK_exclude_unnamedSNPs_1.list --make-bed --out temp1_-1
else
	rename1=1
	mv ${file1}.bed temp1_-1.bed
	mv ${file1}.bim temp1_-1.bim
	mv ${file1}.fam temp1_-1.fam
fi

BIM2unnamedSNPs=0
for ((i = 0; i <= 26; i++))
do
	awk -v chr=${i} '{if ($1 == chr) print;}' ${file2}.bim > temp2_-2.bim
	unnamedSNPsInChr=$(cut -f4 temp2_-2.bim | sort -n | uniq -d | wc -l)
	
	if [ ${unnamedSNPsInChr} -ne 0 ]
	then
		(( BIM2unnamedSNPs += unnamedSNPsInChr ))
		cut -f4 temp2_-2.bim | sort -n | uniq -d > unnamedSNP_position_2.list
		python ${execPath}/removeUnnamedSNPs.py temp2_-2.bim unnamedSNP_position_2.list PLINK_exclude_unnamedSNPs_2.list
		rm unnamedSNP_position_2.list
	fi
	rm temp2_-2.bim
done

printf "...in dataset 2: found "${BIM2unnamedSNPs}" SNPs\n"
if [ ${BIM2unnamedSNPs} -ne 0 ]
then
	rename2=0
	plink --bfile ${file2} --exclude PLINK_exclude_unnamedSNPs_2.list --make-bed --out temp2_-1
else
	rename2=1
	mv ${file2}.bed temp2_-1.bed
	mv ${file2}.bim temp2_-1.bim
	mv ${file2}.fam temp2_-1.fam
fi



printf "\n\n0) unifying SNP names across two datasets...\n\n"

# giving SNPs unified names
mv temp1_-1.bim temp1_-1.bim.backup
awk 'BEGIN{FS="\t";OFS="\t"} {print $1,"chr"$1"at"$4,$3,$4,$5,$6;}' temp1_-1.bim.backup > temp1_-1.bim

mv temp2_-1.bim temp2_-1.bim.backup
awk 'BEGIN{FS="\t";OFS="\t"} {print $1,"chr"$1"at"$4,$3,$4,$5,$6;}' temp2_-1.bim.backup > temp2_-1.bim



printf "\n\n1) checking for missing calls...\n\n"

# find SNPs with missing calls in the first dataset
BIM1A1=$(cut -f5 temp1_-1.bim | grep 0 | wc -l)
BIM1A2=$(cut -f6 temp2_-1.bim | grep 0 | wc -l)

# remove those SNPs
if [ $((${BIM1A1}+${BIM1A2})) -ne 0 ]
then
	awk '{if ($5 == 0 || $6 == 0) print;}' temp1_-1.bim > PLINK_exclude_missing_1.list
	plink --bfile temp1_-1 --exclude PLINK_exclude_missing_1.list --make-bed --out temp1_0
fi

# find SNPs with missing calls in the second dataset
BIM2A1=$(cut -f5 temp2_-1.bim | grep 0 | wc -l)
BIM2A2=$(cut -f6 temp2_-1.bim | grep 0 | wc -l)

# remove those SNPs
if [ $((${BIM2A1}+${BIM2A2})) -ne 0 ]
then
	awk '{if ($5 == 0 || $6 == 0) print;}' temp2_-1.bim > PLINK_exclude_missing_2.list
	plink --bfile temp2_-1 --exclude PLINK_exclude_missing_2.list --make-bed --out temp2_0
fi



printf "\n\n2) finding the intersecting set of SNPs...\n\n"

# the remaining SNPs in the first dataset
if [ -f temp1_0.bim ]
then
	cut -f2 temp1_0.bim > SNPs_1.list
else
	cut -f2 temp1_-1.bim > SNPs_1.list
fi

# the remaining SNPs in the second dataset
if [ -f temp2_0.bim ]
then
	cut -f2 temp2_0.bim > SNPs_2.list
else
	cut -f2 temp2_-1.bim > SNPs_2.list
fi

# find the intersection of SNPs in the two datasets
grep -Fx -f SNPs_1.list SNPs_2.list > PLINK_extract_intersection.list
rm SNPs_*

# keep only those SNPs in the first dataset
if [ -f temp1_0.bim ]
then
	#plink --bfile temp1_0 --autosome --extract PLINK_extract_intersection.list --make-bed --out temp1_1
	plink --bfile temp1_0 --extract PLINK_extract_intersection.list --make-bed --out temp1_1
	rm temp1_0*
else
	#plink --bfile temp1_-1 --autosome --extract PLINK_extract_intersection.list --make-bed --out temp1_1
	plink --bfile temp1_-1 --extract PLINK_extract_intersection.list --make-bed --out temp1_1
fi

# keep only those SNPs in the second dataset
if [ -f temp2_0.bim ]
then
	#plink --bfile temp2_0 --autosome --extract PLINK_extract_intersection.list --make-bed --out temp2_1
	plink --bfile temp2_0 --extract PLINK_extract_intersection.list --make-bed --out temp2_1
	rm temp2_0*
else
	#plink --bfile temp2_-1 --autosome --extract PLINK_extract_intersection.list --make-bed --out temp2_1
	plink --bfile temp2_-1 --extract PLINK_extract_intersection.list --make-bed --out temp2_1
fi

if [ ${rename1} -eq 1 ]
then
	mv temp1_-1.fam ${file1}.fam
	mv temp1_-1.bed ${file1}.bed
	mv temp1_-1.bim.backup ${file1}.bim
	rm temp1_-1*
else
	rm temp1_-1*
fi

if [ ${rename2} -eq 1 ]
then
	mv temp2_-1.fam ${file2}.fam
	mv temp2_-1.bed ${file2}.bed
	mv temp2_-1.bim.backup ${file2}.bim
	rm temp2_-1*
else
	rm temp2_-1*
fi



printf "\n\n3) finding SNPs that could lead to problems...\n\n"

# find other SNPs that could lead to trouble (reference allele flips, strand flips, etc.)
python ${execPath}/problemSNPs_PLINK.py temp1_1.bim temp2_1.bim PLINK_flip_flips.list PLINK_update_flips.list PLINK_exclude_flips.list

# remove the potentially troublesome SNPs from the two datasets
plink --bfile temp1_1 --exclude PLINK_exclude_flips.list --make-bed --out temp1_2

plink --bfile temp2_1 --exclude PLINK_exclude_flips.list --flip PLINK_flip_flips.list --make-bed --out temp2_2
#plink --bfile temp2_1_temp1 --update-alleles PLINK_update_flips.list --make-bed --out temp2_1_temp2	# DO NOT USE! this line changes the data!
#plink --bfile temp2_1_temp2 --flip PLINK_flip_flips.list --make-bed --out temp2_2	# this line is merged with the one above; no need to use it anymore.

rm temp1_1*
rm temp2_1*



printf "\n\n4) merging...\n\n"

# merge the datasets into one
plink --bfile temp1_2 --bmerge temp2_2 --allow-no-sex --make-bed --out temp_3



printf "\n\n5) performing LD-based strand assignment cross-check of the merged data...\n\n"

# use LD to find incorrect strand assignments
plink --bfile temp_3 --make-pheno temp1_2.fam '*' --flip-scan --allow-no-sex --out PLINK_flipscan
awk '{if ($7 == 0 && $9 >= 1) print;}' PLINK_flipscan.flipscan | awk '{print $2;}' > PLINK_exclude_flipscan.list
plink --bfile temp_3 --exclude PLINK_exclude_flipscan.list --make-bed --out temp_4
rm temp1_2*
rm temp2_2*
rm temp_3*
rm PLINK_flipscan.???
rm PLINK_flipscan.?????



printf "\n\n6) performing QC...\n\n"

# clean the merged data
plink --bfile temp_4 --geno 0.0125 --make-bed --out temp_5
#plink --bfile temp_5 --mind 0.0125 --maf 0.01 --hwe 0.0001 --make-bed --out ${output}	# only do HWE test for randomly mating populations (i.e., without substructure)
plink --bfile temp_5 --mind 0.0125 --maf 0.005 --make-bed --out ${output}
rm temp_4*
rm temp_5*



#printf "\n\n7) pruning for LD in preparation for ADMIXTURE...\n\n"
#
# prune for LD
#plink --bfile ${output} --indep-pairwise 50 10 0.1 --out PLINK_extract_LD
#plink --bfile ${output} --extract PLINK_extract_LD.prune.in --make-bed --out ${output}_LD
#rm PLINK_extract_LD.???
#rm PLINK_extract_LD.?????
