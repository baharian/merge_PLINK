import sys

# check for correct number of arguments
if (len(sys.argv) == 1) or (len(sys.argv) != 4):
	print 'ERROR: invalid input parameters\n'
	print 'this script checks for problematic SNPs that have different names but appear at the same position\n'
	print 'input file in .BIM format: {sample.bim}'
	print 'input file with positions of unnamed SNPs: {positions.list}'
	print 'output file with SNPs to exclude using PLINK: {SNPsToExclude.list}\n'
	print 'USAGE: removeUnnamedSNPs.py {sample.bim} {positions.list} {SNPsToExclude.list}'
	quit()

# loading the problematic positions from the .list file into memory
#print 'reading positions...'
pos = [x.replace('\n', '') for x in open(sys.argv[2], 'r')]

# loading the SNPs from the .bim file to memory and checking against the list of positions
#print 'validating SNPs...'
BIM = open(sys.argv[1], 'r')
excludeList = open(sys.argv[3], 'a')
for line in BIM:
	parts = line.replace('\n', '').split()

	if parts[3] in pos:
		excludeList.write(parts[1] + '\n')
excludeList.close()
BIM.close()
