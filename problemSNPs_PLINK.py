import sys

# check for correct number of arguments
if (len(sys.argv) == 1) or (len(sys.argv) != 6):
	print 'ERROR: invalid input parameters\n'
	print 'this script checks for strand inconsistency between two samples and outputs SNPs that need to be removed\n'
	print 'input files in .BIM format: {reference.bim} {sample.bim}'
	print 'output files with SNPs to be flipped and excluded using PLINK: {toFlip.list} {toUpdate.list} {toExclude.list}\n'
	print 'USAGE: fixStrand_PLINK.py {reference.bim} {sample.bim} {toFlip.list} {toUpdate.list} {toExclude.list}'
	quit()

def complement(char):
	if (char == 'A'):
		return 'T'
	elif (char == 'T'):
		return 'A'
	elif (char == 'G'):
		return 'C'
	else:
		return 'G'

# loading the reference SNPs from the first .bim file to memory
BIM1 = open(sys.argv[1], 'r')

print 'reading reference SNPs...'
SNPdict = dict()
for line in BIM1:
	parts = line.replace('\n', '').split()
	SNPdict.update({parts[1]:[parts[4], parts[5]]})

BIM1.close()

# comparing to (and, if necessary, fixing) the SNPs in the second .bim file
BIM2 = open(sys.argv[2], 'r')
toflip = open(sys.argv[3], 'w')
toupdate = open(sys.argv[4], 'w')
toremove = open(sys.argv[5], 'w')

print 'comparing sample SNPs...'
flip_Ref = 0			# e.g., (T/C in file 1 vs C/T in file 2) or (C/A in file 1 vs A/C in file 2)
flip_Ref_ambig = 0		# e.g., (T/A in file 1 vs A/T in file 2) or (C/G in file 1 vs G/C in file 2)
flip_Strand = 0			# e.g., (A/C in file 1 vs T/G in file 2)
flip_RefAndStrand = 0	# e.g., (A/C in file 1 vs G/T in file 2)
triAllelic = 0			# e.g., (T/G in file 1 vs T/A in file 2)
for line in BIM2:
	parts = line.replace('\n', '').split()

	BIM1alleles = SNPdict[parts[1]]
	BIM2alleles = [parts[4], parts[5]]

	if not ((BIM1alleles[0] == BIM2alleles[0]) and (BIM1alleles[1] == BIM2alleles[1])):
		if ((BIM1alleles[0] == BIM2alleles[1]) and (BIM1alleles[1] == BIM2alleles[0])):
			if BIM1alleles[0] == complement(BIM1alleles[1]):	# ambiguous flips; will be removed
				flip_Ref_ambig += 1
				toremove.write(parts[1] + '\n')
			else:	# ref/alt flips
				flip_Ref += 1
				toupdate.write(str.join('\t', [parts[1], BIM2alleles[0], BIM2alleles[1], BIM1alleles[0], BIM1alleles[1]]) + '\n')
		elif ((BIM1alleles[0] == complement(BIM2alleles[0])) and (BIM1alleles[1] == complement(BIM2alleles[1]))):	# unambiguous strand flips; can be fixed with PLINK
			flip_Strand += 1
			toflip.write(parts[1] + '\n')
		elif ((BIM1alleles[0] == complement(BIM2alleles[1])) and (BIM1alleles[1] == complement(BIM2alleles[0]))):	# unambiguous ref/alf + strand flips; will be removed
			flip_RefAndStrand += 1
			toremove.write(parts[1] + '\n')
		elif ((BIM1alleles[0] == BIM2alleles[0]) or (BIM1alleles[1] == BIM2alleles[1])):	# non-bi-allelic SNPs; will be removed
			triAllelic += 1
			toremove.write(parts[1] + '\n')
		elif ((BIM1alleles[0] == BIM2alleles[1]) or (BIM1alleles[1] == BIM2alleles[0])):
			triAllelic += 1
			toremove.write(parts[1] + '\n')
		else:	# e.g., (C/G in file 1 vs T/A in file 2)
			triAllelic += 1
			toremove.write(parts[1] + '\n')

toremove.close()
toupdate.close()
toflip.close()
BIM2.close()

print flip_Ref, "non-ambiguous reference/alternate mismatches detected (e.g., T/C vs C/T)"
print flip_Ref_ambig, "ambiguous reference/alternate or strand mismatches detected (e.g., A/T vs T/A)"
print flip_Strand, "strand mismatches detected (e.g., A/C vs T/G)"
print flip_RefAndStrand, "reference/alternate + strand mismatches detected (e.g., A/C vs G/T)"
print triAllelic, "tri(+)allelic SNPs detected (e.g., T/G vs T/A)\n"
