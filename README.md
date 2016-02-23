# merge_PLINK
## Merge two genotype datasets in PLINK format

### General info
Bash and Python scripts to merge two genotype datasets in binary PLINK format while checking for inconsistencies in naming and missingness of data in any (or all) of the datasets.

  1. The pipeline is constructed based on [PLINK2](https://www.cog-genomics.org/plink2) program.
  2. The main shell script, `mergePLINKdata.sh`, requires the existence of the two Python scripts, `removeUnnamedSNPs.py` and `problemSNPs_PLINK.py`, in the same path as itself.

The pipeline first checks each dataset separately and chromosome-by-chromosome for SNPs that are denoted by multiple/different names (rsID, ...) at the same position and removes them, if any. Then, it unifies the SNP names across the datasets according to the criteria *chr<chr#>at<pos#>* (note that the original SNP names will be lost). Next,
