# merge_PLINK
## Merge two genotype datasets in PLINK format

### General info
Bash and Python scripts to merge two genotype datasets in binary PLINK format while checking for inconsistencies in naming and missingness of data in any (or all) of the datasets.

  1. The pipeline is constructed based on [PLINK2](https://www.cog-genomics.org/plink2) program.
  2. The main shell script, `mergePLINKdata.sh`, requires the existence of the two Python scripts, `removeUnnamedSNPs.py` and `problemSNPs_PLINK.py`, in the same path as itself.
