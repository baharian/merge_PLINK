## Merge two genotype datasets in PLINK format

### General info
Bash and Python scripts to merge two genotype datasets in binary PLINK format while checking for inconsistencies in naming and missingness of data in any (or all) of the datasets.

  1. The pipeline is constructed based on [PLINK2](https://www.cog-genomics.org/plink2) program.
  2. The main shell script, `mergePLINKdata.sh`, requires the existence of the two Python scripts, `removeUnnamedSNPs.py` and `problemSNPs_PLINK.py`, in the same path as itself.

The pipeline first checks each dataset separately and chromosome-by-chromosome for SNPs that are denoted by multiple/different names (rsID, ...) at the same position and removes them, if any. Then, it unifies the SNP names across the datasets according to the criteria `chr<chr#>at<pos#>` (note that the original SNP names will be lost). Next, variants with missing calls are found in each dataset and removed. Then, the intersection of SNPs in the two datasets is found, and the remaining SNPs are excluded from the datasets. Next, two datasets will be cross-examined for SNPs that could be problematic due to strand misassignment, reference/alternate allele misassignment, ambiguous strand or reference/alternate allele misassignment, and non-bi-allelic variants. For each category, an appropriate action (flip or remove) is performed. This step ensures that there will be no merge issues when running PLINK. Now, the two datasets are merged. Next, an in-depth LD-based strand assignment cross-check of the merged data is performed to uncover potential strand flips which will be performed, if necessary. Finally, a small QC step is performed on the data.

### Usage
```
mergePLINKdata.sh {dataSet1} {dataSet2} {outputName}
{dataSet1} : prefix for the name of 1st set of .bed/.bim/.fam files
{dataSet2} : prefix for the name of 2nd set of .bed/.bim/.fam files
{outputName} : prefix for the output filesnames
```

Example: `mergePLINKdata.sh data1 data2 merged`
