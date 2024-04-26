# GSMR2
 
GSMR2 (Generalised Summary-data-based Mendelian Randomisation v2) is an improved version of [GSMR](https://github.com/jianyanglab/gsmr), which uses GWAS summary statistics to test for a putative causal association between two phenotypes (e.g., a modifiable risk factor and a disease) based on a multi-SNP model. This version implements a global heterogeneity test to remove invalid instrumental variables and provides a causal estimation that is more robust to directional pleiotropy.

## Installation

### Packages dependency

Please install the following package first. 

```{r}
# survey
install.packages('survey')
```

### Installation of gsmr2

Please run the following to install the `gsmr2` package:
```
devtools::install_github("jianyanglab/gsmr2")
```


## Example analysis pipeline

Here is an example of the `gsmr2` workflow to get started:

1. Load the GWAS summary data
```
library("gsmr2")
data("gsmr")
head(gsmr_data)
```
2. Estimate the LD correlation matrix
```
# Save the genetic variants and effect alleles in a text file using R
write.table(gsmr_data[,c(1,2)], "gsmr_example_snps.allele", col.names=F, row.names=F, quote=F)
# Extract the genotype data from a GWAS dataset using GCTA
gcta64 --bfile gsmr_example --extract gsmr_example_snps.allele --update-ref-allele gsmr_example_snps.allele --recode --out gsmr_example

# Estimate LD correlation matrix using R
snp_coeff_id = scan("gsmr_example.xmat.gz", what="", nlines=1)
snp_coeff = read.table("gsmr_example.xmat.gz", header=F, skip=2)

# Match the SNP genotype data with the summary data
snp_id = Reduce(intersect, list(gsmr_data$SNP, snp_coeff_id))
gsmr_data = gsmr_data[match(snp_id, gsmr_data$SNP),]
snp_order = match(snp_id, snp_coeff_id)
snp_coeff_id = snp_coeff_id[snp_order]
snp_coeff = snp_coeff[, snp_order]

# Calculate the LD correlation matrix
ldrho = cor(snp_coeff)

# Check the size of the correlation matrix and double-check if the order of the SNPs in the LD correlation matrix is consistent with that in the GWAS summary data
colnames(ldrho) = rownames(ldrho) = snp_coeff_id

```
If you don't have the genotype file, you can estimate the LD matrix from an independent dataset of the same ancestry.
Or, you can simply assume there is no LD between the SNPs by running the following
```
ldrho = diag(nrow(gsmr_data))
colnames(ldrho) = rownames(ldrho) = snp_coeff_id = snpid = as.character(gsmr_data$SNP)
```

3. Set the parameters for GSMR analysis
```
bzx = gsmr_data$bzx             # SNP effects on the risk factor
bzx_se = gsmr_data$bzx_se       # standard errors of bzx
bzx_pval = gsmr_data$bzx_pval   # p-values for bzx
bzy = gsmr_data$bzy             # SNP effects on the disease
bzy_se = gsmr_data$bzy_se       # standard errors of bzy
bzy_pval = gsmr_data$bzy_pval   # p-values for bzy
n_ref = 7703    # Sample size of the reference sample used to calculate the LD matrix
gwas_thresh = 5e-8              # GWAS threshold to select SNPs as the instruments for the GSMR analysis
multi_snps_heidi_thresh = 0.01  # p-value threshold for multi-SNP-based global HEIDI-outlier analysis
nsnps_thresh = 10               # the minimum number of instruments required for the GSMR analysis
heidi_outlier_flag = T          # flag for HEIDI-outlier analysis
ld_r2_thresh = 0.05             # LD r2 threshold to remove SNPs in high LD
ld_fdr_thresh = 0.05            # FDR threshold to remove the chance correlations between the SNP instruments
gsmr2_beta = 1                  # 0 - the original HEIDI-outlier method; 1 - the new HEIDI-outlier method that is currently under development 
```

4. Main analysis
```
gsmr_results = gsmr(bzx, bzx_se, bzx_pval, bzy, bzy_se, bzy_pval, ldrho, snp_coeff_id, n_ref, heidi_outlier_flag, gwas_thresh,
multi_snps_heidi_thresh, nsnps_thresh, ld_r2_thresh, ld_fdr_thresh, gsmr2_beta)    # GSMR analysis 
```

5. Check the results
```
cat("The estimated effect of the exposure on outcome: ", gsmr_results$bxy)
```
The estimated effect of the exposure on outcome:  0.3793622

```
cat("Standard error of bxy: ",gsmr_results$bxy_se)
```
Standard error of bxy:  0.02159656

# GCTA version
We also implement the C++ version of GSMR/GSMR2 within the [GCTA](https://github.com/jianyangqt/gcta). 

# Citation

**Angli Xue**, Zhihong Zhu, Huanwei Wang, Longda Jiang, Peter M. Visscher, Jian Zeng, Jian Yang. Unravelling the complex causal effects of substance use behaviours on common diseases. **Communications Medicine**. [Full Text](https://www.nature.com/articles/s43856-024-00473-3)

For questions, please email us at Jian Yang (jian.yang@westlake.edu.cn)

For bugs, please raise an issue in this repository and assign it to [@anglixue](https://github.com/anglixue).
