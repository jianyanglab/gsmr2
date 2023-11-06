# GSMR2
 
GSMR2 (Generalised Summary-data-based Mendelian Randomisation v2) is an improved version of [GSMR](https://github.com/jianyanglab/gsmr), which uses GWAS summary statistics to test for a putative causal association between two phenotypes (e.g., a modifiable risk factor and a disease) based on a multi-SNP model

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

or download the package [here](https://yanglab.westlake.edu.cn/software/gsmr/static/gsmr_1.1.1.tar.gz) and install it using the following command
```
install.packages("~/Downloads/gsmr_1.1.1.tar.gz", type = "source", repo = NULL)
```

## Example analysis pipeline

Here is an example of the `gsmr2` workflow to get started:

1. Prepare data for GSMR analysis
 1.1 Load the GWAS summary data
```
library("gsmr2")
data("gsmr")
head(gsmr_data)
```
 1.2 Estimate the LD correlation matrix
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

2. Prepare data for GSMR analysis
```
bzx = gsmr_data$std_bzx    # SNP effects on the risk factor
bzx_se = gsmr_data$std_bzx_se    # standard errors of bzx
bzx_pval = gsmr_data$bzx_pval   # p-values for bzx
bzy = gsmr_data$bzy    # SNP effects on the disease
bzy_se = gsmr_data$bzy_se    # standard errors of bzy
bzy_pval = gsmr_data$bzy_pval    # p-values for bzy
n_ref = 7703    # Sample size of the reference sample
gwas_thresh = 5e-8    # GWAS threshold to select SNPs as the instruments for the GSMR analysis
single_snp_heidi_thresh = 0.01    # p-value threshold for single-SNP-based HEIDI-outlier analysis
multi_snps_heidi_thresh = 0.01    # p-value threshold for multi-SNP-based HEIDI-outlier analysis
nsnps_thresh = 10   # the minimum number of instruments required for the GSMR analysis
heidi_outlier_flag = T    # flag for HEIDI-outlier analysis
ld_r2_thresh = 0.05    # LD r2 threshold to remove SNPs in high LD
ld_fdr_thresh = 0.05   # FDR threshold to remove the chance correlations between the SNP instruments

gsmr2_beta = 1     # 0 - the original HEIDI-outlier method; 1 - the new HEIDI-outlier method that is currently under development 

gsmr_results = gsmr(bzx, bzx_se, bzx_pval, bzy, bzy_se, bzy_pval, ldrho, snp_coeff_id, n_ref, heidi_outlier_flag, gwas_thresh, single_snp_heidi_thresh, multi_snps_heidi_thresh, nsnps_thresh, ld_r2_thresh, ld_fdr_thresh, gsmr2_beta)    # GSMR analysis 
filtered_index=gsmr_results$used_index


```

3. Check the results
```
cat("The estimated effect of the exposure on outcome: ",gsmr_results$bxy)
```
The estimated effect of the exposure on outcome:  0.4322395

```
cat("Standard error of bxy: ",gsmr_results$bxy_se)
```
Standard error of bxy:  0.02210985

# GCTA version
We also implement the C++ version of GSMR/GSMR2 within the [GCTA](https://github.com/jianyangqt/gcta). 

# Citation

**Angli Xue**, Zhihong Zhu, Huanwei Wang, Longda Jiang, Peter M. Visscher, Jian Zeng, Jian Yang. Unravelling the complex causal effects of substance use behaviours on common diseases. **Communications Medicine**. *Accepted in Principle*. [Preprint](https://www.researchsquare.com/article/rs-3465061/v1)

For questions, please email us at Jian Yang (jian.yang@westlake.edu.cn)

For bugs, please raise an issue in this repository and assign it to [@anglixue](https://github.com/anglixue).
