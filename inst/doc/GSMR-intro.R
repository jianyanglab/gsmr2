## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  # gsmr requires the R-package(s)
#  install.packages(c('survey'));
#  # install gsmr
#  install.packages("http://cnsgenomics.com/software/gsmr/static/gsmr_1.0.9.tar.gz",repos=NULL,type="source")

## -----------------------------------------------------------------------------
library("gsmr")
data("gsmr")
head(gsmr_data)
dim(gsmr_data)

## ----eval=FALSE---------------------------------------------------------------
#  # Save the genetic variants and effect alleles in a text file using R
#  write.table(gsmr_data[,c(1,2)], "gsmr_example_snps.allele", col.names=F, row.names=F, quote=F)
#  # Extract the genotype data from a GWAS dataset using GCTA
#  gcta64 --bfile gsmr_example --extract gsmr_example_snps.allele --update-ref-allele gsmr_example_snps.allele --recode --out gsmr_example

## ----eval=FALSE---------------------------------------------------------------
#  # Estimate LD correlation matrix using R
#  snp_coeff_id = scan("gsmr_example.xmat.gz", what="", nlines=1)
#  snp_coeff = read.table("gsmr_example.xmat.gz", header=F, skip=2)

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
dim(ldrho)
# Show the first 5 rows and columns of the matrix  
ldrho[1:5,1:5]

## ----message=FALSE, warning=FALSE---------------------------------------------
snpfreq = gsmr_data$a1_freq             # allele frequencies of the SNPs
bzx = gsmr_data$bzx     # effects of the instruments on risk factor
bzx_se = gsmr_data$bzx_se       # standard errors of bzx
bzx_n = gsmr_data$bzx_n          # GWAS sample size for the risk factor
std_zx = std_effect(snpfreq, bzx, bzx_se, bzx_n)    # perform standardisation
gsmr_data$std_bzx = std_zx$b    # standardized bzx
gsmr_data$std_bzx_se = std_zx$se    # standardized bzx_se
head(gsmr_data)

## ----message=FALSE, warning=FALSE---------------------------------------------
bzx = gsmr_data$std_bzx    # SNP effects on the risk factor
bzx_se = gsmr_data$std_bzx_se    # standard errors of bzx
bzx_pval = gsmr_data$bzx_pval   # p-values for bzx
bzy = gsmr_data$bzy    # SNP effects on the disease
bzy_se = gsmr_data$bzy_se    # standard errors of bzy
bzy_pval = gsmr_data$bzy_pval    # p-values for bzy
n_ref = 7703    # Sample size of the reference sample
gwas_thresh = 5e-8    # GWAS threshold to select SNPs as the instruments for the GSMR analysis
single_snp_heidi_thresh = 0.01    # p-value threshold for single-SNP-based HEIDI-outlier analysis
multi_snp_heidi_thresh = 0.01    # p-value threshold for multi-SNP-based HEIDI-outlier analysis
nsnps_thresh = 10   # the minimum number of instruments required for the GSMR analysis
heidi_outlier_flag = T    # flag for HEIDI-outlier analysis
ld_r2_thresh = 0.05    # LD r2 threshold to remove SNPs in high LD
ld_fdr_thresh = 0.05   # FDR threshold to remove the chance correlations between the SNP instruments
gsmr2_beta = 0     # 0 - the original HEIDI-outlier method; 1 - the new HEIDI-outlier method that is currently under development 
gsmr_results = gsmr(bzx, bzx_se, bzx_pval, bzy, bzy_se, bzy_pval, ldrho, snp_coeff_id, n_ref, heidi_outlier_flag, gwas_thresh, single_snp_heidi_thresh, multi_snp_heidi_thresh, nsnps_thresh, ld_r2_thresh, ld_fdr_thresh, gsmr2_beta)    # GSMR analysis 
filtered_index=gsmr_results$used_index
cat("The estimated effect of the exposure on outcome: ",gsmr_results$bxy)
cat("Standard error of bxy: ",gsmr_results$bxy_se)
cat("P-value for bxy: ", gsmr_results$bxy_pval)
cat("Indexes of the SNPs used in the GSMR analysis: ", gsmr_results$used_index[1:5], "...")
cat("Number of SNPs with missing estimates in the summary data: ", length(gsmr_results$na_snps))
cat("Number of non-significant SNPs: ", length(gsmr_results$weak_snps))
cat("Number of SNPs in high LD ( LD rsq >", ld_r2_thresh, "): ", length(gsmr_results$linkage_snps))
cat("Number of pleiotropic outliers: ", length(gsmr_results$pleio_snps))

## ----message=FALSE, warning=FALSE---------------------------------------------
gsmr_results = bi_gsmr(bzx, bzx_se, bzx_pval, bzy, bzy_se, bzy_pval, ldrho, snp_coeff_id, n_ref, heidi_outlier_flag, gwas_thresh, single_snp_heidi_thresh, multi_snp_heidi_thresh, nsnps_thresh, ld_r2_thresh, ld_fdr_thresh, gsmr2_beta)    # GSMR analysis 
cat("Effect of risk factor on disease: ",gsmr_results$forward_bxy)
cat("Standard error of bxy in the forward-GSMR analysis: ",gsmr_results$forward_bxy_se)
cat("P-value of bxy in the forward-GSMR analysis: ", gsmr_results$forward_bxy_pval)
cat("Effect of disease on risk factor: ",gsmr_results$reverse_bxy)
cat("Standard error of bxy in the reverse-GSMR analysis: ",gsmr_results$reverse_bxy_se)
cat("P-value of bxy in the reverse-GSMR analysis: ", gsmr_results$reverse_bxy_pval)

## ---- fig.width=6, fig.height=6.5---------------------------------------------
effect_col = colors()[75]
vals = c(bzx[filtered_index]-bzx_se[filtered_index], bzx[filtered_index]+bzx_se[filtered_index])
xmin = min(vals); xmax = max(vals)
vals = c(bzy[filtered_index]-bzy_se[filtered_index], bzy[filtered_index]+bzy_se[filtered_index])
ymin = min(vals); ymax = max(vals)
par(mar=c(5,5,4,2))
plot(bzx[filtered_index], bzy[filtered_index], pch=20, cex=0.8, bty="n", cex.axis=1.1, cex.lab=1.2,
        col=effect_col, xlim=c(xmin, xmax), ylim=c(ymin, ymax),
        xlab=expression(LDL~cholesterol~(italic(b[zx]))),
        ylab=expression(Coronary~artery~disease~(italic(b[zy]))))
abline(0, gsmr_results$forward_bxy, lwd=1.5, lty=2, col="dim grey")

nsnps = length(bzx[filtered_index])
for( i in 1:nsnps ) {
    # x axis
    xstart = bzx[filtered_index [i]] - bzx_se[filtered_index[i]]; xend = bzx[filtered_index[i]] + bzx_se[filtered_index[i]]
    ystart = bzy[filtered_index[i]]; yend = bzy[filtered_index[i]]
    segments(xstart, ystart, xend, yend, lwd=1.5, col=effect_col)
    # y axis
    xstart = bzx[filtered_index[i]]; xend = bzx[filtered_index[i]] 
    ystart = bzy[filtered_index[i]] - bzy_se[filtered_index[i]]; yend = bzy[filtered_index[i]] + bzy_se[filtered_index[i]]
    segments(xstart, ystart, xend, yend, lwd=1.5, col=effect_col)
}


