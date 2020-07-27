

####################################################################################################################################################
####################################################################################################################################################
## part18_ferrington_combined_differential_expression_bootstrap100_version_1.0.R
####################################################################################################################################################
####################################################################################################################################################
library("edgeR")
library("sva")
library("dplyr")

# Change global default setting so every data frame created will not auto-convert to factors unless explicitly instructed
options(stringsAsFactors = FALSE) 

# Useful functions
"%!in%" <- Negate("%in%")

####################################################################################################################################################
####################################################################################################################################################
## INPUT DATA
####################################################################################################################################################
####################################################################################################################################################
# Read in (1) RSEM gene-level expected counts data including Cufflinks RABT assembly results, (2) annotations 38.85 gene annotations 
# (originaly downloaded from Biomart within R), (3) metasheet, and (4) list of Eisenberg housekeeping genes that have been identified as suitable
# for differential expression analysis. 

input_filepath <- "/data/HumanRNAProject/ferrington_combined/18_amd_differential_expression/common_input"

# (1)
expected_counts <- read.delim(file = paste0(input_filepath, "/ferrington_combined_genelevel_expectedcounts_byrid_nooutliers.counts.matrix.tsv"), sep = "\t", check.names = FALSE)
rownames(expected_counts) <- expected_counts[, 1]
expected_counts <- expected_counts[, -1]
expected_counts <- expected_counts[order(rownames(expected_counts)), ] # order by gene IDs
expected_counts <- expected_counts[, order(colnames(expected_counts))] # order by r_IDs

# (2)
ensembl <- read.delim(file = paste0(input_filepath, "/hs_ensembl_85_genelevel_annotations.tsv"), sep = "\t", header = TRUE)
ensembl <- ensembl[ensembl$ensembl_gene_id %in% rownames(expected_counts), ] # keep only info for gene IDs in expected counts
ensembl <- ensembl[order(ensembl$ensembl_gene_id), ] # order by gene IDs
rownames(ensembl) <- ensembl$ensembl_gene_id

# (3)
metasheet <- read.csv(file = paste0(input_filepath, "/meta_ferrington_combined_inferior_retina_nooutliers_withbatches.csv"))
rownames(metasheet) <- metasheet$r_id
metasheet <- metasheet[, -1]
metasheet <- metasheet[order(rownames(metasheet)), ] # order by r_IDs
metasheet <- metasheet[rownames(metasheet) %in% colnames(expected_counts), ]

# (4) 
hkg <- read.table(file = paste0(input_filepath, "/eisenberg_hkg_for_supervised_sva_1cpm5.txt"))
all_hkg <- hkg$V1

####################################################################################################################################################
####################################################################################################################################################
## APPLY GENE FILTER
####################################################################################################################################################
####################################################################################################################################################
# Form DGEList object using sorted RSEM expected counts matrix, MGS group information, and annotation gene annotation
dge_expected_counts <- DGEList(expected_counts, group = metasheet$mgs_level, genes = ensembl)

# Transform expected counts to CPM values so that SVA/Combat can be used for batch correction (ie Combat works on logCPM values).
# logCPM is used for PCA plots. Note that the manual conversion of counts to CPM is: cpm <- apply(countmatrix,2, function(x) (x/sum(x))*1000000). 
# Expected count of 0 equals CPM of 0.
tmmnorm <- calcNormFactors(dge_expected_counts, method = "TMM") # TMM normalization
cpm_tmmnorm <- cpm(tmmnorm, log = FALSE) # CPM transformation

# Filter. Keep genes with CPM >= 1 in at least 10% of all samples
number_of_samples <- dim(dge_expected_counts$samples)[1]

genes_meeting_filter <- rowSums(cpm_tmmnorm >= 1) >= 0.05 * number_of_samples
#table(genes_meeting_filter) ["TRUE"] 
#table(genes_meeting_filter) ["FALSE"] 

filtered_dge_expected_counts <- dge_expected_counts[genes_meeting_filter, keep.lib.size = FALSE] # make sure to recalculate library size
#dim(filtered_dge_expected_counts) 

# Normalize and transform into CPM
filtered_tmmnorm <- calcNormFactors(filtered_dge_expected_counts)
filtered_cpm_tmmnorm <- cpm(filtered_tmmnorm, log = FALSE)
#dim(filtered_cpm_tmmnorm) 

####################################################################################################################################################
####################################################################################################################################################
## DIFFERENTIAL EXPRESSION WITH BOOTSTRAP
####################################################################################################################################################
####################################################################################################################################################
# Set up model matrices. The null model consists of the known variables and covariates that must be included as adjustment variables. Age is included
# as a continuous covariate variable because our outcomes of interest are differentially expressed genes due to AMD as well as age.
# The full model includes all the variables in the null model, as well as the variable of interest.
mod0 <- model.matrix(~0 + scale(age), data = metasheet)
mod1 <- model.matrix(~0 + factor(metasheet$mgs_level, levels = c("1", "2", "3", "4")) + scale(age), data = metasheet)

controls <- replicate(100, sample(all_hkg, 1000, replace = FALSE))

kruskal21_chisquared <- list()
kruskal21_df <- list()
kruskal21_pvalues <- list()

kruskal_controls21_chisquared <- list()
kruskal_controls21_df <- list()
kruskal_controls21_pvalues <- list()

kruskal31_chisquared <- list()
kruskal31_df <- list()
kruskal31_pvalues <- list()

kruskal_controls31_chisquared <- list()
kruskal_controls31_df <- list()
kruskal_controls31_pvalues <- list()

kruskal41_chisquared <- list()
kruskal41_df <- list()
kruskal41_pvalues <- list()

kruskal_controls41_chisquared <- list()
kruskal_controls41_df <- list()
kruskal_controls41_pvalues <- list()

kruskal32_chisquared <- list()
kruskal32_df <- list()
kruskal32_pvalues <- list()

kruskal_controls32_chisquared <- list()
kruskal_controls32_df <- list()
kruskal_controls32_pvalues <- list()

kruskal43_chisquared <- list()
kruskal43_df <- list()
kruskal43_pvalues <- list()

kruskal_controls43_chisquared <- list()
kruskal_controls43_df <- list()
kruskal_controls43_pvalues <- list()

common_filepath <- "/data/HumanRNAProject/ferrington_combined/18_amd_differential_expression"
version <- "/1cpm5/version1.0/ferrington_combined_amd_1cpm5_de1.0"

for (i in 1:100) {
  sva_controls <- as.data.frame(rownames(filtered_cpm_tmmnorm))
  sva_controls$control <- NA
  sva_controls$control[sva_controls$`rownames(filtered_cpm_tmmnorm)` %in% controls[ , i]] <- TRUE
  sva_controls$control[sva_controls$`rownames(filtered_cpm_tmmnorm)` %!in% controls[ , i]] <- FALSE

  batch_sup_sva <- svaseq(filtered_cpm_tmmnorm, mod = mod1, mod0 = mod0, controls = sva_controls$control, method = "supervised")

  # Set up the new model matrix by using all surrogate variables
  mod_svaseq <- cbind(mod1, batch_sup_sva$sv)

  # Use voom to convert read counts into log2-cpm with associated weights so that linear modeling can then be performed
  voom_after_tmmnorm_with_model <- voom(filtered_tmmnorm, design = mod_svaseq, plot = FALSE) # after TMM normalization and with model_matrix

  # Find genes that are differentially expressed
  fit <- lmFit(voom_after_tmmnorm_with_model, design = mod_svaseq) # fit linear model
  
  filename <- paste0(common_filepath, version, "_bootstrap", i, "_plots.pdf")
  pdf(file = filename)
  fit21 <- contrasts.fit(fit,  contrasts = c(-1, 1, 0, 0, rep(0, ncol(mod_svaseq) - 4)))
  fit21 <- eBayes(fit21, robust = TRUE) # empirical Bayes statistics and differential expression assessment
  hist(fit21$p.value, breaks = 20)
  qqt(fit21$t, df = fit21$df.total, pch = 16, cex = 0.2, xlab = "Theoretical Quantiles of Standard Normal", ylab = "t") # moderated t-statistics
  abline(0,1)
  observed21 <- as.vector(fit21$p.value)
  theoretical21 <- qunif(ppoints(length(observed21)))
  qqplot(-log10(theoretical21), -log10(observed21), main = "P-value Q-Q Plot", xlab = "-log10 Theoretical Uniform", ylab = "-log10 Observed p-values")
  abline(0, 1)
  kruskal21 <- kruskal.test(list(theoretical21, observed21))
  mtext(text = paste0(kruskal21$method,
                      ": chi-squared = ", round(kruskal21$statistic, digits = 3),
                      "  df = ", kruskal21$parameter,
                      "  p-value = ", round(kruskal21$p.value, digits = 7)), cex = 0.8)
  abline(0,1)

  only_controls21 <- fit21$genes$ensembl_gene_id %in% controls[ , i]
  controls_fit21 <- fit21[only_controls21, ]
  hist(controls_fit21$p.value, breaks = 20)
  observed_controls21 <- as.vector(controls_fit21$p.value)
  theoretical_controls21 <- qunif(ppoints(length(observed_controls21)))
  qqplot(-log10(theoretical_controls21), -log10(observed_controls21), main = "P-value Q-Q Plot", xlab = "-log10 Theoretical Uniform", ylab = "-log10 Observed p-values")
  abline(0, 1)
  kruskal_controls21 <- kruskal.test(list(theoretical_controls21, observed_controls21))
  mtext(text = paste0(kruskal_controls21$method,
                      ": chi-squared = ", round(kruskal_controls21$statistic, digits = 3),
                      "  df = ", kruskal_controls21$parameter,
                      "  p-value = ", round(kruskal_controls21$p.value, digits = 7)), cex = 0.8)
  abline(0,1)
  
  kruskal21_chisquared[[i]] <- kruskal21$statistic
  kruskal21_df[[i]] <- kruskal21$parameter
  kruskal21_pvalues[[i]] <- kruskal21$p.value
  
  kruskal_controls21_chisquared[[i]] <- kruskal_controls21$statistic
  kruskal_controls21_df[[i]] <- kruskal_controls21$parameter
  kruskal_controls21_pvalues[[i]] <- kruskal_controls21$p.value

  fit31 <- contrasts.fit(fit, contrasts = c(-1, 0, 1, 0, rep(0, ncol(mod_svaseq) - 4)))
  fit31 <- eBayes(fit31, robust = TRUE) # empirical Bayes statistics and differential expression assessment
  hist(fit31$p.value, breaks = 20)
  qqt(fit31$t, df = fit31$df.total, pch = 16, cex = 0.2, xlab = "Theoretical Quantiles of Standard Normal", ylab = "t") # moderated t-statistics
  abline(0,1)
  observed31 <- as.vector(fit31$p.value)
  theoretical31 <- qunif(ppoints(length(observed31)))
  qqplot(-log10(theoretical31), -log10(observed31), main = "P-value Q-Q Plot", xlab = "-log10 Theoretical Uniform", ylab = "-log10 Observed p-values")
  abline(0, 1)
  kruskal31 <- kruskal.test(list(theoretical31, observed31))
  mtext(text = paste0(kruskal31$method,
                      ": chi-squared = ", round(kruskal31$statistic, digits = 3),
                      "  df = ", kruskal31$parameter,
                      "  p-value = ", round(kruskal31$p.value, digits = 7)), cex = 0.8)
  abline(0,1)
  
  
  only_controls31 <- fit31$genes$ensembl_gene_id %in% controls[ , i]
  controls_fit31 <- fit31[only_controls31, ]
  hist(controls_fit31$p.value, breaks = 20)
  observed_controls31 <- as.vector(controls_fit31$p.value)
  theoretical_controls31 <- qunif(ppoints(length(observed_controls31)))
  qqplot(-log10(theoretical_controls31), -log10(observed_controls31), main = "P-value Q-Q Plot", xlab = "-log10 Theoretical Uniform", ylab = "-log10 Observed p-values")
  abline(0, 1)
  kruskal_controls31 <- kruskal.test(list(theoretical_controls31, observed_controls31))
  mtext(text = paste0(kruskal_controls31$method,
                      ": chi-squared = ", round(kruskal_controls31$statistic, digits = 3),
                      "  df = ", kruskal_controls31$parameter,
                      "  p-value = ", round(kruskal_controls31$p.value, digits = 7)), cex = 0.8)
  abline(0,1)

  kruskal31_chisquared[[i]] <- kruskal31$statistic
  kruskal31_df[[i]] <- kruskal31$parameter
  kruskal31_pvalues[[i]] <- kruskal31$p.value
  
  kruskal_controls31_chisquared[[i]] <- kruskal_controls31$statistic
  kruskal_controls31_df[[i]] <- kruskal_controls31$parameter
  kruskal_controls31_pvalues[[i]] <- kruskal_controls31$p.value

  fit41 <- contrasts.fit(fit, contrasts = c(-1, 0, 0, 1, rep(0, ncol(mod_svaseq) - 4)))
  fit41 <- eBayes(fit41, robust = TRUE) # empirical Bayes statistics and differential expression assessment
  hist(fit41$p.value, breaks = 20)
  qqt(fit41$t, df = fit41$df.total, pch = 16, cex = 0.2, xlab = "Theoretical Quantiles of Standard Normal", ylab = "t") # moderated t-statistics
  abline(0,1)
  observed41 <- as.vector(fit41$p.value)
  theoretical41 <- qunif(ppoints(length(observed41)))
  qqplot(-log10(theoretical41), -log10(observed41), main = "P-value Q-Q Plot", xlab = "-log10 Theoretical Uniform", ylab = "-log10 Observed p-values")
  abline(0, 1)
  kruskal41 <- kruskal.test(list(theoretical41, observed41))
  mtext(text = paste0(kruskal41$method,
                      ": chi-squared = ", round(kruskal41$statistic, digits = 3),
                      "  df = ", kruskal41$parameter,
                      "  p-value = ", round(kruskal41$p.value, digits = 7)), cex = 0.8)
  abline(0,1)
  
  only_controls41 <- fit41$genes$ensembl_gene_id %in% controls[ , i]
  controls_fit41 <- fit41[only_controls41, ]
  hist(controls_fit41$p.value, breaks = 20)
  observed_controls41 <- as.vector(controls_fit41$p.value)
  theoretical_controls41 <- qunif(ppoints(length(observed_controls41)))
  qqplot(-log10(theoretical_controls41), -log10(observed_controls41), main = "P-value Q-Q Plot", xlab = "-log10 Theoretical Uniform", ylab = "-log10 Observed p-values")
  abline(0, 1)
  kruskal_controls41 <- kruskal.test(list(theoretical_controls41, observed_controls41))
  mtext(text = paste0(kruskal_controls41$method,
                      ": chi-squared = ", round(kruskal_controls41$statistic, digits = 3),
                      "  df = ", kruskal_controls41$parameter,
                      "  p-value = ", round(kruskal_controls41$p.value, digits = 7)), cex = 0.8)
  abline(0,1)
  
  kruskal41_chisquared[[i]] <- kruskal41$statistic
  kruskal41_df[[i]] <- kruskal41$parameter
  kruskal41_pvalues[[i]] <- kruskal41$p.value
  
  kruskal_controls41_chisquared[[i]] <- kruskal_controls41$statistic
  kruskal_controls41_df[[i]] <- kruskal_controls41$parameter
  kruskal_controls41_pvalues[[i]] <- kruskal_controls41$p.value
  
  fit32 <- contrasts.fit(fit, contrasts = c(0, -1, 1, 0, rep(0, ncol(mod_svaseq) - 4)))
  fit32 <- eBayes(fit32, robust = TRUE) # empirical Bayes statistics and differential expression assessment
  hist(fit32$p.value, breaks = 20)
  qqt(fit32$t, df = fit32$df.total, pch = 16, cex = 0.2, xlab = "Theoretical Quantiles of Standard Normal", ylab = "t") # moderated t-statistics)
  abline(0,1)
  observed32 <- as.vector(fit32$p.value)
  theoretical32 <- qunif(ppoints(length(observed32)))
  qqplot(-log10(theoretical32), -log10(observed32), main = "P-value Q-Q Plot", xlab = "-log10 Theoretical Uniform", ylab = "-log10 Observed p-values")
  abline(0, 1)
  kruskal32 <- kruskal.test(list(theoretical32, observed32))
  mtext(text = paste0(kruskal32$method,
                      ": chi-squared = ", round(kruskal32$statistic, digits = 3),
                      "  df = ", kruskal32$parameter,
                      "  p-value = ", round(kruskal32$p.value, digits = 7)), cex = 0.8)
  abline(0,1)
  
  only_controls32 <- fit32$genes$ensembl_gene_id %in% controls[ , i]
  controls_fit32 <- fit32[only_controls32, ]
  hist(controls_fit32$p.value, breaks = 20)
  observed_controls32 <- as.vector(controls_fit32$p.value)
  theoretical_controls32 <- qunif(ppoints(length(observed_controls32)))
  qqplot(-log10(theoretical_controls32), -log10(observed_controls32), main = "P-value Q-Q Plot", xlab = "-log10 Theoretical Uniform", ylab = "-log10 Observed p-values")
  abline(0, 1)
  kruskal_controls32 <- kruskal.test(list(theoretical_controls32, observed_controls32))
  mtext(text = paste0(kruskal_controls32$method,
                      ": chi-squared = ", round(kruskal_controls32$statistic, digits = 3),
                      "  df = ", kruskal_controls32$parameter,
                      "  p-value = ", round(kruskal_controls32$p.value, digits = 7)), cex = 0.8)
  abline(0,1)
  
  kruskal32_chisquared[[i]] <- kruskal32$statistic
  kruskal32_df[[i]] <- kruskal32$parameter
  kruskal32_pvalues[[i]] <- kruskal32$p.value
  
  kruskal_controls32_chisquared[[i]] <- kruskal_controls32$statistic
  kruskal_controls32_df[[i]] <- kruskal_controls32$parameter
  kruskal_controls32_pvalues[[i]] <- kruskal_controls32$p.value
  
  fit43 <- contrasts.fit(fit, contrasts = c(0, 0, -1, 1, rep(0, ncol(mod_svaseq) - 4)))
  fit43 <- eBayes(fit43, robust = TRUE) # empirical Bayes statistics and differential expression assessment
  hist(fit43$p.value, breaks = 20)
  qqt(fit43$t, df = fit43$df.total, pch = 16, cex = 0.2, xlab = "Theoretical Quantiles of Standard Normal", ylab = "t") # moderated t-statistics)
  abline(0,1)
  observed43 <- as.vector(fit43$p.value)
  theoretical43 <- qunif(ppoints(length(observed43)))
  qqplot(-log10(theoretical43), -log10(observed43), main = "P-value Q-Q Plot", xlab = "-log10 Theoretical Uniform", ylab = "-log10 Observed p-values")
  abline(0, 1)
  kruskal43 <- kruskal.test(list(theoretical43, observed43))
  mtext(text = paste0(kruskal43$method,
                      ": chi-squared = ", round(kruskal43$statistic, digits = 3),
                      "  df = ", kruskal43$parameter,
                      "  p-value = ", round(kruskal43$p.value, digits = 7)), cex = 0.8)
  abline(0,1)
  
  only_controls43 <- fit43$genes$ensembl_gene_id %in% controls[ , i]
  controls_fit43 <- fit43[only_controls43, ]
  hist(controls_fit43$p.value, breaks = 20)
  observed_controls43 <- as.vector(controls_fit43$p.value)
  theoretical_controls43 <- qunif(ppoints(length(observed_controls43)))
  qqplot(-log10(theoretical_controls43), -log10(observed_controls43), main = "P-value Q-Q Plot", xlab = "-log10 Theoretical Uniform", ylab = "-log10 Observed p-values")
  abline(0, 1)
  kruskal_controls43 <- kruskal.test(list(theoretical_controls43, observed_controls43))
  mtext(text = paste0(kruskal_controls43$method,
                      ": chi-squared = ", round(kruskal_controls43$statistic, digits = 3),
                      "  df = ", kruskal_controls43$parameter,
                      "  p-value = ", round(kruskal_controls43$p.value, digits = 7)), cex = 0.8)
  abline(0,1)

  kruskal43_chisquared[[i]] <- kruskal43$statistic
  kruskal43_df[[i]] <- kruskal43$parameter
  kruskal43_pvalues[[i]] <- kruskal43$p.value
  
  kruskal_controls43_chisquared[[i]] <- kruskal_controls43$statistic
  kruskal_controls43_df[[i]] <- kruskal_controls43$parameter
  kruskal_controls43_pvalues[[i]] <- kruskal_controls43$p.value
  
  dev.off()

####################################################################################################################################################
####################################################################################################################################################
## OUTPUT DIFFERENTIAL EXPRESSION RESULTS AND KRUSKAL WALLIS TESTS
####################################################################################################################################################
####################################################################################################################################################
  # Obtain output for significant genes
  table21 <- topTable(fit21, number = Inf, coef = 1) # p.value is actually the adjusted p-value
  table21$FC <- 2^table21$logFC
  table21 <- table21[ , c(1:8, 14, 9:13)]
  table21 <- table21[order(table21$P.Value), ]
  write.table(table21, file = paste0(common_filepath, version, "_bootstrap", i, "_mgs21.txt"), row.names = FALSE, quote = FALSE)
  
  table31 <- topTable(fit31, number = Inf, coef = 1) # p.value is actually the adjusted p-value
  table31$FC <- 2^table31$logFC
  table31 <- table31[ , c(1:8, 14, 9:13)]
  table31 <- table31[order(table31$P.Value), ]
  write.table(table31, file = paste0(common_filepath, version, "_bootstrap", i, "_mgs31.txt"), row.names = FALSE, quote = FALSE)
  
  table41 <- topTable(fit41, number = Inf, coef = 1) # p.value is actually the adjusted p-value
  table41$FC <- 2^table41$logFC
  table41 <- table41[ , c(1:8, 14, 9:13)]
  table41 <- table41[order(table41$P.Value), ]
  write.table(table41, file = paste0(common_filepath, version, "_bootstrap", i, "_mgs41.txt"), row.names = FALSE, quote = FALSE)
  
  table32 <- topTable(fit32, number = Inf, coef = 1) # p.value is actually the adjusted p-value
  table32$FC <- 2^table32$logFC
  table32 <- table32[ , c(1:8, 14, 9:13)]
  table32 <- table32[order(table32$P.Value), ]
  write.table(table32, file = paste0(common_filepath, version, "_bootstrap", i, "_mgs32.txt"), row.names = FALSE, quote = FALSE)
  
  table43 <- topTable(fit43, number = Inf, coef = 1) # p.value is actually the adjusted p-value
  table43$FC <- 2^table43$logFC
  table43 <- table43[ , c(1:8, 14, 9:13)]
  table43 <- table43[order(table43$P.Value), ]
  write.table(table43, file = paste0(common_filepath, version, "_bootstrap", i, "_mgs43.txt"), row.names = FALSE, quote = FALSE)
}

# Obtain output for Kruskal-Wallis tests
kruskal21_output <- data.frame(
  kruskal_chisquared = unlist(kruskal21_chisquared),
  kruskal_df = unlist(kruskal21_df),
  kruskal_pvalues = unlist(kruskal21_pvalues))
write.table(kruskal21_output, file = paste0(common_filepath, version, "_mgs21_kruskalstats.txt"), row.names = FALSE, quote = FALSE)

kruskal31_output <- data.frame(
  kruskal_chisquared = unlist(kruskal31_chisquared),
  kruskal_df = unlist(kruskal31_df),
  kruskal_pvalues = unlist(kruskal31_pvalues))
write.table(kruskal31_output, file = paste0(common_filepath, version, "_mgs31_kruskalstats.txt"), row.names = FALSE, quote = FALSE)

kruskal41_output <- data.frame(
  kruskal_chisquared = unlist(kruskal41_chisquared),
  kruskal_df = unlist(kruskal41_df),
  kruskal_pvalues = unlist(kruskal41_pvalues))
write.table(kruskal41_output, file = paste0(common_filepath, version, "_mgs41_kruskalstats.txt"), row.names = FALSE, quote = FALSE)

kruskal32_output <- data.frame(
  kruskal_chisquared = unlist(kruskal32_chisquared),
  kruskal_df = unlist(kruskal32_df),
  kruskal_pvalues = unlist(kruskal32_pvalues))
write.table(kruskal32_output, file = paste0(common_filepath, version, "_mgs32_kruskalstats.txt"), row.names = FALSE, quote = FALSE)

kruskal43_output <- data.frame(
  kruskal_chisquared = unlist(kruskal43_chisquared),
  kruskal_df = unlist(kruskal43_df),
  kruskal_pvalues = unlist(kruskal43_pvalues))
write.table(kruskal43_output, file = paste0(common_filepath, version, "_mgs43_kruskalstats.txt"), row.names = FALSE, quote = FALSE)

kruskal_controls21_output <- data.frame(
  kruskal_chisquared = unlist(kruskal_controls21_chisquared),
  kruskal_df = unlist(kruskal_controls21_df),
  kruskal_pvalues = unlist(kruskal_controls21_pvalues))
write.table(kruskal_controls21_output, file = paste0(common_filepath, version, "_mgs21_controlkruskalstats.txt"), row.names = FALSE, quote = FALSE)

kruskal_controls31_output <- data.frame(
  kruskal_chisquared = unlist(kruskal_controls31_chisquared),
  kruskal_df = unlist(kruskal_controls31_df),
  kruskal_pvalues = unlist(kruskal_controls31_pvalues))
write.table(kruskal_controls31_output, file = paste0(common_filepath, version, "_mgs31_controlkruskalstats.txt"), row.names = FALSE, quote = FALSE)

kruskal_controls41_output <- data.frame(
  kruskal_chisquared = unlist(kruskal_controls41_chisquared),
  kruskal_df = unlist(kruskal_controls41_df),
  kruskal_pvalues = unlist(kruskal_controls41_pvalues))
write.table(kruskal_controls41_output, file = paste0(common_filepath, version, "_mgs41_controlkruskalstats.txt"), row.names = FALSE, quote = FALSE)

kruskal_controls32_output <- data.frame(
  kruskal_chisquared = unlist(kruskal_controls32_chisquared),
  kruskal_df = unlist(kruskal_controls32_df),
  kruskal_pvalues = unlist(kruskal_controls32_pvalues))
write.table(kruskal_controls32_output, file = paste0(common_filepath, version, "_mgs32_controlkruskalstats.txt"), row.names = FALSE, quote = FALSE)

kruskal_controls43_output <- data.frame(
  kruskal_chisquared = unlist(kruskal_controls43_chisquared),
  kruskal_df = unlist(kruskal_controls43_df),
  kruskal_pvalues = unlist(kruskal_controls43_pvalues))
write.table(kruskal_controls43_output, file = paste0(common_filepath, version, "_mgs43_controlkruskalstats.txt"), row.names = FALSE, quote = FALSE)

####################################################################################################################################################
####################################################################################################################################################
## SUMMARY STATS
####################################################################################################################################################
####################################################################################################################################################
# SUMMARY STATS FOR MGS 2 VS MGS 1
summarystats_filepath <- "/data/HumanRNAProject/ferrington_combined/18_amd_differential_expression/1cpm5/version1.0/"

files <- list.files(path = summarystats_filepath, pattern = "_mgs21.txt")

data <- do.call("rbind", lapply(files, function(x) read.table(paste0(summarystats_filepath, x), header = TRUE, stringsAsFactors = FALSE)))

summary_stats <- data %>%
  group_by(ensembl_gene_id) %>%
  summarize(
            external_gene_name = unique(external_gene_name),
            strand = unique(strand),
            chromosome_name = unique(chromosome_name),
            start_position = unique(start_position),
            end_position = unique(end_position),
            gene_biotype = unique(gene_biotype),
            mean_logFC = mean(logFC),
            stdev_logFC = sd(logFC),
            mean_FC = mean(FC),
            stdev_FC = sd(FC),
            mean_AveExpr = mean(AveExpr),
            stdev_AveExpr = sd(AveExpr),
            mean_t = mean(t),
            stdev_t = sd(t),
            mean_P.Value = mean(P.Value),
            stdev_P.Value = sd(P.Value),
            mean_adj.P.Value = mean(adj.P.Val),
            stdev_adj.P.Value = sd(adj.P.Val),
            mean_B = mean(B),
            stdev_B = sd(B))

write.table(summary_stats, file = paste0(common_filepath, version, "_mgs21_summarystats.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE)

# SUMMARY STATS FOR MGS 3 VS MGS 1
files <- list.files(path = summarystats_filepath, pattern = "_mgs31.txt")

data <- do.call("rbind", lapply(files, function(x) read.table(paste0(summarystats_filepath, x), header = TRUE, stringsAsFactors = FALSE)))

summary_stats <- data %>%
  group_by(ensembl_gene_id) %>%
  summarize(
            external_gene_name = unique(external_gene_name),
            strand = unique(strand),
            chromosome_name = unique(chromosome_name),
            start_position = unique(start_position),
            end_position = unique(end_position),
            gene_biotype = unique(gene_biotype),
            mean_logFC = mean(logFC),
            stdev_logFC = sd(logFC),
            mean_FC = mean(FC),
            stdev_FC = sd(FC),
            mean_AveExpr = mean(AveExpr),
            stdev_AveExpr = sd(AveExpr),
            mean_t = mean(t),
            stdev_t = sd(t),
            mean_P.Value = mean(P.Value),
            stdev_P.Value = sd(P.Value),
            mean_adj.P.Value = mean(adj.P.Val),
            stdev_adj.P.Value = sd(adj.P.Val),
            mean_B = mean(B),
            stdev_B = sd(B))

write.table(summary_stats, file = paste0(common_filepath, version, "_mgs31_summarystats.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE)

# SUMMARY STATS FOR MGS 4 VS MGS 1
files <- list.files(path = summarystats_filepath, pattern = "_mgs41.txt")

data <- do.call("rbind", lapply(files, function(x) read.table(paste0(summarystats_filepath, x), header = TRUE, stringsAsFactors = FALSE)))

summary_stats <- data %>%
  group_by(ensembl_gene_id) %>%
  summarize(
            external_gene_name = unique(external_gene_name),
            strand = unique(strand),
            chromosome_name = unique(chromosome_name),
            start_position = unique(start_position),
            end_position = unique(end_position),
            gene_biotype = unique(gene_biotype),
            mean_logFC = mean(logFC),
            stdev_logFC = sd(logFC),
            mean_FC = mean(FC),
            stdev_FC = sd(FC),
            mean_AveExpr = mean(AveExpr),
            stdev_AveExpr = sd(AveExpr),
            mean_t = mean(t),
            stdev_t = sd(t),
            mean_P.Value = mean(P.Value),
            stdev_P.Value = sd(P.Value),
            mean_adj.P.Value = mean(adj.P.Val),
            stdev_adj.P.Value = sd(adj.P.Val),
            mean_B = mean(B),
            stdev_B = sd(B))

write.table(summary_stats, file = paste0(common_filepath, version, "_mgs41_summarystats.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE)

# SUMMARY STATS FOR MGS 3 VS MGS 2
files <- list.files(path = summarystats_filepath, pattern = "_mgs32.txt")

data <- do.call("rbind", lapply(files, function(x) read.table(paste0(summarystats_filepath, x), header = TRUE, stringsAsFactors = FALSE)))

summary_stats <- data %>%
  group_by(ensembl_gene_id) %>%
  summarize(
            external_gene_name = unique(external_gene_name),
            strand = unique(strand),
            chromosome_name = unique(chromosome_name),
            start_position = unique(start_position),
            end_position = unique(end_position),
            gene_biotype = unique(gene_biotype),
            mean_logFC = mean(logFC),
            stdev_logFC = sd(logFC),
            mean_FC = mean(FC),
            stdev_FC = sd(FC),
            mean_AveExpr = mean(AveExpr),
            stdev_AveExpr = sd(AveExpr),
            mean_t = mean(t),
            stdev_t = sd(t),
            mean_P.Value = mean(P.Value),
            stdev_P.Value = sd(P.Value),
            mean_adj.P.Value = mean(adj.P.Val),
            stdev_adj.P.Value = sd(adj.P.Val),
            mean_B = mean(B),
            stdev_B = sd(B))

write.table(summary_stats, file = paste0(common_filepath, version, "_mgs32_summarystats.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE)

# SUMMARY STATS FOR MGS 4 VS MGS 3
files <- list.files(path = summarystats_filepath, pattern = "_mgs43.txt")

data <- do.call("rbind", lapply(files, function(x) read.table(paste0(summarystats_filepath, x), header = TRUE, stringsAsFactors = FALSE)))

summary_stats <- data %>%
  group_by(ensembl_gene_id) %>%
  summarize(
            external_gene_name = unique(external_gene_name),
            strand = unique(strand),
            chromosome_name = unique(chromosome_name),
            start_position = unique(start_position),
            end_position = unique(end_position),
            gene_biotype = unique(gene_biotype),
            mean_logFC = mean(logFC),
            stdev_logFC = sd(logFC),
            mean_FC = mean(FC),
            stdev_FC = sd(FC),
            mean_AveExpr = mean(AveExpr),
            stdev_AveExpr = sd(AveExpr),
            mean_t = mean(t),
            stdev_t = sd(t),
            mean_P.Value = mean(P.Value),
            stdev_P.Value = sd(P.Value),
            mean_adj.P.Value = mean(adj.P.Val),
            stdev_adj.P.Value = sd(adj.P.Val),
            mean_B = mean(B),
            stdev_B = sd(B))

write.table(summary_stats, file = paste0(common_filepath, version, "_mgs43_summarystats.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE)
