#Question 1 Fit a linear model and a logistic regression model 
#to the data for the 3rd SNP. What are the coefficients for the SNP variable?
#How are they interpreted?

BiocManager::install("snpStats")

library(snpStats)
library(broom)
data(for.exercise)
use <- seq(1, ncol(snps.10), 10)
sub.10 <- snps.10[,use]
snpdata = sub.10@.Data
status = subject.support$cc
View(snpdata)

snp_3rd = as.numeric(snpdata[, 3])

# Recode 0 values to NA
snp_3rd[snp_3rd == 0] <- NA

# Fit the linear model
linear_model = lm(status ~ snp_3rd, na.action = na.omit)
linear_coefficients = tidy(linear_model)
print(linear_coefficients)

# Fit the logistic regression model
logistic_model = glm(status ~ snp_3rd, family = binomial, na.action = na.omit)
logistic_coefficients = tidy(logistic_model)
print(logistic_coefficients)


#Question2 Fit a logistic regression model on a recessive (need 2 copies of minor allele to confer risk)
#and additive scale for the 10th SNP. Make a table of the fitted values versus the case/control status.
#Does one model fit better than the other?


# Install and load necessary packages
BiocManager::install("snpStats")
library(snpStats)
library(broom)

# Load data
data(for.exercise)
use <- seq(1, ncol(snps.10), 10)
sub.10 <- snps.10[, use]
snpdata = sub.10@.Data
status = subject.support$cc

# Extract the 10th SNP and convert to numeric
snp_10th = as.numeric(snpdata[, 10])

# Recode 0 values to NA
snp_10th[snp_10th == 0] <- NA

# Filter out rows with NA values in snp_10th
valid_indices = !is.na(snp_10th)
snp_10th = snp_10th[valid_indices]
status = status[valid_indices]

# Create recessive variable (1 if 2 copies of minor allele, 0 otherwise)
recessive = ifelse(snp_10th == 2, 1, 0)

# Fit logistic regression models
logistic_model_additive = glm(status ~ snp_10th, family = binomial)
logistic_model_recessive = glm(status ~ recessive, family = binomial)

# Get fitted values
fitted_additive = fitted(logistic_model_additive)
fitted_recessive = fitted(logistic_model_recessive)

# Create a data frame of the fitted values versus the case/control status
results = data.frame(
  status = status,
  fitted_additive = fitted_additive,
  fitted_recessive = fitted_recessive
)

# Create a table of the fitted values versus case/control status
table_additive = table(cut(results$fitted_additive, breaks = c(0, 0.5, 1)), results$status)
table_recessive = table(cut(results$fitted_recessive, breaks = c(0, 0.5, 1)), results$status)

# Print the tables
print("Table of fitted values (Additive) vs. case/control status:")
print(table_additive)

print("Table of fitted values (Recessive) vs. case/control status:")
print(table_recessive)

#Question 3 Fitted all snp and calculate the
#maximum and the minimum effect size 

# Initialize vectors to store effect sizes
effect_sizes <- numeric(ncol(snpdata))

# Fit additive logistic regression model to each SNP and extract effect sizes

# Fit additive logistic regression model to each SNP and extract effect sizes
for (i in 1:ncol(snpdata)) {
  snp = as.numeric(snpdata[, i])
  snp[snp == 0] <- NA
  valid_indices = !is.na(snp)
  snp = snp[valid_indices]
  status_valid = status[valid_indices]
  
  cat("Processing SNP", i, "with", length(snp), "valid entries\n")
  
  if (length(unique(snp)) > 1) {
    tryCatch({
      model = glm(status_valid ~ snp, family = binomial)
      coefficients = tidy(model)
      # Store the effect size (coefficient for the SNP)
      effect_sizes[i] <- coefficients$estimate[2]
      cat("SNP", i, "effect size:", coefficients$estimate[2], "\n")
    }, error = function(e) {
      cat("Error processing SNP", i, ":", e$message, "\n")
      effect_sizes[i] <- NA
    })
  } else {
    cat("SNP", i, "has no variation and is skipped\n")
    effect_sizes[i] <- NA
  }
}

View(snpdata)
# Remove NA values from effect sizes
effect_sizes <- effect_sizes[!is.na(effect_sizes)]

# Calculate summary statistics
average_effect_size = mean(effect_sizes)
max_effect_size = max(effect_sizes)
min_effect_size = min(effect_sizes)

cat("Average effect size:", average_effect_size, "\n")
cat("Max effect size:", max_effect_size, "\n")
cat("Min effect size:", min_effect_size, "\n")

#########################################################################

# Install and load necessary packages
BiocManager::install("snpStats")
library(snpStats)
library(broom)

# Load data
data(for.exercise)
use <- seq(1, ncol(snps.10), 10)
sub.10 <- snps.10[, use]
snpdata = sub.10@.Data
status = subject.support$cc

# Initialize vectors to store squared coefficients
squared_coefficients <- numeric(ncol(snpdata))

# Fit additive logistic regression model to each SNP and extract squared coefficients
for (i in 1:ncol(snpdata)) {
  snp = as.numeric(snpdata[, i])
  snp[snp == 0] <- NA
  valid_indices = !is.na(snp)
  snp = snp[valid_indices]
  status_valid = status[valid_indices]
  
  if (length(unique(snp)) > 1) {
    tryCatch({
      model = glm(status_valid ~ snp, family = binomial)
      coefficients = tidy(model)
      # Store the squared effect size (squared coefficient for the SNP)
      squared_coefficients[i] <- coefficients$estimate[2]^2
    }, error = function(e) {
      squared_coefficients[i] <- NA
    })
  } else {
    squared_coefficients[i] <- NA
  }
}

# Remove NA values from squared coefficients
squared_coefficients <- squared_coefficients[!is.na(squared_coefficients)]

# Perform snp.rhs.tests to get chi-squared statistics
chi_squared_results <- snp.rhs.tests(status ~ 1, snps.10)

# Extract chi-squared statistics
chi_squared_stats <- chi_squared_results$Chi.squared[use]

# Ensure that lengths match
min_length <- min(length(squared_coefficients), length(chi_squared_stats))
squared_coefficients <- squared_coefficients[1:min_length]
chi_squared_stats <- chi_squared_stats[1:min_length]

# Calculate the correlation
correlation <- cor(squared_coefficients, chi_squared_stats)

# Print the correlation
cat("Correlation between squared coefficients and chi-squared statistics:", correlation, "\n")




########################################################################################################


#Question 5:Do the log2(data + 1) transform and fit calculate F-statistics for the difference between studies/populations using genefilter:rowFtests 
#and using genefilter:rowttests. Do you get the same statistic? Do you get the same p-value?

BiocManager::install("genefilter")
library(genefilter)

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)

#log_2 transfrom
edata_transformed <- log2(edata + 1)

# Check the structure of pdata 
str(pdata)

# Assuming the variable representing studies/populations is 'study'
study <- pdata$study

# Calculate F-statistics using genefilter::rowFtests
f_test_results <- rowFtests(edata_transformed, as.factor(study))

# Calculate t-statistics using genefilter::rowttests
t_test_results <- rowttests(edata_transformed, factor(study))

# Extract statistics and p-values
f_stats <- f_test_results$statistic
t_stats <- t_test_results$statistic
f_pvalues <- f_test_results$p.value
t_pvalues <- t_test_results$p.value

#####################################################################################

library(DESeq2)
library(limma)

# Assuming the variable representing studies/populations is 'study'
study <- factor(pdata$study)

# 1. Use DESeq2 to test for differences between studies
dds <- DESeqDataSetFromMatrix(countData = edata, colData = pdata, design = ~ study)
dds <- DESeq(dds)
res_deseq2 <- results(dds)

# Extract DESeq2 statistics and p-values
deseq2_stat <- res_deseq2$stat
deseq2_pvalue <- res_deseq2$pvalue

# 2. Transform the data and use limma
edata_transformed <- log2(edata + 1)

# Create a design matrix for limma
design <- model.matrix(~ study)

# Fit the linear model using limma
fit <- lmFit(edata_transformed, design)
fit <- eBayes(fit)
limma_results <- topTable(fit, coef = "study", number = Inf, sort.by = "none")

# Extract limma statistics and p-values
limma_stat <- limma_results$t
limma_pvalue <- limma_results$P.Value

# 3. Calculate the correlation between the statistics
correlation_stats <- cor(deseq2_stat, limma_stat, use = "complete.obs")
cat("Correlation between DESeq2 and limma statistics:", correlation_stats, "\n")

# 4. Create an MA-plot
plotMA <- function(stat, pvalue, method) {
  plot(log2(rowMeans(edata_transformed)), stat, col = ifelse(pvalue < 0.05, "red", "black"), main = paste("MA-plot for", method))
}

par(mfrow = c(1, 2))
plotMA(deseq2_stat, deseq2_pvalue, "DESeq2")
plotMA(limma_stat, limma_pvalue, "limma")

# Compare the differences for large and small statistics
differences <- abs(deseq2_stat - limma_stat)
correlation_large_stats <- cor(deseq2_stat[deseq2_stat > 2], limma_stat[deseq2_stat > 2], use = "complete.obs")
correlation_small_stats <- cor(deseq2_stat[deseq2_stat <= 2], limma_stat[deseq2_stat <= 2], use = "complete.obs")

cat("Correlation for large statistics (|stat| > 2):", correlation_large_stats, "\n")
cat("Correlation for small statistics (|stat| <= 2):", correlation_small_stats, "\n")