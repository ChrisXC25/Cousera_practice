con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)

edata_log2 = log2(edata + 1)
svd = svd(edata_log2)
plot(svd$d^2/sum(svd$d^2))

edata_centered = edata - rowMeans(edata)
svd1= svd(edata_centered)
plot(svd1$d^2/sum(svd1$d^2))


edata_log2_centered = edata_log2 - rowMeans(edata_log2)

set.seed(333)
kmeans_result <- kmeans(t(edata_log2_centered), centers = 2)
cluster_assignments <- kmeans_result$cluster
svd_result <- svd(edata_log2_centered)
first_singular_vector <- svd_result$v[, 1]

# Convert cluster assignments to a numeric vector
cluster_indicator <- as.numeric(cluster_assignments)

# Calculate the correlation
correlation <- cor(first_singular_vector, cluster_indicator)
print(correlation)


#####################################################################


con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
edata = exprs(bm)
pdata_bm=pData(bm)

# Extract the counts for the first gene
first_gene_counts <- edata[1, ]

# Assume pdata contains a column 'num.tech.reps' for the number of technical replicates
# Here we are assuming pdata has this column, if not, it needs to be created or modified accordingly.
num_tech_reps <- as.factor(pdata$num.tech.reps)

# Fit the linear model
lm_model <- lm(first_gene_counts ~ num_tech_reps)

# Summarize the linear model
summary(lm_model)

# Plot the data
plot(num_tech_reps, first_gene_counts, xlab = "Number of Technical Replicates", ylab = "First Gene Counts", main = "First Gene Counts vs. Number of Technical Replicates")

# Add the linear model fit
abline(lm_model, col = "red")

View(pdata_bm)

lm1 <- lm(edata[1,]~ pdata_bm$age + pdata_bm$gender)
summary(lm1)

#####################################################

population = as.factor(pdata$population)

matrix1 <- model.matrix(~population)

fit = lm.fit(matrix1, t(edata_log2))

names(fit)

coeff <- fit$coefficients
dim(coeff)

res <- fit$residuals
dim(res)

eff <- fit$effects
dim(eff)

####################################################################

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
edata = exprs(bm)
pdata_bm=pData(bm)

# Remove samples with missing values for age
complete_samples <- !is.na(pdata$age)
edata <- edata[, complete_samples]
pdata <- pdata[complete_samples, ]

BiocManager::install("limma")
options(repos = c(CRAN = "http://cran.r-project.org"))
install.packages("statmod")

install.packages("statmod")
library(limma)
# Create the design matrix with age as the outcome
design_matrix <- model.matrix(~ age, data = pdata_bm)

# Fit the linear model
fit_1 <- lmFit(edata, t(design_matrix)
fit_1
# Apply empirical Bayes smoothing to the standard errors
fit1 <- eBayes(fit1)

# Extract the coefficients
coefficients <- fit1$coefficients

coef_age_1000th_gene <- coefficients[1000, "age"]
print(coef_age_1000th_gene)



