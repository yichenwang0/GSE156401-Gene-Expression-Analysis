library(GEOquery)
geo_data <- getGEO("GSE156401", GSEMatrix = TRUE)
df <- exprs(geo_data[[1]])
dim(df)
head(df)
colnames(df)
unique_labels <- c("shControl.1", "shControl.2", "shControl.3", "shMEK5.1", "shMEK5.2", "shMEK5.3")
colnames(df) <- unique_labels

# Check for outliers
# Cluster tree
dat <- t(df)
dat_dist <- dist(dat, method = "euclidean")
dat_clust <- hclust(dat_dist, method = "single")
plot(dat_clust, labels = names(dat), cex = 0.8)

# CV vs. Mean Plot
df_mean <- apply(log2(df), 2, function(x) mean(x, na.rm = TRUE))
df_sd <- sqrt(apply(log2(df), 2, function(x) var(x, na.rm = TRUE)))
df_cv <- df_sd/df_mean
plot(df_mean, df_cv, main = "PC3 Human Prostate Cancer Cell Dataset - Sample CV vs. Mean (log2-transformed)", xlab = "Mean", ylab = "CV", col = "blue", cex = 1.5, type = "n")
points(df_mean, df_cv, bg = "lightblue", col = 1, pch = 21)
text(df_mean, df_cv, label = colnames(df), pos = 1, cex = 0.7)

# Average Correlation Plot
df_cor <- cor(df, use = "pairwise.complete.obs")
df_avg <- apply(df_cor, 1, mean)
plot(c(1, length(df_avg)), range(df_avg), type = "n", xlab = "", ylab = "Avg Correlation", main = "Average Correlation Plot of Control/shMEK5 Samples", axes = F)
points(df_avg, bg = "red", col = 1, pch = 21, cex = 1.25)
axis(1, at = c(1:length(df_avg)), labels = colnames(df), las= 2, cex.lab = 0.8, cex.axis = 1)
axis(2)

# Filter out genes with low expression levels
# Check and remove NAs
any(is.na(df))
any(is.infinite(df))
df_cleaned <- na.omit(df)
dim(df_cleaned)

# Histogram of CVs
gene_sd <- apply(df_cleaned, 1, sd)
gene_avg <- apply(df_cleaned, 1, mean)
gene_cv_percent <- (gene_sd/gene_avg)*100
hist(gene_cv_percent, n = 20, col = "cyan", main = "Distribution of CVs for Genes", xlab = "CV (%)", xlim = c(0, max(gene_cv_percent)))

# Remove genes with low CVs
sum(gene_cv_percent < 4)
retained_genes <- which(gene_cv_percent >= 4)
df_filtered <- df_cleaned[retained_genes, ]
dim(df_filtered)

# Feature Selection with empirical Bayes method
library(multtest)
library(limma)
groups <- factor(rep(c("Control", "shMEK5"), each = 3))
design <- model.matrix(~groups)
fit <- lmFit(df_filtered, design)
eb_fit <- eBayes(fit)
table_cp <- topTable(eb_fit, number=nrow(df_filtered))
table_cp <- table_cp[order(table_cp$P.Value), ]
ebayes_p <- table_cp$P.Value
bon_ebayes_p <- mt.rawp2adjp(ebayes_p, "Bonferroni")

# Non-Adjusted vs. Adjusted P-Values Plot
p_x <- 1:length(bon_ebayes_p$adjp[, 1])
plot(p_x, bon_ebayes_p$adjp[, 1], type = "l", col = "blue", xlab = "Number of Comparison", ylab = "Sorted P-Value", main = "Non-Adjusted vs. Adjusted P-Values")
lines(p_x, bon_ebayes_p$adjp[, 2], type = "l", col = "red")
legend("bottomright", legend = c("Non-Adjusted", "Adjusted (Bonferroni)"), col = c("blue", "red"), lty = 1)

# Extract only the genes with an adjusted p-value < 0.01
adjusted_p <- bon_ebayes_p$adjp[, 2]
names(adjusted_p) <- rownames(table_cp)
sig_gene_p <- adjusted_p[which(adjusted_p < 0.01)]
sig_gene_p
sig_gene <- names(sig_gene_p)
hist(sig_gene_p, n = 10, col = "cyan", main = "Distribution of P-Values for Genes (p < 0.01)", xlab = "P-Value", xlim = c(0, 0.01))

# Sample Visulization with PCA plot
df_subset <- df_filtered[sig_gene, ]
dat.pca <- prcomp(t(df_subset))
loading <- dat.pca$x[, 1:2]
plot(range(loading[, 1]), range(loading[, 2]), type = "n", xlab = "p1", ylab = "p2", main = "PCA Plot of PC3 Human Prostate Cancer Cell Data\np2 vs. p1")
points(loading[, 1][1:3], loading[, 2][1:3], col = 1, bg = "red", pch = 21, cex = 1.5)
points(loading[, 1][4:6], loading[, 2][4:6], col = 1, bg = "blue", pch = 21, cex = 1.5)
legend("bottomright", legend = c("Control", "shMEK5"), col = c("red", "blue"), pch = 15, cex = 1.5)

# Scree plot
dat.pca.var <- round(dat.pca$sdev^2 / sum(dat.pca$sdev^2)*100, 2)
plot(c(1:length(dat.pca.var)), dat.pca.var, type = "b", xlab = "# components", ylab = "% variance", pch = 21, col = 1, bg = 3, cex = 1.5)
title("Scree plot showing % variability explained by each eigenvalue\nControl/shMEK5 dataset")
dat.pca.var[1] + dat.pca.var[2]

# KMeans clustering
kmeans_2 <- kmeans(loading, centers = 2)
plot(loading, col = kmeans_2$cluster, cex = 1, pch = 16, main = "PCA Scatter Plot with KMeans Cluster Membership")
points(kmeans_2$centers, col = 1:2, pch = "*", cex = 2.5)
text(loading[, 1], loading[, 2], labels = colnames(df_subset), pos = 1)

# Classification with Linear Discriminant Analysis
library(MASS)
classes <- rep(c("Control", "shMEK5"), each = 3)
t_df <- t(df_subset)
dm <- data.frame(classes, t_df)
dat.lda <- lda(classes~., dm)
dat.pred <- predict(dat.lda, dm)
table(dat.pred$class, classes)
plot(dat.pred$x, col = as.factor(dat.pred$class), pch = 16, main = "First 2 Discriminant Functions for PC3 Human Prostate Cancer Cell Dataset", xlab = "Discriminant Function 1", ylab = "Discriminant Function 2")
legend("bottomright", legend = levels(as.factor(dat.pred$class)), col = 1:2, pch = 16, title = "Predicted Classes")
text(dat.pred$x, labels = classes, pos = 3)

# Top 5 Discriminant Genes
# Calculate log2 fold changes
ct_mean <- apply(log2(df_subset)[, 1:3], 1, mean)
MEK5_mean <- apply(log2(df_subset)[, 4:6], 1, mean)
fold <- MEK5_mean - ct_mean
neg <- names(fold[which(fold < 0)])
pos <- names(fold[which(fold > 0)])
top_5_pos <- names(sig_gene_p[pos][1:5])
top_5_neg <- names(sig_gene_p[neg][1:5])
probe_map <- getGEO("GPL13497", GSEMatrix = TRUE)
probe_table <- Table(probe_map)
colnames(probe_table)
probe_table[probe_table$ID %in% top_5_pos, c("ID", "GENE_SYMBOL", "GENE_NAME", "DESCRIPTION")]
probe_table[probe_table$ID %in% top_5_neg, c("ID", "GENE_SYMBOL", "GENE_NAME", "DESCRIPTION")]