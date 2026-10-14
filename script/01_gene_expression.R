# ============================================================
# Cancer Gene Expression Analysis and Machine Learning
# ============================================================

# Project: GeneExpression-ML-R
# Author: Raania Andleeb
# Language: R

# ============================================================
# 1. PROJECT GOAL
# ============================================================

# Goal:
# To analyze gene expression data from different cancer types
# and investigate whether gene expression patterns can be used
# to distinguish between cancer types using statistical analysis,
# dimensionality reduction, clustering, and machine learning.


# ============================================================
# 2. WHAT DO WE WANT TO ACHIEVE?
# ============================================================

# We want to:
# - Understand the structure of the biological dataset.
# - Explore differences in gene expression between cancer types.
# - Identify important patterns and variation in the data.
# - Use PCA to visualize high-dimensional gene expression data.
# - Investigate whether cancer types naturally form distinct groups.
# - Use clustering to explore similarities between samples.
# - Identify useful gene-expression features for prediction.
# - Build machine learning models to classify cancer types.
# - Compare model performance using appropriate evaluation metrics.
# - Present the findings clearly using graphs, tables, and
#   simple interpretations.


# ============================================================
# 3. WHY ARE WE DOING THIS ANALYSIS?
# ============================================================

# Gene expression data contains information about the activity
# of thousands of genes in biological samples. Because the data
# has many variables, statistical and machine learning methods
# can help us identify patterns that may not be visible by
# looking at individual genes.


# ============================================================
# 4. QUESTIONS WE WANT TO ANSWER
# ============================================================

# 1. How is the dataset structured?
# 2. How many samples belong to each cancer type?
# 3. How variable are the gene-expression measurements?
# 4. Can we visualize differences between cancer types?
# 5. Can PCA reveal separation between cancer types?
# 6. Can clustering identify groups of similar samples?
# 7. Which features are useful for cancer-type classification?
# 8. How accurately can machine learning classify cancer types?
# 9. Which model performs best?
# 10. What can we learn from the final results?


# ============================================================
# 5. ANALYSIS PLAN
# ============================================================

# Step 1: Load and understand the dataset
# Step 2: Data quality and structure checks
# Step 3: Exploratory data analysis
# Step 4: Explore gene-expression variability
# Step 5: Principal Component Analysis (PCA)
# Step 6: Clustering
# Step 7: Feature selection / dimensionality reduction
# Step 8: Machine learning classification
# Step 9: Model evaluation and comparison
# Step 10: Final predictions and conclusions


# ============================================================
# 6. PROJECT OUTPUTS
# ============================================================

# Important graphs will be saved in the "plots" folder.
# Important results/tables will be saved in the "outputs" folder.
# The README will document the project and display important
# graphs and findings for GitHub visitors.


# ============================================================
# 7. TOOLS
# ============================================================

# R
# RStudio
# tidyverse
# ggplot2
# statistical analysis
# PCA
# clustering
# machine learning

#---------------------------------------------------
# Check dataset contents
# --------------------------------------------------
list.files(path = "~/Downloads",full.names = TRUE)
list.files(path = "C:/Users/AR/Downloads",full.names = TRUE)
list.files("C:/Users/AR/Downloads/gene+expression+cancer+rna+seq",recursive = TRUE,full.names = TRUE)
# Why: Identify the files inside the compressed dataset before extracting it.
untar("C:/Users/AR/Downloads/gene+expression+cancer+rna+seq/TCGA-PANCAN-HiSeq-801x20531.tar.gz",list = TRUE)
untar("C:/Users/AR/Downloads/gene+expression+cancer+rna+seq/TCGA-PANCAN-HiSeq-801x20531.tar.gz",exdir = "data")
list.files("data",recursive = TRUE,full.names = TRUE)

# Load gene expression data
# Why: Load the gene-expression measurements so we can understand the dataset.
gene_data <- read.csv("data/TCGA-PANCAN-HiSeq-801x20531/data.csv",check.names = FALSE)

# Load cancer labels
# Why: Load the cancer-type labels so we can connect each sample with its cancer type.
labels <- read.csv("data/TCGA-PANCAN-HiSeq-801x20531/labels.csv",check.names = FALSE)
dim(gene_data)
dim(labels)
head(gene_data[, 1:6])
head(labels)

# ----------------------------Interpretation-----------------------------------
# The dataset contains 801 cancer samples and 20,531 gene-expression features.
# Each sample has a corresponding cancer-type label, which we will use for
# comparing cancer groups and building classification models.

# Cancer type distribution
# Why: Check how many samples belong to each cancer type before comparing groups.
table(labels$Class)
prop.table(table(labels$Class)) * 100

#------------------------------ Interpretation------------------------------
# BRCA represents the largest group (37.45%), while COAD represents the
# smallest group (9.74%). The class distribution is therefore imbalanced,
# which we should consider when evaluating classification performance.

# Cancer type distribution plot
# Why: Visualize the number of samples in each cancer type and identify class imbalance.
library(ggplot2)
class_plot <- ggplot(labels, aes(x = Class)) +geom_bar() +labs(title = "Number of Samples by Cancer Type",x = "Cancer Type",y = "Number of Samples") +theme_minimal()
class_plot
ggsave("plots/cancer_type_distribution.png",plot = class_plot,width = 8,height = 5,dpi = 300)

#----------------------- Interpretation--------------------------------------
# BRCA has the highest number of samples, while COAD has the lowest.
# The differences in sample counts may affect how well our models
# classify each cancer type, so we will check class-wise performance later.

# Check missing values
# Why: Identify missing data before performing statistical analysis and machine learning.
sum(is.na(gene_data))
sum(is.na(labels))

#--------------------- Interpretation----------------------------------------
# No missing values were found in the gene-expression data or cancer labels.
# This means we do not need to handle missing values before continuing.

# Check duplicate sample IDs
# Why: Identify repeated sample identifiers before analyzing gene expression.
sum(duplicated(rownames(gene_data)))
sum(duplicated(rownames(labels)))

#------------------ Interpretation--------------------------
# No duplicated row names were found in either dataset.

# Check sample IDs
# Why: Confirm that each sample ID is unique and matches between both datasets.
sum(duplicated(gene_data[[1]]))
sum(duplicated(labels[[1]]))
identical(gene_data[[1]], labels[[1]])

#------------------------- Interpretation---------------------------------------
# All sample IDs are unique in both datasets, and their order matches exactly.
# This confirms that each gene-expression profile corresponds to the correct
# cancer label, so we can safely combine the data for further analysis.

# Prepare gene-expression features
# Why: Separate sample identifiers from gene measurements so IDs are not
# mistakenly treated as biological features in statistical models.
sample_ids <- gene_data[[1]]
expression_data <- gene_data[, -1, drop = FALSE]
# Check dimensions and data types
dim(expression_data)
str(expression_data[, 1:6])

# ----------------------Interpretation------------------------------------
# After removing sample IDs, the dataset contains 801 samples and 20,531
# numeric gene-expression features. The data is ready for further checks
# of gene variability and biological patterns.

# Check genes with zero variance
# Why: Identify genes with constant expression that cannot distinguish samples.
gene_variances <- sapply(expression_data, var)
sum(gene_variances == 0)
summary(gene_variances)

#------------------------ Interpretation---------------
# After removing sample IDs, the dataset contains 801 samples and 20,531
# numeric gene-expression features. The data is ready for further checks
# of gene variability and biological patterns.

# Check genes with zero variance
# Why: Identify genes with constant expression that cannot distinguish samples.
gene_variances <- sapply(expression_data, var)
sum(gene_variances == 0)
summary(gene_variances)

#--------------------- Interpretation------------------------------------
# A total of 267 genes have zero variance across all samples.
# These genes cannot help distinguish cancer types, so we will remove
# them before dimensionality reduction and machine learning.

# Remove zero-variance genes
# Why: Keep only genes whose expression varies across samples.
expression_data_filtered <- expression_data[, gene_variances > 0,drop = FALSE]
dim(expression_data_filtered)
sum(sapply(expression_data_filtered, var) == 0)

#----------------------------- Interpretation--------------------------------
# We removed 267 genes with constant expression, leaving 20,264 variable
# genes across 801 samples. This reduces unnecessary features before
# dimensionality reduction and machine learning.

# Summarize gene-expression values
# Why: Understand the overall range and distribution of expression values.
summary(as.vector(as.matrix(expression_data_filtered)))

#------------------------- Interpretation----------------------------------
# Expression values range from 0 to 20.779, with a median of 7.683.
# The mean is lower than the median, suggesting a possible skew toward
# lower values. We will visualize the distribution to examine this further.

# Gene-expression distribution
# Why: Visualize the overall distribution of expression values.
expression_values <- as.vector(as.matrix(expression_data_filtered))
expression_plot <- ggplot(data.frame(expression = expression_values),aes(x = expression)) +geom_histogram(bins = 50) +labs(title = "Distribution of Gene Expression Values",x = "Gene Expression Value",y = "Frequency") +theme_minimal()
expression_plot
ggsave("plots/gene_expression_distribution.png",plot = expression_plot,width = 8,height = 5,dpi = 300)

#------------------------ Interpretation-------------------------------
# The distribution has a large spike at zero and another concentration
# around 8–10. This shows that many gene measurements are zero or low,
# while many others have moderate expression values. The distribution
# is not normally shaped, so we should consider scaling before PCA.

# Standardize gene-expression features
# Why: Put gene features on a comparable scale before PCA.
expression_scaled <- scale(expression_data_filtered)
dim(expression_scaled)

#------------------------------ Interpretation----------------------------
# The standardized expression matrix retains all 801 samples and 20,264
# variable genes. Each gene is now centered and scaled, preparing the data
# for PCA so genes with larger variances do not dominate the analysis.

# Principal Component Analysis (PCA)
# Why: Reduce gene-expression dimensions to explore patterns among cancer types.
pca_result <- prcomp(expression_scaled,center = FALSE,scale. = FALSE)
# Check the proportion of variance explained
pca_variance <- pca_result$sdev^2 /sum(pca_result$sdev^2)
round(pca_variance[1:10] * 100, 2)

#------------------------ Interpretation-----------------------------------
# PC1 explains 10.54% of the variation and PC2 explains 8.75%.
# Together, they explain 19.29% of the total variation.
# The PCA plot will show an overview of the cancer samples, but it will
# not capture all the differences present in the gene-expression data.

# Prepare PCA plot data
# Why: Combine the first two principal components with cancer labels for plotting.
pca_plot_data <- data.frame(PC1 = pca_result$x[, 1],PC2 = pca_result$x[, 2],Class = labels$Class)
dim(pca_plot_data)

# Visualize PCA results
# Why: Check whether samples from different classes form separate clusters
pca_plot <- ggplot(pca_plot_data,aes(x = PC1, y = PC2, color = Class)) +geom_point(size = 2, alpha = 0.7) +labs(title = "PCA of Gene Expression Data",x = paste0("PC1 (", round(pca_variance[1] * 100, 2), "% variance)"),y = paste0("PC2 (", round(pca_variance[2] * 100, 2), "% variance)"),color = "Cancer Type") +theme_minimal()
pca_plot
ggsave("plots/pca_cancer_types.png",plot = pca_plot,width = 9,height = 6,dpi = 300)

#------------------------ Interpretation -------------------------------
# KIRC samples form a relatively distinct cluster on the right side
# of the plot, while the other cancer types show more overlap.
# BRCA, LUAD, PRAD, and COAD are less clearly separated using the
# first two principal components.
# PC1 explains 10.54% and PC2 explains 8.75% of the total variation,
# so the plot shows only part of the overall gene-expression structure.

# K-means clustering
# Why: Group samples with similar gene-expression patterns without using cancer labels.
set.seed(123)
kmeans_result <- kmeans(expression_scaled,centers = 5,nstart = 10,iter.max = 100)
table(kmeans_result$cluster)

#------------------------- Interpretation-----------------------------------
# The five clusters contain different numbers of samples, with Cluster 3
# being the largest (250 samples) and Cluster 1 the smallest (73 samples).
# However, cluster sizes alone do not tell us whether the groups correspond
# to the five cancer types. We need to compare clusters with the actual
# cancer labels to understand how well they align.

# Compare clusters with known cancer type
# Why: Check whether gene-expression clusters align with the actual cancer labels.
cluster_comparison <- table(Cluster = kmeans_result$cluster,Cancer_Type = labels$Class)
cluster_comparison
write.csv(as.data.frame.matrix(cluster_comparison),"outputs/kmeans_cluster_comparison.csv")

# -------------------------Interpretation---------------------------------
# Cluster 1 contains all COAD samples, Cluster 2 contains almost all
# KIRC samples, and Cluster 4 contains nearly all PRAD samples.
# Cluster 3 mainly contains BRCA samples, while Cluster 5 contains
# most LUAD samples along with some BRCA and other cancer samples.
# Overall, K-means separates several cancer types well, although BRCA
# and LUAD samples show some mixing across clusters.

# Visualize K-means clusters using PCA
# Why: See how the clusters formed by K-means are distributed in two dimensions.
pca_cluster_data <- data.frame(PC1 = pca_result$x[, 1],PC2 = pca_result$x[, 2],Cluster = factor(kmeans_result$cluster))
pca_cluster_plot <- ggplot(pca_cluster_data,aes(x = PC1, y = PC2, color = Cluster)) +geom_point(size = 2, alpha = 0.7) +labs(title = "K-means Clusters of Gene Expression Data",x = paste0("PC1 (", round(pca_variance[1] * 100, 2), "% variance)"),y = paste0("PC2 (", round(pca_variance[2] * 100, 2), "% variance)"),color = "Cluster") +theme_minimal()
pca_cluster_plot
ggsave("plots/kmeans_clusters_pca.png",plot = pca_cluster_plot,width = 9,height = 6,dpi = 300)

#--------------------------Interpretation------------------------------------
# Cluster 2 forms a relatively distinct group on the right side of the plot,
# while Cluster 3 is mainly concentrated in the lower region.
# Clusters 1, 4, and 5 show more overlap near the center of the plot.
# This suggests that some clusters are clearly separated in the first
# two principal components, while others have more similar expression patterns.

# Evaluate clustering using Adjusted Rand Index (ARI)
# Why: Measure the agreement between K-means clusters and known cancer labels.
install.packages("mclust")
library(mclust)
ari_score <- adjustedRandIndex(kmeans_result$cluster,labels$Class)
ari_score

#-------------------------- Interpretation----------------------------------------
# The Adjusted Rand Index (ARI) is 0.795, indicating strong agreement
# between the K-means clusters and the actual cancer types.
# This suggests that the gene-expression patterns capture much of the
# structure associated with cancer type, although the clustering is not perfect.

# Split data into training and testing sets
# Why: Train the model on one subset and evaluate its performance on unseen samples.
set.seed(123)
# Select 80% of samples from each cancer type
sample_index <- unlist(lapply(split(seq_len(nrow(labels)), labels$Class),function(i) sample(i, size = floor(0.80 * length(i)))),use.names = FALSE)
# Shuffle the selected training indices
sample_index <- sample(sample_index)
# Create training and testing sets
x_train <- expression_scaled[sample_index, , drop = FALSE]
x_test <- expression_scaled[-sample_index, , drop = FALSE]
y_train <- labels$Class[sample_index]
y_test <- labels$Class[-sample_index]
# Check dimensions
dim(x_train)
dim(x_test)
# Check cancer-type distribution
table(y_train)
table(y_test)

#---------------------- Interpretation-------------------------------
# The dataset is divided into 638 training samples and 163 testing
# samples. All five cancer types are represented in both sets, with
# approximately 80% used for training and 20% for testing.
# This allows us to evaluate model performance on unseen samples.

# Prepare training and testing data without data leakage
# Why: Calculate scaling parameters from training samples only to avoid using information from the test set.
set.seed(123)
# Calculate gene means and standard deviations using training samples
train_means <- colMeans(expression_data_filtered[sample_index, , drop = FALSE])
train_sds <- apply(expression_data_filtered[sample_index, , drop = FALSE], 2, sd)
# Remove genes with zero or invalid training standard deviation
valid_genes <- is.finite(train_sds) & train_sds > 0
train_means <- train_means[valid_genes]
train_sds <- train_sds[valid_genes]
# Scale training and testing data using training parameters
x_train_scaled <- scale(expression_data_filtered[sample_index, valid_genes, drop = FALSE],center = train_means,scale = train_sds)
x_test_scaled <- scale(expression_data_filtered[-sample_index, valid_genes, drop = FALSE],center = train_means,scale = train_sds)
# Check dimensions
dim(x_train_scaled)
dim(x_test_scaled)

# Check class distribution
# Why: Ensure all cancer types are represented in training and testing data.
table(y_train)
prop.table(table(y_train))
table(y_test)
prop.table(table(y_test))

# Train Random Forest model
# Why: Learn patterns in gene expression data to classify cancer types.
library(randomForest)
set.seed(123)
rf_model <- randomForest(x = x_train_scaled,y = as.factor(y_train),ntree = 300,importance = TRUE)
print(rf_model)

#-----------------------------Interpretation---------------------------
#The Random Forest model achieved 99.53% OOB accuracy, with only 3 misclassifications 
#among 638 training samples. This indicates strong separation between the
#five cancer types, but test-set evaluation is needed to confirm 
#performance on unseen data.

# Predict cancer types in the test set
# Why: Evaluate how well the model classifies unseen samples.
rf_predictions <- predict(rf_model,newdata = x_test_scaled)
# Compare predicted and actual cancer types
table(Actual = y_test,Predicted = rf_predictions)


#---------------------- Interpretation--------------------------------
# The model correctly classified all COAD, KIRC, and PRAD samples,
# and most BRCA and LUAD samples. Only one LUAD sample was
# misclassified as BRCA, showing strong classification performance.

# Calculate test accuracy
# Why: Measure the proportion of test samples correctly classified by the model.
accuracy <- mean(rf_predictions == y_test)
accuracy

# -----------------Interpretation---------------------------
# The Random Forest model achieved 99.39% test accuracy,
# correctly classifying 162 out of 163 unseen samples.

# Evaluate classification performance
# Why: Calculate precision, recall, and F1-score for each cancer type.
install.packages("yardstick")
library(yardstick)
library(dplyr)
results <- data.frame(truth = factor(y_test),estimate = factor(rf_predictions, levels = levels(factor(y_test))))
# Calculate metrics for each cancer type
results %>%precision(truth, estimate, estimator = "macro")
results %>%recall(truth, estimate, estimator = "macro")
results %>%f_meas(truth, estimate, estimator = "macro")

#------------------- Interpretation-----------------------------
# The model achieved 99.7% macro precision, 99.3% macro recall,
# and a 99.5% macro F1-score, showing strong classification
# performance across the five cancer types.

# Train an SVM model
# Why: Compare a different classification algorithm with Random Forest.

install.packages("e1071")  # Run once if not installed
library(e1071)
svm_model <- svm(x = x_train_scaled,y = as.factor(y_train),kernel = "radial")
# Predict cancer types in the test set
svm_predictions <- predict(svm_model, x_test_scaled)
# Calculate test accuracy
svm_accuracy <- mean(svm_predictions == y_test)
svm_accuracy

# ----------------Interpretation------------------------------
# The SVM model achieved 98.77% test accuracy, slightly lower
# than Random Forest's 99.39% accuracy on the same test set.

#Compare model performance
# Why: Identify which model achieves higher test accuracy.
model_comparison <- data.frame(Model = c("Random Forest", "SVM"),Test_Accuracy = c(accuracy, svm_accuracy))
model_comparison$Test_Accuracy <- round(model_comparison$Test_Accuracy * 100, 2)
print(model_comparison)

#----------------- Interpretation-------------------------------------
# Random Forest achieved 99.39% test accuracy, outperforming
# SVM (98.77%) by 0.62 percentage points.
# Both models showed excellent cancer classification performance.

# Visualize model performance
# Why: Compare test accuracy visually and save the result for GitHub.
model_plot <- ggplot(model_comparison,aes(x = Model, y = Test_Accuracy, fill = Model)) +geom_col(width = 0.6) +geom_text(aes(label = paste0(Test_Accuracy, "%")),vjust = -0.5) +labs(title = "Model Comparison: Test Accuracy",x = "Machine Learning Model",y = "Test Accuracy (%)") +coord_cartesian(ylim = c(95, 100)) +theme_minimal() +theme(legend.position = "none")
print(model_plot)
ggsave("Plots/model_comparison.png",plot = model_plot,width = 8,height = 6,dpi = 300)
 
# ------------------------Interpretation---------------------------------
# The  plot shows that Random Forest achieved 99.39% test accuracy,
# while SVM achieved 98.77%. Random Forest performed slightly better
# than SVM in classifying the five cancer types.

# Examine Random Forest feature importance
# Why: Identify genes that contribute most to cancer type classification.
gene_importance <- importance(rf_model)
head(gene_importance[order(gene_importance[, "MeanDecreaseGini"],decreasing = TRUE),],10)

#-------------------- Interpretation---------------------------
# Gene_17801 had the highest Mean Decrease Gini importance,
# followed by gene_15589 and gene_5578. These genes contributed
# strongly to the Random Forest model's classification decisions,
# but their biological significance requires further investigation.


# Visualize top 10 important genes
# Why: Show which genes contribute most to Random Forest classification.
top_genes <- data.frame(Gene = rownames(gene_importance),Importance = gene_importance[, "MeanDecreaseGini"])
top_genes <- top_genes[order(top_genes$Importance, decreasing = TRUE),]
top_genes <- head(top_genes, 10)
top_genes$Gene <- factor(top_genes$Gene,levels = rev(top_genes$Gene))
importance_plot <- ggplot(top_genes,aes(x = Gene, y = Importance)) +geom_col() +coord_flip() +labs(title = "Top 10 Important Genes",x = "Gene",y = "Mean Decrease Gini") +theme_minimal()
print(importance_plot)
ggsave("plots/top_10_gene_importance.png",plot = importance_plot,width = 9,height = 6,dpi = 300)

#----------------------- Interpretation--------------------------------
# gene_17801 had the highest Mean Decrease Gini importance (~4.53),
# while the other top genes had lower, relatively similar scores.
# It is a candidate for further investigation, but its biological
# significance must be validated before calling it a biomarker.

# Identify top genes for each cancer type
# Why: Find genes that contribute most to distinguishing each cancer class.
class_importance <- importance(rf_model)
cancer_types <- c("BRCA", "COAD", "KIRC", "LUAD", "PRAD")
for (cancer in cancer_types) {cat("\nTop 10 genes for", cancer, "\n")
class_scores <- class_importance[order(class_importance[, cancer], decreasing = TRUE),]
print(head(class_scores[, cancer, drop = FALSE], 10))}

#------------------------ Interpretation---------------
# Different cancer types have different top-ranked genes.
# gene_5578 ranked highest for BRCA, gene_9175 for COAD,
# gene_5729 for KIRC, gene_11903 for LUAD, and gene_12245
# for PRAD. These genes may help distinguish cancer types,
# but require further biological validation.

# Create gene importance heatmap
# Why: Compare how selected genes contribute to distinguishing cancer types.
# Select the top 10 genes overall
top_gene_names <- names(sort(gene_importance[, "MeanDecreaseGini"], decreasing = TRUE))[1:10]
# Extract class-specific importance scores
heatmap_data <- as.data.frame(class_importance[top_gene_names, cancer_types, drop = FALSE])
heatmap_data$Gene <- rownames(heatmap_data)
# Convert to long format
heatmap_data <- reshape(heatmap_data,varying = cancer_types,v.names = "Importance",timevar = "Cancer_Type",times = cancer_types,direction = "long")
gene_heatmap <- ggplot(heatmap_data,aes(x = Cancer_Type, y = Gene, fill = Importance)) +geom_tile(color = "white") +labs(title = "Top 10 Genes: Class-Specific Importance",x = "Cancer Type",y = "Gene",fill = "Importance") +theme_minimal()
print(gene_heatmap)
ggsave("plots/gene_importance_heatmap.png",plot = gene_heatmap,width = 8,height = 6,dpi = 300)

#----------------------Interpretation----------------------------
  # gene_5578 and gene_17801 show high importance for BRCA,
  # while gene_5729 is important for KIRC. gene_17801 also
  # contributes to LUAD classification. Importance varies across
  # cancer types, suggesting class-specific predictive patterns.

# Examine gene importance using Mean Decrease Accuracy
# Why: Identify genes whose disruption most reduces model accuracy.
accuracy_importance <- importance(rf_model)
top_accuracy_genes <- head(accuracy_importance[order(accuracy_importance[, "MeanDecreaseAccuracy"],decreasing = TRUE),"MeanDecreaseAccuracy",drop = FALSE],10)
top_accuracy_genes

#------------------- Interpretation-----------------------------------
# gene_11903 had the highest Mean Decrease Accuracy, followed by
# gene_5578 and gene_15894. These genes are important for model
# predictions, supporting further investigation of their roles
# in cancer type classification.


# Compare gene importance rankings
# Why: Identify genes that rank highly under both importance measures.
gini_ranking <- sort(gene_importance[, "MeanDecreaseGini"],decreasing = TRUE)
accuracy_ranking <- sort(gene_importance[, "MeanDecreaseAccuracy"],decreasing = TRUE)
top_gini_genes <- names(gini_ranking)[1:10]
top_accuracy_genes <- names(accuracy_ranking)[1:10]
common_genes <- intersect(top_gini_genes, top_accuracy_genes)
common_genes

#-------------------- Interpretation------------------------------
# gene_17801, gene_5578, and gene_5729 ranked in the top 10
# for both importance measures, making them useful candidates
# for further investigation in cancer type classification.

# Why: Compare both models on the same unseen test samples.
model_comparison <- data.frame(Model = c("Random Forest", "SVM"),Test_Accuracy = c(accuracy, svm_accuracy))
model_comparison$Test_Accuracy <- round(model_comparison$Test_Accuracy * 100, 2)
model_comparison


