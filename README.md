# GeneExpression-ML-R

## Project Overview

This project analyzes gene-expression data from 801 samples across five
cancer types using exploratory data analysis, PCA, K-means clustering,
Random Forest, and SVM classification.

## Objectives

-   Examine dataset structure, quality, and cancer-type distribution.
-   Explore gene-expression variability.
-   Visualize high-dimensional expression data using PCA.
-   Compare K-means clusters with known cancer labels.
-   Train and compare classifiers.
-   Explore genes contributing to Random Forest predictions.

## Dataset

The dataset contains 801 samples and 20,531 gene-expression features.

  Cancer type     Samples
  ------------- ---------
  BRCA                300
  COAD                 78
  KIRC                146
  LUAD                141
  PRAD                136
  **Total**       **801**

No missing values were found. Removing sample IDs from the feature
matrix and filtering out 267 zero-variance genes left 20,264 variable
genes.

## Workflow

1.  Load the expression data and cancer labels.
2.  Check dimensions, missing values, sample IDs, and class
    distribution.
3.  Explore expression values and plot their distribution.
4.  Remove zero-variance genes.
5.  Use PCA to visualize sample patterns.
6.  Apply K-means clustering with five clusters.
7.  Compare clusters with known labels using the Adjusted Rand Index.
8.  Train Random Forest and radial-kernel SVM classifiers.
9.  Compare test accuracy and inspect Random Forest feature importance.

## Key Findings

### PCA

PC1 explains **10.54%** of variation and PC2 explains **8.75%**, for
**19.29%** combined. KIRC samples form a relatively distinct group in
the first two components, while the other cancer types show more
overlap.

### K-means clustering

The Adjusted Rand Index was **0.795**, indicating strong agreement
between cluster assignments and known cancer types. Several clusters
aligned closely with COAD, KIRC, and PRAD, while BRCA and LUAD showed
more mixing.

### Classification results

The script reports the following results on the same test split:

  Model                 Test accuracy
  ------------------- ---------------
  Random Forest            **99.39%**
  Radial-kernel SVM        **98.77%**

Random Forest outperformed SVM by **0.62 percentage points**. The Random
Forest confusion table showed 162 correctly classified samples out of
163, with one LUAD sample misclassified as BRCA.

The script reports these macro-averaged Random Forest metrics:

  Metric        Result
  ----------- --------
  Precision      99.7%
  Recall         99.3%
  F1-score       99.5%

These are the results recorded in the supplied script and have not been
independently recomputed for this README.

### Gene importance

Random Forest importance was examined using Mean Decrease Gini and Mean
Decrease Accuracy. `gene_17801`, `gene_5578`, and `gene_5729` appear
among the genes ranking highly under both measures. These are
exploratory model-derived features, not confirmed biomarkers; biological
significance requires independent validation.

## Visualizations

The script saves the following plots:

-   `plots/cancer_type_distribution.png`
-   `plots/gene_expression_distribution.png`
-   `plots/pca_cancer_types.png`
-   `plots/kmeans_clusters_pca.png`
-   `plots/top_10_gene_importance.png`
-   `plots/gene_importance_heatmap.png`
-   `Plots/model_comparison.png`

**Note:** The model-comparison plot uses `Plots/` with a capital P,
while the other plots use `plots/`. Use one consistent folder name,
especially for GitHub.

## Tools and Packages

-   R and RStudio
-   `ggplot2`
-   `mclust`
-   `randomForest`
-   `e1071`
-   `dplyr`
-   `yardstick`

## Conclusion

The analysis found meaningful structure in the gene-expression data.
K-means clustering agreed strongly with the known cancer labels (ARI ≈
0.795), and both supervised models achieved high test accuracy. Random
Forest performed slightly better than SVM on the reported split and
provided gene-importance rankings for further exploration.

These findings suggest that gene-expression patterns can help
distinguish the five cancer types in this dataset. The high accuracy
should be confirmed using repeated or external validation before drawing
conclusions about performance on new datasets. The gene-importance
results are exploratory and should not be interpreted as clinically
validated biomarkers.

## Reproducibility Notes

-   The script uses `set.seed(123)` for the train/test split and
    K-means.
-   The split contains 638 training samples and 163 test samples.
-   For supervised modeling, calculate scaling parameters on the
    training set only and apply the same parameters to the test set. The
    supplied script includes a corrected training-only scaling section;
    use `x_train_scaled` and `x_test_scaled` for model fitting and
    prediction.
-   Remove or make optional machine-specific file-listing and extraction
    commands before sharing the script.
-   Avoid running `install.packages()` every time the script runs.

## Suggested Repository Structure

``` text
GeneExpression-ML-R/
├── data/
│   └── TCGA-PANCAN-HiSeq-801x20531/
│       ├── data.csv
│       └── labels.csv
├── plots/
├── outputs/
├── script/
└── GeneExpression-ML-R.Rproj
```

This is the intended structure; include only folders and files that
exist in the repository. Check the dataset's original source and license
terms before redistributing raw data publicly.
