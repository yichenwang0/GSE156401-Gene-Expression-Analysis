# MEK5-Mediated Gene Expression Changes in PC3 Prostate Cancer Cells

This repository contains code and documentation for analyzing a publicly available dataset from the Gene Expression Omnibus (GEO), which can be accessed [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156401). The dataset includes gene expression measurements in PC3 cells expressing a scrambled shRNA (shControl) or MEK5 shRNA, with three replicates for each condition. 

## Project Overview  

The analysis pipeline includes the following steps:

1. **Data Preprocessing**: Checking for outlier samples.
2. **Gene Filtering**: Removing genes with low expression values.
3. **Feature Selection**: Utilizing the empirical Bayes method for identifying significant genes with multiplicity adjustment
4. **Sample Visualization**: Visualizing samples in two-dimensional space with clustering and dimensionality reduction methods 
5. **Classification**: Predicting sample classes based on gene expression profiles using Linear Discriminant Analysis (LDA)
6. **Functional Annotation of Top Discriminant Genes**: Exploring the biological relevance of top discriminant genes.

## Repository Contents

- **Analysis.R**: This R script contains all the code for the analysis pipeline.
- **README.md**: You are here! This file provides an overview of the project and its contents.
## Usage
