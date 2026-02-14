# RNA-seq Differential Expression Analysis Using DESeq2

## Overview

This repository contains a reproducible RNA-seq differential gene expression workflow implemented in R using the DESeq2 package.

The analysis compares gene expression between two biological timepoints (T1 vs T24) and includes statistical modelling, visualisation, and functional annotation integration.

This project demonstrates a complete end-to-end differential expression pipeline suitable for transcriptomics analysis in bacterial systems.

---

## Objectives

- Perform library size normalisation
- Estimate dispersion parameters
- Conduct Wald test for differential expression
- Identify significantly upregulated and downregulated genes
- Integrate functional annotation
- Analyse plasmid-associated gene expression


---

## Methods

### 1. Data Input
- Raw count matrix (RData format)
- Gene annotation table (RData format)

### 2. DESeq2 Workflow
- Construction of `DESeqDataSet`
- Size factor normalisation
- Dispersion estimation
- Negative binomial Wald testing
- Multiple testing correction (Benjamini–Hochberg)

### 3. Filtering Criteria
Significant genes defined as:

- Adjusted p-value (padj) < 0.05
- |log2 Fold Change| > 1

### 4. Transformations & Visualisation
- rlog transformation
- PCA plot
- Sample distance heatmap
- MA plot
- Plasmid gene expression barplots

### 5. Annotation Integration
- Merge DESeq2 output with functional annotation
- Extraction of top differentially expressed genes
- Focused analysis of plasmid-associated genes

---



