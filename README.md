# TCGA-BRCA-RNAseq-DEG-Analysis-R

## Project Overview

This project performs RNA-seq differential expression analysis on breast cancer data using R and Bioconductor. The dataset was retrieved from the TCGA-BRCA project, which contains breast cancer transcriptomic data.

The main objective of this project is to compare breast cancer tumor samples with normal tissue samples and identify genes that are significantly upregulated or downregulated.

This project demonstrates a complete beginner-to-intermediate bioinformatics workflow including data retrieval, preprocessing, differential gene expression analysis, and visualization.

---

## Project Title

**TCGA-BRCA RNA-seq Differential Expression Analysis Using R**

---

## Dataset

The data used in this project was retrieved from:

**TCGA-BRCA: The Cancer Genome Atlas Breast Invasive Carcinoma Project**

For computational efficiency, a small subset of samples was used:

- 3 Primary Tumor samples
- 3 Solid Tissue Normal samples

This subset was used to demonstrate the complete RNA-seq analysis workflow in a lightweight and reproducible way.

---

## Objective

The main objectives of this project are:

1. Retrieve breast cancer RNA-seq data from TCGA using R.
2. Prepare count matrix and sample metadata.
3. Perform differential expression analysis using DESeq2.
4. Identify upregulated and downregulated genes.
5. Generate visualizations such as PCA plot and volcano plot.
6. Export results as CSV files for further interpretation.

---

## Tools and Libraries Used

### 1. R

R was used as the main programming language for the complete analysis. It is widely used in bioinformatics, genomics, statistics, and data visualization.

### 2. BiocManager

`BiocManager` was used to install Bioconductor packages required for bioinformatics analysis.

Bioconductor packages are commonly used for genomic and transcriptomic data analysis.

### 3. TCGAbiolinks

`TCGAbiolinks` was used to query, download, and prepare TCGA-BRCA RNA-seq data directly from the Genomic Data Commons.

In this project, it was used for:

- Accessing TCGA-BRCA data
- Querying RNA-seq gene expression data
- Downloading selected tumor and normal samples
- Preparing data into a format suitable for downstream analysis

### 4. SummarizedExperiment

`SummarizedExperiment` was used to handle the TCGA data object prepared by TCGAbiolinks.

It stores:

- Count matrix
- Sample metadata
- Gene information

This makes it easier to manage RNA-seq data in a structured format.

### 5. DESeq2

`DESeq2` was used for differential gene expression analysis.

In this project, DESeq2 was used to:

- Create a DESeq2 dataset object
- Normalize raw count data
- Estimate variation between samples
- Compare tumor vs normal conditions
- Identify differentially expressed genes

Genes were classified as significant based on:

- Adjusted p-value < 0.05
- Absolute log2 fold change > 1

### 6. dplyr

`dplyr` was used for data manipulation and filtering.

It helped in:

- Selecting tumor and normal samples
- Filtering data
- Creating sample groups
- Organizing result tables

### 7. ggplot2

`ggplot2` was used for data visualization.

In this project, it was used to generate:

- Volcano plot
- Customized visualization of differential expression results

The volcano plot shows genes based on log2 fold change and adjusted p-value.

### 8. pheatmap

`pheatmap` was used to generate a heatmap of top differentially expressed genes.

The heatmap helps visualize gene expression patterns across tumor and normal samples.

---

## Workflow

The project follows this workflow:

1. Install and load required R packages.
2. Query TCGA-BRCA RNA-seq data using TCGAbiolinks.
3. Select a small subset of tumor and normal samples.
4. Download and prepare the selected samples.
5. Extract raw count matrix.
6. Extract and prepare sample metadata.
7. Create tumor and normal condition labels.
8. Perform differential expression analysis using DESeq2.
9. Identify significant differentially expressed genes.
10. Classify genes as upregulated or downregulated.
11. Generate PCA plot.
12. Generate volcano plot.
13. Save all results and figures.

---

## Repository Structure

```text
Breast-Cancer-RNAseq-DEG-Analysis-R/
│
├── README.md
│
├── scripts/
│   └── TCGA_BRCA_RNAseq_DEG_Analysis.R
│
├── data/
│   ├── TCGA_BRCA_small_raw_counts.csv
│   └── TCGA_BRCA_small_sample_metadata.csv
│
├── results/
│   ├── TCGA_BRCA_small_DESeq2_results.csv
│   ├── TCGA_BRCA_small_significant_DEGs.csv
│   ├── TCGA_BRCA_small_upregulated_genes.csv
│   └── TCGA_BRCA_small_downregulated_genes.csv
│
└── figures/
    ├── TCGA_BRCA_small_pca_plot.png
    ├── TCGA_BRCA_small_volcano_plot.png
    └── TCGA_BRCA_small_heatmap_top_30_genes.png

## Visualizations and Interpretation

### PCA Plot

![PCA Plot](figures/TCGA_BRCA_small_pca_plot.png)

The PCA plot was used to visualize the overall gene expression variation between TCGA-BRCA tumor and normal samples.

**Interpretation:**

- Each point represents one RNA-seq sample.
- The samples are colored based on their biological group.
- **Pink/Salmon points represent Normal samples.**
- **Blue/Turquoise points represent Tumor samples.**
- PC1 explains **77% of the total variance**, while PC2 explains **9% of the variance**.
- The tumor and normal samples show clear separation along PC1, indicating that their overall gene expression profiles are different.
- This separation suggests that breast cancer tumor samples have distinct transcriptomic patterns compared to normal tissue samples.

---

### Volcano Plot

![Volcano Plot](figures/TCGA_BRCA_small_volcano_plot.png)

The volcano plot shows differentially expressed genes between TCGA-BRCA tumor and normal samples.

**Interpretation:**

- Each point represents one gene.
- The X-axis shows **log2 fold change**.
- The Y-axis shows **-log10 adjusted p-value**.
- Genes farther from the center show larger expression differences.
- Genes higher on the plot are more statistically significant.

**Color meaning:**

- **Blue points:** Upregulated genes  
  These genes show higher expression in tumor samples compared to normal samples.

- **Red/Pink points:** Downregulated genes  
  These genes show lower expression in tumor samples compared to normal samples.

- **Green points:** Not significant genes  
  These genes do not pass the selected significance threshold.

The significance criteria used in this project were:

```text
Adjusted p-value < 0.05
Absolute log2 fold change > 1
