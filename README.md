# TCGA-BRCA RNA-seq Differential Expression Analysis Using R

## Project Overview

This project performs RNA-seq differential expression analysis on breast cancer data using R and Bioconductor. The dataset was retrieved from the TCGA-BRCA project, which contains breast cancer transcriptomic data.

The main objective of this project is to compare breast cancer tumor samples with normal tissue samples and identify genes that are significantly upregulated or downregulated.

This project demonstrates a beginner-to-intermediate bioinformatics workflow including data retrieval, preprocessing, differential gene expression analysis, and visualization.

---

## Dataset

The data used in this project was retrieved from:

**TCGA-BRCA: The Cancer Genome Atlas Breast Invasive Carcinoma Project**

For computational efficiency, a small subset of samples was used:

* 3 Primary Tumor samples
* 3 Solid Tissue Normal samples

This subset was used to demonstrate the RNA-seq analysis workflow in a lightweight and reproducible way.

---

## Objective

The main objectives of this project are:

1. Retrieve breast cancer RNA-seq data from TCGA using R.
2. Prepare count matrix and sample metadata.
3. Perform differential expression analysis using DESeq2.
4. Identify upregulated and downregulated genes.
5. Generate PCA plot, volcano plot, and heatmap.
6. Export results as CSV files for further interpretation.

---

## Tools and Libraries Used

### R

R was used as the main programming language for the complete analysis. It is widely used in bioinformatics, genomics, statistics, and data visualization.

### BiocManager

`BiocManager` was used to install Bioconductor packages required for bioinformatics analysis.

### TCGAbiolinks

`TCGAbiolinks` was used to query, download, and prepare TCGA-BRCA RNA-seq data directly from the Genomic Data Commons.

In this project, it was used for:

* Accessing TCGA-BRCA data
* Querying RNA-seq gene expression data
* Downloading selected tumor and normal samples
* Preparing data into a format suitable for downstream analysis

### SummarizedExperiment

`SummarizedExperiment` was used to handle the TCGA data object prepared by TCGAbiolinks.

It stores:

* Count matrix
* Sample metadata
* Gene information

This makes it easier to manage RNA-seq data in a structured format.

### DESeq2

`DESeq2` was used for differential gene expression analysis.

In this project, DESeq2 was used to:

* Create a DESeq2 dataset object
* Normalize raw count data
* Estimate variation between samples
* Compare tumor vs normal conditions
* Identify differentially expressed genes

Genes were classified as significant based on:

* Adjusted p-value < 0.05
* Absolute log2 fold change > 1

### dplyr

`dplyr` was used for data manipulation and filtering.

It helped in:

* Selecting tumor and normal samples
* Filtering data
* Creating sample groups
* Organizing result tables

### ggplot2

`ggplot2` was used for data visualization.

In this project, it was used to generate the volcano plot and customize visualization of differential expression results.

### pheatmap

`pheatmap` was used to generate a heatmap of the top differentially expressed genes.

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
13. Generate heatmap.
14. Save results and figures.

---

## Repository Structure

* README.md
* TCGA_BRCA_RNAseq_DEG_Analysis.R
* TCGA_BRCA_small_raw_counts.csv
* TCGA_BRCA_small_sample_metadata.csv
* TCGA_BRCA_small_DESeq2_results.csv
* TCGA_BRCA_small_significant_DEGs.csv
* TCGA_BRCA_small_upregulated_genes.csv
* TCGA_BRCA_small_downregulated_genes.csv
* TCGA_BRCA_small_pca_plot.png
* TCGA_BRCA_small_volcano_plot.png
* TCGA_BRCA_small_heatmap_top_30_genes.png

---

## Output Files

### DESeq2 Results

**TCGA_BRCA_small_DESeq2_results.csv**

This file contains the complete differential expression result table.

Important columns include:

* baseMean: average normalized gene expression
* log2FoldChange: expression difference between tumor and normal samples
* lfcSE: standard error of log2 fold change
* stat: test statistic
* pvalue: raw p-value
* padj: adjusted p-value
* gene_id: Ensembl gene ID
* regulation: upregulated, downregulated, or not significant

---

### Significant DEGs

**TCGA_BRCA_small_significant_DEGs.csv**

This file contains genes that passed the significance criteria:

* Adjusted p-value < 0.05
* Absolute log2 fold change > 1

---

### Upregulated Genes

**TCGA_BRCA_small_upregulated_genes.csv**

This file contains genes that showed higher expression in tumor samples compared to normal samples.

Criteria:

* padj < 0.05
* log2FoldChange > 1

---

### Downregulated Genes

**TCGA_BRCA_small_downregulated_genes.csv**

This file contains genes that showed lower expression in tumor samples compared to normal samples.

Criteria:

* padj < 0.05
* log2FoldChange < -1

---

## Visualizations and Interpretation

This project generated three main visual outputs to understand the gene expression differences between TCGA-BRCA breast cancer tumor and normal samples.

---

### 1. PCA Plot

The PCA plot was generated to visualize the overall gene expression variation between tumor and normal samples.

**Purpose of PCA Plot:**

Principal Component Analysis helps reduce complex gene expression data into major components so that sample-level variation can be visualized easily. In this project, PCA was used to check whether tumor and normal samples show different overall expression patterns.

**Interpretation:**

* Each point in the PCA plot represents one RNA-seq sample.
* Samples are grouped based on their overall gene expression profiles.
* Normal and tumor samples show separation from each other.
* PC1 explains the largest amount of variation in the dataset.
* Clear separation between tumor and normal samples suggests that their transcriptomic profiles are different.
* This supports the idea that breast cancer tumor samples have distinct gene expression patterns compared to normal tissue samples.

**Color Meaning:**

* Pink/Salmon points represent normal samples.
* Blue/Turquoise points represent tumor samples.

**Biological Meaning:**

The separation between tumor and normal samples indicates that breast cancer tissues have major transcriptomic differences compared to normal breast tissue. This provides a strong basis for performing differential gene expression analysis.

---

### 2. Volcano Plot

The volcano plot was generated to identify and visualize differentially expressed genes between tumor and normal samples.

**Purpose of Volcano Plot:**

A volcano plot combines both statistical significance and fold change in one figure. It helps identify genes that are strongly and significantly different between two biological conditions.

**Interpretation:**

* Each point represents one gene.
* The X-axis represents log2 fold change.
* The Y-axis represents -log10 adjusted p-value.
* Genes on the right side show higher expression in tumor samples.
* Genes on the left side show lower expression in tumor samples.
* Genes higher on the plot are more statistically significant.
* Genes far from the center have stronger expression differences.

**Color Meaning:**

* Blue points represent upregulated genes.
* Red/Pink points represent downregulated genes.
* Green points represent genes that are not statistically significant.

**Significance Criteria Used:**

* Adjusted p-value < 0.05
* Absolute log2 fold change > 1

**Biological Meaning:**

The volcano plot shows that several genes are significantly upregulated or downregulated in breast cancer tumor samples. These genes may be involved in breast cancer development, tumor progression, altered cellular pathways, or cancer-related molecular changes.

---

### 3. Heatmap

The heatmap was generated to visualize the expression patterns of the top 30 differentially expressed genes across tumor and normal samples.

**Purpose of Heatmap:**

A heatmap helps compare gene expression levels across samples. It allows us to observe whether selected genes show different expression patterns between tumor and normal groups.

**Interpretation:**

* Each row represents one gene.
* Each column represents one RNA-seq sample.
* The top annotation bar represents the sample group.
* Tumor and normal samples show distinct expression patterns.
* Samples with similar expression profiles cluster together.
* The top differentially expressed genes help separate tumor samples from normal samples.

**Color Meaning:**

* Red/Orange shades represent higher relative gene expression.
* Blue shades represent lower relative gene expression.
* White/Light shades represent intermediate expression.

**Sample Group Colors:**

* Blue/Turquoise annotation represents normal samples.
* Pink/Salmon annotation represents tumor samples.

**Biological Meaning:**

The heatmap supports the differential expression results by showing that the top genes have clear expression differences between tumor and normal samples. This suggests that these genes may contribute to breast cancer-related molecular changes.

---

## Summary of Visual Findings

Overall, the visualizations show clear differences between breast cancer tumor and normal samples.

* The PCA plot shows separation between tumor and normal samples.
* The volcano plot identifies significantly upregulated and downregulated genes.
* The heatmap shows distinct expression patterns among the top differentially expressed genes.

Together, these visualizations support the conclusion that TCGA-BRCA tumor samples have different gene expression profiles compared to normal breast tissue samples.

---

## How to Run This Project

1. Download or clone this repository.
2. Open the R script file in RStudio.
3. Run the script step by step.
4. The script will download selected TCGA-BRCA samples, prepare count data, run DESeq2 analysis, and generate result files and plots.

Main script:

**TCGA_BRCA_RNAseq_DEG_Analysis.R**

---

## Significance Criteria

Genes were considered significant if they met the following criteria:

* Adjusted p-value < 0.05
* Absolute log2 fold change > 1

Gene regulation was classified as:

**Upregulated genes:**

* padj < 0.05
* log2FoldChange > 1

**Downregulated genes:**

* padj < 0.05
* log2FoldChange < -1

**Not significant genes:**

* Genes that did not meet the above criteria

---

## Skills Demonstrated

This project demonstrates skills in:

* Bioinformatics data analysis
* RNA-seq analysis
* TCGA data retrieval
* R programming
* Bioconductor workflow
* DESeq2 differential expression analysis
* Data preprocessing
* Count matrix handling
* Metadata preparation
* Data visualization
* PCA analysis
* Volcano plot generation
* Heatmap visualization
* GitHub project documentation

---

## Importance of This Project

RNA-seq differential expression analysis is an important method in bioinformatics and cancer genomics. It helps researchers identify genes that may be involved in disease progression, tumor biology, and potential biomarker discovery.

This project provides a reproducible workflow for analyzing breast cancer RNA-seq data using real public cancer genomics data.

---

## Limitations

This project uses a small subset of TCGA-BRCA samples for demonstration purposes.

Only:

* 3 tumor samples
* 3 normal samples

were used to keep the workflow lightweight and easy to run on a personal computer.

Because of the small sample size, the results should be considered as a demonstration of the workflow rather than a final biological conclusion.

---

## Future Work

Future improvements may include:

* Expanding the analysis to the full TCGA-BRCA cohort
* Adding GO and KEGG pathway enrichment analysis
* Converting Ensembl IDs to official gene symbols
* Adding interactive visualizations
* Building a Shiny dashboard
* Including biological interpretation of top genes
* Comparing results with published breast cancer biomarkers

---

## Conclusion

This project demonstrates an end-to-end RNA-seq differential expression analysis workflow using TCGA-BRCA breast cancer data. It includes database-driven data retrieval, differential expression analysis using DESeq2, and visualization of results using PCA, volcano plot, and heatmap.

The project highlights practical skills in R programming, Bioconductor, transcriptomics, and cancer genomics analysis.

---

## Author

**Anushka Singh**
Bioinformatic Student


