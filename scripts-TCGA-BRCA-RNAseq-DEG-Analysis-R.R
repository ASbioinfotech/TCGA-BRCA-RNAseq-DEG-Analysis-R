# ============================================================
# PART 1: SETUP
# TCGA-BRCA RNA-seq Differential Expression Analysis
# ============================================================

# Install BiocManager if not installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install required Bioconductor packages
BiocManager::install(
  c("TCGAbiolinks", "SummarizedExperiment"),
  ask = FALSE,
  update = FALSE
)

# Install required CRAN packages
install.packages(
  c("dplyr", "ggplot2"),
  dependencies = TRUE
)

library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(ggplot2)

getwd()

dir.create("data", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

list.files()

query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)

query_results <- getResults(query)

dim(query_results)
head(query_results)
table(query_results$sample_type)

dim(query_results)
head(query_results)
table(query_results$sample_type)


# ============================================================
# PART 6: SELECT SMALL TEST DATASET
# 3 Primary Tumor samples + 3 Solid Tissue Normal samples
# ============================================================

tumor_samples <- query_results %>%
  filter(sample_type == "Primary Tumor") %>%
  slice(1:3)

normal_samples <- query_results %>%
  filter(sample_type == "Solid Tissue Normal") %>%
  slice(1:3)

selected_samples <- c(
  tumor_samples$cases,
  normal_samples$cases
)

selected_samples

length(selected_samples)

tumor_samples[, c("cases", "sample_type")]
normal_samples[, c("cases", "sample_type")]

# ============================================================
# PART 7: CREATE SMALL QUERY
# ============================================================

query_small <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal"),
  barcode = selected_samples
)

small_results <- getResults(query_small)

dim(small_results)
head(small_results)
table(small_results$sample_type)

# ============================================================
# PART 8: DOWNLOAD SMALL TCGA-BRCA DATA
# ============================================================

GDCdownload(query_small)

# ============================================================
# PART 9: PREPARE DOWNLOADED DATA
# ============================================================

tcga_data <- GDCprepare(query_small)

tcga_data

# ============================================================
# PART 10: EXTRACT COUNT MATRIX
# ============================================================

count_data <- assay(tcga_data)

dim(count_data)
count_data[1:5, 1:5]

# ============================================================
# PART 11: EXTRACT SAMPLE METADATA
# ============================================================

sample_metadata <- as.data.frame(colData(tcga_data))

dim(sample_metadata)
head(sample_metadata)
table(sample_metadata$sample_type)

# ============================================================
# PART 12: CREATE CONDITION COLUMN
# ============================================================

sample_metadata$condition <- ifelse(
  sample_metadata$sample_type == "Primary Tumor",
  "Tumor",
  "Normal"
)

sample_metadata$condition <- factor(
  sample_metadata$condition,
  levels = c("Normal", "Tumor")
)

table(sample_metadata$condition)

# Save simplified metadata only

sample_metadata_simple <- data.frame(
  sample_id = rownames(sample_metadata),
  sample_type = sample_metadata$sample_type,
  condition = sample_metadata$condition
)

write.csv(
  sample_metadata_simple,
  "data/TCGA_BRCA_small_sample_metadata.csv",
  row.names = FALSE
)

list.files("data")

colData = sample_metadata

sample_metadata_simple

# ============================================================
# PART 15: CHECK SAMPLE MATCHING BEFORE DESEQ2
# ============================================================

colnames(count_data)
rownames(sample_metadata)

all(colnames(count_data) == rownames(sample_metadata))

sample_metadata <- sample_metadata[colnames(count_data), ]

all(colnames(count_data) == rownames(sample_metadata))

# ============================================================
# PART 16: LOAD DESEQ2
# ============================================================

BiocManager::install(
  "DESeq2",
  ask = FALSE,
  update = FALSE
)

library(DESeq2)

# ============================================================
# PART 17: CREATE DESEQ2 DATASET
# ============================================================

dds <- DESeqDataSetFromMatrix(
  countData = round(count_data),
  colData = sample_metadata,
  design = ~ condition
)

dds

# ============================================================
# PART 18: FILTER LOW-COUNT GENES
# ============================================================

dds <- dds[rowSums(counts(dds)) >= 10, ]

dim(dds)

# ============================================================
# PART 19: RUN DIFFERENTIAL EXPRESSION ANALYSIS
# ============================================================

dds <- DESeq(dds)

# ============================================================
# PART 20: EXTRACT RESULTS
# ============================================================

res <- results(
  dds,
  contrast = c("condition", "Tumor", "Normal")
)

res <- res[order(res$padj), ]

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

res_df <- na.omit(res_df)

dim(res_df)
head(res_df)

# ============================================================
# PART 21: CLASSIFY GENES
# ============================================================

res_df$regulation <- "Not Significant"

res_df$regulation[
  res_df$padj < 0.05 & res_df$log2FoldChange > 1
] <- "Upregulated"

res_df$regulation[
  res_df$padj < 0.05 & res_df$log2FoldChange < -1
] <- "Downregulated"

table(res_df$regulation)

# ============================================================
# PART 22: SAVE DEG RESULTS
# ============================================================

significant_genes <- res_df[
  res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1,
]

upregulated_genes <- res_df[
  res_df$regulation == "Upregulated",
]

downregulated_genes <- res_df[
  res_df$regulation == "Downregulated",
]

write.csv(
  res_df,
  "results/TCGA_BRCA_small_DESeq2_results.csv",
  row.names = FALSE
)

write.csv(
  significant_genes,
  "results/TCGA_BRCA_small_significant_DEGs.csv",
  row.names = FALSE
)

write.csv(
  upregulated_genes,
  "results/TCGA_BRCA_small_upregulated_genes.csv",
  row.names = FALSE
)

write.csv(
  downregulated_genes,
  "results/TCGA_BRCA_small_downregulated_genes.csv",
  row.names = FALSE
)

list.files("results")

# ============================================================
# PART 23: SAVE ANALYSIS SUMMARY
# ============================================================

summary_table <- data.frame(
  Category = c(
    "Total genes analyzed",
    "Significant DEGs",
    "Upregulated genes",
    "Downregulated genes"
  ),
  Count = c(
    nrow(res_df),
    nrow(significant_genes),
    nrow(upregulated_genes),
    nrow(downregulated_genes)
  )
)

print(summary_table)

write.csv(
  summary_table,
  "results/TCGA_BRCA_small_analysis_summary.csv",
  row.names = FALSE
)

list.files("results")

# ============================================================
# PART 24: PCA PLOT
# ============================================================

vsd <- vst(dds, blind = FALSE)

pca_plot <- plotPCA(
  vsd,
  intgroup = "condition"
) +
  ggtitle("PCA Plot: TCGA-BRCA Tumor vs Normal") +
  theme_minimal()

print(pca_plot)

ggsave(
  "figures/TCGA_BRCA_small_pca_plot.png",
  plot = pca_plot,
  width = 7,
  height = 5,
  dpi = 300
)

list.files("figures")

# ============================================================
# PART 25: VOLCANO PLOT
# ============================================================

volcano_plot <- ggplot(
  res_df,
  aes(
    x = log2FoldChange,
    y = -log10(padj),
    color = regulation
  )
) +
  geom_point(alpha = 0.7, size = 1.5) +
  theme_minimal() +
  labs(
    title = "Volcano Plot: TCGA-BRCA Tumor vs Normal",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Regulation"
  )

print(volcano_plot)

ggsave(
  "figures/TCGA_BRCA_small_volcano_plot.png",
  plot = volcano_plot,
  width = 8,
  height = 6,
  dpi = 300
)

list.files("figures")

# ============================================================
# HEATMAP CODE: TOP 30 GENES
# ============================================================

if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}

library(pheatmap)

dir.create("figures", showWarnings = FALSE)

vsd_matrix <- assay(vsd)

top_genes <- res_df[order(res_df$padj), "gene_id"]
top_genes <- top_genes[!is.na(top_genes)]

common_genes <- intersect(top_genes, rownames(vsd_matrix))

if (length(common_genes) == 0) {
  stop("No matching genes found between res_df and vsd matrix.")
}

top_30_genes <- common_genes[1:min(30, length(common_genes))]

heatmap_data <- vsd_matrix[top_30_genes, ]

annotation_col <- data.frame(
  Condition = sample_metadata$condition
)

rownames(annotation_col) <- rownames(sample_metadata)

png(
  filename = "figures/TCGA_BRCA_small_heatmap_top_30_genes.png",
  width = 1000,
  height = 800
)

pheatmap(
  heatmap_data,
  annotation_col = annotation_col,
  scale = "row",
  show_rownames = FALSE,
  show_colnames = FALSE,
  main = "Top 30 Differentially Expressed Genes"
)

dev.off()

list.files("figures")

# ============================================================
# DOWNLOAD / EXPORT PROJECT RESULTS AS ZIP
# ============================================================

# Check current working directory
getwd()

# Create final output folder
dir.create("final_project_outputs", showWarnings = FALSE)

# Copy data files
if (dir.exists("data")) {
  file.copy(
    from = list.files("data", full.names = TRUE),
    to = "final_project_outputs",
    overwrite = TRUE
  )
}

# Copy result files
if (dir.exists("results")) {
  file.copy(
    from = list.files("results", full.names = TRUE),
    to = "final_project_outputs",
    overwrite = TRUE
  )
}

# Copy figure files
if (dir.exists("figures")) {
  file.copy(
    from = list.files("figures", full.names = TRUE),
    to = "final_project_outputs",
    overwrite = TRUE
  )
}

# Check copied files
list.files("final_project_outputs")

# Create ZIP file
zip(
  zipfile = "TCGA_BRCA_RNAseq_Project_Outputs.zip",
  files = list.files("final_project_outputs", full.names = TRUE)
)

# Check ZIP file created
list.files(pattern = "TCGA_BRCA_RNAseq_Project_Outputs.zip")