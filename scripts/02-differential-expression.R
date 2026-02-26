#############################################
# RNA-seq Differential Expression
#############################################

# ---------- 0) Install packages (run once) ----------
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("sva", "edgeR", "limma", "Biobase", "biomaRt", 
                       "clusterProfiler", "EnhancedVolcano", "org.Hs.eg.db"))

install.packages(c("data.table", "readxl", "stringr", "ggplot2", "ggrepel", 
                   "ggfortify", "ggprism", "pheatmap", "VennDiagram", 
                   "corrplot", "Hmisc", "stats", "tidyverse"))

# ---------- 1) Load libraries ----------
library(data.table)
library(edgeR)
library(limma)
library(ggplot2)
library(ggrepel)
library(ggfortify)
library(stats)
library(sva)
library(EnhancedVolcano)
library(pheatmap)

# ---------- 2) Set working directory ----------
setwd("/Users/vydang/rna-seq-analysis/")

# ---------- 3) Locate the featureCounts output ----------
file_name <- list.files(
  "/Users/vydang/rna-seq-analysis/results/counts",
  pattern = "counts_s2.txt$", 
  recursive = TRUE, 
  full.names = TRUE
  )
# Reads the featureCounts table into R.
fc <- fread(file_name, comment.char = "#", data.table = FALSE)

# ---------- 4) Build a clean counts matrix ----------
counts_raw <- fc[, c("Geneid", colnames(fc)[7:10])]
rownames(counts_raw) <- counts_raw$Geneid # name the rows as the geneid
counts_raw$Geneid <- NULL # removes the geneid column

# ---------- 5) Clean sample names ----------
colnames(counts_raw) <- basename(colnames(counts_raw)) 
colnames(counts_raw) <- gsub("\\.bam$", "", colnames(counts_raw)) # removes the .bam

# ---------- 6) Create sample metadata ----------
# tells R which samples belong to which condition.
sample_metadata <- data.frame(
  sample = colnames(counts_raw),
  group = c("Control", "Control", "NUDT21_KD", "NUDT21_KD")
)

# ---------- 7) Save raw counts matrix ----------
dir.create("results/counts", recursive = TRUE, showWarnings = FALSE)
fwrite(counts_raw, "results/counts/counts_raw_matrix.tsv", sep = "\t")

# ---------- 8) Create edgeR object ----------
# DGEList stores counts + library sizes + normalization factors + group.
dge <- DGEList(counts = counts_raw, group = sample_metadata$group)
str(dge)
dge$samples

# ---------- 9) Filter out low-expression genes ----------
# filterByExpr uses group info to decide what to keep.
keep <- filterByExpr(dge, group = sample_metadata$group)
# Subset to kept genes and recalculate library sizes after filtering
dge <- dge[keep, , keep.lib.sizes = FALSE]

# ---------- 10) Normalize (TMM) ----------
# TMM corrects for library size & composition bias so samples are comparable.
dge <- calcNormFactors(dge, method = "TMM")

# ---------- 11) Create the design matrix ----------
# encodes your experiment for linear modeling.
# (Intercept) = mean expression in Control
# conditionNUDT21_KD = (NUDT21_KD - Control)
sample_metadata$group <- factor(sample_metadata$group, 
                                levels = c("Control", "NUDT21_KD"))
levels(sample_metadata$group)
design <- model.matrix(~ group, data = sample_metadata)
colnames(design)

# ---------- 12) voom transform ----------
# voom converts counts -> log2 CPM and estimates mean-variance relationship.
# computes precision weights for each observation.
dir.create("results/differential-gene-expression", recursive = TRUE, showWarnings = FALSE)
dge_v <- voom(dge, design = design, plot = TRUE)
saveRDS(dge_v, "results/differential-gene-expression/differential-gene-expression-voom.rds")

# ---------- 13) Fit linear model per gene ----------
# fits a weighted linear model for each gene using the design matrix.
fit <- lmFit(dge_v, design)

# ---------- 14) Empirical Bayes moderation ----------
# eBayes stabilizes variance estimates by borrowing information across genes.
fit <- eBayes(fit)

# ---------- 15) Extract differential expression results ----------
res <- topTable(fit, coef = "groupNUDT21_KD", number = Inf)
head(res)
res$gene <- rownames(res)
fwrite(res, "results/differential-gene-expression/differential-gene-expression.tsv", sep = "\t")

# ---------- 16) How many significant genes? ----------
# adj.P.Val is FDR (multiple-testing corrected p-value).
sum(res$adj.P.Val < 0.05)
# ~18% of genes are significantly differentially expressed
sum(res$adj.P.Val < 0.05 & res$logFC > 0) # upregulated
sum(res$adj.P.Val < 0.05 & res$logFC < 0) # downregulated
sum(res$adj.P.Val < 0.05 & abs(res$logFC) > 1) # statistically significant
sum(res$adj.P.Val < 0.05 & res$logFC > 1) # strongly upregulated
sum(res$adj.P.Val < 0.05 & res$logFC < -1) # strongly downregulated

# ---------- 17) PCA preparation ----------
expr_matrix <- dge_v$E
expr_t <- t(expr_matrix)

# ---------- 18) Run PCA ----------
pca_res <- prcomp(expr_t, scale. = TRUE)

# ---------- 19) Plot PCA ----------
pca_plot <- autoplot(pca_res,
                     data = sample_metadata,
                     color = "group",
                     shape = "group"
                     ) +
  geom_text_repel(aes(label = sample), size = 4) +
  theme_minimal()
pca_plot

# ---------- 20) Save PCA plot ----------
ggsave("results/differential-gene-expression/pca-plot.png", 
       plot = pca_plot, 
       width = 16, 
       height = 14, 
       units = "cm")

# ---------- 21) Volcano plot ----------
p_threshold <- 0.05
fc_threshold <- 2

# prepares data
res$adj.P.Val <- pmax(res$adj.P.Val, 1e-300)
res$neg_log10_adj.P.Val <- -log10(res$adj.P.Val)
res$significance <- ifelse(
  res$adj.P.Val < p_threshold & abs(res$logFC) > log2(fc_threshold),
  ifelse(res$logFC > 0, "Upregulated", "Downregulated"),
  "Non-Significant")
res$significance <- factor(
  res$significance,
  levels = c("Downregulated", "Non-Significant", "Upregulated")
)

# label top 20 most significant
genes_to_label <- res$gene[order(res$adj.P.Val)][1:20]
volcano_plot <- ggplot(res, aes(x = logFC,
                           y = neg_log10_adj.P.Val)) +
  geom_point(aes(color = significance), size = 0.6, alpha = 0.6) +
  # threshold lines
  scale_color_manual(values = c(
    "Downregulated"   = "blue",
    "Non-Significant" = "grey50",
    "Upregulated"     = "red"
  )) +
  geom_vline(xintercept = c(-log2(fc_threshold), log2(fc_threshold)),
             linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(p_threshold),
             linetype = "dashed", color = "grey50") +
  # gene labels
  geom_text_repel(data = subset(res, gene %in% genes_to_label),
                  aes(label = gene, color = significance),
                  size = 3,
                  max.overlaps = Inf,
                  box.padding = 0.5) +
  labs(x = expression(log[2]~"Fold Change"),
       y = expression(-log[10]~"adjusted p-value"),
       color = "Differential Expression"
       ) +
  coord_cartesian(ylim = c(0, 20)) +
  theme_minimal() +
  theme(legend.position = "right",
        panel.grid.minor = element_blank(),
        text = element_text(size = 12),
        axis.text = element_text(size = 10))
volcano_plot

ggsave(
  "results/differential-gene-expression/volcano-plot.png",
  volcano_plot,
  width = 14,
  height = 14,
  dpi = 400
)

# ---------- 22) Heatmap ----------
col_annot_df <- data.frame(
  Treatment = sample_metadata$group,
  row.names = sample_metadata$sample
)

deg_sig <- res$gene[abs(res$logFC) > log2(fc_threshold) & 
    res$adj.P.Val < p_threshold]
expr_deg <- dge_v$E[deg_sig, rownames(col_annot_df), drop = FALSE]

heatmap <- pheatmap(
  expr_deg,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = TRUE,
  color = colorRampPalette(c("navy", "white", "red"))(100),
  annotation_col = col_annot_df,
  main = "Differential Expression Heatmap",
  filename = "results/differential-gene-expression/heatmap.png",
  width = 14,
  height = 14
)
