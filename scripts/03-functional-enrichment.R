#############################################
# RNA-seq Enrichment Analysis
#############################################

# ---------- 0) Setup ----------
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
renv::restore()

# ---------- 1) Load libraries ----------
library(ggplot2)
library(enrichplot)
library(org.Hs.eg.db)
library(data.table)
library(clusterProfiler)

# ---------- 2) Load DEG table ----------
deg <- fread(file.path("analysis", "differential_expression","differential-gene-expression.tsv"))

# ---------- 3) Rank genes for Gene Set Enrichment Analysis ----------
gene_rank <- deg$t
names(gene_rank) <- deg$gene
gene_rank <- sort(gene_rank, decreasing = TRUE)
# Break ties due to identical t-value (71.85% ties)
set.seed(42)
gene_rank2 <- gene_rank + rnorm(length(gene_rank), 0, 1e-10)
gene_rank2 <- sort(gene_rank2, decreasing = TRUE)

# ---------- 4) Gene Set Enrichment Analysis ----------
enrich_gsea <- gseGO(
  geneList = gene_rank2, 
  OrgDb = org.Hs.eg.db, 
  ont = "ALL", 
  pvalueCutoff = 0.05, 
  keyType ="ENSEMBL", 
  eps = 0,
  verbose= FALSE)
dotplot_enrich_gsea <- dotplot(enrich_gsea,
                               showCategory = 10,
                               orderBy="NES")
dir.create(file.path("analysis", "gsea"),
           recursive = TRUE,
           showWarnings = FALSE)
ggsave(file.path("analysis", "gsea", "dotplot-gsea.png"),
       dotplot_enrich_gsea,
       device = "png",
       units = "cm",
       width = 16,
       height = 18)

# ---------- 5) Plot GSEA enrichment map ----------
enrich_gsea_sim <- pairwise_termsim(enrich_gsea)
emapplot(enrich_gsea_sim,
         showCategory = 20,
        node_label_size = 3)
ggsave(file.path("analysis", "gsea", "emapplot-gsea.png"), width=18, height=16, units="cm")

# ---------- 7) Top 10 Upregulated and Downregulated ----------
enrich_go_gsea_tbl <- as.data.table(as.data.frame(enrich_gsea))
#  Filter significant pathways
enrich_go_gsea_tbl_sig <- enrich_go_gsea_tbl[p.adjust < 0.05]
up_pathways <- enrich_go_gsea_tbl_sig[NES > 0][order(-NES)]
down_pathways <- enrich_go_gsea_tbl_sig[NES < 0][order(NES)]
fwrite(enrich_go_gsea_tbl, file.path("analysis", "gsea", "gsea_GO_all.tsv"), sep = "\t")
top_up_10 <- up_pathways[1:10]
top_down_10 <- down_pathways[1:10]
fwrite(top_up_10, file.path("analysis", "gsea", "GO_upregulated_top10.tsv"), sep = "\t")
fwrite(top_down_10, file.path("analysis", "gsea", "GO_downregulated_top10.tsv"), sep = "\t")

