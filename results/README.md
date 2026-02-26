# Results Directory

This folder contains all output generated from the RNA-seq analysis workflow.  
## figures/

Publication-style visualizations generated from differential expression and enrichment analyses:
- PCA plot (voom-transformed expression)
- Volcano plot (log2FC vs FDR)
- Heatmap of significant genes
- GSEA dotplot
- GSEA enrichment map
- Mean–variance (voom) plot
- MultiQC alignment summary
These figures summarize sample separation, differential expression patterns, and pathway-level biological interpretation.

---

## docs/

Tabular result files used for interpretation and downstream analysis:

- `top10_upregulated.tsv`
- `top10_downregulated.tsv`
- `gsea_GO_all.tsv`
- `gsea_GO_upregulated_top10.tsv`
- `gsea_GO_downregulated_top10.tsv`
These files contain ranked differential expression results and gene set enrichment outputs for reproducibility and further analysis.
