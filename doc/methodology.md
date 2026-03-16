# RNA-seq Differential Expression (NUDT21 Knockdown) - Methodology

## Dataset
RNA-seq samples retrieved from SRA (BioProject: PRJNA1305742)
 - 2 Control
 - 2 NUDT21 knockdown  
**Genome reference**: GRCh38  
**Annotation**: Ensembl Homo_sapiens.GRCh38

---

## 1. Preprocessing Data (BASH)
*Ensures high-quality reads prior to alignment* 

- Raw read quality assessed with `FastQC` pre and post-trimming to confirm adapter removal and quality improvement
- Adapters trimmed from raw reads using `Trimmomatic` 
- Parameters used: `ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:36`

## 2. Alignment & Annotation
*Pair-end RNA-seq reads aligned to human reference genome, sorted BAM files generated with samtools. Produce count-matrix for downstream DEG analysis by quantifying read counts*  
- Reads aligned to GRCh38 index using `HISAT2`
    - Preferred HISAT2 over STAR due to memory requirements
    - Splice-aware aligner designed for RNA-seq alignment
- `samtools flagstat` assessed alignment rates
- `MultiQC` aggregated fastQC and HISAT2 alignment summaries
- `featureCounts` summarized reads to gene-level counts
    - Strandedness set to reverse stranded `-s 2` to match library prep
    - Extracted sample count from column 7 
    - count matrix was inputed into R for differential expression analysis

## 3. Differential Expression (R)
- Low-count genes were filtered & normalization applied using `edgeR`
    - `filterByExpr` group labels to set CPM threshold (retains genes expressed above a minimum level, prevent noise from low-expressed genes
    - `keep.lib.sizes = FALSE` recalculates library sizes
    - TMM normalization corrects for compositional difference in library size
- PCA conducted on normalized log-CPM (counts per million) for quality control
    - Confirms sample clustering by condition and detects outliers
- `voom` transforms count data into log-CPM and assigns each gene a weight based on its expression level
    - In RNA-seq data, highly expressed genes tends to be more variable. `voom` ensures all genes are modeled fairly
- Linear model fit using `limma-voom` to contrast NUDT21 knockdown vs control
    - Used over DESeq2 due to small sample size
    - `eBayes` stabilized the variance
- Genes ranked by t-statistics from `eBayes`
    - t-statistics balances effect size & confidence
    - 71.85% genes had identical t-values, a negligible random noise (1e-10) added to break ties

### 3.1 Volcano Plot
- Significant thresholds: FDR < 0.05 & fold change > 2 (log2FC > 1)
- Top 20 most significant genes labeled by adjusted p-value

### 3.2 Heatmap
- Genes filtered by FDR < 0.05 & fold-change > 2 to show strongly significant DEGs
- 
## 4. Functional Enrichment (R)
- GSEA performed on the ranked gene list using `clusterProfiler`
    - Preferred over ORA as GSEA use full ranked list without cutoff
    - ran enrichment on all GO ontologies
    - Exported top 10 upregrulated and top 10 downregulated 
- Genes are annotated using `org.Hs.eg.db`
- Pairwise term similarity calculated with `pairwise_termsim()` and visualized as `emapplot`
- Visualizations were made using `enrichplot`
    - Dot plot, GSEA enrichment plot and concept maps
 
---

## Parameter Decisions
| Step | Parameter | Value | Reason |
| :--- | :------------ | :---- | :-------- |
| Trimmomatic | ILLUMINACLIP | TruSeq3-PE.fa:2:30:10 | Removes TruSeq paired-end adapters |
| Trimmomatic | SLIDINGWINDOW | 4:20 | Balances quality stringency with read retention |
| Trimmomatic | MINLEN | 36 | Removes reads too short for alignment
| featureCounts | strandedness | -s 2 (reverse stranded) | Must match library prep to avoid systematic miscounting |
| featureCounts | paired-end | enabled | Matches paired-end library preparation |
| edgeR | filterByExpr | min group size = 2 | Retains genes expressed in at least all samples of the smallest group |
| edgeR | keep.lib.sizes | FALSE | Recalculates library sizes after filtering for accuracy |
| edgeR | normalization | TMM | Corrects for compositional bias in library sizes |
| limma-voom | trend | enabled | Accounts for mean-variance relationship in count data |
| clusterProfiler | ranking metric | t-statistic | Balances effect size, directionality, and confidence |
| clusterProfiler | ont | ALL | Captures BP, MF, and CC ontologies in one run |
| clusterProfiler | eps | 0 | Preserves very significant results below default p-value floor |
| clusterProfiler | pvalueCutoff | 0.05 | Standard significance threshold for enrichment |
| Volcano plot | p-value floor | 1e-300 | Prevents -log10(0) = Inf from breaking the plot |
| Volcano plot | y-axis cap | 20 | Improves readability without losing meaningful data |
| Volcano plot | FC threshold | 2 (log2FC > 1) | Standard cutoff for biologically meaningful fold-change |
| Heatmap | scale | row | Z-score per gene; prevents highly expressed genes dominating color range |
| Heatmap | cluster_cols | FALSE | Preserves sample grouping order (Control vs. KD) |
| Heatmap | cluster_rows | TRUE | Groups co-regulated genes together visually |

