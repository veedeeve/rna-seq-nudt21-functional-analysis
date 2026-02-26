# RNA-seq Functional Analysis of NUDT21 Knockdown

## Project Overview
This project implements a fully reproducible RNA-seq analysis workflow to investigate the functional impact of **NUDT21** knockdown in human samples.

This pipeline includes: 
1. Preprocessing Data
2. Alignment & Annotation
3. Differential Expression
4. Functional Enrichment

## Objectives

- Process raw RNA-seq data from SRA
- Perform alignment and gene-level quantification
- Model differential expression using limma-voom
- Identify coordinated biological programs using GO-based GSEA

## Dataset
RNA-seq samples retrieved from SRA (BioProject: PRJNA1305742)
 - 2 Control
 - 2 NUDT21 knockdown
**Genome reference**: GRCh38
**Annotation**: Ensembl Homo_sapiens.GRCh38

## Project Structure
### 1. Preprocessing Data (BASH)
*Ensures high-quality reads prior to alignment*  
#### 1.1 Data Download - Retrieved from NCBI SRA using ```prefetch```  
#### 1.2 Quality Control - Assessed using ```FastQC```  
#### 1.3 Adapter & Quality Trimming - Performed using ```Trimmomatic```  
#### 1.4 Post-trimming quality control - Re-evaluated using ```FastQC```  

### 2. Alignment & Annotation
*Pair-end RNA-seq reads aligned to human reference genome, sorted BAM files generated with samtools. Produce count-matrix for downstream DEG analysis by quantifying read counts*  
#### 2.1 Alignment - Align to GRCh38 index using ```HISAT2```  
#### 2.2 Alignment Quality Assessment - ```samtools flagstat``` + ```MultiQC```  
#### 2.3 Gene-level Quantification - ```featureCounts```  

### 3. Differential Expression (R)
#### 3.1 Filtering + TMM normalization - ```edgeR```  
#### 3.3 Linear modeling - ```limma-voom```  
#### 3.4 T-statistics ranking  
#### 3.5 PCA for QC

### 4. Functional Enrichment (R)
#### GO Gene Set Enrichment Analysis - ``` clusterProfiler``` 
#### Human gene annotation - ``` org.Hs.eg.db``` 
#### Enrichment visualization - ``` enrichplot ``` 

## File Structure
```
rna-seq-analysis/
├── data/
│   ├── ref/
│   ├── annot/
│   ├── trimmed/
├── results/
│   ├── figures/
├── scripts/
│   ├── 01-rna-seq-processing.sh
│   ├── 02-differential-expression.R
│   └── 03-functional-enrichment.R
└── README.md
```

## Results
### Alignment & Quality Control
All samples demonstrated high alignment efficiency (98.7–98.9% mapped reads) to the GRCh38 reference genome. Sequencing depth ranged from 56–85 million mapped reads per sample, providing sufficient coverage for differential gene expression.  
<p align="center">
  <img src="results/figures/multiqc-alignment-summary.png" width="600">
</p> 

### Mean–Variance Modeling (voom)
The voom transformation modeled the mean–variance relationship of log-counts and generated precision weights for linear modeling, ensuring appropriate stabilization prior to differential expression analysis.
<p align="center">
  <img src="results/figures/voom-mean-variance.png" width="600">
</p> 

### Significant Genes
A total of 18,180 genes were tested for differential expression.  
- **11,450 genes (63%)** were significant at FDR < 0.05  
- 5,965 genes were upregulated  
- 5,485 genes were downregulated  
- 2,953 genes showed |log2FC| > 1  
- 1,701 genes were strongly upregulated (log2FC > 1)  
- 1,252 genes were strongly downregulated (log2FC < -1)  

These results indicate widespread transcriptomic changes following NUDT21 knockdown.

### Principal Component Analysis (PCA)
PCA of voom-transformed expression values revealed strong separation between control and NUDT21 knockdown samples along PC1 (74.04% variance explained). Reducing NUDT21 caused a major change in gene activity.
<p align="center">
  <img src="results/figures/pca-plot.png" width="600">
</p> 

### Volcano Plot
The volcano plot illustrates the magnitude and statistical significance of gene expression changes. A large proportion of genes exhibited significant differential expression, consistent with extensive regulatory disruption.  
<p align="center">
  <img src="results/figures/volcano-plot.png" width="600">
</p> 

### Heatmap of Differentially Expressed Genes
Hierarchical clustering of significant genes shows clear separation between control and knockdown samples, confirming consistent expression patterns within groups and strong differential regulation between conditions. The genes flipped their behavior based on the condition, consistently changed gene programs.  
<p align="center">
  <img src="results/figures/heatmap.png" width="600">
</p>
