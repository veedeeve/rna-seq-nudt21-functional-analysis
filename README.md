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
![MultiQC Alignment Summary](results/figures/multiqc-alignment-summary.png)
