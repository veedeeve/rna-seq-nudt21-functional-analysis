# RNA-seq Functional Analysis of NUDT21 Knockdown

## Project Overview
This project implements a fully reproducible RNA-seq analysis workflow to investigate the functional impact of **NUDT21** knockdown in human samples.

This pipeline includes: 
 - SRA data retrieval
 - Quality control (FastQC, MultiQC)
 - Adapter trimming (Trimmomatic)
 - Alignment to GRCh38 (HISAT2)
 - Gene-level quanitification (featureCounts)
 - Differential expression analysis (limma-voom)
 - Gene Ontology (GO)
 - Gene Set Enrichment Analysis (GSEA)

## Objectives

- Process raw RNA-seq data from SRA
- Perform alignment and gene-level quantification
- Model differential expression using limma-voom
- Identify coordinated biological programs using GO-based GSEA
- Deliver reproducible, publication-ready outputs
- 
## Dataset
RNA-seq samples retrieved from SRA (BioProject: PRJNA1305742)
 - 2 Control
 - 2 NUDT21 knockdown
Genome reference: GRCh38
Annotation: Ensembl Homo_sapiens.GRCh38

## Tools Used
**RNA-seq Processing**

- SRA Toolkit
- FastQC
- Trimmomatic
- HISAT2
- Samtools
- featureCounts
- MultiQC

**Differential Expression**
- edgeR
- limma-voom

**Functional Enrichment**
- clusterProfiler
- org.Hs.eg.db
- enrichplot

## Project Structure
### 1. RNA-seq Processing (Bash/Linux)
- Data retrieval - SRA Toolkit
- Quality control - FastQC + MultiQC
- Adapter trimming - Trimmomatic
- Gene alignment - HISAT2
- BAM Processing - Samtools
- Gene quantification - featureCounts

### 2. Differential Expression (R)
- Filtering + TMM normalization - edgeR
- Linear modeling - limma-voom
- T-statistics ranking
- PCA for QC

### 3. Functional Enrichment
- GO Gene Set Enrichment Analysis - clusterProfiler
- Human gene annotation - org.Hs.eg.db
- Enrichment visualization - enrichplot 

```
myfastqs/
├── sample1
│   ├── sample1_R1.fastq.gz
│   └── sample1_R2.fastq.gz
├── sample2
│   ├── sample2_R1.fastq.gz
│   └── sample2_R2.fastq.gz
└── sample3
    ├── sample3_R1.fastq.gz
    └── sample3_R2.fastq.gz

```
