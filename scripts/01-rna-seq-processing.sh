#!/usr/bin/env bash
set -Eeuo pipefail
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

cd "$(dirname "$0")"
source config.sh 

# ####################################
# 1. Download (SRR: SRR34987610)
######################################

# Control files

log "Downloading reads ..."

: << 'COMMENT'
prefetch "$READS1" --output-directory "$RAW_DATA_DIR" > "$LOGS_DIR/01_download_read1.log" 2>&1
fasterq-dump "$RAW_DATA_DIR/$READS1/" \
 --split-files \
 -O "$RAW_DATA_DIR" \
 --threads "$THREADS" \
 > "$LOGS_DIR/02_split_file.log" 2>&1
gzip "$RAW_DATA_DIR/$READS1"*.fastq

rm -rf "$RAW_DATA_DIR/$READS1/"

log "Read 1 completed"

prefetch "$READS2" --output-directory "$RAW_DATA_DIR" > "$LOGS_DIR/01_download_read2.log" 2>&1
fasterq-dump "$RAW_DATA_DIR/$READS2/" \
 --split-files \
 -O "$RAW_DATA_DIR" \
 --threads "$THREADS" \
 > "$LOGS_DIR/02_read2_split_file.log" 2>&1
gzip "$RAW_DATA_DIR/$READS2"*.fastq

rm -rf "$RAW_DATA_DIR/$READS2"

log "Read 2 complete"
COMMENT

# NUDT21 Knockdown files
prefetch "$READS3" --output-directory "$RAW_DATA_DIR" > "$LOGS_DIR/01_download_read3.log" 2>&1
fasterq-dump "$RAW_DATA_DIR/$READS3/" \
--split-files \
-O "$RAW_DATA_DIR" \
--threads "$THREADS" \
> "$LOGS_DIR/02_read3_split_file.log" 2>&1
gzip "$RAW_DATA_DIR/$READS3"*.fastq

rm -rf "$RAW_DATA_DIR/$READS3"

log "Read 3 complete"

prefetch "$READS4" --output-directory "$RAW_DATA_DIR" > "$LOGS_DIR/01_download_read4.log" 2>&1
fasterq-dump "$RAW_DATA_DIR/$READS4/" \
--split-files \
-O "$RAW_DATA_DIR" \
--threads "$THREADS" \
> "$LOGS_DIR/02_read4_split_file.log" 2>&1
gzip "$RAW_DATA_DIR/$READS4"*.fastq

rm -rf "$RAW_DATA_DIR/$READS4"

log "Read 4 complete"
: << 'COMMENT'
# ---------- 2) Quality Control (before trim)  ----------
mkdir -p ~/rna-seq-analysis/results/fastqc/before-trim
fastqc ~/rna-seq-analysis/data/SRR*_*.fastq.gz -o ~/rna-seq-analysis/data/results/fastqc/before-trim

# ---------- 3) Adapter Trimming  ----------
mkdir -p ~/rna-seq-analysis/data/trimmed
trimmomatic PE -threads 8 -phred33 ~/rna-seq-analysis/data/SRR34987610_1.fastq.gz ~/rna-seq-analysis/data/SRR34987610_2.fastq.gz \
~/rna-seq-analysis/data/trimmed/SRR34987610_1_paired.fastq.gz ~/rna-seq-analysis/data/trimmed/SRR34987610_1_unpaired.fastq.gz \
~/rna-seq-analysis/data/trimmed/SRR34987610_2_paired.fastq.gz ~/rna-seq-analysis/data/trimmed/SRR34987610_2_unpaired.fastq.gz \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
SLIDINGWINDOW:4:20 \
MINLEN:36
echo 'SRR34987610 trimming completed!'
trimmomatic PE -threads 8 -phred33 ~/rna-seq-analysis/data/SRR34987611_1.fastq.gz ~/rna-seq-analysis/data/SRR34987611_2.fastq.gz \
~/rna-seq-analysis/data/trimmed/SRR34987611_1_paired.fastq.gz ~/rna-seq-analysis/data/trimmed/SRR34987611_1_unpaired.fastq.gz \
~/rna-seq-analysis/data/trimmed/SRR34987611_2_paired.fastq.gz ~/rna-seq-analysis/data/trimmed/SRR34987611_2_unpaired.fastq.gz \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
SLIDINGWINDOW:4:20 \
MINLEN:36
echo 'SRR34987611 trimming completed!'
trimmomatic PE -threads 8 -phred33 ~/rna-seq-analysis/data/SRR34987612_1.fastq.gz ~/rna-seq-analysis/data/SRR34987612_2.fastq.gz \
~/rna-seq-analysis/data/trimmed/SRR34987612_1_paired.fastq.gz ~/rna-seq-analysis/data/trimmed/SRR34987612_1_unpaired.fastq.gz \
~/rna-seq-analysis/data/trimmed/SRR34987612_2_paired.fastq.gz ~/rna-seq-analysis/data/trimmed/SRR34987612_2_unpaired.fastq.gz \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
SLIDINGWINDOW:4:20 \
MINLEN:36
echo 'SRR34987612 trimming completed!'
trimmomatic PE -threads 8 -phred33 ~/rna-seq-analysis/data/SRR34987613_1.fastq.gz ~/rna-seq-analysis/data/SRR34987613_2.fastq.gz \
~/rna-seq-analysis/data/trimmed/SRR34987613_1_paired.fastq.gz ~/rna-seq-analysis/data/trimmed/SRR34987613_1_unpaired.fastq.gz \
~/rna-seq-analysis/data/trimmed/SRR34987613_2_paired.fastq.gz ~/rna-seq-analysis/data/trimmed/SRR34987613_2_unpaired.fastq.gz \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
SLIDINGWINDOW:4:20 \
MINLEN:36
echo 'SRR34987613 trimming completed!'

# ---------- 4) Quality Control (after trim)  ----------
mkdir -p ~/rna-seq-analysis/results/fastqc/after-trim
fastqc ~/rna-seq-analysis/data/trimmed/*_paired.fastq.gz -o ~/rna-seq-analysis/results/fastqc/after-trim

# ---------- 5) Alignment (HISAT2)  ----------
# Download HISAT2 index file
mkdir -p ~/rna-seq-analysis/data/ref
wget -P ~/rna-seq-analysis/data/ref https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
tar -xvzf ~/rna-seq-analysis/data/ref/grch38_genome.tar.gz -C ~/rna-seq-analysis/data/ref
mkdir -p ~/rna-seq-analysis/results/alignment
hisat2 -p 8 -x ~/rna-seq-analysis/data/ref/grch38/genome \
-1 ~/rna-seq-analysis/data/trimmed/SRR34987610_1_paired.fastq.gz \
-2 ~/rna-seq-analysis/data/trimmed/SRR34987610_2_paired.fastq.gz  \
| samtools sort -@ 8 -o ~/rna-seq-analysis/results/alignment/SRR34987610.bam
echo 'SRR34987610 alignment complete!'
hisat2 -p 8 -x ~/rna-seq-analysis/data/ref/grch38/genome \
-1 ~/rna-seq-analysis/data/trimmed/SRR34987611_1_paired.fastq.gz \
-2 ~/rna-seq-analysis/data/trimmed/SRR34987611_2_paired.fastq.gz  \
| samtools sort -@ 8 -o ~/rna-seq-analysis/results/alignment/SRR34987611.bam
echo 'SRR34987611 alignment complete!'
hisat2 -p 8 -x ~/rna-seq-analysis/data/ref/grch38/genome \
-1 ~/rna-seq-analysis/data/trimmed/SRR34987612_1_paired.fastq.gz \
-2 ~/rna-seq-analysis/data/trimmed/SRR34987612_2_paired.fastq.gz  \
| samtools sort -@ 8 -o ~/rna-seq-analysis/results/alignment/SRR34987612.bam
echo 'SRR34987612 alignment complete!'
hisat2 -p 8 -x ~/rna-seq-analysis/data/ref/grch38/genome \
-1 ~/rna-seq-analysis/data/trimmed/SRR34987613_1_paired.fastq.gz \
-2 ~/rna-seq-analysis/data/trimmed/SRR34987613_2_paired.fastq.gz  \
| samtools sort -@ 8 -o ~/rna-seq-analysis/results/alignment/SRR34987613.bam
echo 'SRR34987613 alignment complete!'

# ---------- 6) Quality Control (MultiQC)  ----------
# Create logs for MultiQC
samtools flagstat ~/rna-seq-analysis/results/alignment/SRR34987610.bam > ~/rna-seq-analysis/results/alignment/SRR34987610.flagstat.txt
samtools flagstat ~/rna-seq-analysis/results/alignment/SRR34987611.bam > ~/rna-seq-analysis/results/alignment/SRR34987611.flagstat.txt
samtools flagstat ~/rna-seq-analysis/results/alignment/SRR34987612.bam > ~/rna-seq-analysis/results/alignment/SRR34987612.flagstat.txt
samtools flagstat ~/rna-seq-analysis/results/alignment/SRR34987613.bam > ~/rna-seq-analysis/results/alignment/SRR34987613.flagstat.txt
mkdir -p ~/rna-seq-analysis/results/multiqc/
multiqc ~/rna-seq-analysis/ \
  -o ~/rna-seq-analysis/results/multiqc/

# ---------- 7) Gene Count Matrix (featureCounts)  ----------
# Download gene annotation
mkdir -p ~/rna-seq-analysis/results/counts
mkdir -p ~/rna-seq-analysis/data/annot
wget -P ~/rna-seq-analysis/data/annot https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz
gunzip -c ~/rna-seq-analysis/data/annot/Homo_sapiens.GRCh38.115.gtf.gz > ~/rna-seq-analysis/data/annot/Homo_sapiens.GRCh38.115.gtf
featureCounts -T 8 -p --countReadPairs -s 2 \
-t exon -g gene_id \
-a ~/rna-seq-analysis/data/annot/Homo_sapiens.GRCh38.115.gtf \
-o ~/rna-seq-analysis/results/counts/counts_s2.txt \
~/rna-seq-analysis/results/alignment/*.bam

echo "$((SECONDS/60)) minutes and $((SECONDS%60)) seconds elapsed."
COMMENT