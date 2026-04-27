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

prefetch "$READ1_ACC" --output-directory "$RAW_DATA_DIR" > "$LOGS_DIR/01_download_read1.log" 2>&1
fasterq-dump "$RAW_DATA_DIR/$READ1_ACC/" \
 --split-files \
 -O "$RAW_DATA_DIR" \
 --threads "$THREADS" \
 > "$LOGS_DIR/02_split_file.log" 2>&1
gzip "$RAW_DATA_DIR/$READ1_ACC"*.fastq

rm -rf "$RAW_DATA_DIR/$READ1_ACC/"

log "Read 1 completed"

prefetch "$READ2_ACC" --output-directory "$RAW_DATA_DIR" > "$LOGS_DIR/01_download_read2.log" 2>&1
fasterq-dump "$RAW_DATA_DIR/$READ2_ACC/" \
 --split-files \
 -O "$RAW_DATA_DIR" \
 --threads "$THREADS" \
 > "$LOGS_DIR/02_read2_split_file.log" 2>&1
gzip "$RAW_DATA_DIR/$READ2_ACC"*.fastq

rm -rf "$RAW_DATA_DIR/$READ2_ACC"

log "Read 2 complete"

# NUDT21 Knockdown files
prefetch "$READ3_ACC" --output-directory "$RAW_DATA_DIR" > "$LOGS_DIR/01_download_read3.log" 2>&1
fasterq-dump "$RAW_DATA_DIR/$READ3_ACC/" \
--split-files \
-O "$RAW_DATA_DIR" \
--threads "$THREADS" \
> "$LOGS_DIR/02_read3_split_file.log" 2>&1
gzip "$RAW_DATA_DIR/$READ3_ACC"*.fastq

rm -rf "$RAW_DATA_DIR/$READ3_ACC"

log "Read 3 complete"

prefetch "$READ4_ACC" --output-directory "$RAW_DATA_DIR" > "$LOGS_DIR/01_download_read4.log" 2>&1
fasterq-dump "$RAW_DATA_DIR/$READ4_ACC/" \
--split-files \
-O "$RAW_DATA_DIR" \
--threads "$THREADS" \
> "$LOGS_DIR/02_read4_split_file.log" 2>&1
gzip "$RAW_DATA_DIR/$READ4_ACC"*.fastq

rm -rf "$RAW_DATA_DIR/$READ4_ACC"

log "Read 4 complete"

####################################
# 2. FastQC (raw_data)
####################################

log "Starting FastQC ..."

fastqc "$READ1_R1" "$READ1_R2" "$READ2_R1" "$READ2_R2" "$READ3_R1" "$READ3_R2" "$READ4_R1" "$READ4_R2" \
-o "$FASTQC_RAW_DIR" \
> "$LOGS_DIR/03_fastqc_raw.log" 2>&1

log "FastQC completed"

####################################
# 3. Trimmomatic
####################################

log "Starting trimming ..."

trimmomatic PE \
-threads "$THREADS" \
-phred33 "$READ1_R1" "$READ1_R2" \
"$TRIMMED_DIR/$READ1_ACC"_1_paired.fastq.gz "$TRIMMED_DIR/$READ1_ACC"_1_unpaired.fastq.gz \
"$TRIMMED_DIR/$READ1_ACC"_2_paired.fastq.gz "$TRIMMED_DIR/$READ1_ACC"_2_unpaired.fastq.gz \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
SLIDINGWINDOW:4:20 \
MINLEN:36 \
> "$LOGS_DIR/04_trimming_read1.log" 2>&1

log "Completed trimming on read 1"

trimmomatic PE \
-threads "$THREADS" \
-phred33 "$READ2_R1" "$READ2_R2" \
"$TRIMMED_DIR/$READ2_ACC"_1_paired.fastq.gz "$TRIMMED_DIR/$READ2_ACC"_1_unpaired.fastq.gz \
"$TRIMMED_DIR/$READ2_ACC"_2_paired.fastq.gz "$TRIMMED_DIR/$READ2_ACC"_2_unpaired.fastq.gz \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
SLIDINGWINDOW:4:20 \
MINLEN:36 \
> "$LOGS_DIR/04_trimming_read2.log" 2>&1

log "Completed trimming on read 2"

trimmomatic PE \
-threads "$THREADS" \
-phred33 "$READ3_R1" "$READ3_R2" \
"$TRIMMED_DIR/$READ3_ACC"_1_paired.fastq.gz "$TRIMMED_DIR/$READ3_ACC"_1_unpaired.fastq.gz \
"$TRIMMED_DIR/$READ3_ACC"_2_paired.fastq.gz "$TRIMMED_DIR/$READ3_ACC"_2_unpaired.fastq.gz \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
SLIDINGWINDOW:4:20 \
MINLEN:36 \
> "$LOGS_DIR/04_trimming_read3.log" 2>&1

log "Completed trimming on read 3"

trimmomatic PE \
-threads "$THREADS" \
-phred33 "$READ4_R1" "$READ4_R2" \
"$TRIMMED_DIR/$READ4_ACC"_1_paired.fastq.gz "$TRIMMED_DIR/$READ4_ACC"_1_unpaired.fastq.gz \
"$TRIMMED_DIR/$READ4_ACC"_2_paired.fastq.gz "$TRIMMED_DIR/$READ4_ACC"_2_unpaired.fastq.gz \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
SLIDINGWINDOW:4:20 \
MINLEN:36 \
> "$LOGS_DIR/04_trimming_read4.log" 2>&1

log "Completed trimming on read 4"


####################################
# 4. FastQC (trimmed_data)
####################################

log "Starting FastQC on trimmed reads ..."

fastqc "$TRIMMED_DIR/"*_paired.fastq.gz -o "$FASTQC_TRIMMED_DIR" > "$LOGS_DIR/05_fastqc_trimmed.log" 2>&1

log "FastQC on trimmed reads completed"


####################################
# 5. Alignment (HISAT2)
####################################

# Download HISAT2 index file

log "Downloading HISAT2 index file ..."

wget -O "$HISAT2_INDEX" https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz > "$LOGS_DIR/06_hisat_index.log" 2>&1
tar -xvzf "$HISAT2_INDEX" -C "$REF_DIR"

rm -rf "$HISAT2_INDEX"

log "Completed index file"

log "Starting alignment ..."

PAIRS=( 
  "$READ1_R1_PAIRED|$READ1_R2_PAIRED|READ1" 
  "$READ2_R1_PAIRED|$READ2_R2_PAIRED|READ2" 
  "$READ3_R1_PAIRED|$READ3_R2_PAIRED|READ3"
  "$READ4_R1_PAIRED|$READ4_R2_PAIRED|READ4"
  )

for PAIR in "${PAIRS[@]}"
do
    IFS="|" read -r R1 R2 SAMPLE <<< "$PAIR"

    {
    hisat2 -x "$HISAT2_INDEX_DIR" \
        -1 "$R1" \
        -2 "$R2" \
        | samtools sort -@ 8 -o "$HISAT_ALIGN_DIR/${SAMPLE}.bam"
        
        samtools index "$HISAT_ALIGN_DIR/${SAMPLE}.bam" \
        
        echo "Completed $SAMPLE alignment"
    } > "$LOGS_DIR/07_${SAMPLE}_alignment.log" 2>&1

done

log "Completed alignment"

####################################
# 6. Quality Control (MultiQC)
####################################

log "Starting MultiQC ..."

multiqc "$PROJECT_ROOT" \
  -o "$MULTIQC_DIR" \
  > "$LOGS_DIR/08_multiqc.log" 2>&1

log "Completed MultiQC"

####################################
# 7. Gene Count Matrix (featureCounts)
####################################
# Download gene annotation
log "Downloading gene annotation ..."

wget -O "$COUNT_DIR/Homo_sapiens.GRCh38.115.gtf.gz" https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz \
> "$LOGS_DIR/09_gene_annotation.log" 2>&1
gunzip -c "$COUNT_DIR/"Homo_sapiens.GRCh38.115.gtf.gz > "$COUNT_DIR/"Homo_sapiens.GRCh38.115.gtf

log "Completed download"


log "Starting count matrix ..."

featureCounts -T "$THREADS" \
-p --countReadPairs \
-s 2 \
-t exon \
-g gene_id \
-a "$ANNOTATION_FILE" \
-o "$COUNT_DIR"/counts_s2.txt \
"$HISAT_ALIGN_DIR"/*.bam \
> "$LOGS_DIR/10_gene_count.log" 2>&1

log "Completed count matrix"