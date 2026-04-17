#!/usr/bin/env bash

# Gets project root automatically 
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

############################
# DIRECTORIES
############################

RAW_DATA_DIR="$PROJECT_ROOT/raw_data"
ANALYSIS_DIR="$PROJECT_ROOT/analysis"
LOGS_DIR="$PROJECT_ROOT/logs"


READS1="SRR34987610"
READS2="SRR34987611"
READS3="SRR34987612"
READS4="SRR34987613"

THREADS="8"