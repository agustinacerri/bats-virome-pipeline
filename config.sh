#!/usr/bin/env bash
# =============================================================================
# config.sh - Central configuration for the bat virome pipeline.
#Edit the paths below; the numbered scripts load this file automatically.
#=============================================================================

CONFIG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

#------------------------------ Input / output --------------------------------
#Directory with raw paired-end reads named <sample>_R1.fastq.gz / <sample>_R2.fastq.gz
RAW_DIR="/path/to/raw_reads"
#Output directory (all results are written here)
OUTDIR="/path/to/output"
#Text file with one sample name per line (blank lines and lines starting with # are ignored)
SAMPLES_FILE="$CONFIG_DIR/bat_samples"

#-------------------------------- Databases -----------------------------------
#Folder containing all databases. If the environment variable $db is set, it is used.
DB_DIR="${db:-/path/to/databases}"
NR_DMND="$DB_DIR/nr.dmnd"                                   #DIAMOND-formatted NCBI nr
CORE_NT_DB="$DB_DIR/core_nt/core_nt"                        #BLAST core_nt database (prefix)
TAXONKIT_DATA_DIR="${TAXONKIT_DATA_DIR:-$HOME/.taxonkit}"   #NCBI taxdump files for TaxonKit

#------------------------------- Output folders -------------------------------
TRIMMED_DIR="$OUTDIR/01_trimmed"
BBDUK_STATS_DIR="$OUTDIR/01_bbduk_stats"
ASSEMBLY_DIR="$OUTDIR/02_assemblies"
POOLED_DIR="$OUTDIR/03_contigs_pooled"
CDHIT_DIR="$OUTDIR/03_cdhit"
DIAMOND_DIR="$OUTDIR/04a_diamond_blastx"
BLASTN_DIR="$OUTDIR/04b_blastn"
MAP_DIR="$OUTDIR/05_mapping"
CLASSIFY_DIR="$OUTDIR/06_classification"

#--------------------------------- Safeguards ---------------------------------
if [[ "$OUTDIR" == "/path/to/output" || "$RAW_DIR" == "/path/to/raw_reads" ]]; then
    echo "ERROR: please edit config.sh (RAW_DIR, OUTDIR, database paths) before running." >&2
    exit 1
fi

#Read sample names (ignores blank lines, comments and Windows carriage returns)
mapfile -t SAMPLES < <(grep -Ev '^[[:space:]]*(#|$)' "$SAMPLES_FILE" | tr -d '\r')
if [[ ${#SAMPLES[@]} -eq 0 ]]; then
    echo "ERROR: no samples found in $SAMPLES_FILE" >&2
    exit 1
fi
