#!/usr/bin/env bash
# Step 1 - Quality trimming and filtering of paired-end reads with BBDuk (BBTools).
#   qtrim=rl : trim low-quality bases from both read ends (Phred < 15)
#   maq      : discard reads with average quality below 20
#   minlen   : discard reads shorter than 120 bp after trimming
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/config.sh"

mkdir -p "$TRIMMED_DIR" "$BBDUK_STATS_DIR"

for sample in "${SAMPLES[@]}"; do
    echo "[bbduk] $sample"
    bbduk.sh -Xmx1g threads=16 \
        in1="$RAW_DIR/${sample}_R1.fastq.gz" in2="$RAW_DIR/${sample}_R2.fastq.gz" \
        out1="$TRIMMED_DIR/${sample}_R1.fastq.gz" out2="$TRIMMED_DIR/${sample}_R2.fastq.gz" \
        tpe qtrim=rl trimq=15 maq=20 minlen=120 \
        stats="$BBDUK_STATS_DIR/stats_${sample}.txt"
done
