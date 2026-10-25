#!/usr/bin/env bash
# Step 5 - Map the (trimmed) reads of each sample back to its non-redundant contigs
# with Bowtie2 and compute per-contig read counts (samtools idxstats).
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/config.sh"

for sample in "${SAMPLES[@]}"; do
    echo "[mapping] $sample"
    sample_dir="$MAP_DIR/$sample"
    mkdir -p "$sample_dir"
    index="$sample_dir/${sample}_index"
    bam="$sample_dir/${sample}.bowtie2.sorted.bam"

    bowtie2-build "$CDHIT_DIR/${sample}.fasta" "$index" > /dev/null

    # Align and sort in one stream (avoids writing a large intermediate SAM file).
    # Bowtie2's alignment summary is saved in a log file.
    bowtie2 -p 20 -x "$index" --very-sensitive \
        -1 "$TRIMMED_DIR/${sample}_R1.fastq.gz" -2 "$TRIMMED_DIR/${sample}_R2.fastq.gz" \
        2> "$sample_dir/${sample}.bowtie2.log" \
        | samtools sort -@ 20 -o "$bam" -

    samtools index "$bam"
    samtools idxstats "$bam" > "$sample_dir/${sample}.idxstats.txt"

    # Remove the Bowtie2 index files
    rm -f "${index}".*.bt2 "${index}".*.bt2l
done
