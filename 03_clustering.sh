#!/usr/bin/env bash
# Step 3 - Pool all assemblies of each sample and remove redundancy with CD-HIT-EST
# -T 0 = use all available cores
# 95% sequence identity, default word size; one representative sequence per cluster
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/config.sh"

mkdir -p "$POOLED_DIR" "$CDHIT_DIR"

for sample in "${SAMPLES[@]}"; do
    echo "[pool + cd-hit-est] $sample"

    # Concatenate the contigs of all assemblers; seqkit rename makes duplicated
    # contig IDs unique (different assemblers can produce identical headers).
    cat "$ASSEMBLY_DIR/$sample"/*.fasta | seqkit rename - > "$POOLED_DIR/${sample}.fasta"

    cd-hit-est \
        -i "$POOLED_DIR/${sample}.fasta" \
        -o "$CDHIT_DIR/${sample}.fasta" \
        -c 0.95 -T 0
done
