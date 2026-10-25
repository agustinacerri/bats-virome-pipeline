#!/usr/bin/env bash
# Step 2 - De novo assembly of each sample with several assemblers/parameter sets:
#   4 x MEGAHIT (different k-mer ranges) + metaSPAdes + metaviralSPAdes + rnaviralSPAdes.
# Contigs shorter than 500 bp are discarded with seqtk.
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/config.sh"

# MEGAHIT parameter sets: run n uses k-min = KMIN[n], k-step = KSTEP[n]
MEGAHIT_KMIN=(21 23 25 27)
MEGAHIT_KSTEP=(18 20 22 24)
MEGAHIT_KMAX=149

# SPAdes runs: output name and the corresponding mode flag
SPADES_NAMES=(metaspades metaviralspades rnaviralspades)
SPADES_FLAGS=(--meta --metaviral --rnaviral)
SPADES_KMERS="21,33,55,77,89,99"
# Resources used: MEGAHIT -t 125 threads, -m in bytes (~250 GB); SPAdes -t 60 threads, -m in GB. Adjust to your machine.

for sample in "${SAMPLES[@]}"; do
    r1="$TRIMMED_DIR/${sample}_R1.fastq.gz"
    r2="$TRIMMED_DIR/${sample}_R2.fastq.gz"
    sample_dir="$ASSEMBLY_DIR/$sample"
    tmp="$sample_dir/tmp_assembly"
    mkdir -p "$sample_dir"

    # ---- MEGAHIT ----
    for idx in "${!MEGAHIT_KMIN[@]}"; do
        n=$((idx + 1))
        out="$sample_dir/megahit${n}.fasta"
        if [[ -s "$out" ]]; then echo "[skip] $out exists"; continue; fi

        echo "[MEGAHIT $n] $sample"
        rm -rf "$tmp"   # MEGAHIT refuses to write into an existing directory
        megahit \
            --k-min "${MEGAHIT_KMIN[$idx]}" --k-max "$MEGAHIT_KMAX" --k-step "${MEGAHIT_KSTEP[$idx]}" \
            --min-count 1 -t 125 -m 250000000000 \
            -1 "$r1" -2 "$r2" -o "$tmp"
        seqtk seq -L 500 "$tmp/final.contigs.fa" > "$out"
        rm -rf "$tmp"
    done

    # ---- SPAdes (meta / metaviral / rnaviral) ----
    for idx in "${!SPADES_NAMES[@]}"; do
        name="${SPADES_NAMES[$idx]}"
        out="$sample_dir/${name}.fasta"
        if [[ -s "$out" ]]; then echo "[skip] $out exists"; continue; fi

        echo "[$name] $sample"
        rm -rf "$tmp"
        spades.py -1 "$r1" -2 "$r2" -o "$tmp" \
            -k "$SPADES_KMERS" -t 60 -m 120 \
            "${SPADES_FLAGS[$idx]}"
        seqtk seq -L 500 "$tmp/contigs.fasta" > "$out"
        rm -rf "$tmp"
    done
done
