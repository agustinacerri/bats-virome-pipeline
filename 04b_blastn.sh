#!/usr/bin/env bash
# Step 4b - Nucleotide-level classification of the non-redundant contigs:
# BLASTn against NCBI core_nt (tabular output), followed by lineage annotation with TaxonKit.
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/config.sh"

mkdir -p "$BLASTN_DIR"

# Column names of the tabular output (must match the -outfmt 6 fields below)
HEADER_BASE="qseqid\tsseqid\tpident\tbitscore\tevalue\tqstart\tqend\tlength\tslen\tqlen\tqcovhsp\tqcovus\tstitle\tstaxid"
HEADER_TAXON="${HEADER_BASE}\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies"

for sample in "${SAMPLES[@]}"; do
    echo "[BLASTn] $sample"
    tab="$BLASTN_DIR/${sample}_blastn_tab"

    blastn -task blastn -evalue 0.001 -max_target_seqs 20 -num_threads 20 \
        -query "$CDHIT_DIR/${sample}.fasta" -db "$CORE_NT_DB" \
        -outfmt '6 qseqid sseqid pident bitscore evalue qstart qend length slen qlen qcovhsp qcovus stitle staxid' \
        -out "$tab"

    # Add header
    sed -i "1i ${HEADER_BASE}" "$tab"

    # Add lineage (kingdom -> species) from the subject taxid (column 14)
    echo "[TaxonKit] $sample"
    taxonkit reformat "$tab" --data-dir "$TAXONKIT_DATA_DIR" \
        -f "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}" -I 14 \
        > "${tab}_taxonomy"

    # Replace TaxonKit's header line with the full column names
    sed -i "1s|.*|${HEADER_TAXON}|" "${tab}_taxonomy"
done
