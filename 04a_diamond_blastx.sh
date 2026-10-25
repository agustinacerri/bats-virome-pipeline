#!/usr/bin/env bash
# Step 4a - Protein-level classification of the non-redundant contigs:
# DIAMOND blastx against nr (tabular output), followed by lineage annotation with TaxonKit.
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/config.sh"

mkdir -p "$DIAMOND_DIR"

# Column names of the tabular output (must match the --outfmt 6 fields below)
HEADER_BASE="qseqid\tsseqid\tpident\tbitscore\tevalue\tqframe\tqstart\tqend\tsstart\tsend\tqlen\tslen\tlength\tqcovhsp\tscovhsp\tstitle\tsstaxids\tsscinames"
HEADER_TAXON="${HEADER_BASE}\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies"

for sample in "${SAMPLES[@]}"; do
    echo "[DIAMOND blastx] $sample"
    tab="$DIAMOND_DIR/${sample}_diamond_tab.tsv"

    diamond blastx --very-sensitive -b 4 -d "$NR_DMND" --threads 20 \
        -q "$CDHIT_DIR/${sample}.fasta" -o "$tab" \
        -f 6 qseqid sseqid pident bitscore evalue qframe qstart qend sstart send \
             qlen slen length qcovhsp scovhsp stitle staxids sscinames

    # Add header
    sed -i "1i ${HEADER_BASE}" "$tab"

    # Add lineage (kingdom -> species) from the subject taxid (column 17)
    echo "[TaxonKit] $sample"
    taxonkit reformat "$tab" --data-dir "$TAXONKIT_DATA_DIR" \
        -f "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}" -I 17 \
        > "$DIAMOND_DIR/${sample}_diamond_tab_taxonomy.tsv"

    # Replace TaxonKit's header line with the full column names
    sed -i "1s|.*|${HEADER_TAXON}|" "$DIAMOND_DIR/${sample}_diamond_tab_taxonomy.tsv"
done
