#!/usr/bin/env bash
# Step 6 - Filter the BLASTn / DIAMOND tables, keep viral contigs, assign taxonomy and
# write the final classification table and contig counts (see 06_filter_and_classify.py).
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

python3 "$SCRIPT_DIR/06_filter_and_classify.py" \
    --blastn-dir "$BLASTN_DIR" \
    --blastx-dir "$DIAMOND_DIR" \
    --outdir "$CLASSIFY_DIR" \
    --samples "${SAMPLES[@]}"
