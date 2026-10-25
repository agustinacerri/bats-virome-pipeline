#!/usr/bin/env python3
"""
Step 6 - Filter BLASTn and DIAMOND blastx tables, keep viral contigs and assign taxonomy.

For every sample:
  1. BLASTn (nucleotide). Hits with qcovus >= 50 are kept and the best hit per contig is
     selected (highest qcovus, then highest pident). Viral best hits are then classified:
       - hits of the contig (any taxon) belong to more than one phylum -> lineage left empty
       - pident >= 80 and qcovus >= 80                                 -> full lineage kept
       - otherwise                                                     -> lineage kept down to phylum
  2. DIAMOND blastx (protein). Contigs already assigned by BLASTn (any kingdom) are discarded.
     Viral hits with qcovhsp > 80 or scovhsp > 80 are kept, one hit per contig (highest
     coverage). Protein identities are low, so only the phylum is reported, and only when all
     the viral hits of the contig agree on it.
  3. Both tables are combined (column `source` tells which method classified each contig) and
     a per-sample summary of contig counts is written.

The viral realm is inferred from the phylum (see PHYLA_BY_REALM). Contigs left without a realm
should be reviewed manually.
"""

import argparse
from pathlib import Path

import pandas as pd

# ------------------------------- Thresholds ----------------------------------
BLASTN_MIN_QCOVUS = 50      # minimum query coverage (qcovus) to keep a BLASTn hit
BLASTN_HIGH_PIDENT = 80     # "high confidence" BLASTn hit: pident >= 80 ...
BLASTN_HIGH_QCOVUS = 80     # ... and qcovus >= 80
BLASTX_MIN_COVERAGE = 80    # DIAMOND hit kept if qcovhsp > 80 or scovhsp > 80

# ------------------------------ Realm assignment ------------------------------
PHYLA_BY_REALM = {
    "Adnaviria": ["taleaviricota"],
    "Duplodnaviria": ["peploviricota", "uroviricota"],
    "Monodnaviria": ["hofneiviricota", "phixviricota", "commensaviricota", "cossaviricota",
                     "cressdnaviricota", "calorviricota", "saleviricota"],
    "Riboviria": ["ambiviricota", "artimaviricota", "duplornaviricota", "kitrinoviricota",
                  "lenarviricota", "negarnaviricota", "pisuviricota", "artverviricota"],
    "Singelaviria": ["dividoviricota"],
    "Varidnaviria": ["produgelaviricota", "nucleocytoviricota", "preplasmiviricota"],
}
REALM_BY_PHYLUM = {phylum: realm for realm, phyla in PHYLA_BY_REALM.items() for phylum in phyla}

# RNA virus families that lack a phylum in the NCBI taxonomy (used for DIAMOND hits only)
RIBOVIRIA_FAMILIES_WITHOUT_PHYLUM = ["permutotetraviridae", "birnaviridae",
                                     "polymicoviridae", "tonesaviridae"]

# ------------------------------- Output columns -------------------------------
RANKS = ["phylum", "class", "order", "family", "genus", "species"]
BLASTN_COLUMNS = ["qseqid", "sseqid", "pident", "bitscore", "evalue", "qstart", "qend", "length",
                  "slen", "qlen", "qcovhsp", "qcovus", "stitle", "staxid", "kingdom", "realm"] + RANKS
BLASTX_COLUMNS = ["qseqid", "sseqid", "pident", "bitscore", "evalue", "qstart", "qend", "length",
                  "slen", "qlen", "qcovhsp", "scovhsp", "stitle", "staxid", "kingdom", "realm",
                  "phylum"]
COMBINED_COLUMNS = ["qseqid", "sseqid", "pident", "bitscore", "evalue", "qstart", "qend", "length",
                    "slen", "qlen", "qcovhsp", "qcovus", "scovhsp", "stitle", "staxid", "kingdom",
                    "realm", "phylum", "class", "order", "family", "genus", "species", "source"]


def any_column_contains(df, text):
    """Boolean mask of the rows in which any column contains `text` (case-insensitive)."""
    as_text = df.astype(str)
    return as_text.apply(lambda col: col.str.contains(text, case=False, regex=False)).any(axis=1)


def collapse_unique(df, key, column):
    """Per `key`, the sorted unique non-null values of `column` joined with ';'."""
    return df.groupby(key)[column].apply(lambda x: ";".join(sorted(set(x.dropna().astype(str)))))


def single_value(joined):
    """Return `joined` if it holds exactly one value, otherwise an empty string."""
    return joined if isinstance(joined, str) and joined and ";" not in joined else ""


def assign_realm(species, phylum, all_phyla="", family="", use_fallback=False):
    """Infer the viral realm from a lineage. Rules are applied in priority order.

    With `use_fallback`, contigs without a single phylum are still assigned to Riboviria when
    any of their phyla (or their family) is a known RNA virus taxon.
    """
    species = str(species).lower()
    phylum = "" if pd.isna(phylum) else str(phylum).lower()

    if "riboviria sp" in species:
        return "Riboviria"
    if phylum in REALM_BY_PHYLUM:
        return REALM_BY_PHYLUM[phylum]
    if use_fallback and phylum == "":
        phyla = str(all_phyla).lower().split(";")
        if any(p in phyla for p in PHYLA_BY_REALM["Riboviria"]):
            return "Riboviria"
        if any(f in str(family).lower() for f in RIBOVIRIA_FAMILIES_WITHOUT_PHYLUM):
            return "Riboviria"
    return ""


def classify_blastn(blastn):
    """Best BLASTn hit per contig, restricted to viral contigs and classified by confidence.

    Returns (viral table, IDs of all contigs assigned by BLASTn, summary counts).
    """
    # Contigs whose best hit (highest qcovus, no coverage filter) mentions "Virus" in any column
    best_unfiltered = blastn.loc[blastn.groupby("qseqid")["qcovus"].idxmax()]
    n_viral_unfiltered = int(any_column_contains(best_unfiltered, "Virus").sum())

    # Best hit per contig among hits with enough coverage, keeping the original contig order
    passing = blastn[blastn["qcovus"] >= BLASTN_MIN_QCOVUS]
    best = (passing.sort_values(["qcovus", "pident"], ascending=False)
                   .drop_duplicates("qseqid")
                   .sort_index()
                   .copy())

    # All phyla/classes/orders found among the hits of each contig (used to flag ambiguity)
    for rank in ("phylum", "class", "order"):
        best[f"all_{rank}"] = best["qseqid"].map(collapse_unique(passing, "qseqid", rank))

    viral = best[any_column_contains(best, "Viruses")].copy()
    viral["realm"] = [assign_realm(sp, ph) for sp, ph in zip(viral["species"], viral["phylum"])]
    viral = viral.sort_values("qcovus", ascending=False)

    # Report only the taxonomic ranks that the evidence supports
    ambiguous = viral["all_phylum"].astype(str).str.contains(";", regex=False)
    high_confidence = ((viral["pident"] >= BLASTN_HIGH_PIDENT)
                       & (viral["qcovus"] >= BLASTN_HIGH_QCOVUS))
    viral[RANKS] = viral[RANKS].astype(object)
    viral.loc[ambiguous, RANKS] = ""
    viral.loc[~(ambiguous | high_confidence), RANKS[1:]] = ""
    viral["source"] = "blastn"

    counts = {
        "blastn_best_hit_viral_unfiltered": n_viral_unfiltered,
        "blastn_contigs_assigned": len(best),
        "blastn_viral_contigs": len(viral),
    }
    return viral[BLASTN_COLUMNS + ["source"]], set(best["qseqid"]), counts


def classify_blastx(blastx, excluded_ids):
    """Viral DIAMOND hits with enough coverage, for contigs not assigned by BLASTn.

    Returns (viral table, summary counts).
    """
    remaining = blastx[~blastx["qseqid"].isin(excluded_ids)]

    covered = (remaining["qcovhsp"] > BLASTX_MIN_COVERAGE) | (remaining["scovhsp"] > BLASTX_MIN_COVERAGE)
    is_viral = remaining["kingdom"].str.contains("Viruses", case=False, na=False)
    viral_hits = remaining[covered & is_viral].copy()

    # One hit per contig: the one with the highest coverage (protein fragment or whole protein)
    viral_hits["max_cov"] = viral_hits[["qcovhsp", "scovhsp"]].max(axis=1)
    best = viral_hits.loc[viral_hits.groupby("qseqid")["max_cov"].idxmax()].copy()

    # Taxonomy is trusted only down to phylum, and only if all viral hits agree on it
    for rank in ("phylum", "class", "order"):
        best[f"all_{rank}"] = best["qseqid"].map(collapse_unique(viral_hits, "qseqid", rank))
    best["final_phylum"] = best["all_phylum"].map(single_value)
    best["realm"] = [
        assign_realm(sp, ph, all_ph, fam, use_fallback=True)
        for sp, ph, all_ph, fam in zip(best["species"], best["final_phylum"],
                                       best["all_phylum"], best["family"])
    ]

    final = best.rename(columns={"sstaxids": "staxid"}).assign(phylum=best["final_phylum"])
    final = final[BLASTX_COLUMNS].assign(source="blastx")

    counts = {
        "blastx_contigs_with_hits": remaining["qseqid"].nunique(),
        "blastx_contigs_covered": remaining.loc[covered, "qseqid"].nunique(),
        "blastx_viral_contigs": len(final),
    }
    return final, counts


def combine(blastn_viral, blastx_viral):
    """Stack the BLASTn and DIAMOND tables (a contig is never present in both)."""
    combined = pd.concat([blastn_viral, blastx_viral], ignore_index=True)
    return combined[COMBINED_COLUMNS].sort_values("qseqid", kind="stable").reset_index(drop=True)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    parser.add_argument("--blastn-dir", required=True, type=Path,
                        help="folder with <sample>_blastn_tab_taxonomy files (step 4b)")
    parser.add_argument("--blastx-dir", required=True, type=Path,
                        help="folder with <sample>_diamond_tab_taxonomy.tsv files (step 4a)")
    parser.add_argument("--outdir", required=True, type=Path, help="output folder")
    parser.add_argument("--samples", required=True, nargs="+", help="sample names")
    return parser.parse_args()


def main():
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    summary = []
    for sample in args.samples:
        print(f"[classify] {sample}")
        blastn = pd.read_csv(args.blastn_dir / f"{sample}_blastn_tab_taxonomy", sep="\t",
                             low_memory=False)
        blastx = pd.read_csv(args.blastx_dir / f"{sample}_diamond_tab_taxonomy.tsv", sep="\t",
                             low_memory=False)

        blastn_viral, assigned_ids, blastn_counts = classify_blastn(blastn)
        blastx_viral, blastx_counts = classify_blastx(blastx, assigned_ids)
        combined = combine(blastn_viral, blastx_viral)

        combined.to_csv(args.outdir / f"{sample}_viral_classification.tsv", sep="\t", index=False)
        summary.append({"sample": sample, **blastn_counts, **blastx_counts,
                        "total_viral_contigs": len(combined)})

    summary = pd.DataFrame(summary)
    summary.to_csv(args.outdir / "summary_counts.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
