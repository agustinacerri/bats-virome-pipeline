# Bat virome pipeline

Scripts used for the viral metagenomic analysis of bat samples: read trimming, de novo assembly,
redundancy removal, taxonomic classification of contigs (DIAMOND blastx, BLASTn) with TaxonKit lineages, read mapping,
and filtering of the hits to obtain a final table of viral contigs.

> Raw data and results are **not** included in this repository.

## Associated publication

This pipeline was used in the following preprint:

> Cerri A., Morro L., Zorec T. M., Montani M. E., Di Domenica V., Poljak M., Bolatti E. M., Giri A. A.
> *Virome Analysis of South American Arthropodophagous Bat Species With Sequence- and Profile-Based Homology Approaches Reveals Hidden RNA Virus Diversity.*
> Preprints.org 2026, 202610.1900 (v1). https://www.preprints.org/manuscript/202610.1900
> DOI: [10.20944/preprints202610.1900.v1](https://doi.org/10.20944/preprints202610.1900.v1)

This repository contains the sequence-similarity workflow (trimming, assembly, BLASTn/DIAMOND classification and filtering). The profile-based (HMMER) searches, the viral dark matter analysis and the community-composition analyses described in the preprint are not included.

## Workflow

| Step | Script | Tools | Output |
|------|--------|-------|--------|
| 1 | `01_qc_trimming.sh` | BBDuk (BBTools) | trimmed reads + trimming stats |
| 2 | `02_assembly.sh` | MEGAHIT (4 k-mer sets), metaSPAdes, metaviralSPAdes, rnaviralSPAdes, seqtk | contigs >= 500 bp per assembler |
| 3 | `03_clustering.sh` | seqkit, CD-HIT-EST | non-redundant contigs per sample |
| 4a | `04a_diamond_blastx.sh` | DIAMOND, TaxonKit | tabular blastx hits vs nr + lineage |
| 4b | `04b_blastn.sh` | BLAST+, TaxonKit | tabular blastn hits vs core_nt + lineage |
| 5 | `05_mapping.sh` | Bowtie2, samtools | sorted BAM, per-contig read counts |
| 6 | `06_filter_and_classify.sh` (runs `06_filter_and_classify.py`) | Python, pandas | final table of viral contigs per sample + summary counts |

## Usage

1. Install dependencies: `conda env create -f environment.yml`.
2. Edit **`config.sh`**: input/output paths and database locations. You do not need to run this file: each numbered script loads it automatically (`source config.sh`), so just keep it in the same folder as the scripts.
3. List your samples in `bat_samples` (one name per line; reads named `<sample>_R1.fastq.gz` / `<sample>_R2.fastq.gz`).
4. Run the scripts in order:

```bash
bash 01_qc_trimming.sh
bash 02_assembly.sh
bash 03_clustering.sh
bash 04a_diamond_blastx.sh
bash 04b_blastn.sh
bash 05_mapping.sh
bash 06_filter_and_classify.sh
```

## Step 6: how viral contigs are selected and classified

BLASTn is used first; DIAMOND blastx is used only for the contigs that BLASTn did not assign.

**BLASTn (core_nt)**

1. Keep hits with query coverage `qcovus >= 50`.
2. Select one hit per contig: highest `qcovus`, then highest `pident`.
3. Keep contigs whose best hit is viral and report the lineage according to the evidence:
   - hits of the contig (any taxon) belong to more than one phylum: lineage left empty;
   - `pident >= 80` and `qcovus >= 80`: full lineage (phylum to species);
   - otherwise: lineage down to phylum only.

**DIAMOND blastx (nr)**

1. Discard contigs already assigned by BLASTn.
2. Keep viral hits with `qcovhsp > 80` or `scovhsp > 80` (protein fragment or whole protein), and select the hit with the highest coverage per contig.
3. Report the phylum only (protein identities are low), and only if all the viral hits of the contig agree on it.

**Realm.** It is inferred from the phylum with the rules in `PHYLA_BY_REALM`. For DIAMOND hits without a single phylum, contigs are assigned to *Riboviria* when any of their phyla, or their family, is a known RNA virus taxon. Contigs left without a realm should be reviewed manually.

**Outputs** (`06_classification/`): `<sample>_viral_classification.tsv` (BLASTn and DIAMOND results together; the `source` column tells which method classified each contig) and `summary_counts.tsv` (contig counts at each filtering step). Thresholds are defined at the top of `06_filter_and_classify.py`.

## Databases required

- NCBI nr formatted for DIAMOND (`diamond makedb`, with taxonomy files)
- NCBI core_nt (BLAST database, downloaded with `update_blastdb.pl`)
- NCBI taxdump for TaxonKit (`~/.taxonkit` by default)
