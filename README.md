# Schizochytrium sp. HS6 — Comparative Genomics

Comparative genomics, transcriptomics, and synthetic-biology part discovery across three *Schizochytrium* strains:

- **HS6** — EMS-mutagenized, carotenoid-inhibitor-selected mutant bred for higher DHA/lipid yield and faster growth
- **#4 (WT)** — HS6's unmutagenized parental strain
- **#7** — an independently isolated *Schizochytrium* strain, used as an outgroup comparator

**Live site:** https://jinho1987.github.io/schizochytrium-hs6-genomics/

## Contents

| Page | Description |
|---|---|
| [`figure1.html`](figure1.html) | HS6 vs. #7 — lipid, isoprenoid & respiratory-chain metabolic map (journal-style figure) |
| [`figure2.html`](figure2.html) | #4 (WT) vs. HS6 — expression-integrated regulatory & metabolic map |
| [`report1.html`](report1.html) | HS6 vs. #7 comparative genomics report (genome relatedness, gene content, the MVK duplication finding) |
| [`report2.html`](report2.html) | HS6 vs. #4 RNA-seq findings report (the original mutation-hunt + regulatory-rewiring story) |
| [`report3.html`](report3.html) | HS6 vs. #4 analytical reproducibility log (methods, corrections, QC) |
| [`report4.html`](report4.html) | **Extended functional genomics**: promoters, introns/UTRs for synthetic biology, transposable elements, RNA-seq scaffolding, independent re-annotation, expression-constrained flux model |

## Data & reproducibility

This repo contains the **figures, reports, analysis scripts, and small result tables** — not the raw sequencing data or full genome assemblies (too large for git; available on request). Every result in `report4.html` was produced by a script in [`scripts/`](scripts/) against:

- `HS6_genome.fasta` / `strain7_genome.fasta` — PacBio HiFi de novo assemblies
- `braker.gff3` — BRAKER2/AUGUSTUS gene models
- `merged.stringtie.gtf` — RNA-seq-guided UTR/isoform refinement (StringTie, merged across all 6 libraries)
- 6 RNA-seq libraries (3 timepoints × 2 genotypes, n=1 each — **no biological replicates**, see the caveats in `report3.html`)
- WT4-vs-HS6 short-read resequencing variant calls (bcftools, haploid model)

Small result tables live in [`results/`](results/) (ranked candidate lists, classifications, overlay tables, sequences & primers). The combined differential-expression workbook — all 13,529 genes, genotype-effect stats + per-timepoint log2FC + normalized counts + annotation in one sheet, color-coded, with a sorted significant-only sheet and GO enrichment tables — is [`results/HS6_vs_WT4_DEG_combined.xlsx`](results/HS6_vs_WT4_DEG_combined.xlsx). Conda/pip environment setup for every tool used is documented at the bottom of `report4.html`, including several real tool-install bugs hit and fixed along the way (wrong CLI flags, dead download URLs, gene-ID format mismatches) — kept in the writeup so they don't cost the next person the same time.

## Site structure

The homepage (`index.md`, built with a custom Jekyll layout in [`_layouts/default.html`](_layouts/default.html)) groups everything by the two comparisons this repo actually contains — HS6 vs. #7 (genome) and HS6 vs. #4 (RNA-seq) — rather than a flat report list. A shared navigation bar ([`assets/css/site-nav.css`](assets/css/site-nav.css)) with the same grouping appears on every page so you can move between reports and figures without going back through the homepage.

## Statistical caveats that apply throughout

- **n=1 per genotype/timepoint** — no biological replicates in the RNA-seq. DESeq2 results are internally consistent but a single divergent culture could produce the same signal.
- Findings are explicitly tiered by evidence strength: genome-sequence-confirmed (e.g. the MVK copy-number duplication, independently verified from raw HiFi reads on **both** strains) vs. transcriptome-correlational (e.g. the MYB3R timing argument) vs. purely computational/untested (e.g. the flux model, the IME candidates).
- Where an earlier claim in this project turned out to be wrong (a buggy `--max-target-seqs` diamond search, a blank KEGG pathway map from a pathview rendering bug, an "HS6-specific" gene cluster that was actually a KO-assignment artifact), the correction is documented in place rather than silently fixed — see the correction banners in `report1.html`/`report2.html`.
