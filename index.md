---
layout: default
title: Home
---

# Schizochytrium HS6 — comparative genomics & transcriptomics

**HS6 (EMS mutant) vs. #4 (WT parent) vs. Schizochytrium sp. #7 (independent isolate) — differential expression, genome comparison, and synthetic-biology part discovery**

> ⚠️ **Statistical caveat: n=1 per genotype/timepoint.** The RNA-seq design is 2 genotype (WT #4, HS6) × 3 timepoint (20h/44h/68h fermentation), **no biological replicates**. DESeq2 results are internally consistent, but a single divergent culture could in principle produce the same signal. Findings below are explicitly tiered by evidence strength — see [Limitations](#limitations--evidence-tiers).

## Project summary

| | |
|---|---|
| Organism | *Schizochytrium* sp. — HS6 (EMS-mutagenized, carotenoid-inhibitor-selected), #4 (its unmutagenized WT parent), #7 (independently isolated outgroup strain) |
| Phenotype of interest | HS6 was bred for **higher DHA/lipid yield and faster growth** than #4, via EMS mutagenesis + carotenoid-synthesis inhibitor + Nile Red/OD screening |
| Growth medium | GYP medium |
| Design | 2 genotype (WT #4, HS6) × 3 timepoint (20h, 44h, 68h) × **n=1** (no biological replicates) |
| Sequencing | RNA-seq: Illumina, 6 libraries. Genome: PacBio HiFi de novo assembly (HS6, #7); Illumina short-read resequencing (#4/WT) |
| Reference genome | HS6 PacBio HiFi assembly, 34 scaffolds, 61.6 Mb, BUSCO 96.1% (stramenopiles_odb12, lineage-specific) |
| Gene set | 16,088 genes (BRAKER2/AUGUSTUS), 13,529 with RNA-seq count data, StringTie-refined UTR/isoform structure |

## Headline findings

- **A confirmed, genome-level structural difference between HS6 and #7**: HS6 carries an extra copy of **mevalonate kinase (*MVK*)** at a locus (Scaffolds_28) that sits at the head of the isoprenoid precursor pathway. This is independently confirmed from raw long-read depth on **both** strains — HS6's own HiFi reads show ~79.5× depth at that locus vs. a ~41× genome-wide baseline (a clean 2×); strain #7's own HiFi reads, mapped the same way, show 42.9× vs. a ~43× baseline (no duplication at all). This is the single most solid genomic finding in the whole project, and a leading candidate for HS6's growth edge over #7.
- **A coordinated regulatory shift, not one broken gene, separates HS6 from its own parent #4**: multiple paralogs of the **MYB3R cell-cycle transcription-factor family** are dysregulated (some up, some down) in HS6 vs. #4. This lines up with two visible phenotypic programs — HS6 puts more resources into DHA/lipid synthesis and less into cytoskeleton/growth machinery. The cytoskeleton signature is backed by strong formal statistics: GO term **actin filament binding** is enriched among genes down in HS6 at **p.adjust = 1.7×10⁻¹³** (30/291 genes).
- **Time-course evidence (suggestive, not proof) that MYB3R acts early**: most MYB3R-family genes are already substantially dysregulated at the earliest sampled timepoint (20h) and stay flat/sustained from there, while core lipid-turnover genes (*DGAT*, fatty-acyl-CoA ligases) show an escalating pattern that keeps building through 68h — an ordering more consistent with MYB3R as an upstream driver than a downstream consequence. See [Extended Analyses](report4.html).
- **Candidate driver mutations**: three point mutations sit near significantly DE genes — *g6423* (RanGAP-like nuclear-transport gene), *g1760* (G-type lectin receptor kinase), and *g6124* (ribosomal S6 kinase, a growth-signaling kinase, found via a genome-wide extension of the original mutation search). All three are leads, not confirmed causes.
- **A corrected gene-content comparison**: HS6 vs. #7 differ by 102 and 47 genes respectively (HS6-only / #7-only) — corrected down from an initially-reported 472/3,983 after catching a `--max-target-seqs` bug in the verification search. See the [correction notice](report1.html) in the main comparative genomics report.
- **A retraction, on the record**: an apparently HS6-specific respiratory-chain/ubiquinone gene cluster was flagged by an automated ortholog comparison, then **directly verified against both genomes and found to be present in both** — a KO-assignment pipeline artifact, not real biology. Reported as tested-and-excluded rather than quietly dropped. Same treatment for a KEGG pathway map (`ko01040`) that rendered blank due to a `pathview` bug — traced, fixed, and the real (non-empty) map is now in the report.

Full detail in the [Comparative Genomics Report](report1.html) (HS6 vs. #7), the [RNA-seq Findings Report](report2.html) (HS6 vs. #4), and [Extended Analyses](report4.html) (synthetic-biology parts, transposable elements, flux modeling, and more).

## Differential expression counts

| Contrast | Genes tested | Significant (padj&lt;0.05) | Up in HS6 | Down in HS6 |
|---|---|---|---|---|
| Genotype effect (HS6 vs. #4, adjusted for time) | 13,529 | **1,223** | 398 | 825 |
| Genotype × time interaction (exploratory) | 13,529 | 43 | — | — |

The lopsided up/down split (398 vs. 825) is itself part of the story — it's the same direction as the strong cytoskeleton-down enrichment above, not an independent observation.

## Navigate

- [**HS6 vs. #7 Comparative Genomics Report**](report1.html) — genome relatedness, gene content, the MVK duplication, telomere/chromosome completeness
- [**HS6 vs. #4 RNA-seq Findings Report**](report2.html) — the mutation hunt, MYB3R regulatory rewiring, candidate driver genes
- [Analytical Reproducibility Log](report3.html) — methods, corrections made mid-project, QC, statistical rigor assessment
- [Extended Functional Genomics & Synthetic-Biology Parts](report4.html) — endogenous promoters, introns/UTRs for intron-mediated enhancement, transposable elements, RNA-seq-guided scaffolding, independent re-annotation, expression-constrained flux model
- [Figure 1](figure1.html) — HS6 vs. #7 lipid/isoprenoid/respiratory-chain metabolic map
- [Figure 2](figure2.html) — #4 (WT) vs. HS6 expression-integrated regulatory & metabolic map
- [Scripts](scripts/) — every analysis script used to produce the results in Extended Analyses
- [Data downloads](results/) — ranked candidate lists, classification tables, overlay results

## Limitations & evidence tiers

Findings across this project sit at genuinely different confidence levels, and are labeled as such rather than presented uniformly:

- **Genome-sequence-confirmed** (strongest): the MVK copy-number duplication — verified from raw HiFi reads on both strains independently.
- **Statistically robust, correlational**: the cytoskeleton-down GO enrichment (p.adjust=1.7×10⁻¹³) — a real, strong pattern in the expression data, but expression correlation, not a causal mechanism.
- **Suggestive, single-timepoint-design-limited**: the MYB3R-timing argument, the candidate driver mutations (g6423/g1760/g6124) — plausible leads that would need replicated RNA-seq or functional validation (knockout/complementation) to move beyond "candidate."
- **Exploratory / crude by design**: the expression-constrained flux model (built on a generic, not-organism-specific biomass equation) and computational IME-candidate predictions (architecture is right; actual expression boost unconfirmed without a reporter assay).

No biological replicates exist anywhere in this project's RNA-seq data. That single fact is the largest lever on how much weight any expression-based finding here can bear, and it applies uniformly across every report in this repo.
