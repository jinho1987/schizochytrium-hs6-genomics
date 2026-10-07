# HS6 vs. #4 primer design -- Sanger genotyping + qPCR expression validation

Prepared 2026-09-23. Covers all 54 genes carrying a genotype-level SNP shift
between HS6 and #4 (WT), per the SKKU no-filter analysis, cross-checked
independently before any primer was designed.

## Candidate selection and verification (why this list is "exclusive, error-free")

Source: `4.SNP/SNP정리_250801.xlsx` (SKKU's curated mild+no-filter SNP summary,
per instruction to prioritize SKKU results over the vendor 통로바이오 data, and
to use the relaxed/no-filter set). That sheet lists 54 distinct genes.

Every one of the 54 was independently re-verified before being trusted for
primer design -- nothing was taken on the spreadsheet's word alone:

1. Looked up each gene's **exact transcript-level** (mRNA) coordinates in
   `HS6_Genome/braker.gff3` -- not the gene-level span, which can differ when
   a locus has multiple predicted isoforms (this caught a real bug: g11453's
   listed sequence matched its `g11453.t2` GeneMark isoform, 2071bp, not the
   2339bp AUGUSTUS gene-level span -- using the wrong span would have shifted
   every downstream coordinate for that gene).
2. Confirmed the HS6 sequence in the spreadsheet **matches `HS6_scaf.fasta`**
   at those exact coordinates, byte for byte (strand-corrected).
3. Diffed the HS6 vs. #4 sequences to find the exact SNP position(s), and
   confirmed the number found matches the spreadsheet's own reported count.
4. Cross-checked every individual mutated base against the reference genome
   independently a second time.

Result: **54/54 genes passed all four checks cleanly** -- no gene needed to
be dropped, no manual patch was required. (One gene, g11453, required fixing
the coordinate source in step 1 before it passed; it's included normally
below.)

## Files in this folder

- **`sanger_primers.tsv`** -- 74 primer pairs (some genes needed >1 amplicon;
  see below) for confirming each mutation by genomic-DNA Sanger sequencing.
  Columns: `gene`, `transcript_id`, `cluster` (amplicon number within the
  gene), `locus` (genomic span the mutation(s) sit in), `strand`,
  `n_mutations`, `mutations` (exact genomic position + `HS6>#4` base change,
  **always given on the plus strand** -- matches `amplicon_seq` directly,
  regardless of the gene's own strand), `DE_padj`/`DE_log2FC` (expression
  context from the RNA-seq DE table), `fwd_primer`/`rev_primer` + their Tm,
  `product_size_bp` (501-900bp for every pair -- a single clean Sanger read
  comfortably covers each one), `amplicon_seq`, and `purpose`.
- **`qpcr_primers.tsv`** -- 54 primer pairs (one per gene), for RT-qPCR
  expression validation -- i.e. independently confirming the RNA-seq DE call
  by qPCR, not genotyping. Designed on the spliced CDS, product size
  80-200bp (standard qPCR range). `DE_significance` flags the 7 genes that
  are actually significant DEGs (padj<0.05) -- prioritize these first if
  bench time is limited (see below).
- **`flanking_500bp_context.tsv`** / **`flanking_500bp_context.fasta`** --
  the raw 500bp-upstream + 500bp-downstream genomic context around every
  mutation cluster (the sequence primers were designed from), for anyone who
  wants to design alternative primers by hand or just inspect the sequence
  directly. FASTA headers encode gene, transcript, cluster number, exact
  genomic span, and the mutation list.

Both TSVs now carry five extra columns from the BLAST verification below:
`fwd_true_sites`, `rev_true_sites` (how many places in `HS6_scaf.fasta` each
primer matches at >=95% identity over its full length), `severity`
(CLEAN/MODEST/ELEVATED/SEVERE, see below), and `fwd_other_loci`/`rev_other_loci`
(up to 5 of the other genomic locations, for anyone who wants to inspect them).

## Why some genes have more than one Sanger amplicon

A handful of genes (e.g. g11453, which carries 11 SNPs) have their mutations
spread across more than ~1.5kb of the transcript -- too wide for one Sanger
read to reliably cover start to finish. Rather than force one oversized,
partially-unreadable amplicon, mutations were grouped into clusters (nearby
SNPs sharing one amplicon; distant ones getting their own), so **every
amplicon in this file is a realistic, single-read-covered product (501-900bp)**.
g11453 accordingly has 4 amplicons covering its 11 SNPs, not 1.

## The 7 genes that are both mutated AND significantly differentially expressed

These are the strongest candidates for explaining the HS6 growth/glucose-
consumption phenotype, and are the highest-priority qPCR targets:

| Gene | padj | log2FC | Note |
|---|---|---|---|
| g2636 | 0.0018 | +1.25 | most significant single hit; CDS-level SNP |
| g6055 (RanGAP2) | 0.0094 | -1.55 | intronic mutation only |
| g1664 | 0.0144 | -1.11 | p.H205Q missense, Armadillo-repeat domain |
| g6423 (RanGAP1) | 0.0189 | -1.11 | intronic mutation only |
| g14345 | 0.0330 | +0.76 | intronic; prior RNA-fold/coiled-coil analysis exists |
| g15200 | 0.0426 | -0.99 | intronic mutation only |
| g845 (RanGAP3) | 0.0490 | -1.29 | intronic mutation only |

Note the RanGAP trio (g6055/g6423/g845) all move together -- same direction,
similar magnitude, same mutation class -- a family-wide effect, not an
isolated gene. See report4.html's RanGAP docking section for the structural
side of this story.

## BLAST verification -- an important, unresolved finding

Every primer (all 128 pairs) was BLASTed against `HS6_scaf.fasta`
(blastn-short, then a strict full-length/>=95%-identity check for what counts
as a genuine binding site -- a loose "90% over 18-20bp" filter turns out to
be nearly useless at this length, since that much of the sequence can match
by chance; full-length near-identity is the only meaningful bar). For each
design, Primer3 was asked for 5 ranked candidate pairs, and the
highest-ranked one that BLASTed clean (<=2 true sites for both primers) was
selected automatically where one existed.

**Result: most of the panel is not unique in this genome, and it isn't
fixable by picking a different Primer3 candidate.**

| Tier | Sanger (74) | qPCR (54) | Meaning |
|---|---|---|---|
| CLEAN | 6 | 11 | single genomic match -- safe to order as-is |
| MODEST | 10 | 8 | 2-3 matches, likely a close paralog -- probably workable |
| ELEVATED | 12 | 10 | 4-10 matches -- confirm by gel size / melt-curve before trusting results |
| SEVERE | 46 | 25 | >10 matches (some in the hundreds) -- will not act as a specific primer, needs a different design strategy entirely |

The SEVERE cases are not a bug in this pipeline -- they were traced to actual
repeated sequence in `HS6_scaf.fasta` itself. One example: the 20-mer
`AGCGACCGGTTCTCTTTCAG` was independently chosen by Primer3 for three
unrelated genes (g6423, g1302, g845) because it's a genuinely strong-looking
primer candidate (good Tm/GC) -- and it also appears ~740 times across
multiple scaffolds in a regular, tandem-repeat-like spacing pattern. Tried
masking with this genome's own EarlGrey-derived repeat family library
(`working_files/4.repeat/HS6_earlgray_summary/HS6-families.fa.strained.clstrd.fa`)
as a Primer3 mispriming library, but the specific repeat responsible isn't
softmasked by that annotation (spot-checked directly), and using the library
as a hard filter rejected nearly every candidate genome-wide, so it wasn't a
usable fix either. Distinct primers landing on a **small, recurring set of
site counts** (740, 560, 508, 223...) across otherwise-unrelated genes is the
signature of a modest number of highly abundant repeat/duplication families
in this draft assembly, not scattered coincidences -- worth investigating as
a genome-quality question in its own right, separate from this primer task.

**What this means practically:** the 6 Sanger + 11 qPCR CLEAN-tier primers
are safe to order now. Everything else needs a decision -- accept the risk
for MODEST/ELEVATED tiers with wet-lab confirmation (gel size check for
Sanger; melt-curve for qPCR), or invest in a proper fix for the SEVERE tier
(e.g. genome self-alignment to map out duplicated/repetitive blocks first,
then design specifically outside them, or fall back to longer/nested primers
anchored further from the repeat-dense regions). Not decided here --
flagging for your call before anything gets ordered.
