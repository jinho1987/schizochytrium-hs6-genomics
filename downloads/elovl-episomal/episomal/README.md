# Genome screen for an episomal plasmid in Schizochytrium (HS6 genome, checked against #7)

Genome: 61.6 Mb, 34 scaffolds, GC 45.1%. An episomal plasmid needs (1) something that makes it
replicate and be passed on, (2) a way to express the marker/cargo, (3) codons the host reads
well. This screen looks for candidates for each. Everything in (1) is a hypothesis from
sequence alone; none of it is known to work until tested by transformation.

## 1. Replication and maintenance candidates

**Telomere repeat (strongest result).** A 5-base unit, AACCT on the strand shown (the TTAGG
type), forms 320-830 bp tandem arrays at **41 of the 68 scaffold ends**. This identifies the
host telomere repeat. Possible use: telomere seeds at the ends of a linear plasmid.
Two ~1 kb end sequences (each with an 830 bp array) are in `telomere_seed_candidates.fa`.

**ARS-like (replication-origin-like) regions.** No replication origins have been mapped in
this organism, so these are guesses based on yeast-type features: AT-rich (>=78.5% AT, top
0.5% of intergenic windows; genome-wide intergenic mean is 58.5%), >=300 bp, between genes,
single-copy in HS6, outside repeat families, and conserved in #7. 254 regions pass.
- Warning: 81 of 254 sit between divergently transcribed genes, which means they are probably
  shared promoter regions, not origins. 117 lie between same-direction genes and 56 between
  convergent genes (3' ends facing each other).
- `ARS_like_candidates_classified.tsv` lists all 254 with orientation class, AT%, yeast-ARS-like
  motif count, #7 identity and neighboring genes. The 10 convergent-gene regions with the most
  motifs are listed in `screen1.log` / the summary. 254 is too many to test; plan on picking
  ~6-10 spread across scaffolds.

**Centromere-like regions (weak).** Satellite-repeat clusters of 9-13 kb exist, but most sit
near scaffold ends and all have genes within 20 kb, so there is no clear gene-poor centromere
candidate. We haven't identified a centromere. Not pursued further.

## 2. Native expression elements (`expression_element_candidates.tsv`, `promoters_terminators.fa`)

From 6 RNA-seq samples (WT and MT at 20, 44, 68 h): genes in the top 5% by expression with
CV <= 0.25 across samples, >=600 bp free upstream and >=250 bp free downstream, no repeat
overlap; 202 genes pass; the 15 listed also have promoter and terminator conserved in #7 (94-100%
identity). Promoter = up to 1,000 bp upstream of the start codon (includes any 5' UTR);
terminator = up to 500 bp downstream of the stop. Sequences are in gene orientation.

Most conventional housekeeping picks among them: **g227 (elongation factor 2, mean 18,408,
CV 0.16)** and **g12164 (ADP-ribosylation factor, CV 0.23)**. g827 (plasma-membrane H+-ATPase,
highest expression, CV 0.09) is also steady. Several others are transporters or proteases whose
promoters may be condition-dependent.

Limit: "steady" means steady across our three time points in one culture condition. Test 3-4
promoters side by side with a reporter before committing.

## 3. Codon usage (`codon_usage_HS6.tsv`)

15,960 coding sequences, 8.8 M codons. GC3 = 47.2%, close to the genome average, so there is
no strong GC bias; codon adaptation of a marker gene should be modest. The table gives
per-codon fractions and RSCU for adapting a marker.

## Not covered

Mitochondrial/rDNA origins, and any published Schizochytrium replication element, were not
searched (not available from our data). Yeast-type origin features are an assumption for this
organism.
