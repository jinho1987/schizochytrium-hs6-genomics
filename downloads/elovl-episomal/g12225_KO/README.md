# g12225 (ELOVL3/6-type elongase) knockout design -- HS6 homology arms, no Cas9

Strategy: replace the whole g12225 coding sequence with a selectable-marker cassette flanked
by HS6 homology arms (homologous recombination only). The same construct is intended for
both HS6 and #7. Marker/cassette not chosen yet.

Target: HS6 `g12225` (Scaffolds_31:162,976-164,202, minus strand) and #7 `LOCUS_004835`
(contig7:1,768,160-1,769,386). Single exon, 1,227 bp including the stop, 408 aa. The coding
sequence is identical in both strains (both translate cleanly to the same 408 aa protein).
Sequences are in gene orientation, not plus-strand.

## Construct layout

`[5' arm 800 bp] -- [marker cassette] -- [3' arm 800 bp]`
replacing the entire 1,227 bp CDS. Arms are directly adjacent to the ATG (5' arm) and the
stop codon (3' arm). File: `HS6_arms_final.fa`.
HS6 plus-strand coordinates: 5' arm Scaffolds_31:164,203-165,002; 3' arm Scaffolds_31:162,176-162,975.

## Arm-amplification primers (`arm_amplification_primers.tsv`)

Gene-specific part only (add overlap tails / restriction sites once the marker and cloning
method are chosen). Each has one site in the HS6 genome and one in #7.

| Primer | Sequence 5'->3' | Tm |
|---|---|---|
| 5arm_F | CGAACAGGCCCACGAGCTGTTT | 65.8 |
| 5arm_R | GGTCCCAAGAGTTCCGTCTACT | 61.1 |
| 3arm_F | ATCTTCGACGAGGTCGATCAAC | 60.2 |
| 3arm_R | TTCTATTCCCGGGTCTTCTTGC | 60.1 |

5arm_F runs 4-5 C hotter than the rest; use a touchdown or annealing ~60 C.

## How well do the HS6 arms match #7? (`arm_differences_vs_S7.tsv`)

| Arm | Identity | Differences | Longest identical stretch |
|---|---|---|---|
| 5' arm | 95.1% | 11 substitutions + 6 indels (29 bp total; the largest is a 16 bp HS6-only block at -514) | 160 bp (-739 to -578) |
| 3' arm | 98.9% | 9 substitutions, no indels | 193 bp (+1 to +195) |

In #7 the 3' arm side is well matched right up to the stop codon (194 identical bases).
The 5' arm side has a mismatch 14 bp upstream of the ATG (13 identical bases at the junction)
and several indels, so it is the likely weak point for recombination in #7. If knockout
efficiency in #7 is poor, switch to a #7-specific 5' arm (`homology_arms_800bp.fa` still
contains the #7 versions of both arms).

## Genotyping primers (`genotyping_primers.tsv`)

Each primer matches exactly in both strains and has one site per genome. Outer primers lie
outside both arms, so they can't amplify the donor construct.

| Purpose | Forward | Reverse | HS6 | #7 |
|---|---|---|---|---|
| Outer pair (diagnostic) | GAAGAGCTCGTGGTACGTTGA | GCATACGAGCAGCCTTATCCT | 3,903 bp | 3,890 bp |
| Internal CDS pair (diagnostic): band = gene present, none = deleted | TTTGGTGTACGCTGCCCTTAT | CAGTGCCAACAACCATCTGTG | 665 bp | 665 bp |
| 5' control (band in WT and KO) | GAAGAGCTCGTGGTACGTTGA | GGCTTGCTATCACCCCACATA | 722 bp | 724 bp |
| 3' control (band in WT and KO) | TGACAAGTGTCTCCAATCCCG | GCATACGAGCAGCCTTATCCT | 739 bp | 739 bp |

- Expected KO outer-pair product = WT product - 1,227 bp + cassette length.
- 5' junction PCR: outer forward primer + a reverse primer inside the marker.
  3' junction PCR: a forward primer inside the marker + outer reverse primer.
  Marker-side primers can't be designed until the cassette is chosen.

## Not used

Cas9 guide candidates were designed earlier and are kept in `not_used_Cas9/` in case the
plan changes. They are not part of this design.

## Checks and limits

- Not tested in the lab. Specificity was checked against the two assemblies we have.
- An earlier run of the BLAST specificity check failed silently and gave false zeros; this was
  caught (every primer must match its own site at least once), fixed and re-run. All numbers
  here come from the corrected run.
- Without Cas9, recombination efficiency depends on arm identity; the 5' arm in #7 is the
  main risk (see above).
