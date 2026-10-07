# Schizochytrium episomal test plasmids -- first-pass design

Purpose: test whether any of our candidate replication/maintenance elements lets a plasmid
persist as an episome in HS6 / #7. Nothing here is known to work in Schizochytrium; the
design is built so that a positive result can't be confused with random genome integration.

## Architecture (all circular maps are GenBank files with annotated parts)

Backbone: pUC19c (NCBI L09137.2), kept byte-for-byte identical except the polylinker, so the
*E. coli* origin and the ampicillin-resistance gene (bla, verified as the intact 286-aa ORF)
are untouched. The polylinker between EcoRI and HindIII (45 bp) is replaced by:

`PacI - [replication candidate slot] - AscI - promoter g227 - ble - terminator g227 - SbfI`

- **Marker cassette:** native g227 (elongation factor 2) promoter (1,000 bp) and terminator
  (500 bp) from HS6, driving *Sh ble* (zeocin resistance, UniProt P17493, 124 aa). The ble
  coding sequence (375 bp including stop) was back-translated with HS6 codon frequencies
  (GC 54%, stop TAA, the most common HS6 stop) and checked: it translates exactly to P17493, and
  has no homopolymer of 6+, no cloning-enzyme sites, and no telomere-like repeats.
- **PacI/AscI slot:** swap replication candidates by cutting with PacI + AscI. PacI, AscI,
  SbfI and PmeI each occur once (PmeI twice in the telomere version) in every map; checked.
- The g227 promoter contains one BsaI and one BsmBI site, and the terminator one BsaI site. Irrelevant for this
  restriction-cloning design, but they would block a later Golden Gate switch.

## The plasmids

| File | Size | What it tests |
|---|---|---|
| pSCHZ_noARS_control | 4,540 bp | **Essential control.** Same plasmid, no replication candidate. Any zeocin-resistant colonies = background integration/random events. |
| pSCHZ_ARS_Scaffolds_18_2316401-2316900 | 5,040 bp | ARS-like candidate 1 (500 bp, 85% AT, 2 yeast-ARS-like motifs, single copy in both genomes) |
| pSCHZ_ARS_Scaffolds_19_846151-846700 | 5,090 bp | candidate 2 (550 bp) |
| pSCHZ_ARS_Scaffolds_10_453601-454200 | 5,140 bp | candidate 3 (600 bp) |
| pSCHZ_TEL | 5,276 bp | Linear-plasmid test. Cut with PmeI (releases an 8 bp spacer) to get a 5,268 bp linear molecule with a 360 bp native (CCTAA)n telomere block at the left end and its reverse complement (TTAGG)n at the right end. |

`ARS_candidate_library.tsv` has 6 candidates (all single-copy in both genomes, no cloning-site
conflicts; a 7th was rejected for containing a PmeI site) with PCR primers carrying PacI/AscI
tails (leading ATCG clamp), so the other 3 can be added by PCR from HS6 genomic DNA, or ordered
as fragments (the full sequences are in the table; the PCR product is slightly shorter than
the full candidate). `parts.fa` has the cassette parts and telomere blocks.

## Suggested test logic

1. First, zeocin kill curves on HS6 and #7 in your normal medium. **Risk:** zeocin is less
   active at high salt and in rich/marine-type media, which is how many Schizochytrium media
   are made. If wild-type cells tolerate your medium's zeocin, change the marker (the
   promoter/terminator carry over; only the ble coding sequence would be swapped).
2. Transform each plasmid in parallel by the same method. Compare colony counts with the
   no-ARS control; an element only counts if it clearly beats the control.
3. For colonies that grow: grow without selection for several generations, then check whether
   resistance is lost (episomes usually are) versus kept (integration). Recover plasmid into
   *E. coli* from total DNA and compare restriction patterns with the original; PCR across
   the ble cassette and PacI/AscI junction.
4. A positive result for the telomere version (circular vs linear) tells you whether linear
   telomere-ended maintenance is possible.

## Limits (what is not established)

- The replication candidates and the telomere block are sequence-based hypotheses.
- The PmeI cut leaves 4 non-telomeric bases at each end (AAAC at the left, GTTT at the right).
  If telomere addition seems to fail, generate the linear molecule by PCR with primers that
  begin directly with the repeat instead.
- No transformation or growth data for this host were available to check marker choice, promoter
  strength, or the transformation method; those need to be established experimentally.
- Constructs were verified in silico only (sizes, unique sites, ORF translations, unchanged
  backbone), not sequenced.

## Addendum: yeast CEN6 test plasmids

Two plasmids were added to test whether a centromere element helps. We have no confirmed
Schizochytrium centromere (the gene-poor regions found by the screen are only candidates), so
the elements used are borrowed from yeast; whether they do anything in Schizochytrium is unknown.

- **pSCHZ_CEN6-ARSH4_yeast (5,051 bp):** the standard 511 bp yeast CEN6/ARSH4 fragment in the slot, inserted contiguously (CEN6 + the 11 bp native linker + ARSH4, exactly as in pRS416).
  Taken from pRS416 (U03450) and located by matching: CEN6 = pRS416 4,320-4,444 (100% identical to
  yeast chromosome VI 148,504-148,628, 125 bp); ARSH4 = pRS416 4,456-4,830 (98.4% identical to
  yeast chromosome II 254,788-255,162, 375 bp). It has no meaningful similarity to the HS6 or #7
  genomes (one short 72-79 bp local match, probably chance).
- **pSCHZ_ARS18_plus_CEN6 (5,165 bp):** native ARS candidate 1 followed by CEN6, to ask whether a
  centromere element adds stability to a native candidate.
- Both pass the same checks as the other plasmids (unique PacI/AscI/SbfI sites, ble translates
  correctly, bla intact, backbone unchanged). Compare against the no-ARS control exactly as above.
- Diagnostic primers: see the "CEN plasmids" sheet in plasmid_diagnostic_primers.xlsx.

- **No CEN6-like sequence in HS6 or #7.** A sensitive search with CEN6 alone gave only chance-level hits (best about 80 bp at ~80% identity, comparable to a shuffled control), and a scan for the point-centromere architecture (CDEI motif, 70-100 bp AT-rich spacer at least 88% AT, CDEIII core) found 0 matches in either genome. Point centromeres are largely a budding-yeast feature, so the CEN6 plasmids may show no centromere effect in this host; the replication-element part (ARSH4 or the native candidate) is more likely to matter. The scan is strict and could miss a highly diverged element.
