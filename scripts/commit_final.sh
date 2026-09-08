#!/bin/bash
set -e
cd /mnt/c/Users/USER/Desktop/schizochytrium-hs6-genomics

git add results/lysis_enzymes/ scripts/*.py scripts/*.sh
git status --short scripts/ results/ | grep -v '^??.*figures' || true
git commit -m "$(cat <<'EOF'
Expand lysis panel to 18 enzymes and add pH-dependent re-docking

Adds 5 enzymes matching a real wet-lab cell-wall-lysis screen (pectinase,
beta-galactosidase, a Savinase-type protease, an M4 neutral protease, a
B. subtilis xylanase host variant), each verified against UniProt before
folding (2 of 5 initial accession guesses were wrong-organism and were
corrected during fetch). Diagnosed and fixed a real GPU issue: beta-
galactosidase A (1005 aa, far the longest sequence in the panel) was
silently thrashing against the 8GB VRAM ceiling via WSL GPU-passthrough
driver-level memory failures that never raised a catchable OOM, so the
pipeline's own retry-with-smaller-chunk-size logic never fired -- folded
on CPU instead once diagnosed.

Also adds pH-dependent re-docking: each of the 18 structures re-protonated
at pH 4.0 and 9.0 via PDB2PQR/PROPKA, re-docked with identical Vina
settings against both ligands (72 new runs), using the existing pH~7
docking as baseline. Derived per-enzyme gap-vs-pH summary in
ph_gap_summary.json.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
EOF
)"

git add figures/lysis_enzymes/
git commit -m "$(cat <<'EOF'
Full image/GIF coverage for all 18 enzymes

Every enzyme now has a structure thumbnail, a top-pose PNG vs. the
galactan probe, and a 360-degree rotation GIF (previously only 2 of 13
had GIFs). The 4 mechanistically-flagged enzymes (agarase A, pectinase,
beta-galactosidase, neutral protease) also get a canonical-substrate
pose. Fixed the GIF-frame race condition from the prior session by
giving each enzyme its own frame subdirectory instead of a shared one.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
EOF
)"

git add lysis-enzymes.html
git commit -m "$(cat <<'EOF'
Integrate 18-enzyme results, pH sensitivity, and a Vina-choice rationale

- Full 18-row results table (previously only had the original 13);
  pectinase and beta-galactosidase -- added purely to mirror a wet-lab
  reagent list, not for any chemistry reason -- turn out to show the
  strongest and near-strongest galactan-preferring gaps in the whole
  panel, refining the earlier "galactose chemistry predicts affinity"
  idea into a more specific binding-pocket-promiscuity hypothesis (beta-
  agarase A, by contrast, binds its own real substrate with the highest
  specificity in the panel and shows one of the largest gaps).
- New pH & temperature section: the computed pH4/7/9 re-docking results
  (all 3 reversed enzymes stay reversed at every pH tested; lysozyme M1
  reaches -2.08 kcal/mol at pH4, near its own literature-documented
  activity optimum of pH 5.3; two enzymes flip sign only outside their
  documented comfort zone) plus a literature-sourced pH/temperature
  optimum table for all 18 enzymes, cross-referenced against both
  Schizochytrium's native culture conditions and the ~50-60C/mild-
  alkaline condition typically recommended for enzymatic lysis.
- New callout explaining why AutoDock Vina was used instead of a newer
  ML docking method (DiffDock/DiffDock-L): those are trained almost
  entirely on drug-like PDBbind complexes, a poor match for this panel's
  mostly-oligosaccharide ligand set.
- Fold-time table and summary text corrected to actually include the 5
  new enzymes and accurately describe the CPU fallback for one of them.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
EOF
)"

echo "--- log ---"
git log --oneline -8
echo "--- status ---"
git status --short
