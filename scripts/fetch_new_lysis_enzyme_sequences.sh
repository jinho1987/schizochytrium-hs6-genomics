#!/bin/bash
# Fetch canonical Swiss-Prot sequences for the 5 NEW lysis-panel enzymes added
# to mirror the wet-lab Experiment 4 screen. Appends to the existing
# results/lysis_enzymes/sequences/length_check.csv (does not touch the
# original 13 rows already in that file).
set -euo pipefail

OUT=results/lysis_enzymes/sequences
mkdir -p "$OUT"

# accession label note
declare -a NEW_PANEL=(
  "P26213 pectinase_An Endopolygalacturonase_I_Aspergillus_niger"
  "Q2UCU3 betagalactosidase_Ao Beta-galactosidase_A_Aspergillus_oryzae"
  "P29600 subtilisin_savinase_Bsp Subtilisin_Savinase_Lederbergia_lenta_fka_Bacillus_lentus"
  "P68736 neutralprotease_Bs Bacillolysin_NprE_Bacillus_subtilis"
  "P18429 xylanaseC_Bs Endo-1_4-beta-xylanase_A_Bacillus_subtilis"
)

for entry in "${NEW_PANEL[@]}"; do
  read -r acc label note <<< "$entry"
  fasta="$OUT/${label}.fasta"
  echo "=== Fetching $acc ($label) : $note ==="
  curl -sS "https://rest.uniprot.org/uniprotkb/${acc}.fasta" -o "$fasta"
  if [ ! -s "$fasta" ]; then
    echo "  ERROR: empty response for $acc"
    exit 1
  fi
  seq=$(grep -v '^>' "$fasta" | tr -d '\n')
  actlen=${#seq}
  echo "  header: $(head -1 "$fasta")"
  echo "  length: $actlen"
  echo "$acc,$label,,$actlen,NEW" >> "$OUT/length_check.csv"
  sleep 0.3
done

echo
echo "=== Summary (new rows) ==="
tail -5 "$OUT/length_check.csv"
