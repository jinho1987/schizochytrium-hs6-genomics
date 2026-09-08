#!/bin/bash
# Fetch canonical Swiss-Prot sequences for the 13-enzyme lysis panel from UniProt.
# Saves one FASTA per enzyme to results/lysis_enzymes/sequences/<label>.fasta
# and checks length against the published lysis-enzymes.html table as a sanity check.
set -euo pipefail

OUT=results/lysis_enzymes/sequences
mkdir -p "$OUT"

# accession label expected_length
declare -a PANEL=(
  "P00780 subtilisin_carlsberg 379"
  "P00782 subtilisin_BPN 382"
  "P04189 subtilisin_E 381"
  "P07518 subtilisin_pumilus 275"
  "P00698 lysozyme_C_chicken 147"
  "P00720 endolysin_T4 164"
  "P25310 lysozyme_M1_strepto 294"
  "P07981 endoglucanase_I_Tr 459"
  "O74706 endoglucanase_B_An 331"
  "G0RUP7 xylanase2_Tr 223"
  "P55330 xylanaseB_An 225"
  "G0L322 betaagaraseA_Zg 539"
  "Q9RGX8 betaagaraseB_Zg 353"
)

echo "accession,label,expected_len,actual_len,match" > "$OUT/length_check.csv"

for entry in "${PANEL[@]}"; do
  read -r acc label explen <<< "$entry"
  fasta="$OUT/${label}.fasta"
  echo "=== Fetching $acc ($label) ==="
  curl -sS "https://rest.uniprot.org/uniprotkb/${acc}.fasta" -o "$fasta"
  if [ ! -s "$fasta" ]; then
    echo "  ERROR: empty response for $acc"
    exit 1
  fi
  seq=$(grep -v '^>' "$fasta" | tr -d '\n')
  actlen=${#seq}
  match="OK"
  if [ "$actlen" != "$explen" ]; then
    match="MISMATCH"
  fi
  echo "  length: $actlen (expected $explen) -> $match"
  echo "$acc,$label,$explen,$actlen,$match" >> "$OUT/length_check.csv"
  sleep 0.3
done

echo
echo "=== Summary ==="
cat "$OUT/length_check.csv"
