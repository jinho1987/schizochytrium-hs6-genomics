#!/bin/bash
# Try PubChem PUG-REST name lookup for each of the 5 canonical ligands.
# Writes results to results/lysis_enzymes/ligands/pubchem_lookup.tsv
set -uo pipefail

OUT=results/lysis_enzymes/ligands
mkdir -p "$OUT"

declare -a NAMES=(
  "succinyl-Ala-Ala-Pro-Phe-4-nitroanilide"
  "chitohexaose"
  "cellopentaose"
  "xylopentaose"
  "neoagarohexaose"
)
declare -a KEYS=(
  "protease_substrate"
  "lysozyme_substrate"
  "cellulase_substrate"
  "xylanase_substrate"
  "agarase_substrate"
)

echo -e "key\tname\tsmiles\tstatus" > "$OUT/pubchem_lookup.tsv"

for i in "${!NAMES[@]}"; do
  name="${NAMES[$i]}"
  key="${KEYS[$i]}"
  enc=$(python3 -c "import urllib.parse,sys; print(urllib.parse.quote(sys.argv[1]))" "$name")
  url="https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/${enc}/property/CanonicalSMILES/TXT"
  echo "=== $key : $name ==="
  resp=$(curl -sS "$url")
  echo "$resp"
  if echo "$resp" | grep -qi "NotFound\|Fault\|<html"; then
    echo -e "${key}\t${name}\t\tNOT_FOUND" >> "$OUT/pubchem_lookup.tsv"
  else
    smiles=$(echo "$resp" | tr -d '\r' | sed -n '1p')
    echo -e "${key}\t${name}\t${smiles}\tFOUND" >> "$OUT/pubchem_lookup.tsv"
  fi
  sleep 0.5
done

echo
cat "$OUT/pubchem_lookup.tsv"
