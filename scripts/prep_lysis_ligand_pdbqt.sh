#!/bin/bash
set -e
cd /mnt/c/Users/USER/Desktop/schizochytrium-hs6-genomics/results/lysis_enzymes/ligands
for f in *.sdf; do
  key="${f%.sdf}"
  echo "=== $key ==="
  mk_prepare_ligand.py -i "$f" -o "${key}.pdbqt"
done
echo "--- resulting pdbqt files ---"
ls -la *.pdbqt
