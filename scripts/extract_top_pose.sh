#!/bin/bash
# Extract the top (mode 1) pose from a multi-model Vina PDBQT output as a clean PDB
# for rendering with PyMOL.
set -e
DOCK_DIR=/mnt/c/Users/USER/Desktop/schizochytrium-hs6-genomics/results/lysis_enzymes/docking
cd "$DOCK_DIR"
for f in *.pdbqt; do
  case "$f" in
    *_receptor.pdbqt) continue ;;
  esac
  base="${f%.pdbqt}"
  out="${base}_top_pose.pdb"
  obabel "$f" -O "$out" -f 1 -l 1 2>&1 | tail -1
done
echo "done"
ls -la *_top_pose.pdb 2>/dev/null | wc -l
