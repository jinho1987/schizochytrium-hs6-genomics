#!/bin/bash
# Wait until all 26 docking runs are recorded, then extract top poses and render figures.
set -e
cd /mnt/c/Users/USER/Desktop/schizochytrium-hs6-genomics
RESULTS=results/lysis_enzymes/docking_results.json

echo "Waiting for all 26 docking results..."
while true; do
  count=$(python3 -c "import json,os; p='$RESULTS'; print(len(json.load(open(p))) if os.path.exists(p) else 0)")
  echo "  $(date +%T) : $count/26 docking runs recorded"
  if [ "$count" -ge 26 ]; then
    break
  fi
  sleep 30
done

echo "All docking runs done. Extracting top poses..."
bash scripts/extract_top_pose.sh

echo "Rendering structure thumbnails (13)..."
source ~/miniconda3/etc/profile.d/conda.sh
conda activate pymol_render
pymol -cq scripts/render_lysis_structures.py

echo "Rendering docking pose figures (13 static + 2 gifs)..."
pymol -cq scripts/render_lysis_docking.py

echo "RENDER PIPELINE COMPLETE"
ls -la figures/lysis_enzymes/
