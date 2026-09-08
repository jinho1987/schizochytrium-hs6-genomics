#!/bin/bash
# Wait for all 18 ESMFold structures, then run the extended docking batch
# (10 new runs -> 36 total) and the full render pass (all 18 enzymes get a
# static structure PNG, a galactan-probe pose PNG+GIF, and the 3 new special
# enzymes + already-existing agarase A also get a canonical-substrate pose).
#
# IMPORTANT: only run ONE instance of this script (or of render_lysis_docking.py
# directly) at a time -- see the race-condition note in render_lysis_docking.py.
set -e
cd /mnt/c/Users/USER/Desktop/schizochytrium-hs6-genomics
SUMMARY=results/lysis_enzymes/esmfold_summary.json

echo "Waiting for all 18 ESMFold structures..."
while true; do
  count=$(grep -c '"label"' "$SUMMARY" 2>/dev/null || echo 0)
  echo "  $(date +%T) : $count/18 folds recorded"
  if [ "$count" -ge 18 ]; then
    break
  fi
  sleep 30
done

echo "All 18 folds present. Running extended docking batch (docking_env)..."
source /home/jin/miniconda3/etc/profile.d/conda.sh
conda activate docking_env
python3 scripts/dock_lysis_panel.py

echo "Extracting top poses..."
bash scripts/extract_top_pose.sh

echo "Rendering structure thumbnails (pymol_render)..."
conda activate pymol_render
pymol -cq scripts/render_lysis_structures.py

echo "Rendering docking pose figures (18 static+gif vs galactan, plus 3 canonical poses)..."
pymol -cq scripts/render_lysis_docking.py

echo "NEW_PIPELINE_COMPLETE"
ls -la figures/lysis_enzymes/ | tail -40
