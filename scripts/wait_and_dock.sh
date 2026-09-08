#!/bin/bash
# Wait until all 13 ESMFold structures exist, then run the full docking batch.
set -e
cd /mnt/c/Users/USER/Desktop/schizochytrium-hs6-genomics
STRUCT_DIR=results/lysis_enzymes/structures

LABELS=(subtilisin_carlsberg subtilisin_BPN subtilisin_E subtilisin_pumilus \
        lysozyme_C_chicken endolysin_T4 lysozyme_M1_strepto \
        endoglucanase_I_Tr endoglucanase_B_An xylanase2_Tr xylanaseB_An \
        betaagaraseA_Zg betaagaraseB_Zg)

echo "Waiting for all 13 ESMFold structures..."
while true; do
  count=0
  for l in "${LABELS[@]}"; do
    [ -f "$STRUCT_DIR/${l}.pdb" ] && count=$((count+1))
  done
  echo "  $(date +%T) : $count/13 structures present"
  if [ "$count" -eq 13 ]; then
    break
  fi
  sleep 30
done

echo "All 13 structures present. Starting docking batch."
source ~/miniconda3/etc/profile.d/conda.sh
conda activate docking_env
python3 scripts/dock_lysis_panel.py
echo "DOCKING BATCH COMPLETE"
