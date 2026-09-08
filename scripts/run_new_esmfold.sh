#!/bin/bash
set -e
cd /mnt/c/Users/USER/Desktop/schizochytrium-hs6-genomics
source /home/jin/miniconda3/etc/profile.d/conda.sh
conda activate esm_scoring
python3 scripts/run_esmfold_lysis_panel.py
echo "ESMFOLD_NEW_BATCH_COMPLETE"
