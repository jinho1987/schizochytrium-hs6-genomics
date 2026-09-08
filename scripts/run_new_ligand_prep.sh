#!/bin/bash
set -e
cd /mnt/c/Users/USER/Desktop/schizochytrium-hs6-genomics
source /home/jin/miniconda3/etc/profile.d/conda.sh
conda activate docking_env
python3 scripts/prep_new_lysis_ligands.py
bash scripts/prep_lysis_ligand_pdbqt.sh
echo "NEW_LIGAND_PREP_COMPLETE"
