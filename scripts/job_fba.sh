#!/bin/bash
set -e
source /home/user/miniforge3/etc/profile.d/conda.sh

echo "=== installing diamond via conda (only the compiled binary needs conda) ==="
conda create -y -n fba -c bioconda -c conda-forge python=3.10 diamond pip 2>&1 | tail -40

conda activate fba
echo "=== installing carveme + cobrapy via pip (conda solver kept conflicting on old python pins) ==="
pip install --quiet cobra carveme
pip install --quiet swiglpk 2>&1 | tail -10

OUT=/mnt/d/HS6_vs_7_comparison/fba
mkdir -p "$OUT"
cd "$OUT"

echo "=== building draft genome-scale model with CarveMe (HS6 proteome) ==="
carve /mnt/d/HS6_vs_7_comparison/genomes/HS6_protein.fasta \
  --fbc2 -o hs6_draft_model.xml -v 2>&1 | tee carveme_hs6.log

echo "=== CarveMe done, model summary ==="
python3 -c "
import cobra
m = cobra.io.read_sbml_model('hs6_draft_model.xml')
print('Reactions:', len(m.reactions))
print('Metabolites:', len(m.metabolites))
print('Genes:', len(m.genes))
sol = m.optimize()
print('Default FBA growth objective:', sol.objective_value)
"

echo "=== DONE building draft model, ready for expression-constrained FBA in next step ==="
