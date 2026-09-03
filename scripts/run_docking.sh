#!/bin/bash
set -e
source /home/user/miniforge3/etc/profile.d/conda.sh
conda activate docking

RESULTS=/mnt/d/esmfold_results
OUT=/mnt/d/docking_g6124_g3557
mkdir -p "$OUT"
cd "$OUT"

cp "$RESULTS/g6124_HS6.pdb" receptor.pdb
cp "$RESULTS/g3557_HS6_RPS6.pdb" ligand.pdb

echo "=== setup (auto swarm count from surface density) ==="
lightdock3_setup.py receptor.pdb ligand.pdb --noxt --noh 2>&1

echo "=== running docking (100 steps, dfire scoring) ==="
lightdock3.py setup.json 100 -s dfire -c 4 2>&1

NSWARMS=$(ls -d swarm_* 2>/dev/null | wc -l)
echo "=== generating final structures for $NSWARMS swarms ==="
for swarm_dir in swarm_*; do
  if [ -d "$swarm_dir" ]; then
    cd "$swarm_dir"
    lgd_generate_conformations.py ../receptor.pdb ../ligand.pdb gso_100.out 10 2>&1 | tail -3
    cd ..
  fi
done

echo "=== ranking across all $NSWARMS swarms ==="
lgd_rank.py "$NSWARMS" 100 2>&1

echo "=== DONE ==="
ls -la "$OUT"
