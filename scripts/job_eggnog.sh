#!/bin/bash
set -e
source /home/user/miniforge3/etc/profile.d/conda.sh

echo "=== reusing eggnog env if it already has diamond/hmmer from the prior attempt ==="
conda activate eggnog 2>/dev/null && which diamond && which hmmsearch || {
  conda create -y -n eggnog -c bioconda -c conda-forge "python=3.10" pip diamond hmmer 2>&1 | tail -30
  conda activate eggnog
}
pip install --quiet "eggnog-mapper==2.1.13"
emapper.py --version

OUT=/mnt/d/HS6_vs_7_comparison/eggnog
rm -rf "$OUT/data"
mkdir -p "$OUT/data"
cd "$OUT"

echo "=== patching dead eggnogdb.embl.de domain (known issue with this release) ==="
sed -i 's|http://eggnogdb.embl.de|http://eggnog5.embl.de|g' "$(dirname "$(which emapper.py)")/download_eggnog_data.py"

echo "=== downloading eggNOG database (this is the big download, ~50GB) ==="
download_eggnog_data.py -y --data_dir "$OUT/data"

echo "=== running emapper.py on HS6 proteome ==="
emapper.py -i /mnt/d/HS6_vs_7_comparison/genomes/HS6_protein.fasta \
  --output hs6_eggnog --output_dir "$OUT" \
  --data_dir "$OUT/data" --cpu 8 -m diamond --itype proteins

echo "=== DONE ==="
wc -l "$OUT/hs6_eggnog.emapper.annotations"
