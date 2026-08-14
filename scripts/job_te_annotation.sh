#!/bin/bash
set -e
source /home/user/miniforge3/etc/profile.d/conda.sh

echo "=== conda env repeatmod already installed, reusing ==="
conda activate repeatmod

OUT=/mnt/d/HS6_vs_7_comparison/te_annotation
rm -rf "$OUT/hs6" "$OUT/strain7"
mkdir -p "$OUT/hs6" "$OUT/strain7"

echo "=== HS6: building RepeatModeler database ==="
cd "$OUT/hs6"
BuildDatabase -name hs6_db /mnt/d/HS6_vs_7_comparison/genomes/HS6_genome.fasta
echo "=== HS6: running RepeatModeler2 (de novo TE family discovery) ==="
RepeatModeler -database hs6_db -pa 2
if [ ! -f hs6_db-families.fa ]; then echo "FATAL: RepeatModeler produced no families file for HS6"; exit 1; fi

echo "=== HS6: running RepeatMasker with the de novo library ==="
RepeatMasker -pa 8 -gff -lib hs6_db-families.fa -dir hs6_masked /mnt/d/HS6_vs_7_comparison/genomes/HS6_genome.fasta

echo "=== strain7: building RepeatModeler database ==="
cd "$OUT/strain7"
BuildDatabase -name s7_db /mnt/d/HS6_vs_7_comparison/genomes/strain7_genome.fasta
echo "=== strain7: running RepeatModeler2 ==="
RepeatModeler -database s7_db -pa 2
if [ ! -f s7_db-families.fa ]; then echo "FATAL: RepeatModeler produced no families file for strain7"; exit 1; fi

echo "=== strain7: running RepeatMasker with the de novo library ==="
RepeatMasker -pa 8 -gff -lib s7_db-families.fa -dir s7_masked /mnt/d/HS6_vs_7_comparison/genomes/strain7_genome.fasta

echo "=== DONE ==="
echo "HS6 TE summary:"
cat "$OUT/hs6/hs6_masked"/*.tbl 2>/dev/null
echo "strain7 TE summary:"
cat "$OUT/strain7/s7_masked"/*.tbl 2>/dev/null
