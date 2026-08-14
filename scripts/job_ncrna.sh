#!/bin/bash
set -e
source /home/user/miniforge3/etc/profile.d/conda.sh

echo "=== installing ncRNA tools ==="
conda create -y -n ncrna -c bioconda -c conda-forge trnascan-se barrnap infernal 2>&1 | tail -30

conda activate ncrna

OUT=/mnt/d/HS6_vs_7_comparison/ncrna
mkdir -p "$OUT"
cd "$OUT"

HS6=/mnt/d/HS6_vs_7_comparison/genomes/HS6_genome.fasta
S7=/mnt/d/HS6_vs_7_comparison/genomes/strain7_genome.fasta

echo "=== tRNAscan-SE: HS6 ==="
tRNAscan-SE -E -o hs6_trnas.txt -f hs6_trnas.struct -s hs6_trnas.iso --thread 8 "$HS6" 2>&1 | tail -20

echo "=== tRNAscan-SE: strain7 ==="
tRNAscan-SE -E -o s7_trnas.txt -f s7_trnas.struct -s s7_trnas.iso --thread 8 "$S7" 2>&1 | tail -20

echo "=== barrnap (rRNA): HS6 ==="
barrnap --kingdom euk --threads 8 "$HS6" > hs6_rrna.gff 2>hs6_barrnap.log

echo "=== barrnap (rRNA): strain7 ==="
barrnap --kingdom euk --threads 8 "$S7" > s7_rrna.gff 2>s7_barrnap.log

echo "=== DONE ==="
echo "HS6 tRNA count:"; grep -c '^Scaffolds' hs6_trnas.txt 2>/dev/null || wc -l hs6_trnas.txt
echo "strain7 tRNA count:"; grep -c '^contig' s7_trnas.txt 2>/dev/null || wc -l s7_trnas.txt
echo "HS6 rRNA hits:"; wc -l hs6_rrna.gff
echo "strain7 rRNA hits:"; wc -l s7_rrna.gff
