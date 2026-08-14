#!/bin/bash
set -e
source /home/user/miniforge3/etc/profile.d/conda.sh

OUT=/mnt/d/HS6_vs_7_comparison/rascaf_hs6
mkdir -p "$OUT"
cd "$OUT"

echo "=== installing rascaf ==="
conda create -y -n rascaf -c bioconda -c conda-forge rascaf samtools 2>&1 | tail -20

conda activate rascaf

GENOME=/mnt/d/HS6_vs_7_comparison/genomes/HS6_genome.fasta
BAMDIR=/home/user/rnaseq_hs6/align

echo "=== merging all 6 RNA-seq BAMs ==="
samtools merge -f -@ 4 merged_rnaseq.bam \
  "$BAMDIR"/WT_20.sorted.bam "$BAMDIR"/WT_44.sorted.bam "$BAMDIR"/WT_68.sorted.bam \
  "$BAMDIR"/MT_20.sorted.bam "$BAMDIR"/MT_44.sorted.bam "$BAMDIR"/MT_68.sorted.bam
samtools sort -@ 4 -o merged_rnaseq.sorted.bam merged_rnaseq.bam
samtools index merged_rnaseq.sorted.bam

echo "=== running rascaf ==="
rascaf -b merged_rnaseq.sorted.bam -f "$GENOME" -o hs6_rascaf

echo "=== running rascaf-join ==="
rascaf-join -r hs6_rascaf.out -o hs6_rascaf_scaffolded

echo "=== DONE ==="
grep -c '^>' "$GENOME"
grep -c '^>' hs6_rascaf_scaffolded.fa
