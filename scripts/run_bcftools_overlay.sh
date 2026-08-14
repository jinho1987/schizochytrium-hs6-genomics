#!/bin/bash
set -e
source /home/user/miniforge3/etc/profile.d/conda.sh
conda activate varcall

REGIONS=$(cat /mnt/c/Users/user/AppData/Local/Temp/claude/C--Users-user-Desktop-HS6-vs--4-RNA-seq-260720/6e32b67b-d05a-4791-b9c3-db6271bce82e/scratchpad/de_gene_regions.txt)
VCF="/home/user/rnaseq_hs6/varcall/WT_vs_HS6.haploid.vcf.gz"

bcftools view -H -r "$REGIONS" "$VCF" > /mnt/c/Users/user/AppData/Local/Temp/claude/C--Users-user-Desktop-HS6-vs--4-RNA-seq-260720/6e32b67b-d05a-4791-b9c3-db6271bce82e/scratchpad/bcftools_hits.txt

wc -l /mnt/c/Users/user/AppData/Local/Temp/claude/C--Users-user-Desktop-HS6-vs--4-RNA-seq-260720/6e32b67b-d05a-4791-b9c3-db6271bce82e/scratchpad/bcftools_hits.txt
