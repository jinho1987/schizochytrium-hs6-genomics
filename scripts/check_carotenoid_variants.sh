#!/bin/bash
set -e
VCF=/home/user/rnaseq_hs6/varcall/WT_vs_HS6.haploid.vcf.gz
BCFTOOLS=/home/user/miniforge3/envs/varcall/bin/bcftools

declare -A GENES
GENES[g1472]="Scaffolds_11:748701-749720"
GENES[g3852]="Scaffolds_17:297998-298864"
GENES[g4282]="Scaffolds_18:176565-177422"
GENES[g5428]="Scaffolds_19:411279-411989"
GENES[g6151]="Scaffolds_2:425025-429011"
GENES[g7205]="Scaffolds_21:705171-706814"
GENES[g7815]="Scaffolds_22:190239-191855"
GENES[g14287]="Scaffolds_4:2617466-2619016"

echo "gene	region_checked	variants_found"
for gid in "${!GENES[@]}"; do
  region="${GENES[$gid]}"
  scaf="${region%%:*}"
  range="${region##*:}"
  start="${range%%-*}"
  end="${range##*-}"
  win_start=$((start - 2000))
  if [ "$win_start" -lt 1 ]; then win_start=1; fi
  win_end=$((end + 2000))
  hits=$("$BCFTOOLS" view -H "$VCF" "${scaf}:${win_start}-${win_end}" 2>/dev/null | wc -l)
  echo "$gid	${scaf}:${win_start}-${win_end}	$hits"
  if [ "$hits" -gt 0 ]; then
    "$BCFTOOLS" view -H "$VCF" "${scaf}:${win_start}-${win_end}" 2>/dev/null
  fi
done
