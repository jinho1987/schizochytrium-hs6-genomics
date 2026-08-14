import re, csv, subprocess

GFF = r"C:\Users\user\Desktop\HS6 vs. 4_RNA-seq_260720\HS6_Genome\braker.gff3"
TABLES = r"C:\Users\user\Desktop\HS6 vs. 4_RNA-seq_260720\RNA-seq_analysis_results\tables"
VCF = "/mnt/c/Users/user/Desktop/HS6 vs. 4_RNA-seq_260720/RNA-seq_analysis_results/variant_calling/WT_vs_HS6.haploid.vcf.gz"
UPSTREAM = 2000

# 1) load significant DE genes
sig = {}
with open(f"{TABLES}\\DE_genotype_MTvsWT_adj_for_time.tsv", encoding="utf-8") as f:
    reader = csv.reader(f, delimiter="\t")
    next(reader)
    for row in reader:
        if not row or not row[0]:
            continue
        gid = row[0].replace(".t1", "")
        try:
            lfc = float(row[2]); padj = float(row[6]) if row[6] not in ("", "NA") else None
        except (ValueError, IndexError):
            continue
        if padj is not None and padj < 0.05:
            sig[gid] = (lfc, padj)

print(f"Significant DE genes: {len(sig)}")

# 2) gene coordinates
coords = {}
with open(GFF, encoding="utf-8") as f:
    for line in f:
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 9 or p[2] != "gene":
            continue
        m = re.search(r"ID=([^;]+)", p[8])
        if m and m.group(1) in sig:
            coords[m.group(1)] = (p[0], int(p[3]), int(p[4]), p[6])

print(f"Coordinates found for: {len(coords)} / {len(sig)}")

# 3) build a BED-like region list: gene body + upstream promoter window (strand-aware)
regions = []  # (gid, scaf, region_start, region_end)
for gid, (scaf, start, end, strand) in coords.items():
    if strand == "+":
        rstart = max(1, start - UPSTREAM)
        rend = end
    else:
        rstart = start
        rend = end + UPSTREAM
    regions.append((gid, scaf, rstart, rend, start, end, strand))

# 4) write the region string to a file for a separate bcftools query step (run via Bash tool)
region_str = ",".join(f"{scaf}:{rstart}-{rend}" for gid, scaf, rstart, rend, s, e, strand in regions)
with open(r"C:\Users\user\AppData\Local\Temp\claude\C--Users-user-Desktop-HS6-vs--4-RNA-seq-260720\6e32b67b-d05a-4791-b9c3-db6271bce82e\scratchpad\de_gene_regions.txt", "w") as f:
    f.write(region_str)
print(f"Wrote region string ({len(region_str)} chars) to de_gene_regions.txt -- run bcftools separately, then rerun this script's second half.")

# read back bcftools output if it already exists
import os
BCFOUT = r"C:\Users\user\AppData\Local\Temp\claude\C--Users-user-Desktop-HS6-vs--4-RNA-seq-260720\6e32b67b-d05a-4791-b9c3-db6271bce82e\scratchpad\bcftools_hits.txt"
if not os.path.exists(BCFOUT):
    print("bcftools_hits.txt not found yet -- stopping here.")
    raise SystemExit(0)
with open(BCFOUT, encoding="utf-8") as f:
    lines = [l for l in f.read().strip().split("\n") if l]
print(f"Raw variant lines returned: {len(lines)}")

variants = []
for line in lines:
    parts = line.split("\t")
    if len(parts) < 5:
        continue
    scaf, pos, vid, ref, alt = parts[0], int(parts[1]), parts[2], parts[3], parts[4]
    variants.append((scaf, pos, ref, alt))

# 5) assign each variant to overlapping gene region(s), classify as CDS-body vs promoter-window
results = []
for gid, scaf, rstart, rend, start, end, strand in regions:
    for vscaf, vpos, ref, alt in variants:
        if vscaf == scaf and rstart <= vpos <= rend:
            in_body = start <= vpos <= end
            loc = "gene_body" if in_body else "promoter_window"
            lfc, padj = sig[gid]
            results.append((gid, scaf, vpos, ref, alt, loc, lfc, padj))

print(f"\nSignificant DE genes with a variant in body or {UPSTREAM}bp promoter window: {len(set(r[0] for r in results))}")
print(f"{'gene':10s} {'scaffold':14s} {'pos':>10s} {'ref>alt':12s} {'location':16s} {'log2FC':>8s} {'padj':>10s}")
for gid, scaf, pos, ref, alt, loc, lfc, padj in sorted(results, key=lambda r: r[7]):
    print(f"{gid:10s} {scaf:14s} {pos:10d} {ref+'>'+alt:12s} {loc:16s} {lfc:8.2f} {padj:10.2e}")

with open(r"C:\Users\user\AppData\Local\Temp\claude\C--Users-user-Desktop-HS6-vs--4-RNA-seq-260720\6e32b67b-d05a-4791-b9c3-db6271bce82e\scratchpad\genomewide_mutation_overlay.tsv", "w", encoding="utf-8") as out:
    out.write("gene\tscaffold\tpos\tref_alt\tlocation\tlog2FC\tpadj\n")
    for gid, scaf, pos, ref, alt, loc, lfc, padj in results:
        out.write(f"{gid}\t{scaf}\t{pos}\t{ref}>{alt}\t{loc}\t{lfc:.3f}\t{padj:.4e}\n")
print("\nwrote full results to genomewide_mutation_overlay.tsv")
