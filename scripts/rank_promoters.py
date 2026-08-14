import csv, statistics

TABLES = r"C:\Users\user\Desktop\HS6 vs. 4_RNA-seq_260720\RNA-seq_analysis_results\tables"

# 1) load normalized counts
counts = {}
with open(f"{TABLES}\\normalized_counts.tsv", encoding="utf-8") as f:
    reader = csv.reader(f, delimiter="\t")
    header = next(reader)
    for row in reader:
        gid = row[0]
        vals = [float(x) for x in row[1:]]
        counts[gid] = vals

# 2) load DE results (genotype effect) to flag/exclude genes that differ between HS6 and WT
de = {}
with open(f"{TABLES}\\DE_genotype_MTvsWT_adj_for_time.tsv", encoding="utf-8") as f:
    reader = csv.reader(f, delimiter="\t")
    header = next(reader)
    # find log2FC and padj columns
    lfc_idx = header.index("log2FoldChange") if "log2FoldChange" in header else 2
    padj_idx = header.index("padj") if "padj" in header else 6
    for row in reader:
        if not row or not row[0]:
            continue
        gid = row[0].replace(".t1", "")
        try:
            lfc = float(row[lfc_idx]) if row[lfc_idx] not in ("", "NA") else None
            padj = float(row[padj_idx]) if row[padj_idx] not in ("", "NA") else None
        except (ValueError, IndexError):
            lfc, padj = None, None
        de[gid] = (lfc, padj)

# 3) load annotation for descriptions (swissprot hits)
desc = {}
try:
    with open(f"{TABLES}\\hs6_swissprot_hits.tsv", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        hdr = next(reader)
        for row in reader:
            if row:
                gid = row[0].replace(".t1", "")
                desc[gid] = row[-1] if len(row) > 1 else ""
except FileNotFoundError:
    pass

# 4) rank: require expressed in all 6 libraries (>0), compute mean and CV
candidates = []
for gid, vals in counts.items():
    if min(vals) <= 5:  # require robust expression in every library, not just on average
        continue
    mean = statistics.mean(vals)
    sd = statistics.pstdev(vals)
    cv = sd / mean if mean > 0 else float("inf")
    lfc, padj = de.get(gid, (None, None))
    is_de = (padj is not None and padj < 0.05)
    candidates.append((gid, mean, cv, lfc, padj, is_de))

# sort: prioritize NOT significantly DE (stable across genotype), then high mean, then low CV
non_de = [c for c in candidates if not c[5]]
non_de.sort(key=lambda c: (-c[1] / (1 + c[2]*3)))  # composite: high mean, penalize high CV

print(f"Total genes with count data: {len(counts)}")
print(f"Robustly expressed (min>5 across all 6 libs): {len(candidates)}")
print(f"...of which NOT significantly genotype-DE (stable candidates): {len(non_de)}")
print()
print("=== TOP 25 STRONG, STABLE, ENDOGENOUS PROMOTER CANDIDATES ===")
print(f"{'gene':10s} {'mean_expr':>10s} {'CV':>6s} {'log2FC':>8s} {'padj':>10s}  desc")
for gid, mean, cv, lfc, padj, is_de in non_de[:25]:
    lfc_s = f"{lfc:.2f}" if lfc is not None else "NA"
    padj_s = f"{padj:.3f}" if padj is not None else "NA"
    d = desc.get(gid, "")[:60]
    print(f"{gid:10s} {mean:10.1f} {cv:6.3f} {lfc_s:>8s} {padj_s:>10s}  {d}")

# save full ranked list
with open(r"C:\Users\user\AppData\Local\Temp\claude\C--Users-user-Desktop-HS6-vs--4-RNA-seq-260720\6e32b67b-d05a-4791-b9c3-db6271bce82e\scratchpad\promoter_candidates_ranked.tsv", "w", encoding="utf-8") as out:
    out.write("gene\tmean_expr\tCV\tlog2FC\tpadj\tdescription\n")
    for gid, mean, cv, lfc, padj, is_de in non_de:
        lfc_s = f"{lfc:.3f}" if lfc is not None else "NA"
        padj_s = f"{padj:.4f}" if padj is not None else "NA"
        out.write(f"{gid}\t{mean:.1f}\t{cv:.4f}\t{lfc_s}\t{padj_s}\t{desc.get(gid,'')}\n")
print("\nwrote full ranked list to promoter_candidates_ranked.tsv")
