import csv, statistics, math

TABLES = r"C:\Users\user\Desktop\HS6 vs. 4_RNA-seq_260720\RNA-seq_analysis_results\tables"

# load full DE table
rows = []
with open(f"{TABLES}\\DE_genotype_MTvsWT_adj_for_time.tsv", encoding="utf-8") as f:
    reader = csv.reader(f, delimiter="\t")
    header = next(reader)
    for r in reader:
        if not r or not r[0]:
            continue
        gid = r[0].replace(".t1", "")
        try:
            lfc = float(r[2]) if r[2] not in ("", "NA") else None
            pval = float(r[5]) if r[5] not in ("", "NA") else None
            padj = float(r[6]) if r[6] not in ("", "NA") else None
        except (ValueError, IndexError):
            lfc, pval, padj = None, None, None
        protein = r[7] if len(r) > 7 else ""
        rows.append({"gid": gid, "lfc": lfc, "pval": pval, "padj": padj, "protein": protein})

print(f"Total genes tested: {len(rows)}")

# ============ 1) p-value calibration check ============
pvals = [r["pval"] for r in rows if r["pval"] is not None]
print(f"\n=== P-value calibration check ===")
print(f"Genes with a p-value: {len(pvals)}")
bins = [0]*10
for p in pvals:
    idx = min(int(p*10), 9)
    bins[idx] += 1
print("Histogram (bin = p-value decile, 0.0-0.1 ... 0.9-1.0):")
for i, c in enumerate(bins):
    frac = c/len(pvals)
    bar = "#" * int(frac*100)
    print(f"  [{i/10:.1f}-{(i+1)/10:.1f}) {c:5d} ({frac*100:5.1f}%) {bar}")
uniform_expected = len(pvals)/10
excess_near_zero = bins[0] / uniform_expected
print(f"Enrichment in [0.0-0.1) bin vs. uniform-null expectation: {excess_near_zero:.2f}x")
# flat tail check (bins 5-9 should be close to uniform if test is well-calibrated)
tail_bins = bins[5:]
tail_mean = statistics.mean(tail_bins)
tail_cv = statistics.pstdev(tail_bins)/tail_mean if tail_mean else float('inf')
print(f"Tail bins [0.5-1.0) mean={tail_mean:.0f}, CV={tail_cv:.3f} (low CV = well-calibrated null)")

# ============ 2) annotation completeness audit ============
sig = [r for r in rows if r["padj"] is not None and r["padj"] < 0.05]
print(f"\n=== Annotation completeness among {len(sig)} significant DE genes ===")
def is_annotated(protein):
    if not protein or protein.strip() in ("", "NA"):
        return False
    low = protein.lower()
    if "hypothetical" in low or "uncharacterized" in low or "unknown" in low:
        return False
    return True

n_annotated = sum(1 for r in sig if is_annotated(r["protein"]))
print(f"Annotated (real protein name, not hypothetical/uncharacterized): {n_annotated} / {len(sig)} ({100*n_annotated/len(sig):.1f}%)")
print(f"Unannotated / hypothetical / uncharacterized: {len(sig)-n_annotated} / {len(sig)} ({100*(len(sig)-n_annotated)/len(sig):.1f}%)")

# same audit for ALL tested genes (background rate) for comparison
n_annotated_all = sum(1 for r in rows if is_annotated(r["protein"]))
print(f"(Background rate, all {len(rows)} tested genes: {100*n_annotated_all/len(rows):.1f}% annotated)")

# ============ 3) effect size distribution ============
print(f"\n=== Effect size (log2FC) distribution among {len(sig)} significant DE genes ===")
lfcs = [r["lfc"] for r in sig if r["lfc"] is not None]
abs_lfcs = [abs(x) for x in lfcs]
bins_lfc = {"0.5-1": 0, "1-2": 0, "2-3": 0, "3-5": 0, ">5": 0}
for a in abs_lfcs:
    if a < 1: bins_lfc["0.5-1"] += 1
    elif a < 2: bins_lfc["1-2"] += 1
    elif a < 3: bins_lfc["2-3"] += 1
    elif a < 5: bins_lfc["3-5"] += 1
    else: bins_lfc[">5"] += 1
for k, v in bins_lfc.items():
    print(f"  |log2FC| {k:6s}: {v:5d} ({100*v/len(abs_lfcs):5.1f}%)")
print(f"Median |log2FC|: {statistics.median(abs_lfcs):.2f}")
print(f"Max |log2FC|: {max(abs_lfcs):.2f} (gene: {sig[abs_lfcs.index(max(abs_lfcs))]['gid']})")

# save volcano plot data
with open(r"C:\Users\user\AppData\Local\Temp\claude\C--Users-user-Desktop-HS6-vs--4-RNA-seq-260720\6e32b67b-d05a-4791-b9c3-db6271bce82e\scratchpad\volcano_data.tsv", "w", encoding="utf-8") as out:
    out.write("gene\tlog2FC\tpadj\tneg_log10_padj\tsignificant\tprotein\n")
    for r in rows:
        if r["lfc"] is None or r["padj"] is None:
            continue
        nlp = -math.log10(r["padj"]) if r["padj"] > 0 else 300
        sig_flag = "yes" if r["padj"] < 0.05 else "no"
        out.write(f"{r['gid']}\t{r['lfc']:.4f}\t{r['padj']:.6e}\t{nlp:.3f}\t{sig_flag}\t{r['protein']}\n")
print("\nWrote volcano_data.tsv for plotting")
