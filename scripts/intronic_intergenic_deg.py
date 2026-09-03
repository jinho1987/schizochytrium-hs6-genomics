import csv

VCF_DUMP = r"C:\Users\user\AppData\Local\Temp\claude\C--Users-user-Desktop-HS6-vs--4-RNA-seq-260720\6e32b67b-d05a-4791-b9c3-db6271bce82e\scratchpad\all_84_variants.tsv"
GFF = r"C:\Users\user\Desktop\HS6 vs. 4_RNA-seq_260720\HS6_Genome\braker.gff3"
DE_TABLE = r"C:\Users\user\Desktop\HS6 vs. 4_RNA-seq_260720\RNA-seq_analysis_results\tables\DE_genotype_MTvsWT_adj_for_time.tsv"

# load gene body (min/max) and CDS/exon intervals per gene
genes = {}   # gid -> [chrom, start, end, strand]
exons = {}   # gid -> list of (start,end)  (CDS features, proxy for exonic space)
with open(GFF, encoding="utf-8") as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) < 9:
            continue
        chrom, feat, start, end, strand, attrs = p[0], p[2], int(p[3]), int(p[4]), p[6], p[8]
        if feat == "gene" and attrs.startswith("ID="):
            gid = attrs[3:].rstrip(";")
            genes[gid] = [chrom, start, end, strand]
        elif feat == "CDS":
            # Parent=g1.t1 or similar
            parent = None
            for tok in attrs.split(";"):
                if tok.startswith("Parent="):
                    parent = tok[len("Parent="):]
            if parent:
                gid = parent.split(".")[0]
                exons.setdefault(gid, []).append((chrom, start, end))

de = {}
with open(DE_TABLE, encoding="utf-8") as f:
    reader = csv.reader(f, delimiter="\t")
    next(reader)
    for r in reader:
        if not r or not r[0]:
            continue
        gid = r[0]
        try:
            padj = float(r[6]) if r[6] not in ("", "NA") else None
        except (ValueError, IndexError):
            padj = None
        try:
            lfc = float(r[2])
        except (ValueError, IndexError):
            lfc = None
        de[gid] = (padj, lfc, r[7] if len(r) > 7 else "")

# sort genes by chrom+start for nearest-gene lookup
by_chrom = {}
for gid, (chrom, start, end, strand) in genes.items():
    by_chrom.setdefault(chrom, []).append((start, end, gid))
for chrom in by_chrom:
    by_chrom[chrom].sort()

def classify_and_find_genes(chrom, pos):
    hits = []
    for gid, (c, s, e, strand) in genes.items():
        if c == chrom and s <= pos <= e:
            # inside gene body -- is it in an exon/CDS?
            in_cds = any(cc == chrom and cs <= pos <= ce for cc, cs, ce in exons.get(gid, []))
            hits.append((gid, "CDS/exon" if in_cds else "intron", 0))
    if hits:
        return hits
    # intergenic: find nearest gene(s) up/downstream
    candidates = by_chrom.get(chrom, [])
    nearest = []
    best_dist = None
    for s, e, gid in candidates:
        if pos < s:
            d = s - pos
        elif pos > e:
            d = pos - e
        else:
            d = 0
        if best_dist is None or d < best_dist:
            best_dist = d
            nearest = [(gid, "intergenic", d)]
        elif d == best_dist:
            nearest.append((gid, "intergenic", d))
    return nearest

variants = []
with open(VCF_DUMP, encoding="utf-8") as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) < 4 or not p[1].strip().isdigit():
            continue
        chrom, pos, ref, alt = p[0], int(p[1]), p[2], p[3]
        variants.append((chrom, pos, ref, alt))

print(f"loaded {len(variants)} validated variants\n")

rows = []
for chrom, pos, ref, alt in variants:
    hits = classify_and_find_genes(chrom, pos)
    for gid, loc, dist in hits:
        padj, lfc, ident = de.get(gid, (None, None, ""))
        rows.append({
            "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
            "gene": gid, "location": loc, "distance_bp": dist,
            "padj": padj, "log2fc": lfc, "identity": ident,
        })

rows.sort(key=lambda r: (r["padj"] is None, r["padj"] if r["padj"] is not None else 2))

print(f"{'chrom':12s} {'pos':>10s} {'gene':10s} {'location':12s} {'dist':>6s} {'padj':>10s} {'log2FC':>8s} identity")
for r in rows:
    padj_s = f"{r['padj']:.4f}" if r["padj"] is not None else "NA"
    lfc_s = f"{r['log2fc']:.2f}" if r["log2fc"] is not None else "NA"
    print(f"{r['chrom']:12s} {r['pos']:>10d} {r['gene']:10s} {r['location']:12s} {r['distance_bp']:>6d} {padj_s:>10s} {lfc_s:>8s} {r['identity']}")

sig = [r for r in rows if r["padj"] is not None and r["padj"] < 0.05]
print(f"\n=== {len(sig)} variant-gene links where the gene is DE-significant (padj<0.05) ===")
for r in sig:
    print(f"  {r['gene']} ({r['location']}, {r['distance_bp']}bp): padj={r['padj']:.4f} log2FC={r['log2fc']:.2f} {r['identity']}")

OUT = r"C:\Users\user\AppData\Local\Temp\claude\C--Users-user-Desktop-HS6-vs--4-RNA-seq-260720\6e32b67b-d05a-4791-b9c3-db6271bce82e\scratchpad\intronic_intergenic_deg_links.tsv"
with open(OUT, "w") as f:
    f.write("chrom\tpos\tref\talt\tgene\tlocation\tdistance_bp\tpadj\tlog2FC\tidentity\n")
    for r in rows:
        f.write(f"{r['chrom']}\t{r['pos']}\t{r['ref']}\t{r['alt']}\t{r['gene']}\t{r['location']}\t{r['distance_bp']}\t{r['padj']}\t{r['log2fc']}\t{r['identity']}\n")
print(f"\nwrote {len(rows)} rows to {OUT}")
