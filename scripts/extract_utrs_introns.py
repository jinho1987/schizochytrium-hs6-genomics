import re, csv

GFF = r"C:\Users\user\Desktop\HS6 vs. 4_RNA-seq_260720\HS6_Genome\braker.gff3"
STRINGTIE_GTF = "/mnt/d/HS6_vs_7_comparison/stringtie/merged.stringtie.gtf"  # will read via wsl-mounted path from Windows too
STRINGTIE_GTF_WIN = r"D:\HS6_vs_7_comparison\stringtie\merged.stringtie.gtf"
TABLES = r"C:\Users\user\Desktop\HS6 vs. 4_RNA-seq_260720\RNA-seq_analysis_results\tables"
FASTA = r"C:\Users\user\Desktop\HS6 vs. 4_RNA-seq_260720\HS6_Genome\HS6_scaf.fasta"

# 1) CDS extents per transcript from braker (gives us the start-codon / stop-codon boundary)
cds = {}  # tx_id -> (scaf, cds_min, cds_max, strand)
with open(GFF, encoding="utf-8") as f:
    for line in f:
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 9 or p[2] != "CDS":
            continue
        m = re.search(r"Parent=([^;]+)", p[8])
        if not m:
            continue
        tx = m.group(1)
        scaf, s, e, strand = p[0], int(p[3]), int(p[4]), p[6]
        if tx not in cds:
            cds[tx] = [scaf, s, e, strand]
        else:
            cds[tx][1] = min(cds[tx][1], s)
            cds[tx][2] = max(cds[tx][2], e)

print(f"Transcripts with CDS coords: {len(cds)}")

# 2) StringTie exon structure, grouped by ref_gene_id (the braker gene it matched to) and transcript_id
tx_exons = {}  # stringtie transcript_id -> list of (start,end)
tx_meta = {}   # stringtie transcript_id -> (scaf, strand, ref_gene_id)
with open(STRINGTIE_GTF_WIN, encoding="utf-8") as f:
    for line in f:
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 9 or p[2] != "exon":
            continue
        m_tx = re.search(r'transcript_id "([^"]+)"', p[8])
        m_ref = re.search(r'ref_gene_id "([^"]+)"', p[8])
        if not m_tx or not m_ref:
            continue
        tx_id = m_tx.group(1)
        ref_gene = m_ref.group(1)
        scaf, s, e, strand = p[0], int(p[3]), int(p[4]), p[6]
        tx_exons.setdefault(tx_id, []).append((s, e))
        tx_meta[tx_id] = (scaf, strand, ref_gene)

print(f"StringTie transcripts with ref_gene_id match: {len(tx_meta)}")

# 3) load promoter candidate ranking (already computed) to prioritize
top_genes = []
with open(r"C:\Users\user\AppData\Local\Temp\claude\C--Users-user-Desktop-HS6-vs--4-RNA-seq-260720\6e32b67b-d05a-4791-b9c3-db6271bce82e\scratchpad\promoter_candidates_ranked.tsv", encoding="utf-8") as f:
    next(f)
    for line in f:
        top_genes.append(line.split("\t")[0])

# 4) compute UTR/intron structure for each candidate (using its matching StringTie transcript, best = ref_gene_id match with a CDS)
def analyze(gid):
    # find a stringtie transcript matching this gene as ref_gene_id AND matching braker tx for CDS coords
    braker_tx = f"{gid}.t1"
    if braker_tx not in cds:
        return None
    scaf_c, cds_min, cds_max, strand = cds[braker_tx]

    candidates = [tx for tx, (scaf, s, rg) in tx_meta.items() if rg == gid]
    if not candidates:
        return None
    # pick the stringtie transcript with the most exons (most informative) among matches
    best_tx = max(candidates, key=lambda t: len(tx_exons[t]))
    exons = sorted(tx_exons[best_tx])
    tx_start, tx_end = exons[0][0], exons[-1][1]

    introns = []
    for i in range(len(exons) - 1):
        istart = exons[i][1] + 1
        iend = exons[i + 1][0] - 1
        if iend >= istart:
            introns.append((istart, iend))

    if strand == "+":
        utr5 = (tx_start, cds_min - 1) if tx_start < cds_min else None
        utr3 = (cds_max + 1, tx_end) if tx_end > cds_max else None
    else:
        utr5 = (cds_max + 1, tx_end) if tx_end > cds_max else None
        utr3 = (tx_start, cds_min - 1) if tx_start < cds_min else None

    utr5_introns = []
    if utr5:
        for istart, iend in introns:
            if istart >= utr5[0] and iend <= utr5[1]:
                utr5_introns.append((istart, iend))

    return {
        "gid": gid, "tx": best_tx, "scaf": scaf_c, "strand": strand,
        "n_exons": len(exons), "n_introns": len(introns),
        "utr5": utr5, "utr5_len": (utr5[1]-utr5[0]+1) if utr5 else 0,
        "utr3": utr3, "utr3_len": (utr3[1]-utr3[0]+1) if utr3 else 0,
        "utr5_introns": utr5_introns,
        "first_intron": introns[0] if introns else None,
    }

results = []
for gid in top_genes:  # scan the FULL stable/robustly-expressed candidate pool (~9,236 genes)
    r = analyze(gid)
    if r:
        results.append(r)

print(f"\nAnalyzed {len(results)} / {len(top_genes)} stable-expression candidates with usable StringTie structure")
print(f"\nTop 20 by expression (for reference):")
print(f"{'gene':10s} {'exons':>6s} {'introns':>8s} {'5UTR_len':>9s} {'3UTR_len':>9s} {'5UTR_introns':>13s}")
for r in results[:20]:
    print(f"{r['gid']:10s} {r['n_exons']:6d} {r['n_introns']:8d} {r['utr5_len']:9d} {r['utr3_len']:9d} {len(r['utr5_introns']):13d}")

# highlight IME candidates: genes with an intron INSIDE the 5'UTR (classic intron-mediated enhancement architecture)
ime_candidates = [r for r in results if r["utr5_introns"]]
print(f"\n=== Intron-mediated-enhancement (IME) candidates: intron sitting inside the 5'UTR ===")
for r in ime_candidates:
    print(f"  {r['gid']}: 5'UTR={r['utr5_len']}bp, {len(r['utr5_introns'])} intron(s) inside it, first={r['utr5_introns'][0]}")

with open(r"C:\Users\user\AppData\Local\Temp\claude\C--Users-user-Desktop-HS6-vs--4-RNA-seq-260720\6e32b67b-d05a-4791-b9c3-db6271bce82e\scratchpad\utr_intron_candidates.tsv", "w", encoding="utf-8") as out:
    out.write("gene\tscaf\tstrand\tn_exons\tn_introns\tutr5_start\tutr5_end\tutr5_len\tutr3_start\tutr3_end\tutr3_len\tn_utr5_introns\tfirst_intron_start\tfirst_intron_end\n")
    for r in results:
        u5 = r["utr5"] or ("", "")
        u3 = r["utr3"] or ("", "")
        fi = r["first_intron"] or ("", "")
        out.write(f"{r['gid']}\t{r['scaf']}\t{r['strand']}\t{r['n_exons']}\t{r['n_introns']}\t{u5[0]}\t{u5[1]}\t{r['utr5_len']}\t{u3[0]}\t{u3[1]}\t{r['utr3_len']}\t{len(r['utr5_introns'])}\t{fi[0]}\t{fi[1]}\n")
print("\nwrote utr_intron_candidates.tsv")
