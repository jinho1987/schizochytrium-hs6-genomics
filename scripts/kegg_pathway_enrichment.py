import csv
from math import comb, lgamma, log, exp
from collections import defaultdict

TABLES = r"C:\Users\user\Desktop\HS6 vs. 4_RNA-seq_260720\RNA-seq_analysis_results\tables"
KO_MAP = "/mnt/d/HS6_vs_7_comparison/kegg/hs6_gene_to_ko.tsv"  # gene_id, kegg_gene_id, ko

# load gene->KO
gene_ko = {}
with open(KO_MAP.replace("/mnt/d/", "D:\\").replace("/", "\\"), encoding="utf-8") as f:
    reader = csv.reader(f, delimiter="\t")
    next(reader)
    for r in reader:
        if len(r) >= 3:
            gid = r[0].replace(".t1", "")
            gene_ko[gid] = r[2]

print(f"Genes with a KO assignment: {len(gene_ko)}")

# load DE results
sig = set()
tested = set()
with open(f"{TABLES}\\DE_genotype_MTvsWT_adj_for_time.tsv", encoding="utf-8") as f:
    reader = csv.reader(f, delimiter="\t")
    next(reader)
    for r in reader:
        if not r or not r[0]:
            continue
        gid = r[0].replace(".t1", "")
        try:
            padj = float(r[6]) if r[6] not in ("", "NA") else None
        except (ValueError, IndexError):
            padj = None
        if padj is not None:
            tested.add(gid)
            if padj < 0.05:
                sig.add(gid)

print(f"Tested genes: {len(tested)}, significant: {len(sig)}")

# map KOs to pathway maps used throughout this project (the 8 curated maps)
# use a lightweight static KO->pathway membership fetched via known map IDs already used
PATHWAY_KOS = {}  # will populate from local xml files already downloaded
import re
KEGG_DIR = r"D:\HS6_vs_7_comparison\kegg"
import os
pathway_names = {
    "ko00061": "Fatty acid biosynthesis", "ko00062": "Fatty acid elongation",
    "ko00071": "Fatty acid degradation", "ko00561": "Glycerolipid metabolism",
    "ko00564": "Glycerophospholipid metabolism", "ko00900": "Terpenoid backbone biosynthesis (MVA)",
    "ko00906": "Carotenoid biosynthesis", "ko01040": "Biosynthesis of unsaturated fatty acids",
}
for pid, pname in pathway_names.items():
    xmlpath = os.path.join(KEGG_DIR, pid + ".xml")
    if not os.path.exists(xmlpath):
        continue
    with open(xmlpath, encoding="utf-8") as f:
        xml = f.read()
    kos = set()
    for m in re.finditer(r'<entry\b[^>]*>', xml):
        tag = m.group(0)
        if 'type="ortholog"' not in tag:
            continue
        nm = re.search(r'name="([^"]+)"', tag)
        if nm:
            for tok in nm.group(1).split():
                if tok.startswith("ko:"):
                    kos.add(tok.replace("ko:", ""))
    PATHWAY_KOS[pid] = (pname, kos)

print(f"\nLoaded {len(PATHWAY_KOS)} pathway KO sets")

# hypergeometric test: for each pathway, how many tested/significant genes map into it
def log_comb(a, b):
    if b < 0 or b > a:
        return float("-inf")
    return lgamma(a + 1) - lgamma(b + 1) - lgamma(a - b + 1)

def hypergeom_sf(k, N, K, n):
    # P(X >= k) for hypergeometric, computed in log-space to avoid overflow with large N
    log_denom = log_comb(N, n)
    log_terms = []
    for i in range(k, min(K, n) + 1):
        lt = log_comb(K, i) + log_comb(N - K, n - i) - log_denom
        log_terms.append(lt)
    if not log_terms:
        return 1.0
    m = max(log_terms)
    s = sum(exp(t - m) for t in log_terms)
    return exp(m) * s if m > -700 else 0.0

N = len(tested)
n = len(sig)
print(f"\n=== KEGG pathway enrichment (hypergeometric, N={N} tested, n={n} significant) ===")
results = []
for pid, (pname, kos) in PATHWAY_KOS.items():
    pathway_genes_tested = set(g for g in tested if gene_ko.get(g) in kos)
    pathway_genes_sig = set(g for g in sig if gene_ko.get(g) in kos)
    K = len(pathway_genes_tested)
    k = len(pathway_genes_sig)
    if K == 0:
        continue
    pval = hypergeom_sf(k, N, K, n) if k > 0 else 1.0
    results.append((pid, pname, k, K, pval, pathway_genes_sig))

results.sort(key=lambda x: x[4])
print(f"{'Pathway':10s} {'Name':45s} {'sig/tested':>12s} {'p-value':>12s}")
for pid, pname, k, K, pval, genes in results:
    print(f"{pid:10s} {pname:45s} {f'{k}/{K}':>12s} {pval:12.4f}")
    if k > 0:
        print(f"           genes: {sorted(genes)}")
