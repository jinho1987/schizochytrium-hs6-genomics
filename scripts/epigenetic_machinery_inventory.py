"""Inventory of DNA-methylation and chromatin-machinery genes in the HS6 gene set,
joined to the HS6-vs-#4 differential-expression table.

Run from the repository root:
    python3 scripts/epigenetic_machinery_inventory.py

Inputs  (both already in this repo)
    results/eggnog/hs6_eggnog.emapper.annotations.tsv   eggNOG-mapper output
    results/volcano_data.tsv                            HS6 vs #4, adjusted for time (+ve = higher in HS6)
Output
    results/epigenetic_machinery_inventory.tsv          one row per gene in the set

Matching is by exact Pfam family name (comma-separated PFAMs column), never by free-text
description: description matching picked up DnaJ_CXXCXGXG as "CXXC" and similar false hits.

Limit: eggNOG-mapper annotated only part of the 16,088 gene models. A family that is
"not detected" here is absent from the annotated subset, not proven absent from the genome.
"""
import math
import sys

EGGNOG = sys.argv[1] if len(sys.argv) > 1 else "results/eggnog/hs6_eggnog.emapper.annotations.tsv"
VOLCANO = sys.argv[2] if len(sys.argv) > 2 else "results/volcano_data.tsv"
OUT = sys.argv[3] if len(sys.argv) > 3 else "results/epigenetic_machinery_inventory.tsv"
N_GENE_MODELS = 16088  # BRAKER2 gene models (see index.md)

# category -> Pfam family names. Chosen before looking at expression.
CATEGORIES = {
    "DNA C5-MTase fold": {"DNA_methylase"},
    "MT-A70 (METTL3/14, MTA1 family)": {"MT-A70"},
    "N6/N4 adenine-cytosine MTase": {"N6_N4_Mtase", "N6-adenineMlase"},
    "PWWP": {"PWWP"},
    "MBD": {"MBD"},
    "SAD_SRA (UHRF/VIM reader)": {"SAD_SRA"},
    "zf-CXXC": {"zf-CXXC"},
    "TET/JBP oxidase": {"Tet_JBP"},
    "HhH-GPD glycosylase": {"HhH-GPD"},
    "AlkB/ALKBH": {"2OG-FeII_Oxy_2"},
    "SET (histone KMT)": {"SET"},
    "DOT1": {"DOT1"},
    "JmjC": {"JmjC"},
    "HDAC": {"Hist_deacetyl"},
    "Sirtuin": {"SIR2"},
    "MYST/p300 HAT": {"MOZ_SAS", "HAT_KAT11"},
    "Chromo": {"Chromo"},
    "Bromodomain": {"Bromodomain"},
    "SNF2 remodeler": {"SNF2_N"},
    "Core histone fold": {"Histone"},
}


def gene_of(query):
    return query.split(".t")[0]


rows, header = [], None
with open(EGGNOG, encoding="utf-8", errors="replace") as f:
    for line in f:
        if line.startswith("##"):
            continue
        if line.startswith("query\t"):
            header = line.rstrip("\n").split("\t")
        elif header:
            rows.append(dict(zip(header, line.rstrip("\n").split("\t"))))
annotated = {gene_of(r["query"]) for r in rows}

de = {}
with open(VOLCANO, encoding="utf-8") as f:
    next(f)
    for line in f:
        g, lfc, padj, _, sig, _prot = line.rstrip("\n").split("\t")[:6]
        de[g] = (float(lfc), float(padj), sig == "yes")
n_de, n_sig = len(de), sum(1 for v in de.values() if v[2])

members = {}  # gene -> (categories, eggNOG row)
for r in rows:
    pf = {t.strip() for t in r["PFAMs"].split(",") if t.strip() not in ("", "-")}
    cats = [c for c, names in CATEGORIES.items() if pf & names]
    if cats:
        members[gene_of(r["query"])] = (cats, r)


def gnum(g):
    return int(g[1:]) if g[1:].isdigit() else 0


with open(OUT, "w", encoding="utf-8") as out:
    out.write("gene\tcategories\tpfams\tdescription\tpreferred_name\tko\tbest_hit_level\tevalue\t"
              "log2FC_HS6_vs_4\tpadj\tsignificant_padj05\n")
    for g in sorted(members, key=gnum):
        cats, r = members[g]
        lfc, padj, sig = de.get(g, (None, None, None))
        out.write("\t".join([
            g, "; ".join(cats), r["PFAMs"], r["Description"], r["Preferred_name"], r["KEGG_ko"],
            r["max_annot_lvl"], r["evalue"],
            "" if lfc is None else f"{lfc:.4f}", "" if padj is None else f"{padj:.3e}",
            "" if sig is None else ("yes" if sig else "no")]) + "\n")

print(f"gene models: {N_GENE_MODELS};  eggNOG-annotated genes: {len(annotated)} "
      f"({100 * len(annotated) / N_GENE_MODELS:.1f}%)")
print(f"genes with HS6-vs-#4 statistics: {n_de};  padj<0.05: {n_sig} ({100 * n_sig / n_de:.1f}%)")
print("\ngenes per category (unique genes):")
for c in CATEGORIES:
    gl = sorted((g for g, (cs, _) in members.items() if c in cs), key=gnum)
    print(f"  {len(gl):>3}  {c}" + (f"  -> {', '.join(gl)}" if 0 < len(gl) <= 8 else ""))

tested = [g for g in members if g in de]
sig = [g for g in tested if de[g][2]]
k, n = len(sig), len(tested)
p_ge = sum(math.comb(n_sig, i) * math.comb(n_de - n_sig, n - i)
           for i in range(k, min(n_sig, n) + 1)) / math.comb(n_de, n)
print(f"\nmachinery set: {len(members)} genes, {n} with statistics, {k} significant "
      f"({100 * k / n:.1f}% vs {100 * n_sig / n_de:.1f}% genome-wide; expected {n * n_sig / n_de:.1f})")
print(f"hypergeometric P(>= {k} significant) = {p_ge:.3f}   (high P = no enrichment)")
print("\nsignificant members (+ve log2FC = higher in HS6):")
for g in sorted(sig, key=lambda x: de[x][1]):
    print(f"  {g:<7} log2FC={de[g][0]:+.2f}  padj={de[g][1]:.1e}  [{'; '.join(members[g][0])}]  "
          f"{members[g][1]['Description'][:55]}")
print(f"\nwrote {OUT}")
