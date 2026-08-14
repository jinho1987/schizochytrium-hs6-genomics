import cobra
import csv

MODEL = "/mnt/d/HS6_vs_7_comparison/fba/hs6_draft_model.xml"
COUNTS = "/mnt/c/Users/user/Desktop/HS6 vs. 4_RNA-seq_260720/RNA-seq_analysis_results/tables/normalized_counts.tsv"

model = cobra.io.read_sbml_model(MODEL)
print(f"Model: {len(model.reactions)} reactions, {len(model.genes)} genes")

# load normalized expression, average WT vs MT separately
wt_expr, mt_expr = {}, {}
with open(COUNTS, encoding="utf-8") as f:
    reader = csv.reader(f, delimiter="\t")
    header = next(reader)
    for row in reader:
        gid = row[0]
        vals = [float(x) for x in row[1:]]
        wt_expr[gid] = sum(vals[0:3]) / 3
        mt_expr[gid] = sum(vals[3:6]) / 3

import re
def strip_suffix(gid):
    return re.sub(r"_t\d+$", "", gid)

def eflux_bound(model, expr_dict, max_bound=1000):
    """Simple E-Flux: scale each reaction's upper bound by the max expression
    of its associated genes, normalized to the global max."""
    gene_vals = {g.id: expr_dict.get(strip_suffix(g.id), 0) for g in model.genes}
    n_matched = sum(1 for g in model.genes if strip_suffix(g.id) in expr_dict)
    print(f"  matched {n_matched} / {len(model.genes)} model genes to expression table")
    global_max = max(gene_vals.values()) if gene_vals else 1
    bounds = {}
    for rxn in model.reactions:
        if not rxn.genes:
            continue
        vals = [gene_vals.get(g.id, 0) for g in rxn.genes]
        if not vals:
            continue
        scale = max(vals) / global_max if global_max > 0 else 0
        bounds[rxn.id] = scale * max_bound
    return bounds

for label, expr in [("WT4", wt_expr), ("HS6", mt_expr)]:
    m = model.copy()
    bounds = eflux_bound(m, expr)
    for rxn in m.reactions:
        if rxn.id in bounds:
            b = max(bounds[rxn.id], 1e-3)  # avoid fully zeroing reactions (E-Flux convention: small min flux allowed)
            if rxn.upper_bound > 0:
                rxn.upper_bound = b
            if rxn.lower_bound < 0:
                rxn.lower_bound = -b
    sol = m.optimize()
    print(f"\n{label}: expression-constrained FBA growth objective = {sol.objective_value:.3f}")
    # report top 15 most active lipid/isoprenoid-relevant reactions by flux magnitude
    fluxes = sol.fluxes.abs().sort_values(ascending=False)
    print(f"  Top 10 highest-flux reactions:")
    for rxn_id, flux in fluxes.head(10).items():
        try:
            name = m.reactions.get_by_id(rxn_id).name
        except Exception:
            name = ""
        print(f"    {rxn_id:20s} flux={flux:10.3f}  {name}")

print("\n=== DONE: crude expression-constrained flux comparison ===")
