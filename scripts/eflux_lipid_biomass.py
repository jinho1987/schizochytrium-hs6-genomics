import cobra, csv, re

MODEL = "/mnt/d/HS6_vs_7_comparison/fba/hs6_draft_model.xml"
COUNTS = "/mnt/c/Users/user/Desktop/HS6 vs. 4_RNA-seq_260720/RNA-seq_analysis_results/tables/normalized_counts.tsv"

base_model = cobra.io.read_sbml_model(MODEL)
growth = base_model.reactions.get_by_id("Growth")

# 1) Remove the bacterial peptidoglycan precursor -- biologically absent in a eukaryotic alga
if "uaagmda_c" in [m.id for m in growth.metabolites]:
    met = base_model.metabolites.get_by_id("uaagmda_c")
    coef = growth.metabolites[met]
    growth.add_metabolites({met: -coef})  # zero it out
    print(f"Removed uaagmda_c (bacterial peptidoglycan precursor), was coefficient {coef}")

# 2) Reduce amino acid coefficients by 45% (making room for a lipid fraction --
#    Schizochytrium is reported at 50-77% dry-weight lipid under lipid-accumulation
#    conditions in the literature; the untouched model has ~0% lipid weighting)
aa_ids = ["ala__L_c","arg__L_c","asn__L_c","asp__L_c","cys__L_c","gln__L_c","glu__L_c",
          "gly_c","his__L_c","ile__L_c","leu__L_c","lys__L_c","met__L_c","phe__L_c",
          "pro__L_c","ser__L_c","thr__L_c","trp__L_c","tyr__L_c","val__L_c"]
reduced_mass = 0
for aid in aa_ids:
    try:
        m = base_model.metabolites.get_by_id(aid)
    except KeyError:
        continue
    if m in growth.metabolites:
        coef = growth.metabolites[m]
        delta = -coef * 0.45  # remove 45% of the (negative) coefficient -> less negative
        growth.add_metabolites({m: delta})
        reduced_mass += abs(coef * 0.45)
print(f"Reduced amino-acid coefficients by 45% (freed magnitude ~{reduced_mass:.3f})")

# 3) Add fatty-acid biomass terms using the freed magnitude, weighted toward
#    the saturated/monounsaturated species typically dominant in thraustochytrid
#    storage lipid (palmitate C16:0 is usually the largest single SFA fraction;
#    true DHA (C22:6) has no BiGG-universal metabolite ID available in this draft
#    model, a real limitation of using a generic reconstruction database for this
#    organism -- noted explicitly in the writeup).
fa_weights = {"hdca_c": 0.40, "ocdcea_c": 0.30, "ocdca_c": 0.15, "hdcea_c": 0.10, "ttdca_c": 0.05}
for fid, w in fa_weights.items():
    m = base_model.metabolites.get_by_id(fid)
    growth.add_metabolites({m: -reduced_mass * w})
print("Added fatty-acid biomass terms:", fa_weights)

cobra.io.write_sbml_model(base_model, "/mnt/d/HS6_vs_7_comparison/fba/hs6_lipid_model.xml")
print("\nSaved corrected model to hs6_lipid_model.xml")

# sanity check: default FBA on corrected model
sol = base_model.optimize()
print(f"Default (unconstrained) FBA growth objective on corrected model: {sol.objective_value:.3f}")

# ==== now rerun E-Flux on the corrected model ====
def strip_suffix(gid):
    return re.sub(r"_t\d+$", "", gid)

wt_expr, mt_expr = {}, {}
with open(COUNTS, encoding="utf-8") as f:
    reader = csv.reader(f, delimiter="\t")
    next(reader)
    for row in reader:
        gid = row[0]
        vals = [float(x) for x in row[1:]]
        wt_expr[gid] = sum(vals[0:3]) / 3
        mt_expr[gid] = sum(vals[3:6]) / 3

def eflux_bound(model, expr_dict, max_bound=1000):
    gene_vals = {g.id: expr_dict.get(strip_suffix(g.id), 0) for g in model.genes}
    global_max = max(gene_vals.values()) if gene_vals else 1
    bounds = {}
    for rxn in model.reactions:
        if not rxn.genes:
            continue
        vals = [gene_vals.get(g.id, 0) for g in rxn.genes]
        scale = max(vals) / global_max if global_max > 0 else 0
        bounds[rxn.id] = scale * max_bound
    return bounds

print("\n=== Expression-constrained FBA on the lipid-corrected model ===")
for label, expr in [("WT4", wt_expr), ("HS6", mt_expr)]:
    m = base_model.copy()
    bounds = eflux_bound(m, expr)
    for rxn in m.reactions:
        if rxn.id in bounds:
            b = max(bounds[rxn.id], 1e-3)
            if rxn.upper_bound > 0:
                rxn.upper_bound = b
            if rxn.lower_bound < 0:
                rxn.lower_bound = -b
    sol = m.optimize()
    print(f"\n{label}: growth objective (lipid-corrected biomass) = {sol.objective_value:.4f}")
    fluxes = sol.fluxes.abs().sort_values(ascending=False)
    print("  Top 8 highest-flux reactions:")
    for rxn_id, flux in fluxes.head(8).items():
        try:
            name = m.reactions.get_by_id(rxn_id).name
        except Exception:
            name = ""
        print(f"    {rxn_id:20s} flux={flux:10.4f}  {name}")
