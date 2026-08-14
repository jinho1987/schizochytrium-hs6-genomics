import cobra
model = cobra.io.read_sbml_model("/mnt/d/HS6_vs_7_comparison/fba/hs6_draft_model.xml")
ids = [g.id for g in model.genes]
print("Sample gene IDs from model:", ids[:10])
print("Total genes:", len(ids))
