import json

d = json.load(open("results/lysis_enzymes/docking_results.json"))
by_label = {}
for r in d:
    by_label.setdefault(r["label"], {})[r["ligand"]] = r

rows = []
for label, ligs in by_label.items():
    canon = [v for k, v in ligs.items() if k != "galactan_probe"][0]
    gal = ligs["galactan_probe"]
    gap = gal["top_affinity_kcal_mol"] - canon["top_affinity_kcal_mol"]
    rows.append((gap, label, canon["accession"], canon["class"], canon["host"],
                 canon["top_affinity_kcal_mol"], gal["top_affinity_kcal_mol"]))

rows.sort()
print(f"{'gap':>7} {'label':25} {'acc':8} {'canon':>8} {'galactan':>9}  class / host")
for gap, label, acc, klass, host, canon_aff, gal_aff in rows:
    flag = " <-- REVERSED" if gap < 0 else ""
    print(f"{gap:7.3f} {label:25} {acc:8} {canon_aff:8.3f} {gal_aff:9.3f}  {klass} / {host}{flag}")
