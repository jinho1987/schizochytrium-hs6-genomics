import json

# pH~7 baseline from the main docking run
baseline = json.load(open("results/lysis_enzymes/docking_results.json"))
ph_data = json.load(open("results/lysis_enzymes/ph_docking_results.json"))

by_label_base = {}
for r in baseline:
    by_label_base.setdefault(r["label"], {})[r["ligand"]] = r["top_affinity_kcal_mol"]

by_label_ph = {}
for r in ph_data:
    by_label_ph.setdefault(r["label"], {}).setdefault(r["ph"], {})[r["ligand"]] = r["top_affinity_kcal_mol"]

# canonical ligand key per label
canon_key = {}
for r in baseline:
    if r["ligand"] != "galactan_probe":
        canon_key[r["label"]] = r["ligand"]

rows = []
for label in by_label_base:
    ck = canon_key[label]
    canon7 = by_label_base[label][ck]
    gal7 = by_label_base[label]["galactan_probe"]
    gap7 = gal7 - canon7

    canon4 = by_label_ph[label][4.0][ck]
    gal4 = by_label_ph[label][4.0]["galactan_probe"]
    gap4 = gal4 - canon4

    canon9 = by_label_ph[label][9.0][ck]
    gal9 = by_label_ph[label][9.0]["galactan_probe"]
    gap9 = gal9 - canon9

    gap_range = max(gap4, gap7, gap9) - min(gap4, gap7, gap9)
    canon_range = max(canon4, canon7, canon9) - min(canon4, canon7, canon9)
    gal_range = max(gal4, gal7, gal9) - min(gal4, gal7, gal9)

    rows.append({
        "label": label, "gap4": gap4, "gap7": gap7, "gap9": gap9,
        "gap_range": gap_range, "canon_range": canon_range, "gal_range": gal_range,
        "canon4": canon4, "canon7": canon7, "canon9": canon9,
        "gal4": gal4, "gal7": gal7, "gal9": gal9,
    })

rows.sort(key=lambda r: r["gap7"])
print(f"{'label':25} {'gap@4':>7} {'gap@7':>7} {'gap@9':>7} {'gap_rng':>8} {'canon_rng':>10} {'gal_rng':>8}")
for r in rows:
    print(f"{r['label']:25} {r['gap4']:7.3f} {r['gap7']:7.3f} {r['gap9']:7.3f} {r['gap_range']:8.3f} {r['canon_range']:10.3f} {r['gal_range']:8.3f}")

print("\n--- sign changes (reversed at some pH but not others) ---")
for r in rows:
    signs = set(1 if g >= 0 else -1 for g in (r["gap4"], r["gap7"], r["gap9"]))
    if len(signs) > 1:
        print(f"{r['label']}: gap4={r['gap4']:.3f} gap7={r['gap7']:.3f} gap9={r['gap9']:.3f}")

print("\n--- biggest divergence between canon_range and gal_range (mode-of-binding hint) ---")
rows2 = sorted(rows, key=lambda r: abs(r["canon_range"] - r["gal_range"]), reverse=True)
for r in rows2[:6]:
    print(f"{r['label']}: canon_range={r['canon_range']:.3f} gal_range={r['gal_range']:.3f} diff={abs(r['canon_range']-r['gal_range']):.3f}")

json.dump(rows, open("results/lysis_enzymes/ph_gap_summary.json", "w"), indent=2)
