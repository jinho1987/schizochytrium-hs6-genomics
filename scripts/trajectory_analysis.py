import csv

TABLES = r"C:\Users\user\Desktop\HS6 vs. 4_RNA-seq_260720\RNA-seq_analysis_results\tables"

# load per-timepoint log2FC
traj = {}
with open(f"{TABLES}\\per_timepoint_log2FC_descriptive.tsv", encoding="utf-8") as f:
    reader = csv.reader(f, delimiter="\t")
    header = next(reader)
    for row in reader:
        gid = row[0]
        try:
            lfc20, lfc44, lfc68 = float(row[1]), float(row[2]), float(row[3])
        except ValueError:
            continue
        traj[gid] = (lfc20, lfc44, lfc68, row[5] if len(row) > 5 else "")

# load significant DE genes (genotype effect)
sig_genes = {}
with open(f"{TABLES}\\DE_genotype_MTvsWT_adj_for_time.tsv", encoding="utf-8") as f:
    reader = csv.reader(f, delimiter="\t")
    header = next(reader)
    for row in reader:
        if not row or not row[0]:
            continue
        gid = row[0].replace(".t1", "")
        try:
            lfc = float(row[2])
            padj = float(row[6]) if row[6] not in ("", "NA") else None
        except (ValueError, IndexError):
            continue
        if padj is not None and padj < 0.05:
            sig_genes[gid] = (lfc, padj)

print(f"Significant DE genes (padj<0.05): {len(sig_genes)}")

# classify trajectory shape for each sig gene
def classify(lfc20, lfc44, lfc68):
    mags = [abs(lfc20), abs(lfc44), abs(lfc68)]
    # require consistent sign for a clean call
    signs = [1 if x > 0 else -1 for x in (lfc20, lfc44, lfc68)]
    consistent_sign = signs[0] == signs[1] == signs[2]
    if not consistent_sign:
        return "mixed_direction"
    if mags[0] >= mags[1] >= mags[2] * 0.8:
        return "early_responder"  # strong early, fading or flat
    if mags[2] >= mags[1] >= mags[0] * 0.8 and mags[2] > mags[0] * 1.5:
        return "late_responder"  # builds up over time
    if mags[1] > mags[0] * 1.3 and mags[1] > mags[2] * 1.3:
        return "transient_peak_44h"
    return "roughly_flat"

counts = {}
classified = {}
for gid, (lfc, padj) in sig_genes.items():
    if gid not in traj:
        continue
    lfc20, lfc44, lfc68, desc = traj[gid]
    cls = classify(lfc20, lfc44, lfc68)
    counts[cls] = counts.get(cls, 0) + 1
    classified[gid] = (cls, lfc20, lfc44, lfc68, desc)

print("\n=== Trajectory shape distribution among significant DE genes ===")
for cls, n in sorted(counts.items(), key=lambda x: -x[1]):
    print(f"  {cls:20s}: {n}")

# Now specifically check MYB3R-family and known lipid/cytoskeleton pathway genes
MYB_GENES = ["g13677", "g8263", "g10428", "g15706", "g16063", "g8479", "g9086", "g11711"]  # from earlier project's MYB3R candidate list
PATHWAY_GENES = {
    "g12095": "FabD", "g15357": "Pks13", "g6387": "FA-CoA ligase", "g7952": "FA-CoA ligase",
    "g7955": "FA-CoA ligase", "g7896": "DGAT", "g3889": "alkane monooxygenase", "g8236": "SDP1",
}

print("\n=== MYB3R-family gene trajectories ===")
for gid in MYB_GENES:
    if gid in classified:
        cls, lfc20, lfc44, lfc68, desc = classified[gid]
        print(f"  {gid:10s} {cls:20s} 20h={lfc20:+.2f} 44h={lfc44:+.2f} 68h={lfc68:+.2f}  {desc}")
    elif gid in traj:
        lfc20, lfc44, lfc68, desc = traj[gid]
        print(f"  {gid:10s} {'(not significant)':20s} 20h={lfc20:+.2f} 44h={lfc44:+.2f} 68h={lfc68:+.2f}  {desc}")
    else:
        print(f"  {gid:10s} NOT FOUND in trajectory table")

print("\n=== Core lipid-pathway gene trajectories ===")
for gid, name in PATHWAY_GENES.items():
    if gid in classified:
        cls, lfc20, lfc44, lfc68, desc = classified[gid]
        print(f"  {gid:10s} {name:20s} {cls:20s} 20h={lfc20:+.2f} 44h={lfc44:+.2f} 68h={lfc68:+.2f}")
    elif gid in traj:
        lfc20, lfc44, lfc68, desc = traj[gid]
        print(f"  {gid:10s} {name:20s} {'(not significant)':20s} 20h={lfc20:+.2f} 44h={lfc44:+.2f} 68h={lfc68:+.2f}")
    else:
        print(f"  {gid:10s} {name:20s} NOT FOUND")

# save full classification
with open(r"C:\Users\user\AppData\Local\Temp\claude\C--Users-user-Desktop-HS6-vs--4-RNA-seq-260720\6e32b67b-d05a-4791-b9c3-db6271bce82e\scratchpad\trajectory_classification.tsv", "w", encoding="utf-8") as out:
    out.write("gene\tclass\tlfc20\tlfc44\tlfc68\tdescription\n")
    for gid, (cls, lfc20, lfc44, lfc68, desc) in classified.items():
        out.write(f"{gid}\t{cls}\t{lfc20:.3f}\t{lfc44:.3f}\t{lfc68:.3f}\t{desc}\n")
print("\nwrote full classification to trajectory_classification.tsv")
