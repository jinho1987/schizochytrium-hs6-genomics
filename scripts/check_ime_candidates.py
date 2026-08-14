import csv

TABLES = r"C:\Users\user\Desktop\HS6 vs. 4_RNA-seq_260720\RNA-seq_analysis_results\tables"

clean_candidates = ["g6200", "g13671", "g11549", "g13643", "g5935", "g5747"]
suspicious_candidates = ["g14644", "g7492", "g4134", "g9752", "g13578", "g9813"]

# pull description + expression from the promoter ranking file (already has mean expr + description)
info = {}
with open(r"C:\Users\user\AppData\Local\Temp\claude\C--Users-user-Desktop-HS6-vs--4-RNA-seq-260720\6e32b67b-d05a-4791-b9c3-db6271bce82e\scratchpad\promoter_candidates_ranked.tsv", encoding="utf-8") as f:
    next(f)
    for line in f:
        p = line.rstrip("\n").split("\t")
        info[p[0]] = {"mean": p[1], "cv": p[2], "desc": p[5] if len(p) > 5 else ""}

print("=== CLEAN IME candidates (plausible 5'UTR length, <1.2kb) ===")
for g in clean_candidates:
    d = info.get(g, {})
    print(f"  {g}: mean_expr={d.get('mean','?')}, CV={d.get('cv','?')}, {d.get('desc','no description')}")

print("\n=== SUSPICIOUS (likely StringTie gene-fusion artifact, huge pseudo-UTR) ===")
for g in suspicious_candidates:
    d = info.get(g, {})
    print(f"  {g}: mean_expr={d.get('mean','?')}, {d.get('desc','no description')}")
