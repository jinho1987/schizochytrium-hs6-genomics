# -*- coding: utf-8 -*-
import re

BLOSUM62_RAW = """
A R N D C Q E G H I L K M F P S T W Y V
A 4 -1 -2 -2 0 -1 -1 0 -2 -1 -1 -1 -1 -2 -1 1 0 -3 -2 0
R -1 5 0 -2 -3 1 0 -2 0 -3 -2 2 -1 -3 -2 -1 -1 -3 -2 -3
N -2 0 6 1 -3 0 0 0 1 -3 -3 0 -2 -3 -2 1 0 -4 -2 -3
D -2 -2 1 6 -3 0 2 -1 -1 -3 -4 -1 -3 -3 -1 0 -1 -4 -3 -3
C 0 -3 -3 -3 9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
Q -1 1 0 0 -3 5 2 -2 0 -3 -2 1 0 -3 -1 0 -1 -2 -1 -2
E -1 0 0 2 -4 2 5 -2 0 -3 -3 1 -2 -3 -1 0 -1 -3 -2 -2
G 0 -2 0 -1 -3 -2 -2 6 -2 -4 -4 -2 -3 -3 -2 0 -2 -2 -3 -3
H -2 0 1 -1 -3 0 0 -2 8 -3 -3 -1 -2 -1 -2 -1 -2 -2 2 -3
I -1 -3 -3 -3 -1 -3 -3 -4 -3 4 2 -3 1 0 -3 -2 -1 -3 -1 3
L -1 -2 -3 -4 -1 -2 -3 -4 -3 2 4 -2 2 0 -3 -2 -1 -2 -1 1
K -1 2 0 -1 -3 1 1 -2 -1 -3 -2 5 -1 -3 -1 0 -1 -3 -2 -2
M -1 -1 -2 -3 -1 0 -2 -3 -2 1 2 -1 5 0 -2 -1 -1 -1 -1 1
F -2 -3 -3 -3 -2 -3 -3 -3 -1 0 0 -3 0 6 -4 -2 -2 1 3 -1
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4 7 -1 -1 -4 -3 -2
S 1 -1 1 0 -1 0 0 0 -1 -2 -2 0 -1 -2 -1 4 1 -3 -2 -2
T 0 -1 0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1 1 5 -2 -2 0
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1 1 -4 -3 -2 11 2 -3
Y -2 -2 -2 -3 -2 -1 -2 -3 2 -1 -1 -2 -1 3 -3 -2 -2 2 7 -1
V 0 -3 -3 -3 -1 -2 -2 -3 -3 3 1 -2 1 -1 -2 -2 0 -3 -1 4
"""
lines = [l.split() for l in BLOSUM62_RAW.strip().split("\n")]
aa_order = lines[0]
BLOSUM62 = {}
for row in lines[1:]:
    a = row[0]
    for j, val in enumerate(row[1:]):
        BLOSUM62[(a, aa_order[j])] = int(val)

GRANTHAM = {
    ('S','R'):110,('S','L'):145,('S','P'):74,('S','T'):58,('S','A'):99,('S','V'):124,('S','G'):56,('S','I'):142,
    ('S','F'):155,('S','Y'):144,('S','C'):112,('S','H'):89,('S','Q'):68,('S','N'):46,('S','K'):121,('S','D'):65,
    ('S','E'):80,('S','M'):135,('S','W'):177,('P','R'):103,('P','L'):98,('P','T'):38,('P','A'):27,('P','V'):68,
    ('P','G'):42,('P','I'):95,('P','F'):114,('P','Y'):110,('P','C'):169,('P','H'):77,('P','Q'):76,('P','N'):91,
    ('P','K'):103,('P','D'):108,('P','E'):93,('P','M'):87,('P','W'):147,('T','R'):71,('T','L'):92,('T','A'):58,
    ('T','V'):69,('T','G'):59,('T','I'):89,('T','F'):103,('T','Y'):92,('T','C'):149,('T','H'):47,('T','Q'):42,
    ('T','N'):65,('T','K'):78,('T','D'):85,('T','E'):65,('T','M'):81,('T','W'):128,('A','R'):112,('A','L'):96,
    ('A','V'):64,('A','G'):60,('A','I'):94,('A','F'):113,('A','Y'):112,('A','C'):195,('A','H'):86,('A','Q'):91,
    ('A','N'):111,('A','K'):106,('A','D'):126,('A','E'):107,('A','M'):84,('A','W'):148,('G','R'):125,('G','L'):138,
    ('G','V'):109,('G','I'):135,('G','F'):153,('G','Y'):147,('G','C'):159,('G','H'):98,('G','Q'):87,('G','N'):80,
    ('G','K'):127,('G','D'):94,('G','E'):98,('G','M'):127,('G','W'):184,('C','R'):180,('C','L'):198,('C','V'):192,
    ('C','I'):198,('C','F'):205,('C','Y'):194,('C','H'):174,('C','Q'):154,('C','N'):139,('C','K'):202,('C','D'):154,
    ('C','E'):170,('C','M'):196,('C','W'):215,('D','R'):96,('D','L'):172,('D','V'):152,('D','I'):168,('D','F'):177,
    ('D','Y'):160,('D','H'):81,('D','Q'):61,('D','N'):23,('D','K'):101,('D','E'):45,('D','M'):160,('D','W'):181,
    ('E','R'):54,('E','L'):138,('E','V'):121,('E','I'):134,('E','F'):140,('E','Y'):122,('E','H'):40,('E','Q'):29,
    ('E','N'):42,('E','K'):56,('E','M'):126,('E','W'):152,('N','R'):86,('N','L'):149,('N','V'):133,('N','I'):149,
    ('N','F'):158,('N','Y'):143,('N','H'):68,('N','Q'):46,('N','K'):94,('N','M'):142,('N','W'):174,('Q','R'):43,
    ('Q','L'):113,('Q','V'):96,('Q','I'):109,('Q','F'):116,('Q','Y'):99,('Q','H'):24,('Q','K'):53,('Q','M'):101,
    ('Q','W'):130,('H','R'):29,('H','L'):99,('H','V'):84,('H','I'):94,('H','F'):100,('H','Y'):83,('H','K'):32,
    ('H','M'):87,('H','W'):115,('K','R'):26,('K','L'):107,('K','V'):97,('K','I'):102,('K','F'):102,('K','Y'):85,
    ('K','M'):95,('K','W'):110,('M','R'):91,('M','L'):15,('M','V'):21,('M','I'):10,('M','F'):28,('M','Y'):36,
    ('M','W'):67,('I','R'):97,('I','L'):5,('I','V'):29,('I','F'):21,('I','Y'):33,('I','W'):61,('L','R'):102,
    ('L','V'):32,('L','F'):22,('L','Y'):36,('L','W'):61,('V','F'):50,('V','Y'):55,('V','W'):88,('F','Y'):22,
    ('F','W'):40,('Y','W'):37,
    ('R','F'):205,('R','Y'):194,('R','W'):101,('R','C'):180,('R','G'):125,('R','A'):112,
    ('R','V'):96,('R','P'):103,('R','T'):71,('R','S'):110,
}
def grantham(a, b):
    if a == b:
        return 0
    return GRANTHAM.get((a, b)) or GRANTHAM.get((b, a)) or None

HYDROPHOBICITY = {
    'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,'I':4.5,
    'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2,
}
CHARGE = {'D':-1,'E':-1,'K':1,'R':1,'H':0.5}

def blosum(a, b):
    return BLOSUM62.get((a, b), BLOSUM62.get((b, a)))

def classify(bl, gr):
    if gr is None:
        return "unscored"
    if gr >= 150 or (bl is not None and bl <= -2):
        return "RADICAL"
    if gr >= 100 or (bl is not None and bl <= 0):
        return "MODERATE"
    return "CONSERVATIVE"

MUTATIONS = [
    ("g1664", "H205Q", 0.0144, -1.109, "unannotated"),
    ("g1760", "T894A", 0.2571, -0.965, "G-type lectin receptor kinase"),
    ("g853",  "N83K",  0.5544, -0.508, "unannotated -- LRR/RanGAP-like (InterProScan)"),
    ("g11453","M15I",  0.7264, 0.276, "unannotated -- weak ssDNA-ATPase hit"),
    ("g11453","R298H", 0.7264, 0.276, "unannotated"),
    ("g11453","M379I", 0.7264, 0.276, "unannotated"),
    ("g11453","L506F", 0.7264, 0.276, "unannotated"),
    ("g11453","D511A", 0.7264, 0.276, "unannotated"),
    ("g11453","N547S", 0.7264, 0.276, "unannotated"),
    ("g11453","T570A", 0.7264, 0.276, "unannotated"),
    ("g10614","E51A",  0.9262, 0.101, "unannotated -- LRR-family (InterProScan)"),
    ("g10298","T13R",  None,   0.661, "unannotated -- LRR/Ribonuclease-Inhibitor fold"),
    ("g12826","T624K", None,  -0.570, "unannotated -- LRR/RanGAP-like (InterProScan)"),
    ("g6055", "N236K", 0.0094, -1.554, "unannotated -- LRR-family (InterProScan)"),
    ("g6055", "T280N", 0.0094, -1.554, "unannotated"),
    ("g6055", "K337Q", 0.0094, -1.554, "unannotated"),
    ("g6423", "R27F",  0.0189, -1.109, "RanGAP-like nuclear transport"),
    ("g6423", "D30S",  0.0189, -1.109, "RanGAP-like nuclear transport"),
    ("g6423", "A280V", 0.0189, -1.109, "RanGAP-like nuclear transport"),
    ("g15200","K439Q", 0.0426, -0.987, "unannotated -- LRR/RanGAP-adjacent (InterProScan)"),
    ("g1759", "E251R", 0.0777, -0.896, "unannotated -- LRR-family (InterProScan)"),
    ("g1759", "T263A", 0.0777, -0.896, "unannotated"),
]

print(f"{'Gene':9s} {'Change':8s} {'BLOSUM62':>8s} {'Grantham':>9s} {'Class':>12s} {'padj':>10s} {'Identity'}")
rows = []
for gid, change, padj, lfc, identity in MUTATIONS:
    m = re.match(r"([A-Z])(\d+)([A-Z])", change)
    wt4_aa, pos, hs6_aa = m.group(1), m.group(2), m.group(3)
    bl = blosum(wt4_aa, hs6_aa)
    gr = grantham(wt4_aa, hs6_aa)
    cls = classify(bl, gr)
    hphob_delta = round(HYDROPHOBICITY[hs6_aa] - HYDROPHOBICITY[wt4_aa], 2)
    charge_wt4 = CHARGE.get(wt4_aa, 0)
    charge_hs6 = CHARGE.get(hs6_aa, 0)
    charge_flip = "YES" if (charge_wt4 * charge_hs6 < 0) else ("neutral->charged" if (charge_wt4==0)!=(charge_hs6==0) else "no")
    pro_flag = ""
    if wt4_aa == 'P' and hs6_aa != 'P':
        pro_flag = "PROLINE_REMOVED"
    elif hs6_aa == 'P' and wt4_aa != 'P':
        pro_flag = "PROLINE_INTRODUCED"
    padj_str = f"{padj:.4f}" if padj is not None else "NA"
    print(f"{gid:9s} {change:8s} {bl!s:>8s} {gr!s:>9s} {cls:>12s} {padj_str:>10s} {identity}")
    rows.append([gid, change, wt4_aa, hs6_aa, bl, gr, cls, hphob_delta, charge_flip, pro_flag, padj, lfc, identity])

OUT = r"C:\Users\user\AppData\Local\Temp\claude\C--Users-user-Desktop-HS6-vs--4-RNA-seq-260720\6e32b67b-d05a-4791-b9c3-db6271bce82e\scratchpad\missense_impact_scores.tsv"
with open(OUT, "w", encoding="utf-8") as f:
    cols = ["gene","change","wt4_aa","hs6_aa","blosum62","grantham","impact_class",
            "hydrophobicity_delta_KD","charge_flip","proline_flag","padj","log2FC","identity"]
    f.write("\t".join(cols) + "\n")
    for r in rows:
        f.write("\t".join(str(x) for x in r) + "\n")
print(f"\nsaved {len(rows)} rows to missense_impact_scores.tsv")

radical = [r for r in rows if r[6] == "RADICAL"]
print(f"\n=== {len(radical)} RADICAL-class substitutions ===")
for r in radical:
    print(f"  {r[0]} {r[1]}: BLOSUM62={r[4]}, Grantham={r[5]}")
