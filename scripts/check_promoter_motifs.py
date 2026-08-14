import re

FASTA = r"C:\Users\user\AppData\Local\Temp\claude\C--Users-user-Desktop-HS6-vs--4-RNA-seq-260720\6e32b67b-d05a-4791-b9c3-db6271bce82e\scratchpad\promoter_candidates_top10.fasta"

# IUPAC-aware regex for the significant discovered MYB-associated motif: WGCACGTGCW
iupac = {"A":"A","C":"C","G":"G","T":"T","W":"[AT]","S":"[GC]","R":"[AG]","Y":"[CT]","K":"[GT]","M":"[AC]"}
def iupac_to_regex(s):
    return "".join(iupac.get(c, c) for c in s)

myb_motif = iupac_to_regex("WGCACGTGCW")
myb_motif_rc = iupac_to_regex("WGCACGTGCW"[::-1].translate(str.maketrans("ACGT","TGCA")))

tata_regex = re.compile(r"TATA[AT][AT][AT][GA]")

seqs = {}
cur = None
buf = []
with open(FASTA, encoding="utf-8") as f:
    for line in f:
        if line.startswith(">"):
            if cur:
                seqs[cur] = "".join(buf)
            cur = line[1:].strip()
            buf = []
        else:
            buf.append(line.strip())
    if cur:
        seqs[cur] = "".join(buf)

print(f"{'gene_header':55s} {'len':>5s} {'MYB-motif hits':>15s} {'TATA-like hits':>15s}")
for hdr, seq in seqs.items():
    myb_hits = len(re.findall(myb_motif, seq)) + len(re.findall(myb_motif_rc, seq))
    tata_hits = len(tata_regex.findall(seq))
    print(f"{hdr:55s} {len(seq):5d} {myb_hits:15d} {tata_hits:15d}")
