import re

GFF = r"C:\Users\user\Desktop\HS6 vs. 4_RNA-seq_260720\HS6_Genome\braker.gff3"
FASTA = r"C:\Users\user\Desktop\HS6 vs. 4_RNA-seq_260720\HS6_Genome\HS6_scaf.fasta"
CANDS_FILE = r"C:\Users\user\AppData\Local\Temp\claude\C--Users-user-Desktop-HS6-vs--4-RNA-seq-260720\6e32b67b-d05a-4791-b9c3-db6271bce82e\scratchpad\promoter_candidates_ranked.tsv"

TOP_N = 10
UPSTREAM = 1500

# load top N candidate gene IDs
top_genes = []
with open(CANDS_FILE, encoding="utf-8") as f:
    next(f)
    for line in f:
        gid = line.split("\t")[0]
        top_genes.append(gid)
        if len(top_genes) >= TOP_N:
            break

# parse gff3 for gene coordinates
coords = {}
with open(GFF, encoding="utf-8") as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9 or parts[2] != "gene":
            continue
        m = re.search(r"ID=([^;]+)", parts[8])
        if not m:
            continue
        gid = m.group(1)
        if gid in top_genes:
            coords[gid] = (parts[0], int(parts[3]), int(parts[4]), parts[6])

print("Found coordinates for:", list(coords.keys()))

# load genome fasta into memory per-scaffold (only scaffolds we need)
needed_scaffolds = set(c[0] for c in coords.values())
seqs = {}
cur_id = None
cur_seq = []
with open(FASTA, encoding="utf-8") as f:
    for line in f:
        if line.startswith(">"):
            if cur_id in needed_scaffolds:
                seqs[cur_id] = "".join(cur_seq)
            cur_id = line[1:].split()[0]
            cur_seq = []
        else:
            if cur_id in needed_scaffolds:
                cur_seq.append(line.strip())
    if cur_id in needed_scaffolds:
        seqs[cur_id] = "".join(cur_seq)

print("Loaded scaffolds:", list(seqs.keys()))

COMP = str.maketrans("ACGTacgt", "TGCAtgca")
def revcomp(s):
    return s.translate(COMP)[::-1]

out_lines = []
for gid in top_genes:
    if gid not in coords:
        continue
    scaf, start, end, strand = coords[gid]
    seq = seqs[scaf]
    if strand == "+":
        up_start = max(0, start - 1 - UPSTREAM)
        up_end = start - 1
        promoter = seq[up_start:up_end]
    else:
        up_start = end
        up_end = min(len(seq), end + UPSTREAM)
        promoter = revcomp(seq[up_start:up_end])
    out_lines.append(f">{gid}_promoter_{UPSTREAM}bp_{scaf}:{start}-{end}({strand})\n{promoter}\n")
    print(f"{gid}: {scaf}:{start}-{end} ({strand}) -> promoter length {len(promoter)}")

with open(r"C:\Users\user\AppData\Local\Temp\claude\C--Users-user-Desktop-HS6-vs--4-RNA-seq-260720\6e32b67b-d05a-4791-b9c3-db6271bce82e\scratchpad\promoter_candidates_top10.fasta", "w", encoding="utf-8") as f:
    f.writelines(out_lines)
print("\nwrote promoter_candidates_top10.fasta")
