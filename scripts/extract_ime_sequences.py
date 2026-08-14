import csv

FASTA = r"C:\Users\user\Desktop\HS6 vs. 4_RNA-seq_260720\HS6_Genome\HS6_scaf.fasta"
CANDS = r"C:\Users\user\AppData\Local\Temp\claude\C--Users-user-Desktop-HS6-vs--4-RNA-seq-260720\6e32b67b-d05a-4791-b9c3-db6271bce82e\scratchpad\utr_intron_candidates.tsv"

targets = ["g13643", "g5935", "g6200", "g13671", "g5747", "g11549"]

# load candidate rows
rows = {}
with open(CANDS, encoding="utf-8") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for r in reader:
        if r["gene"] in targets:
            rows[r["gene"]] = r

print("Found rows for:", list(rows.keys()))

# load genome sequences for needed scaffolds
needed = set(r["scaf"] for r in rows.values())
seqs = {}
cur_id, cur_seq = None, []
with open(FASTA, encoding="utf-8") as f:
    for line in f:
        if line.startswith(">"):
            if cur_id in needed:
                seqs[cur_id] = "".join(cur_seq)
            cur_id = line[1:].split()[0]
            cur_seq = []
        else:
            if cur_id in needed:
                cur_seq.append(line.strip())
    if cur_id in needed:
        seqs[cur_id] = "".join(cur_seq)

print("Loaded scaffolds:", list(seqs.keys()))

COMP = str.maketrans("ACGTacgt", "TGCAtgca")
def revcomp(s):
    return s.translate(COMP)[::-1]

utr5_fasta = []
intron_fasta = []
cassette_fasta = []  # 500bp upstream promoter stub + 5'UTR-with-intron, ready for a reporter cassette

PROMOTER_STUB = 500

for gid, r in rows.items():
    scaf = r["scaf"]
    strand = r["strand"]
    seq = seqs[scaf]
    u5s, u5e = int(r["utr5_start"]), int(r["utr5_end"])
    utr5_len = int(r["utr5_len"])
    n_utr5_introns = int(r["n_utr5_introns"])
    fi_s, fi_e = r["first_intron_start"], r["first_intron_end"]

    if strand == "+":
        utr5_seq = seq[u5s-1:u5e]
        prom_stub = seq[max(0,u5s-1-PROMOTER_STUB):u5s-1]
        cassette = prom_stub + utr5_seq
    else:
        utr5_seq = revcomp(seq[u5s-1:u5e])
        prom_stub = revcomp(seq[u5e:u5e+PROMOTER_STUB])
        cassette = prom_stub + utr5_seq

    utr5_fasta.append(f">{gid}_5UTR_{utr5_len}bp_{scaf}:{u5s}-{u5e}({strand})_with_{n_utr5_introns}_intron\n{utr5_seq}\n")
    cassette_fasta.append(f">{gid}_promoter{PROMOTER_STUB}bp_plus_5UTR_intron_cassette_{scaf}:{u5s}-{u5e}({strand})\n{cassette}\n")

    if fi_s and fi_e:
        fi_s, fi_e = int(fi_s), int(fi_e)
        if strand == "+":
            intron_seq = seq[fi_s-1:fi_e]
        else:
            intron_seq = revcomp(seq[fi_s-1:fi_e])
        intron_fasta.append(f">{gid}_intron_{fi_e-fi_s+1}bp_{scaf}:{fi_s}-{fi_e}({strand})\n{intron_seq}\n")

    print(f"{gid}: 5'UTR {utr5_len}bp, cassette {len(cassette)}bp, intron {(fi_e-fi_s+1) if fi_s else 'NA'}bp")

OUT = r"C:\Users\user\AppData\Local\Temp\claude\C--Users-user-Desktop-HS6-vs--4-RNA-seq-260720\6e32b67b-d05a-4791-b9c3-db6271bce82e\scratchpad"
with open(f"{OUT}\\ime_5utr_sequences.fasta", "w") as f:
    f.writelines(utr5_fasta)
with open(f"{OUT}\\ime_intron_only_sequences.fasta", "w") as f:
    f.writelines(intron_fasta)
with open(f"{OUT}\\ime_promoter_plus_5utr_cassette.fasta", "w") as f:
    f.writelines(cassette_fasta)
print("\nwrote 3 FASTA files")
