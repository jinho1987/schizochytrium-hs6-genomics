# -*- coding: utf-8 -*-
import re
import csv

GFF = r"C:\Users\user\Desktop\HS6 vs. 4_RNA-seq_260720\HS6_Genome\braker.gff3"
FASTA = r"C:\Users\user\Desktop\HS6 vs. 4_RNA-seq_260720\HS6_Genome\HS6_scaf.fasta"
OUT = r"C:\Users\user\AppData\Local\Temp\claude\C--Users-user-Desktop-HS6-vs--4-RNA-seq-260720\6e32b67b-d05a-4791-b9c3-db6271bce82e\scratchpad"
PANEL_TSV = f"{OUT}\\qpcr_primer_panel.tsv"

panel_genes = []
with open(PANEL_TSV, encoding="utf-8") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for r in reader:
        if r.get("gene_id"):
            panel_genes.append(r["gene_id"])

print("Panel genes:", len(panel_genes))

gene_coords = {}  # gid -> (scaf, start, end, strand)
with open(GFF, encoding="utf-8") as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9 or parts[2] != "gene":
            continue
        scaf, source, feat, start, end, score, strand, frame, attrs = parts
        m = re.search(r"ID=([^;]+)", attrs)
        if m and m.group(1) in panel_genes:
            gene_coords[m.group(1)] = (scaf, int(start), int(end), strand)

missing = set(panel_genes) - set(gene_coords)
if missing:
    print("MISSING gene coordinates for:", missing)

needed_scafs = set(v[0] for v in gene_coords.values())
print("Scaffolds needed:", needed_scafs)

genome = {}
cur_id, buf = None, []
with open(FASTA, encoding="utf-8") as f:
    for line in f:
        if line.startswith(">"):
            if cur_id in needed_scafs:
                genome[cur_id] = "".join(buf)
            cur_id = line[1:].split()[0]
            buf = []
        else:
            if cur_id in needed_scafs:
                buf.append(line.strip())
    if cur_id in needed_scafs:
        genome[cur_id] = "".join(buf)

COMP = str.maketrans("ACGTacgt", "TGCAtgca")
def revcomp(s):
    return s.translate(COMP)[::-1]

fasta_out = []
coord_rows = []
for gid in panel_genes:
    if gid not in gene_coords:
        coord_rows.append([gid, "", "", "", "", ""])
        continue
    scaf, start, end, strand = gene_coords[gid]
    seq = genome[scaf][start - 1:end]
    if strand == "-":
        seq = revcomp(seq)
    fasta_out.append(f">{gid}.t1 genomic_locus {scaf}:{start}-{end}({strand}) length={len(seq)}bp includes_introns_and_UTRs\n")
    for i in range(0, len(seq), 70):
        fasta_out.append(seq[i:i+70] + "\n")
    coord_rows.append([gid, scaf, start, end, strand, len(seq)])

with open(f"{OUT}\\qpcr_gene_genomic_coords.tsv", "w", encoding="utf-8") as f:
    f.write("gene_id\tscaffold\tgenomic_start\tgenomic_end\tstrand\tgenomic_span_bp\n")
    for row in coord_rows:
        f.write("\t".join(str(x) for x in row) + "\n")

with open(f"{OUT}\\qpcr_gene_genomic_sequences.fasta", "w", encoding="utf-8") as f:
    f.writelines(fasta_out)

print(f"wrote coordinates for {len(coord_rows)} genes and {len(fasta_out)} FASTA blocks")
for row in coord_rows:
    print(row)
