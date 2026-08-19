# -*- coding: utf-8 -*-
import re
import csv
import primer3

GFF = "/mnt/c/Users/user/Desktop/HS6 vs. 4_RNA-seq_260720/HS6_Genome/braker.gff3"
FASTA = "/mnt/c/Users/user/Desktop/HS6 vs. 4_RNA-seq_260720/HS6_Genome/HS6_scaf.fasta"
TABLES = "/mnt/c/Users/user/Desktop/HS6 vs. 4_RNA-seq_260720/RNA-seq_analysis_results/tables"
OUT = "/mnt/c/Users/user/AppData/Local/Temp/claude/C--Users-user-Desktop-HS6-vs--4-RNA-seq-260720/6e32b67b-d05a-4791-b9c3-db6271bce82e/scratchpad"

PANEL = [
    # gene_id, category, identity (fallback if not in DE table), log2FC, padj
    ("g9086",  "MYB regulatory network", "Transcriptional activator Myb",            2.78,  9.7e-17),
    ("g13536", "MYB regulatory network", "CDK-like kinase",                          3.36,  1.6e-18),
    ("g15706", "MYB regulatory network", "Transcription factor MYB98",               1.41,  8.4e-4),
    ("g16063", "MYB regulatory network", "Transcription factor MYB3R-1",             1.32,  1.2e-3),
    ("g8479",  "MYB regulatory network", "Transcription factor MYB3R-5",             1.50,  2.8e-2),
    ("g11711", "MYB regulatory network", "Transcription factor MYB3R-3",            -1.58,  1.1e-5),
    ("g13677", "MYB regulatory network", "Myb-related protein B",                   -1.10,  6.5e-4),
    ("g8263",  "MYB regulatory network", "Transcription factor MYB3R-3 (2nd locus)",-0.96,  1.8e-3),
    ("g10428", "MYB regulatory network", "Myb-related protein A",                   -0.78,  1.2e-2),
    ("g12095", "Lipid/PUFA pathway",     "Malonyl CoA-ACP transacylase (FabD)",     -1.86,  2.4e-6),
    ("g15357", "Lipid/PUFA pathway",     "Polyketide synthase Pks13",                1.54,  1.4e-14),
    ("g6387",  "Lipid/PUFA pathway",     "Long-chain-fatty-acid--CoA ligase 5",      0.59,  3.0e-2),
    ("g7952",  "Lipid/PUFA pathway",     "Long-chain-fatty-acid--CoA ligase ACSBG2",-1.71,  2.4e-4),
    ("g7955",  "Lipid/PUFA pathway",     "Long-chain-fatty-acid--CoA ligase ACSBG2",-2.85,  6.2e-5),
    ("g7896",  "Lipid/PUFA pathway",     "Diacylglycerol O-acyltransferase (DGAT)", -2.63,  2.4e-4),
    ("g3889",  "Lipid/PUFA pathway",     "Alkane 1-monooxygenase 1",                -1.72,  8.7e-2),
    ("g8236",  "Lipid/PUFA pathway",     "Triacylglycerol lipase SDP1",             -1.33,  2.7e-9),
    ("g11722", "Lipid/PUFA pathway",     "Diacylglycerol lipase-alpha",             -1.25,  8.1e-3),
    ("g8276",  "Lipid/PUFA pathway",     "Lipid droplet-associated hydrolase",      -0.84,  5.3e-3),
    ("g6423",  "Candidate driver locus", "RanGAP-like nuclear transport",           -1.11,  1.9e-2),
    ("g3906",  "Candidate driver locus", "Carnitine O-acetyltransferase",            1.91,  1.2e-3),
    ("g14345", "Candidate driver locus", "unannotated",                              0.76,  3.3e-2),
    ("g6124",  "Candidate driver locus", "Ribosomal protein S6 kinase alpha-2",     -0.59,  1.7e-2),
    ("g342",   "Largest-effect gene",    "unannotated",                              6.03,  1.5e-28),
    ("g827",   "qPCR reference/normalizer", "Plasma-membrane H+-ATPase 4",           0.04,  0.969),
    ("g227",   "qPCR reference/normalizer", "Elongation factor 2",                  -0.01,  0.990),
    ("g10543", "qPCR reference/normalizer", "ATP synthase F1 catalytic subunit beta",0.03,  0.978),
]

def load_fasta(path):
    seqs = {}
    cur_id, buf = None, []
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                if cur_id:
                    seqs[cur_id] = "".join(buf)
                cur_id = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if cur_id:
            seqs[cur_id] = "".join(buf)
    return seqs

COMP = str.maketrans("ACGTacgt", "TGCAtgca")
def revcomp(s):
    return s.translate(COMP)[::-1]

# parse GFF3: get scaffold+strand per gene, and all CDS blocks per transcript
gene_scaf_strand = {}
cds_blocks = {}  # gene_id -> list of (start,end) 1-based inclusive, genome order
needed_scafs = set()

with open(GFF, encoding="utf-8") as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9:
            continue
        scaf, source, feat, start, end, score, strand, frame, attrs = parts
        if feat == "gene":
            m = re.search(r"ID=([^;]+)", attrs)
            if m:
                gid = m.group(1)
                gene_scaf_strand[gid] = (scaf, strand)
        elif feat == "CDS":
            m = re.search(r"Parent=([^;.]+)\.t1", attrs)
            if m:
                gid = m.group(1)
                cds_blocks.setdefault(gid, []).append((int(start), int(end)))

panel_genes = [row[0] for row in PANEL]
for gid in panel_genes:
    if gid in gene_scaf_strand:
        needed_scafs.add(gene_scaf_strand[gid][0])

print("Loading genome scaffolds needed:", needed_scafs)
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
print("Loaded:", list(genome.keys()))

def build_spliced_cds(gid):
    scaf, strand = gene_scaf_strand[gid]
    blocks = sorted(cds_blocks[gid], key=lambda b: b[0])
    seq_full = genome[scaf]
    exon_seqs = [seq_full[s - 1:e] for s, e in blocks]
    if strand == "-":
        exon_seqs = [revcomp(s) for s in reversed(exon_seqs)]
    spliced = "".join(exon_seqs)
    # junction positions: cumulative length after each exon except the last
    junctions = []
    running = 0
    for s in exon_seqs[:-1]:
        running += len(s)
        junctions.append(running)
    return spliced, junctions, len(blocks)

def design_qpcr(seq, junctions):
    global_args = {
        'PRIMER_TASK': 'generic',
        'PRIMER_PICK_LEFT_PRIMER': 1,
        'PRIMER_PICK_RIGHT_PRIMER': 1,
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 25,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 58.0,
        'PRIMER_MAX_TM': 62.0,
        'PRIMER_PAIR_MAX_DIFF_TM': 1.5,
        'PRIMER_MIN_GC': 40.0,
        'PRIMER_MAX_GC': 60.0,
        'PRIMER_PRODUCT_SIZE_RANGE': [[70, 200]],
        'PRIMER_NUM_RETURN': 1,
        'PRIMER_MAX_POLY_X': 3,
        'PRIMER_SALT_MONOVALENT': 50.0,
        'PRIMER_DNA_CONC': 50.0,
        'PRIMER_MAX_SELF_ANY': 6.0,
        'PRIMER_MAX_SELF_END': 2.0,
        'PRIMER_PAIR_MAX_COMPL_ANY': 6.0,
        'PRIMER_PAIR_MAX_COMPL_END': 2.0,
    }
    seq_args = {'SEQUENCE_ID': 'target', 'SEQUENCE_TEMPLATE': seq}
    if junctions:
        seq_args['SEQUENCE_OVERLAP_JUNCTION_LIST'] = junctions
        global_args['PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION'] = 7
        global_args['PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION'] = 4
    res = primer3.bindings.design_primers(seq_args=seq_args, global_args=global_args)
    return res

rows = []
fasta_out = []
for gid, category, identity, log2fc, padj in PANEL:
    if gid not in gene_scaf_strand or gid not in cds_blocks:
        print(f"{gid}: MISSING from GFF3, skipping")
        continue
    scaf, strand = gene_scaf_strand[gid]
    spliced, junctions, n_exons = build_spliced_cds(gid)
    res = design_qpcr(spliced, junctions)
    n_pairs = res.get('PRIMER_PAIR_NUM_RETURNED', 0)

    if n_pairs == 0 and junctions:
        # retry without the junction constraint if it over-constrained the search
        res = design_qpcr(spliced, [])
        n_pairs = res.get('PRIMER_PAIR_NUM_RETURNED', 0)
        junction_note = "junction-span FAILED, fallback used"
    elif junctions:
        junction_note = "spans exon junction"
    else:
        junction_note = "single-exon CDS -- no junction possible"

    if n_pairs == 0:
        print(f"{gid}: NO PRIMERS FOUND at all (len={len(spliced)})")
        rows.append([gid, f"{gid}.t1", category, identity, log2fc, padj, n_exons, junction_note,
                     "NONE_FOUND", "", "", "", "", ""])
        continue

    fwd = res['PRIMER_LEFT_0_SEQUENCE']; fwd_tm = res['PRIMER_LEFT_0_TM']
    rev = res['PRIMER_RIGHT_0_SEQUENCE']; rev_tm = res['PRIMER_RIGHT_0_TM']
    prod = res['PRIMER_PAIR_0_PRODUCT_SIZE']
    fwd_pos = res['PRIMER_LEFT_0'][0]
    rev_start = res['PRIMER_RIGHT_0'][0]

    print(f"{gid} ({identity[:30]}): FWD {fwd} REV {rev} product={prod}bp exons={n_exons} [{junction_note}]")
    rows.append([gid, f"{gid}.t1", category, identity, log2fc, padj, n_exons, junction_note,
                 fwd, round(fwd_tm,1), rev, round(rev_tm,1), prod, len(spliced)])

    fasta_out.append(f">{gid}.t1 spliced_CDS {len(spliced)}bp {scaf}({strand}) n_exons={n_exons}\n{spliced}\n")

with open(f"{OUT}/qpcr_primer_panel.tsv", "w", encoding="utf-8") as f:
    cols = ["gene_id","transcript_id","category","identity","log2FC","padj","n_exons",
            "primer_design_note","fwd_primer","fwd_tm","rev_primer","rev_tm","product_size_bp","spliced_cds_len_bp"]
    f.write("\t".join(cols) + "\n")
    for r in rows:
        f.write("\t".join(str(x) for x in r) + "\n")

with open(f"{OUT}/qpcr_target_spliced_cds.fasta", "w", encoding="utf-8") as f:
    f.writelines(fasta_out)

print(f"\nwrote {len(rows)} rows to qpcr_primer_panel.tsv and {len(fasta_out)} sequences to qpcr_target_spliced_cds.fasta")
