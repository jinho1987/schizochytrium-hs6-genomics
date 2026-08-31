# -*- coding: utf-8 -*-
import csv

BASE = r"C:\Users\user\Desktop\HS6 vs. 4_RNA-seq_260720\4.SNP\0.input_files"
HS6_PROT = BASE + r"\HS6\HS6_scaf_prot.fasta"
MILD_PROT = BASE + r"\HS4_mild_filter\HS4_mild_filt_prot.fasta"
NOFILT_PROT = BASE + r"\HS4_no_filter\HS4_no_filt_prot.fasta"
DE_TABLE = r"C:\Users\user\Desktop\HS6 vs. 4_RNA-seq_260720\RNA-seq_analysis_results\tables\DE_genotype_MTvsWT_adj_for_time.tsv"

def load_fasta(path):
    seqs = {}
    cur_id, buf = None, []
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                if cur_id:
                    seqs[cur_id] = "".join(buf)
                cur_id = line[1:].strip().split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if cur_id:
            seqs[cur_id] = "".join(buf)
    return seqs

print("Loading protein FASTAs (genome-wide, 16,088 genes each)...")
hs6 = load_fasta(HS6_PROT)
mild = load_fasta(MILD_PROT)
nofilt = load_fasta(NOFILT_PROT)
print(f"HS6: {len(hs6)}  mild(#4): {len(mild)}  nofilt(#4): {len(nofilt)}")

de_padj = {}
de_lfc = {}
de_name = {}
with open(DE_TABLE, encoding="utf-8") as f:
    reader = csv.reader(f, delimiter="\t")
    next(reader)
    for r in reader:
        if not r or not r[0]:
            continue
        gid = r[0]
        try:
            padj = float(r[6]) if r[6] not in ("", "NA") else None
        except (ValueError, IndexError):
            padj = None
        de_padj[gid] = padj
        try:
            de_lfc[gid] = float(r[2])
        except (ValueError, IndexError):
            de_lfc[gid] = None
        de_name[gid] = r[7] if len(r) > 7 else ""

def scan(hs6_seqs, wt4_seqs, label):
    missense_genes = []
    length_diff_genes = []
    identical = 0
    missing = 0
    for gid, hs6_seq in hs6_seqs.items():
        wt4_seq = wt4_seqs.get(gid)
        if wt4_seq is None:
            missing += 1
            continue
        if hs6_seq == wt4_seq:
            identical += 1
            continue
        if len(hs6_seq) != len(wt4_seq):
            length_diff_genes.append((gid, len(wt4_seq), len(hs6_seq)))
            continue
        diffs = [(i + 1, wt4_seq[i], hs6_seq[i]) for i in range(len(hs6_seq)) if wt4_seq[i] != hs6_seq[i]]
        real_missense = [d for d in diffs if d[1] != '*' and d[2] != '*']
        stop_changes = [d for d in diffs if d[1] == '*' or d[2] == '*']
        if real_missense or stop_changes:
            base_gid = gid.replace(".t1", "").replace(".t2", "")
            missense_genes.append({
                "gid": gid, "base_gid": base_gid,
                "n_missense": len(real_missense), "n_stop_changes": len(stop_changes),
                "changes": diffs,
                "padj": de_padj.get(base_gid), "log2fc": de_lfc.get(base_gid),
                "identity": de_name.get(base_gid, ""),
            })
    print(f"\n=== {label}: identical={identical}, missing_from_wt4_set={missing}, "
          f"length-different(indel/frameshift)={len(length_diff_genes)}, "
          f"same-length-with-aa-diffs={len(missense_genes)} ===")
    return missense_genes, length_diff_genes

mild_missense, mild_lendiff = scan(hs6, mild, "Mild_filter reconstruction")
nofilt_missense, nofilt_lendiff = scan(hs6, nofilt, "No_filter reconstruction")

OUT = r"C:\Users\user\AppData\Local\Temp\claude\C--Users-user-Desktop-HS6-vs--4-RNA-seq-260720\6e32b67b-d05a-4791-b9c3-db6271bce82e\scratchpad"
for missense_genes, fn in [(mild_missense, "genomewide_missense_mild.tsv"), (nofilt_missense, "genomewide_missense_nofilt.tsv")]:
    with open(f"{OUT}\\{fn}", "w", encoding="utf-8") as f:
        f.write("gene_id\tn_missense\tn_stop_changes\tpadj\tlog2FC\tidentity\taa_changes_WT4_to_HS6\n")
        for g in sorted(missense_genes, key=lambda x: (x["padj"] if x["padj"] is not None else 2)):
            changes = ",".join(f"p.{wt4a}{pos}{hs6a}" for pos, wt4a, hs6a in g["changes"])
            f.write(f"{g['gid']}\t{g['n_missense']}\t{g['n_stop_changes']}\t{g['padj']}\t{g['log2fc']}\t{g['identity']}\t{changes}\n")

print("\nSaved genomewide_missense_mild.tsv and genomewide_missense_nofilt.tsv")
