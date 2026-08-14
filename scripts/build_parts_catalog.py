import csv

SCRATCH = "/mnt/c/Users/user/AppData/Local/Temp/claude/C--Users-user-Desktop-HS6-vs--4-RNA-seq-260720/6e32b67b-d05a-4791-b9c3-db6271bce82e/scratchpad"

def load_fasta(path):
    seqs = {}
    cur_id, buf = None, []
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                if cur_id:
                    seqs[cur_id] = "".join(buf)
                cur_id = line[1:].strip()
                buf = []
            else:
                buf.append(line.strip())
        if cur_id:
            seqs[cur_id] = "".join(buf)
    return seqs

# rationale for promoters: mean_expr, CV, description
promoter_rationale = {}
with open(f"{SCRATCH}/promoter_candidates_ranked.tsv", encoding="utf-8") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for r in reader:
        promoter_rationale[r["gene"]] = r

# rationale for IME candidates
ime_rationale = {}
with open(f"{SCRATCH}/utr_intron_candidates.tsv", encoding="utf-8") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for r in reader:
        ime_rationale[r["gene"]] = r

# primer results, keyed by fasta header (part id)
primers = {}
with open(f"{SCRATCH}/part_primers.tsv", encoding="utf-8") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for r in reader:
        primers[r["part"]] = r

promoters = load_fasta(f"{SCRATCH}/promoter_candidates_top10.fasta")
cassettes = load_fasta(f"{SCRATCH}/ime_promoter_plus_5utr_cassette.fasta")
utr5s = load_fasta(f"{SCRATCH}/ime_5utr_sequences.fasta")
introns = load_fasta(f"{SCRATCH}/ime_intron_only_sequences.fasta")

rows = []

# 1) promoters
for header, seq in promoters.items():
    gid = header.split("_")[0]
    r = promoter_rationale.get(gid, {})
    p = primers.get(header, {})
    rows.append({
        "part_id": header,
        "gene_id": gid,
        "category": "Promoter (1.5kb upstream)",
        "rationale": f"mean_expr={r.get('mean_expr','NA')}, CV={r.get('CV','NA')}, stable/non-DE (padj={r.get('padj','NA')}), {r.get('description','')}",
        "length_bp": len(seq),
        "fwd_primer": p.get("fwd_primer",""), "fwd_tm": p.get("fwd_tm",""),
        "rev_primer": p.get("rev_primer",""), "rev_tm": p.get("rev_tm",""),
        "product_size": p.get("product_size",""),
        "sequence": seq,
    })

# 2) IME cassettes (promoter stub + 5'UTR + intron) -- primary lab-ready reporter part
for header, seq in cassettes.items():
    gid = header.split("_")[0]
    r = ime_rationale.get(gid, {})
    p = primers.get(header, {})
    rows.append({
        "part_id": header,
        "gene_id": gid,
        "category": "IME cassette (500bp promoter stub + 5'UTR + native intron)",
        "rationale": f"5'UTR={r.get('utr5_len','NA')}bp with {r.get('n_utr5_introns','NA')} intron(s) inside it; first_intron={r.get('first_intron_start','NA')}-{r.get('first_intron_end','NA')}",
        "length_bp": len(seq),
        "fwd_primer": p.get("fwd_primer",""), "fwd_tm": p.get("fwd_tm",""),
        "rev_primer": p.get("rev_primer",""), "rev_tm": p.get("rev_tm",""),
        "product_size": p.get("product_size",""),
        "sequence": seq,
    })

# 3) 5'UTR-only (no separate primers were designed for these; note "use cassette primers, trim in silico" OR flag as reference-only)
for header, seq in utr5s.items():
    gid = header.split("_")[0]
    rows.append({
        "part_id": header,
        "gene_id": gid,
        "category": "5'UTR only (reference sequence, contained within IME cassette above)",
        "rationale": "Subsequence of the IME cassette for this gene; amplify via the cassette primers then subclone, or order as a synthetic fragment",
        "length_bp": len(seq),
        "fwd_primer": "", "fwd_tm": "", "rev_primer": "", "rev_tm": "", "product_size": "",
        "sequence": seq,
    })

# 4) intron-only (reference; for intron-swap constructs)
for header, seq in introns.items():
    gid = header.split("_")[0]
    rows.append({
        "part_id": header,
        "gene_id": gid,
        "category": "Intron only (reference sequence, contained within IME cassette above)",
        "rationale": "Subsequence of the IME cassette for this gene; for direct intron-swap into a heterologous 5'UTR",
        "length_bp": len(seq),
        "fwd_primer": "", "fwd_tm": "", "rev_primer": "", "rev_tm": "", "product_size": "",
        "sequence": seq,
    })

# write master TSV (no full sequence column -- too wide for a table; separate FASTA holds sequence)
with open(f"{SCRATCH}/sequence_parts_catalog.tsv", "w", encoding="utf-8") as f:
    cols = ["part_id","gene_id","category","rationale","length_bp","fwd_primer","fwd_tm","rev_primer","rev_tm","product_size"]
    f.write("\t".join(cols) + "\n")
    for row in rows:
        f.write("\t".join(str(row[c]) for c in cols) + "\n")

# write master FASTA (all sequences, all categories, in one file)
with open(f"{SCRATCH}/all_parts_master.fasta", "w", encoding="utf-8") as f:
    for row in rows:
        f.write(f">{row['part_id']} | {row['category']}\n")
        seq = row["sequence"]
        for i in range(0, len(seq), 70):
            f.write(seq[i:i+70] + "\n")

print(f"wrote {len(rows)} rows to sequence_parts_catalog.tsv and all_parts_master.fasta")
for row in rows:
    print(f"  {row['part_id'][:50]:50s} {row['category'][:45]:45s} {row['length_bp']}bp primers={'yes' if row['fwd_primer'] else 'no'}")
