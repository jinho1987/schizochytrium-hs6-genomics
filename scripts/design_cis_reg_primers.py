import subprocess

def load_fasta(path):
    seqs = {}
    cur_id, buf = None, []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if cur_id: seqs[cur_id] = "".join(buf)
                cur_id = line[1:].strip().split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if cur_id: seqs[cur_id] = "".join(buf)
    return seqs

REF = "/home/user/rnaseq_hs6/genome/HS6_scaf.fasta"
genome = load_fasta(REF)

def run_primer3(seq_id, seq, excl_start, excl_len, product_min, product_max):
    cfg = f"""SEQUENCE_ID={seq_id}
SEQUENCE_TEMPLATE={seq}
SEQUENCE_EXCLUDED_REGION={excl_start},{excl_len}
PRIMER_TASK=generic
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_OPT_SIZE=20
PRIMER_MIN_SIZE=18
PRIMER_MAX_SIZE=23
PRIMER_OPT_TM=60.0
PRIMER_MIN_TM=57.0
PRIMER_MAX_TM=63.0
PRIMER_PRODUCT_SIZE_RANGE={product_min}-{product_max}
=
"""
    p = subprocess.run(["primer3_core"], input=cfg, capture_output=True, text=True)
    out = {}
    for line in p.stdout.split("\n"):
        if "=" in line:
            k, v = line.split("=", 1)
            out[k] = v
    return out

# target, chrom, variant_min_pos, variant_max_pos (1-based), flank_needed
TARGETS = [
    ("g8342",     "Scaffolds_22", 2083254, 2083281),
    ("g568",      "Scaffolds_1",  369524,  369525),
    ("g6124_site1", "Scaffolds_2", 336298,  336322),
    ("g6124_site2", "Scaffolds_2", 337644,  337644),
    ("g8869",     "Scaffolds_23", 2124233, 2124233),
    ("g14345",    "Scaffolds_4",  2843502, 2843502),
]

FLANK = 700  # >600bp requested, use 700 for margin
results = []
for name, chrom, vmin, vmax in TARGETS:
    win_start = vmin - FLANK
    win_end = vmax + FLANK
    seq = genome[chrom][win_start-1:win_end]  # 0-based slice
    variant_span_len = vmax - vmin + 1
    excl_start_in_template = FLANK - 50  # small buffer so excluded region comfortably covers variant span
    excl_len = variant_span_len + 100
    template_len = len(seq)
    r = run_primer3(name, seq, excl_start_in_template, excl_len, template_len - 300, template_len + 50)
    ok = "PRIMER_LEFT_0_SEQUENCE" in r
    print(f"\n--- {name} {chrom}:{vmin}-{vmax} (window {win_start}-{win_end}) ---")
    if not ok:
        print("  FAILED:", r.get("PRIMER_ERROR", r))
        continue
    fwd, rev = r["PRIMER_LEFT_0_SEQUENCE"], r["PRIMER_RIGHT_0_SEQUENCE"]
    ftm, rtm = r["PRIMER_LEFT_0_TM"], r["PRIMER_RIGHT_0_TM"]
    size = r["PRIMER_PAIR_0_PRODUCT_SIZE"]
    fwd_start_in_template = int(r["PRIMER_LEFT_0"].split(",")[0])
    rev_start_in_template = int(r["PRIMER_RIGHT_0"].split(",")[0])
    fwd_genome_pos = win_start + fwd_start_in_template
    rev_genome_pos = win_start + rev_start_in_template
    dist_fwd_to_variant = vmin - fwd_genome_pos
    dist_rev_to_variant = rev_genome_pos - vmax
    print(f"  FWD: {fwd} (Tm {ftm}) at genome pos ~{fwd_genome_pos}, {dist_fwd_to_variant}bp upstream of variant")
    print(f"  REV: {rev} (Tm {rtm}) at genome pos ~{rev_genome_pos}, {dist_rev_to_variant}bp downstream of variant")
    print(f"  product: {size} bp")
    results.append({
        "gene": name, "chrom": chrom, "variant_range": f"{vmin}-{vmax}",
        "fwd": fwd, "ftm": ftm, "rev": rev, "rtm": rtm, "size": size,
        "flank_fwd_bp": dist_fwd_to_variant, "flank_rev_bp": dist_rev_to_variant,
    })

OUT = "/mnt/c/Users/user/AppData/Local/Temp/claude/C--Users-user-Desktop-HS6-vs--4-RNA-seq-260720/6e32b67b-d05a-4791-b9c3-db6271bce82e/scratchpad/cis_regulatory_sanger_primers.tsv"
with open(OUT, "w") as f:
    f.write("gene\tchrom\tvariant_range\tfwd_primer\tfwd_tm\trev_primer\trev_tm\tproduct_bp\tflank_fwd_bp\tflank_rev_bp\n")
    for r in results:
        f.write(f"{r['gene']}\t{r['chrom']}\t{r['variant_range']}\t{r['fwd']}\t{r['ftm']}\t{r['rev']}\t{r['rtm']}\t{r['size']}\t{r['flank_fwd_bp']}\t{r['flank_rev_bp']}\n")
print(f"\nwrote {len(results)} primer sets")
