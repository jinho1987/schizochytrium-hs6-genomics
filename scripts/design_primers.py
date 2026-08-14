import primer3

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

def design(seq, product_size_range):
    result = primer3.bindings.design_primers(
        seq_args={
            'SEQUENCE_ID': 'target',
            'SEQUENCE_TEMPLATE': seq,
        },
        global_args={
            'PRIMER_TASK': 'generic',
            'PRIMER_PICK_LEFT_PRIMER': 1,
            'PRIMER_PICK_RIGHT_PRIMER': 1,
            'PRIMER_OPT_SIZE': 22,
            'PRIMER_MIN_SIZE': 18,
            'PRIMER_MAX_SIZE': 27,
            'PRIMER_OPT_TM': 60.0,
            'PRIMER_MIN_TM': 57.0,
            'PRIMER_MAX_TM': 63.0,
            'PRIMER_MIN_GC': 35.0,
            'PRIMER_MAX_GC': 65.0,
            'PRIMER_PRODUCT_SIZE_RANGE': [product_size_range],
            'PRIMER_NUM_RETURN': 1,
            'PRIMER_MAX_POLY_X': 4,
            'PRIMER_SALT_MONOVALENT': 50.0,
            'PRIMER_DNA_CONC': 50.0,
        }
    )
    return result

all_results = []

# 1) promoters -- amplify (nearly) the full 1500bp upstream region
proms = load_fasta(f"{SCRATCH}/promoter_candidates_top10.fasta")
for name, seq in proms.items():
    L = len(seq)
    res = design(seq, [max(L-100, 200), L])
    all_results.append((name, "promoter", seq, res))

# 2) IME cassettes -- amplify the full promoter-stub+5'UTR+intron cassette
casettes = load_fasta(f"{SCRATCH}/ime_promoter_plus_5utr_cassette.fasta")
for name, seq in casettes.items():
    L = len(seq)
    res = design(seq, [max(L-150, 200), L])
    all_results.append((name, "IME cassette", seq, res))

# report
out_lines = ["part\ttype\ttemplate_len\tfwd_primer\tfwd_tm\tfwd_gc\trev_primer\trev_tm\trev_gc\tproduct_size\n"]
for name, ptype, seq, res in all_results:
    n = res.get('PRIMER_PAIR_NUM_RETURNED', 0)
    if n == 0:
        print(f"{name}: NO PRIMERS FOUND (len={len(seq)})")
        out_lines.append(f"{name}\t{ptype}\t{len(seq)}\tNONE_FOUND\t\t\t\t\t\t\n")
        continue
    fwd = res['PRIMER_LEFT_0_SEQUENCE']
    fwd_tm = res['PRIMER_LEFT_0_TM']
    fwd_gc = res['PRIMER_LEFT_0_GC_PERCENT']
    rev = res['PRIMER_RIGHT_0_SEQUENCE']
    rev_tm = res['PRIMER_RIGHT_0_TM']
    rev_gc = res['PRIMER_RIGHT_0_GC_PERCENT']
    prod = res['PRIMER_PAIR_0_PRODUCT_SIZE']
    print(f"{name}: FWD {fwd} (Tm={fwd_tm:.1f}) REV {rev} (Tm={rev_tm:.1f}) product={prod}bp")
    out_lines.append(f"{name}\t{ptype}\t{len(seq)}\t{fwd}\t{fwd_tm:.1f}\t{fwd_gc:.1f}\t{rev}\t{rev_tm:.1f}\t{rev_gc:.1f}\t{prod}\n")

with open(f"{SCRATCH}/part_primers.tsv", "w", encoding="utf-8") as f:
    f.writelines(out_lines)
print(f"\nwrote part_primers.tsv ({len(all_results)} parts)")
