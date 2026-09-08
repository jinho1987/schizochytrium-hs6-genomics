import os, sys, time, json, glob

import torch
from transformers import AutoTokenizer, EsmForProteinFolding

REPO = "/mnt/c/Users/USER/Desktop/schizochytrium-hs6-genomics"
SEQ_DIR = f"{REPO}/results/lysis_enzymes/sequences"
OUT_DIR = f"{REPO}/results/lysis_enzymes/structures"
SUMMARY_PATH = f"{REPO}/results/lysis_enzymes/esmfold_summary.json"
os.makedirs(OUT_DIR, exist_ok=True)


def load_fasta(path):
    seqs = {}
    cur_id, buf = None, []
    with open(path) as f:
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


def convert_output_to_pdb(output):
    from transformers.models.esm.openfold_utils.protein import to_pdb, Protein as OFProtein
    from transformers.models.esm.openfold_utils.feats import atom14_to_atom37

    final_atom_positions = atom14_to_atom37(output["positions"][-1], output)
    out = {k: v.to("cpu").numpy() for k, v in output.items()}
    final_atom_positions = final_atom_positions.cpu().numpy()
    final_atom_mask = out["atom37_atom_exists"]
    aa = out["aatype"][0]
    pred_pos = final_atom_positions[0]
    mask = final_atom_mask[0]
    resid = out["residue_index"][0] + 1
    pred = OFProtein(
        aatype=aa, atom_positions=pred_pos, atom_mask=mask,
        residue_index=resid, b_factors=out["plddt"][0],
        chain_index=out["chain_index"][0] if "chain_index" in out else None,
    )
    return to_pdb(pred)


def main():
    only = sys.argv[1:] if len(sys.argv) > 1 else None

    labels_order = [
        "subtilisin_carlsberg", "subtilisin_BPN", "subtilisin_E", "subtilisin_pumilus",
        "lysozyme_C_chicken", "endolysin_T4", "lysozyme_M1_strepto",
        "endoglucanase_I_Tr", "endoglucanase_B_An",
        "xylanase2_Tr", "xylanaseB_An",
        "betaagaraseA_Zg", "betaagaraseB_Zg",
    ]
    # sort shortest-first so cheap ones finish fast and any OOM shows up early
    fastas = {}
    for label in labels_order:
        path = f"{SEQ_DIR}/{label}.fasta"
        seqs = load_fasta(path)
        (_, seq), = seqs.items()
        fastas[label] = seq
    order = sorted(fastas.keys(), key=lambda l: len(fastas[l]))
    if only:
        order = [l for l in order if l in only]

    print("CUDA available:", torch.cuda.is_available())
    if torch.cuda.is_available():
        print("GPU:", torch.cuda.get_device_name(0))

    print("Loading tokenizer + model (facebook/esmfold_v1)...")
    tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
    model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True)
    device = "cuda" if torch.cuda.is_available() else "cpu"
    model = model.to(device)
    model.esm = model.esm.half()
    model.trunk.set_chunk_size(64)
    model.eval()
    print("Model loaded on", device)

    # load existing summary if resuming
    summary = []
    if os.path.exists(SUMMARY_PATH):
        with open(SUMMARY_PATH) as f:
            summary = json.load(f)
    done = {s["label"] for s in summary}

    for label in order:
        if label in done:
            print(f"skip {label} (already in summary)")
            continue
        seq = fastas[label]
        print(f"\n=== Folding {label} ({len(seq)} aa) ===")
        t0 = time.time()
        chunk_size = 64
        pdb_str = None
        mean_plddt = None
        last_err = None
        # retry with smaller chunk size on OOM
        for cs in (64, 32, 16, 8):
            try:
                model.trunk.set_chunk_size(cs)
                tokenized = tokenizer([seq], return_tensors="pt", add_special_tokens=False)
                tokenized = {k: v.to(device) for k, v in tokenized.items()}
                with torch.no_grad():
                    output = model(tokenized["input_ids"])
                per_res_plddt = output["plddt"][0].mean(dim=-1) if output["plddt"].dim() == 3 else output["plddt"][0]
                mean_plddt = per_res_plddt.mean().item()
                pdb_str = convert_output_to_pdb(output)
                chunk_size = cs
                del output
                break
            except torch.cuda.OutOfMemoryError as e:
                last_err = e
                print(f"  OOM at chunk_size={cs}, retrying smaller...")
                torch.cuda.empty_cache()
                continue
        elapsed = time.time() - t0

        if pdb_str is None:
            print(f"  FAILED to fold {label}: {last_err}")
            summary.append({
                "label": label, "length": len(seq), "mean_plddt": None,
                "seconds": round(elapsed, 1), "chunk_size": None, "status": "FAILED",
            })
        else:
            pdb_path = f"{OUT_DIR}/{label}.pdb"
            with open(pdb_path, "w") as f:
                f.write(pdb_str)
            print(f"  mean pLDDT: {mean_plddt:.4f} (took {elapsed:.1f}s, chunk_size={chunk_size})")
            print(f"  saved {pdb_path}")
            summary.append({
                "label": label, "length": len(seq), "mean_plddt": round(mean_plddt, 4),
                "seconds": round(elapsed, 1), "chunk_size": chunk_size, "status": "OK",
            })

        with open(SUMMARY_PATH, "w") as f:
            json.dump(summary, f, indent=2)
        torch.cuda.empty_cache() if torch.cuda.is_available() else None

    print("\n=== DONE ===")
    for s in summary:
        print(s)


if __name__ == "__main__":
    main()
