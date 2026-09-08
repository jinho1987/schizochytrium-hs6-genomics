import os, sys, time, json

import torch
from transformers import AutoTokenizer, EsmForProteinFolding

REPO = "/mnt/c/Users/USER/Desktop/schizochytrium-hs6-genomics"
SEQ_DIR = f"{REPO}/results/lysis_enzymes/sequences"
OUT_DIR = f"{REPO}/results/lysis_enzymes/structures"
SUMMARY_PATH = f"{REPO}/results/lysis_enzymes/esmfold_summary.json"

LABEL = "betagalactosidase_Ao"


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
    (_, seq), = load_fasta(f"{SEQ_DIR}/{LABEL}.fasta").items()
    print(f"Folding {LABEL} ({len(seq)} aa) on CPU -- GPU (8GB) can't hold this sequence's "
          f"O(L^2) pair representation regardless of trunk chunk_size; falling back to CPU "
          f"with ~11GB free system RAM and 16 cores instead of fighting VRAM.", flush=True)

    torch.set_num_threads(16)

    print("Loading tokenizer + model (facebook/esmfold_v1) on CPU...", flush=True)
    tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
    model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True)
    model = model.float()  # fp16 is not well-supported on CPU
    model.trunk.set_chunk_size(64)
    model.eval()
    print("Model loaded on CPU", flush=True)

    t0 = time.time()
    tokenized = tokenizer([seq], return_tensors="pt", add_special_tokens=False)
    print("Running fold (this will take a while on CPU) ...", flush=True)
    with torch.no_grad():
        output = model(tokenized["input_ids"])
    per_res_plddt = output["plddt"][0].mean(dim=-1) if output["plddt"].dim() == 3 else output["plddt"][0]
    mean_plddt = per_res_plddt.mean().item()
    pdb_str = convert_output_to_pdb(output)
    elapsed = time.time() - t0

    summary = []
    if os.path.exists(SUMMARY_PATH):
        with open(SUMMARY_PATH) as f:
            summary = json.load(f)
    summary = [s for s in summary if s.get("label") != LABEL]

    pdb_path = f"{OUT_DIR}/{LABEL}.pdb"
    with open(pdb_path, "w") as f:
        f.write(pdb_str)
    print(f"mean pLDDT: {mean_plddt:.4f} (took {elapsed:.1f}s, device=cpu)", flush=True)
    print(f"saved {pdb_path}", flush=True)
    summary.append({
        "label": LABEL, "length": len(seq), "mean_plddt": round(mean_plddt, 4),
        "seconds": round(elapsed, 1), "chunk_size": 64, "status": "OK", "device": "cpu",
    })
    with open(SUMMARY_PATH, "w") as f:
        json.dump(summary, f, indent=2)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
