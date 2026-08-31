import os, time, json
os.environ.setdefault("HF_HOME", "/mnt/d/hf_cache")
import torch
from transformers import AutoTokenizer, EsmForProteinFolding

def load_fasta(path):
    seqs = {}
    cur_id, buf = None, []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if cur_id: seqs[cur_id] = "".join(buf)
                cur_id = line[1:].strip()
                buf = []
            else:
                buf.append(line.strip())
        if cur_id: seqs[cur_id] = "".join(buf)
    return seqs

OUT = "/mnt/d/esmfold_results"
os.makedirs(OUT, exist_ok=True)

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

seqs = load_fasta("/mnt/c/Users/user/AppData/Local/Temp/claude/C--Users-user-Desktop-HS6-vs--4-RNA-seq-260720/6e32b67b-d05a-4791-b9c3-db6271bce82e/scratchpad/esmfold_targets.fasta")

summary = []
for name, seq in seqs.items():
    print(f"\n=== Folding {name} ({len(seq)} aa) ===")
    t0 = time.time()
    tokenized = tokenizer([seq], return_tensors="pt", add_special_tokens=False)
    tokenized = {k: v.to(device) for k, v in tokenized.items()}
    with torch.no_grad():
        output = model(tokenized["input_ids"])
    plddt = output["plddt"][0, :, 1].mean().item() * 100 if output["plddt"].dim() == 3 else output["plddt"][0].mean().item()
    # per-residue mean pLDDT (average over atoms) if needed
    per_res_plddt = output["plddt"][0].mean(dim=-1) if output["plddt"].dim() == 3 else output["plddt"][0]
    mean_plddt = per_res_plddt.mean().item()
    elapsed = time.time() - t0
    print(f"  mean pLDDT: {mean_plddt:.2f}  (took {elapsed:.1f}s)")

    from transformers.models.esm.openfold_utils.protein import to_pdb, Protein as OFProtein
    from transformers.models.esm.openfold_utils.feats import atom14_to_atom37

    def convert_outputs_to_pdb(outputs):
        final_atom_positions = atom14_to_atom37(outputs["positions"][-1], outputs)
        outputs = {k: v.to("cpu").numpy() for k, v in outputs.items()}
        final_atom_positions = final_atom_positions.cpu().numpy()
        final_atom_mask = outputs["atom37_atom_exists"]
        pdbs = []
        for i in range(outputs["aatype"].shape[0]):
            aa = outputs["aatype"][i]
            pred_pos = final_atom_positions[i]
            mask = final_atom_mask[i]
            resid = outputs["residue_index"][i] + 1
            pred = OFProtein(
                aatype=aa, atom_positions=pred_pos, atom_mask=mask,
                residue_index=resid, b_factors=outputs["plddt"][i],
                chain_index=outputs["chain_index"][i] if "chain_index" in outputs else None,
            )
            pdbs.append(to_pdb(pred))
        return pdbs

    pdb_str = convert_outputs_to_pdb(output)[0]
    pdb_path = f"{OUT}/{name}.pdb"
    with open(pdb_path, "w") as f:
        f.write(pdb_str)
    print(f"  saved {pdb_path}")

    summary.append({"name": name, "length": len(seq), "mean_plddt": round(mean_plddt, 2), "seconds": round(elapsed, 1)})
    del output
    torch.cuda.empty_cache() if torch.cuda.is_available() else None

with open(f"{OUT}/summary.json", "w") as f:
    json.dump(summary, f, indent=2)
print("\n=== DONE ===")
for s in summary:
    print(s)
