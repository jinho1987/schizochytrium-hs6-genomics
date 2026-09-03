import os, time, json
os.environ.setdefault("HF_HOME", "/mnt/d/hf_cache")
import torch
from transformers import AutoTokenizer, EsmForProteinFolding

OUT = "/mnt/d/esmfold_results"
os.makedirs(OUT, exist_ok=True)

print("CUDA available:", torch.cuda.is_available(), flush=True)
tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True)
device = "cuda" if torch.cuda.is_available() else "cpu"
model = model.to(device)
model.esm = model.esm.half()
model.trunk.set_chunk_size(64)
model.eval()
print("Model loaded on", device, flush=True)

seq = open("/mnt/c/Users/user/AppData/Local/Temp/claude/C--Users-user-Desktop-HS6-vs--4-RNA-seq-260720/6e32b67b-d05a-4791-b9c3-db6271bce82e/scratchpad/g3557.fasta").read().split("\n")[1].strip()
name = "g3557_HS6_RPS6"

print(f"=== Folding {name} ({len(seq)} aa) ===", flush=True)
t0 = time.time()
tokenized = tokenizer([seq], return_tensors="pt", add_special_tokens=False)
tokenized = {k: v.to(device) for k, v in tokenized.items()}
with torch.no_grad():
    output = model(tokenized["input_ids"])
per_res_plddt = output["plddt"][0].mean(dim=-1) if output["plddt"].dim() == 3 else output["plddt"][0]
mean_plddt = per_res_plddt.mean().item()
elapsed = time.time() - t0
print(f"  mean pLDDT: {mean_plddt:.2f}  (took {elapsed:.1f}s)", flush=True)

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
print(f"  saved {pdb_path}", flush=True)

with open(f"{OUT}/summary_g3557.json", "w") as f:
    json.dump({"name": name, "length": len(seq), "mean_plddt": round(mean_plddt, 2), "seconds": round(elapsed, 1)}, f, indent=2)
print("=== DONE ===", flush=True)
