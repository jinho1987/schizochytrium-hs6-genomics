"""
Prepare the 2 NEW docking ligands needed for the 5-enzyme lysis-panel expansion
(mirroring the wet-lab Experiment 4 screen):

  - pectinase_substrate: digalacturonic acid (PubChem CID 6857565) for the new
    endopolygalacturonase (pectinase_An).
  - galactosidase_substrate: lactose (PubChem CID 6134) for the new
    beta-galactosidase (betagalactosidase_Ao).

Both resolved directly by PubChem name/CID lookup (no hand-built RDKit
structure was needed, unlike the original galactan_probe).

The 3 other new enzymes (subtilisin_savinase_Bsp, neutralprotease_Bs,
xylanaseC_Bs) REUSE the existing protease_substrate / xylanase_substrate
ligands already prepped under results/lysis_enzymes/ligands/ -- no new ligand
files for those, consistent with the panel's "one substrate per functional
class" convention.

Same embed/optimize/save approach as prep_lysis_ligands.py. Appends to (does
not overwrite) the existing ligand_manifest.json. Run inside `docking_env`.
"""
import os
import json
from rdkit import Chem
from rdkit.Chem import AllChem

REPO = "/mnt/c/Users/USER/Desktop/schizochytrium-hs6-genomics"
OUT = f"{REPO}/results/lysis_enzymes/ligands"
MANIFEST = f"{OUT}/ligand_manifest.json"
os.makedirs(OUT, exist_ok=True)

NEW_LIGANDS = {
    "pectinase_substrate": (
        "digalacturonic acid (alpha-D-galacturonopyranosyl-(1->4)-D-galacturonic acid)",
        "C1(C(C(OC(C1O)OC2C(C(C(OC2C(=O)O)O)O)O)C(=O)O)O)O",
        "PubChem CID 6857565 (name lookup: 'digalacturonic acid')",
    ),
    "galactosidase_substrate": (
        "lactose",
        "C(C1C(C(C(C(O1)OC2C(OC(C(C2O)O)O)CO)O)O)O)O",
        "PubChem CID 6134 (name lookup: 'lactose')",
    ),
}


def prep_one(key, desc, smiles, source):
    print(f"\n=== {key}: {desc} ===")
    print(f"  source: {source}")
    print(f"  SMILES: {smiles}")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit failed to parse SMILES for {key}")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    cid = AllChem.EmbedMolecule(mol, params)
    if cid < 0:
        params.useRandomCoords = True
        cid = AllChem.EmbedMolecule(mol, params)
    if cid < 0:
        raise RuntimeError(f"3D embedding failed for {key}")
    ff_ok = AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
    print(f"  embedded conformer, MMFF optimize return code: {ff_ok} (0=converged)")

    mol.SetProp("_Name", key)
    sdf_path = f"{OUT}/{key}.sdf"
    pdb_path = f"{OUT}/{key}.pdb"
    with Chem.SDWriter(sdf_path) as w:
        w.write(mol)
    Chem.MolToPDBFile(mol, pdb_path)
    print(f"  saved {sdf_path}")
    print(f"  saved {pdb_path}")
    return sdf_path


def main():
    manifest = []
    if os.path.exists(MANIFEST):
        with open(MANIFEST) as f:
            manifest = json.load(f)
    existing_keys = {m["key"] for m in manifest}

    for key, (desc, smiles, source) in NEW_LIGANDS.items():
        if key in existing_keys:
            print(f"skip {key} (already in manifest)")
            continue
        sdf_path = prep_one(key, desc, smiles, source)
        manifest.append({"key": key, "description": desc, "smiles": smiles, "source": source, "sdf": sdf_path})

    with open(MANIFEST, "w") as f:
        json.dump(manifest, f, indent=2)
    print("\n=== DONE: new ligands embedded, manifest updated (old entries preserved) ===")


if __name__ == "__main__":
    main()
