"""
Prepare the 6 docking ligands for the lysis-enzyme panel:
  5 canonical substrates (one per enzyme class, fetched from PubChem by name)
  + 1 modeled "wall probe" (sulfated 3,6-anhydro-L-galactose-linked disaccharide,
    built by hand since no solved/deposited structure of the native
    Schizochytrium wall polysaccharide exists).

For each ligand: embed a 3D conformer (RDKit ETKDGv3 + MMFF94 optimization),
save SDF + PDB, then convert to a Vina-ready PDBQT via Meeko.

Run inside the `docking_env` conda environment.
"""
import os
from rdkit import Chem
from rdkit.Chem import AllChem

REPO = "/mnt/c/Users/USER/Desktop/schizochytrium-hs6-genomics"
OUT = f"{REPO}/results/lysis_enzymes/ligands"
os.makedirs(OUT, exist_ok=True)

LIGANDS = {
    "protease_substrate": (
        "Suc-AAPF-pNA (succinyl-Ala-Ala-Pro-Phe-4-nitroanilide)",
        "CC(C(=O)NC(C)C(=O)N1CCCC1C(=O)NC(CC2=CC=CC=C2)C(=O)NC3=CC=C(C=C3)[N+](=O)[O-])NC(=O)CCC(=O)O",
        "PubChem name lookup",
    ),
    "lysozyme_substrate": (
        "chitohexaose",
        "C(C1C(C(C(C(O1)OCC2C(C(C(C(O2)OC3C(C(C(C(O3)CO)O)OC4C(C(C(C(O4)CO)O)O)N)N)N)OC5C(C(C(C(O5)CO)O)O)N)OC6C(C(C(C(O6)CO)O)O)N)N)O)O)O",
        "PubChem name lookup",
    ),
    "cellulase_substrate": (
        "cellopentaose",
        "C(C1C(C(C(C(O1)OC2C(OC(C(C2O)O)OC3C(OC(C(C3O)O)OC4C(OC(C(C4O)O)OC5C(OC(C(C5O)O)O)CO)CO)CO)CO)O)O)O)O",
        "PubChem name lookup",
    ),
    "xylanase_substrate": (
        "xylopentaose",
        "C1C(C(C(C(O1)OC2COC(C(C2O)O)OC3COC(C(C3O)O)OC4COC(C(C4O)O)OC(CO)C(C(C=O)O)O)O)O)O",
        "PubChem name lookup",
    ),
    "agarase_substrate": (
        "neoagarohexaose",
        "C1C2C(C(O1)C(C(O2)OC3C(C(OC(C3O)OC4C5COC4C(C(O5)OC6C(C(OC(C6O)OC7C8COC7C(C(O8)OC9C(C(OC(C9O)O)CO)O)O)CO)O)O)CO)O)O)O",
        "PubChem name lookup",
    ),
    "galactan_probe": (
        "modeled sulfated 3,6-anhydro-L-galactose-(1->4)-L-galactose-6-sulfate disaccharide "
        "(no PubChem/PDB entry exists; hand-built from the agarobiose scaffold "
        "[PubChem CanonicalSMILES for 'agarobiose'] with a sulfate half-ester added "
        "at C6 of the non-anhydro ring, approximating the >95% L-galactose, "
        "sulfated-galactan composition reported for the native Schizochytrium wall)",
        "C1C2C(C(O1)C(C(O2)O)O)OC3C(C(C(C(O3)COS(=O)(=O)O)O)O)O",
        "RDKit-constructed (no PubChem entry for the real wall polymer)",
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
        # retry with random coords as fallback
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
    for key, (desc, smiles, source) in LIGANDS.items():
        sdf_path = prep_one(key, desc, smiles, source)
        manifest.append({"key": key, "description": desc, "smiles": smiles, "source": source, "sdf": sdf_path})

    import json
    with open(f"{OUT}/ligand_manifest.json", "w") as f:
        json.dump(manifest, f, indent=2)
    print("\n=== DONE: all ligands embedded ===")


if __name__ == "__main__":
    main()
