"""
Blind whole-protein AutoDock Vina docking for the 13-enzyme lysis panel.

For each of the 13 ESMFold structures, dock against:
  (a) its own canonical substrate
  (b) the shared galactan wall-probe
-> 26 docking runs total.

Methodology matches lysis-enzymes.html's documented methods section:
  - receptor prep: Open Babel (obabel), rigid, Gasteiger partial charges
  - ligand prep: RDKit + Meeko (already done by prep_lysis_ligands.py / prep_lysis_ligand_pdbqt.sh)
  - blind search box: covers the entire folded protein (no active site pre-mapped)
  - exhaustiveness=8, 9 poses per run
  - AutoDock Vina via the `vina` python package

Run inside the `docking_env` conda environment.
"""
import os, sys, json, subprocess, time

REPO = "/mnt/c/Users/USER/Desktop/schizochytrium-hs6-genomics"
STRUCT_DIR = f"{REPO}/results/lysis_enzymes/structures"
LIGAND_DIR = f"{REPO}/results/lysis_enzymes/ligands"
DOCK_DIR = f"{REPO}/results/lysis_enzymes/docking"
RESULTS_JSON = f"{REPO}/results/lysis_enzymes/docking_results.json"
os.makedirs(DOCK_DIR, exist_ok=True)

PANEL = [
    ("P00780", "subtilisin_carlsberg", "protease", "B. licheniformis", "protease_substrate"),
    ("P00782", "subtilisin_BPN", "protease", "B. amyloliquefaciens", "protease_substrate"),
    ("P04189", "subtilisin_E", "protease", "B. subtilis", "protease_substrate"),
    ("P07518", "subtilisin_pumilus", "protease", "B. pumilus", "protease_substrate"),
    ("P00698", "lysozyme_C_chicken", "lysozyme (GH22)", "Gallus gallus", "lysozyme_substrate"),
    ("P00720", "endolysin_T4", "lysozyme (GH24)", "phage T4", "lysozyme_substrate"),
    ("P25310", "lysozyme_M1_strepto", "lysozyme (GH25)", "Streptomyces globisporus", "lysozyme_substrate"),
    ("P07981", "endoglucanase_I_Tr", "cellulase", "Trichoderma reesei", "cellulase_substrate"),
    ("O74706", "endoglucanase_B_An", "cellulase", "Aspergillus niger", "cellulase_substrate"),
    ("G0RUP7", "xylanase2_Tr", "xylanase", "Trichoderma reesei", "xylanase_substrate"),
    ("P55330", "xylanaseB_An", "xylanase", "Aspergillus niger", "xylanase_substrate"),
    ("G0L322", "betaagaraseA_Zg", "agarase", "Zobellia galactanivorans", "agarase_substrate"),
    ("Q9RGX8", "betaagaraseB_Zg", "agarase", "Zobellia galactanivorans", "agarase_substrate"),
    # 5 new enzymes added to mirror the wet-lab Experiment 4 cell-wall-lysis screen
    ("P26213", "pectinase_An", "pectinase", "Aspergillus niger", "pectinase_substrate"),
    ("Q2UCU3", "betagalactosidase_Ao", "beta-galactosidase", "Aspergillus oryzae", "galactosidase_substrate"),
    ("P29600", "subtilisin_savinase_Bsp", "protease", "Bacillus sp. (Lederbergia lenta)", "protease_substrate"),
    ("P68736", "neutralprotease_Bs", "neutral protease (M4/bacillolysin)", "Bacillus subtilis", "protease_substrate"),
    ("P18429", "xylanaseC_Bs", "xylanase", "Bacillus subtilis", "xylanase_substrate"),
]

EXHAUSTIVENESS = 8
N_POSES = 9


def pdb_bbox_and_center(pdb_path):
    xs, ys, zs = [], [], []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                xs.append(x); ys.append(y); zs.append(z)
    cx, cy, cz = (max(xs) + min(xs)) / 2, (max(ys) + min(ys)) / 2, (max(zs) + min(zs)) / 2
    sx, sy, sz = (max(xs) - min(xs)), (max(ys) - min(ys)), (max(zs) - min(zs))
    pad = 8.0  # padding so ligand can approach from outside the protein envelope
    return (cx, cy, cz), (sx + pad, sy + pad, sz + pad)


def prep_receptor(label):
    pdb_path = f"{STRUCT_DIR}/{label}.pdb"
    pdbqt_path = f"{DOCK_DIR}/{label}_receptor.pdbqt"
    if os.path.exists(pdbqt_path):
        return pdbqt_path
    cmd = ["obabel", pdb_path, "-O", pdbqt_path, "-xr", "--partialcharge", "gasteiger"]
    print("  ", " ".join(cmd))
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    return pdbqt_path


def run_vina(receptor_pdbqt, ligand_pdbqt, center, box_size, out_pdbqt, log_path):
    from vina import Vina
    v = Vina(sf_name="vina")
    v.set_receptor(receptor_pdbqt)
    v.set_ligand_from_file(ligand_pdbqt)
    v.compute_vina_maps(center=list(center), box_size=list(box_size))
    v.dock(exhaustiveness=EXHAUSTIVENESS, n_poses=N_POSES)
    energies = v.energies(n_poses=N_POSES)
    v.write_poses(out_pdbqt, n_poses=N_POSES, overwrite=True)
    with open(log_path, "w") as f:
        f.write(f"receptor: {receptor_pdbqt}\n")
        f.write(f"ligand: {ligand_pdbqt}\n")
        f.write(f"center: {center}\n")
        f.write(f"box_size: {box_size}\n")
        f.write(f"exhaustiveness: {EXHAUSTIVENESS}\n")
        f.write(f"n_poses: {N_POSES}\n")
        f.write("mode  affinity(kcal/mol)  rmsd_lb  rmsd_ub\n")
        for i, e in enumerate(energies):
            f.write(f"{i+1}\t{e[0]:.4f}\t{e[1]:.4f}\t{e[2]:.4f}\n")
    top_affinity = float(energies[0][0])
    return top_affinity


def main():
    only = sys.argv[1:] if len(sys.argv) > 1 else None
    results = []
    if os.path.exists(RESULTS_JSON):
        with open(RESULTS_JSON) as f:
            results = json.load(f)
    done_keys = {(r["accession"], r["ligand"]) for r in results}

    for accession, label, klass, host, canonical_key in PANEL:
        if only and label not in only:
            continue
        pdb_path = f"{STRUCT_DIR}/{label}.pdb"
        if not os.path.exists(pdb_path):
            print(f"SKIP {label}: no folded structure yet")
            continue

        print(f"\n########## {label} ({accession}) ##########")
        receptor_pdbqt = prep_receptor(label)
        center, box_size = pdb_bbox_and_center(pdb_path)
        print(f"  blind box center={tuple(round(c,1) for c in center)} size={tuple(round(s,1) for s in box_size)}")

        for ligand_key in (canonical_key, "galactan_probe"):
            if (accession, ligand_key) in done_keys:
                print(f"  skip {ligand_key} (already docked)")
                continue
            ligand_pdbqt = f"{LIGAND_DIR}/{ligand_key}.pdbqt"
            out_pdbqt = f"{DOCK_DIR}/{label}_{ligand_key}.pdbqt"
            log_path = f"{DOCK_DIR}/{label}_{ligand_key}.log"
            print(f"  --- docking vs {ligand_key} ---")
            t0 = time.time()
            try:
                affinity = run_vina(receptor_pdbqt, ligand_pdbqt, center, box_size, out_pdbqt, log_path)
                elapsed = time.time() - t0
                print(f"    top affinity: {affinity:.3f} kcal/mol ({elapsed:.1f}s)")
                results.append({
                    "accession": accession, "label": label, "class": klass, "host": host,
                    "ligand": ligand_key, "top_affinity_kcal_mol": round(affinity, 3),
                    "seconds": round(elapsed, 1), "status": "OK",
                })
            except Exception as e:
                elapsed = time.time() - t0
                print(f"    FAILED: {e}")
                results.append({
                    "accession": accession, "label": label, "class": klass, "host": host,
                    "ligand": ligand_key, "top_affinity_kcal_mol": None,
                    "seconds": round(elapsed, 1), "status": f"FAILED: {e}",
                })
            with open(RESULTS_JSON, "w") as f:
                json.dump(results, f, indent=2)

    print("\n=== DONE ===")
    for r in results:
        print(r)


if __name__ == "__main__":
    main()
