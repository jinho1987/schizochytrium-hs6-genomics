"""
pH-dependent re-docking for the 18-enzyme lysis panel.

For each enzyme, at pH 4.0 and pH 9.0 (the existing docking_results.json runs
serve as the pH~7 baseline -- receptor prep there used Open Babel's default
Gasteiger charges with no explicit protonation-state assignment, which is a
reasonable proxy for near-neutral pH):

  1. Protonate the ESMFold structure at that pH with PDB2PQR (PROPKA pKa
     prediction), producing a pH-specific PQR.
  2. Convert to a Vina receptor PDBQT via Open Babel (same convention as the
     pH-7 baseline receptor prep in dock_lysis_panel.py).
  3. Dock against BOTH the enzyme's canonical substrate and the shared
     galactan probe -- same blind whole-protein box, same exhaustiveness=8/
     9 poses as every other run in this project.

Scope note: only the RECEPTOR's protonation state is varied by pH here; the
ligands are reused as-is from the pH~7 prep (their own protonation states are
not separately re-assigned per pH). This isolates the more tractable half of
the question (does the enzyme's own surface/pocket charge state change its
apparent affinity for the wall probe?) without doubling the ligand-prep
complexity -- flagged explicitly in the write-up as a real methodological
simplification, not hidden.

Run inside the `docking_env` conda environment (needs pdb2pqr30 + obabel +
the `vina` python package, all already installed there this session).
"""
import os, subprocess, json, time

REPO = "/mnt/c/Users/USER/Desktop/schizochytrium-hs6-genomics"
STRUCT_DIR = f"{REPO}/results/lysis_enzymes/structures"
LIGAND_DIR = f"{REPO}/results/lysis_enzymes/ligands"
DOCK_DIR = f"{REPO}/results/lysis_enzymes/docking_ph"
RESULTS_JSON = f"{REPO}/results/lysis_enzymes/ph_docking_results.json"
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
    ("P26213", "pectinase_An", "pectinase", "Aspergillus niger", "pectinase_substrate"),
    ("Q2UCU3", "betagalactosidase_Ao", "beta-galactosidase", "Aspergillus oryzae", "galactosidase_substrate"),
    ("P29600", "subtilisin_savinase_Bsp", "protease", "Bacillus sp. (Lederbergia lenta)", "protease_substrate"),
    ("P68736", "neutralprotease_Bs", "neutral protease (M4/bacillolysin)", "Bacillus subtilis", "protease_substrate"),
    ("P18429", "xylanaseC_Bs", "xylanase", "Bacillus subtilis", "xylanase_substrate"),
]

PH_VALUES = [4.0, 9.0]
EXHAUSTIVENESS = 8
N_POSES = 9


def pdb_bbox_and_center(pdb_path):
    xs, ys, zs = [], [], []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                xs.append(x); ys.append(y); zs.append(z)
    cx, cy, cz = (max(xs) + min(xs)) / 2, (max(ys) + min(ys)) / 2, (max(zs) + min(zs)) / 2
    sx, sy, sz = (max(xs) - min(xs)), (max(ys) - min(ys)), (max(zs) - min(zs))
    pad = 8.0
    return (cx, cy, cz), (sx + pad, sy + pad, sz + pad)


def prep_receptor_at_ph(label, ph):
    pdb_in = f"{STRUCT_DIR}/{label}.pdb"
    tag = f"{label}_ph{ph:g}"
    pqr_path = f"{DOCK_DIR}/{tag}.pqr"
    pdbqt_path = f"{DOCK_DIR}/{tag}_receptor.pdbqt"
    if os.path.exists(pdbqt_path):
        return pdbqt_path
    if not os.path.exists(pqr_path):
        cmd = ["pdb2pqr30", "--ff=AMBER", f"--with-ph={ph}", "--titration-state-method=propka",
               pdb_in, pqr_path]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if not os.path.exists(pqr_path):
            raise RuntimeError(f"pdb2pqr30 failed for {label} at pH {ph}: {r.stderr[-2000:]}")
    cmd2 = ["obabel", pqr_path, "-O", pdbqt_path, "-xr", "--partialcharge", "gasteiger"]
    subprocess.run(cmd2, check=True, capture_output=True, text=True)
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
        f.write(f"receptor: {receptor_pdbqt}\nligand: {ligand_pdbqt}\n")
        f.write(f"center: {center}\nbox_size: {box_size}\n")
        f.write(f"exhaustiveness: {EXHAUSTIVENESS}\nn_poses: {N_POSES}\n")
        for i, e in enumerate(energies):
            f.write(f"{i+1}\t{e[0]:.4f}\t{e[1]:.4f}\t{e[2]:.4f}\n")
    return float(energies[0][0])


def main():
    results = []
    if os.path.exists(RESULTS_JSON):
        with open(RESULTS_JSON) as f:
            results = json.load(f)
    done_keys = {(r["accession"], r["ligand"], r["ph"]) for r in results}

    for accession, label, klass, host, canonical_key in PANEL:
        print(f"\n########## {label} ({accession}) ##########")
        for ph in PH_VALUES:
            print(f"  --- pH {ph} ---")
            try:
                receptor_pdbqt = prep_receptor_at_ph(label, ph)
            except Exception as e:
                print(f"    PREP FAILED: {e}")
                for ligand_key in (canonical_key, "galactan_probe"):
                    if (accession, ligand_key, ph) in done_keys:
                        continue
                    results.append({
                        "accession": accession, "label": label, "class": klass, "host": host,
                        "ligand": ligand_key, "ph": ph, "top_affinity_kcal_mol": None,
                        "status": "PREP_FAILED",
                    })
                    with open(RESULTS_JSON, "w") as f:
                        json.dump(results, f, indent=2)
                continue

            center, box_size = pdb_bbox_and_center(f"{STRUCT_DIR}/{label}.pdb")
            for ligand_key in (canonical_key, "galactan_probe"):
                if (accession, ligand_key, ph) in done_keys:
                    print(f"    skip {ligand_key} @ pH{ph} (already docked)")
                    continue
                ligand_pdbqt = f"{LIGAND_DIR}/{ligand_key}.pdbqt"
                tag = f"{label}_ph{ph:g}_{ligand_key}"
                out_pdbqt = f"{DOCK_DIR}/{tag}.pdbqt"
                log_path = f"{DOCK_DIR}/{tag}.log"
                t0 = time.time()
                try:
                    affinity = run_vina(receptor_pdbqt, ligand_pdbqt, center, box_size, out_pdbqt, log_path)
                    elapsed = time.time() - t0
                    print(f"    {ligand_key} @ pH{ph}: {affinity:.3f} kcal/mol ({elapsed:.1f}s)")
                    results.append({
                        "accession": accession, "label": label, "class": klass, "host": host,
                        "ligand": ligand_key, "ph": ph,
                        "top_affinity_kcal_mol": round(affinity, 3),
                        "seconds": round(elapsed, 1), "status": "OK",
                    })
                except Exception as e:
                    print(f"    {ligand_key} @ pH{ph}: FAILED ({e})")
                    results.append({
                        "accession": accession, "label": label, "class": klass, "host": host,
                        "ligand": ligand_key, "ph": ph, "top_affinity_kcal_mol": None,
                        "status": "FAILED",
                    })
                with open(RESULTS_JSON, "w") as f:
                    json.dump(results, f, indent=2)

    print("\n=== DONE ===")
    print(f"{len(results)} total pH-docking entries")


if __name__ == "__main__":
    main()
