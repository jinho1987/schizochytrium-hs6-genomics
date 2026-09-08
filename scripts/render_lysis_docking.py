"""
Render static top-docked-pose PNGs for the lysis-enzyme panel: receptor cartoon
(colored by pLDDT, same palette as render_lysis_structures.py) + top Vina pose
ligand as sticks (orange carbons), white ray-traced background.

Also renders a 360-degree orbiting GIF for the two "single data point" headline
call-outs (subtilisin E reversal, cellulase near-parity), per the task spec.

Run with: pymol -cq render_lysis_docking.py
"""
import os
from pymol import cmd
from PIL import Image

REPO = "/mnt/c/Users/USER/Desktop/schizochytrium-hs6-genomics"
STRUCT_DIR = f"{REPO}/results/lysis_enzymes/structures"
DOCK_DIR = f"{REPO}/results/lysis_enzymes/docking"
OUT_DIR = f"{REPO}/figures/lysis_enzymes"
TMP_DIR = f"{REPO}/results/lysis_enzymes/docking/_gif_frames"
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(TMP_DIR, exist_ok=True)

# label -> (ligand_key, make_gif)
JOBS = {
    "subtilisin_carlsberg":  ("galactan_probe", False),
    "subtilisin_BPN":        ("galactan_probe", False),
    "subtilisin_E":          ("galactan_probe", True),   # headline: reversal
    "subtilisin_pumilus":    ("galactan_probe", False),
    "lysozyme_C_chicken":    ("galactan_probe", False),
    "endolysin_T4":          ("galactan_probe", False),
    "lysozyme_M1_strepto":   ("galactan_probe", False),
    "endoglucanase_I_Tr":    ("galactan_probe", False),
    "endoglucanase_B_An":    ("galactan_probe", True),   # headline: near-parity
    "xylanase2_Tr":          ("galactan_probe", False),
    "xylanaseB_An":          ("galactan_probe", False),
    "betaagaraseA_Zg":       ("agarase_substrate", False),  # headline: best mechanistic match (canonical)
    "betaagaraseB_Zg":       ("galactan_probe", False),
}


def setup_scene(label, ligand_key):
    top_pose_pdb = f"{DOCK_DIR}/{label}_{ligand_key}_top_pose.pdb"
    receptor_pdb = f"{STRUCT_DIR}/{label}.pdb"
    if not (os.path.exists(top_pose_pdb) and os.path.exists(receptor_pdb)):
        print(f"SKIP {label}/{ligand_key}: missing {top_pose_pdb if not os.path.exists(top_pose_pdb) else receptor_pdb}")
        return False
    cmd.reinitialize()
    cmd.load(receptor_pdb, "receptor")
    cmd.load(top_pose_pdb, "ligand")
    cmd.hide("everything")
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.set("antialias", 2)
    cmd.set("ray_shadows", 0)
    cmd.set("cartoon_fancy_helices", 1)

    cmd.show("cartoon", "receptor")
    cmd.spectrum("b", "orange_yellow_cyan_blue", "receptor", minimum=0.0, maximum=1.0)
    cmd.set("cartoon_transparency", 0.15, "receptor")

    cmd.show("sticks", "ligand")
    cmd.color("gray30", "ligand and elem C")
    cmd.util.cnc("ligand")
    cmd.set("stick_radius", 0.25, "ligand")

    cmd.orient("receptor")
    cmd.zoom("all", buffer=6)
    return True


def render_static(label, ligand_key):
    if not setup_scene(label, ligand_key):
        return
    out_png = f"{OUT_DIR}/{label}_docking_top_pose.png"
    cmd.ray(1400, 1000)
    cmd.png(out_png, dpi=150)
    print(f"saved {out_png}")


def render_gif(label, ligand_key, n_frames=36):
    if not setup_scene(label, ligand_key):
        return
    frame_paths = []
    for i in range(n_frames):
        cmd.turn("y", 360.0 / n_frames)
        cmd.ray(700, 500)
        fp = f"{TMP_DIR}/{label}_{i:03d}.png"
        cmd.png(fp, dpi=100)
        frame_paths.append(fp)
    frames = [Image.open(p).convert("RGB") for p in frame_paths]
    out_gif = f"{OUT_DIR}/{label}_docking_top_pose.gif"
    frames[0].save(
        out_gif, save_all=True, append_images=frames[1:],
        duration=90, loop=0, optimize=True,
    )
    print(f"saved {out_gif} ({len(frames)} frames)")
    for p in frame_paths:
        os.remove(p)
    # also save a static PNG of the starting orientation for this enzyme
    static_first = frame_paths[0]
    out_png = f"{OUT_DIR}/{label}_docking_top_pose.png"
    # re-render the static at full resolution from the original (pre-rotation) orientation
    setup_scene(label, ligand_key)
    cmd.ray(1400, 1000)
    cmd.png(out_png, dpi=150)
    print(f"saved {out_png}")


for label, (ligand_key, make_gif) in JOBS.items():
    print(f"\n=== {label} vs {ligand_key} (gif={make_gif}) ===")
    if make_gif:
        render_gif(label, ligand_key)
    else:
        render_static(label, ligand_key)

print("\nDONE rendering docking pose figures")
