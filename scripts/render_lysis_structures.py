"""
Render one static cartoon PNG per folded enzyme in the lysis panel, colored by
per-residue pLDDT (stored as B-factor by run_esmfold_lysis_panel.py), using the
standard AlphaFold/ESMFold confidence palette:
  >90 dark blue, 70-90 light blue/cyan, 50-70 yellow, <50 orange.
White background, ray-traced, matching the site's existing figures/esmfold/ style
(white bg, cartoon, ray-traced) but with a pLDDT spectrum instead of flat color
since these are single fresh predictions (not a mutation-comparison figure).

Run with: pymol -cq render_lysis_structures.py
"""
from pymol import cmd

REPO = "/mnt/c/Users/USER/Desktop/schizochytrium-hs6-genomics"
STRUCT_DIR = f"{REPO}/results/lysis_enzymes/structures"
OUT_DIR = f"{REPO}/figures/lysis_enzymes"

LABELS = [
    "subtilisin_carlsberg", "subtilisin_BPN", "subtilisin_E", "subtilisin_pumilus",
    "lysozyme_C_chicken", "endolysin_T4", "lysozyme_M1_strepto",
    "endoglucanase_I_Tr", "endoglucanase_B_An",
    "xylanase2_Tr", "xylanaseB_An",
    "betaagaraseA_Zg", "betaagaraseB_Zg",
]

import os
os.makedirs(OUT_DIR, exist_ok=True)


def render_one(label):
    pdb = f"{STRUCT_DIR}/{label}.pdb"
    if not os.path.exists(pdb):
        print(f"SKIP {label}: no pdb")
        return
    cmd.reinitialize()
    cmd.load(pdb, "prot")
    cmd.hide("everything")
    cmd.show("cartoon", "prot")
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.set("cartoon_transparency", 0)
    cmd.set("antialias", 2)
    cmd.set("ray_trace_mode", 1)
    cmd.set("ray_trace_color", "black")
    cmd.set("cartoon_fancy_helices", 1)
    # AlphaFold/ESMFold-style pLDDT coloring (b-factor stores plddt as 0-1 fraction)
    cmd.spectrum("b", "orange_yellow_cyan_blue", "prot", minimum=0.0, maximum=1.0)
    cmd.orient("prot")
    cmd.zoom("prot", buffer=4)
    cmd.set("ray_shadows", 0)
    png_path = f"{OUT_DIR}/{label}_structure.png"
    cmd.ray(1200, 900)
    cmd.png(png_path, dpi=150)
    print(f"saved {png_path}")


for label in LABELS:
    render_one(label)

print("DONE rendering structure thumbnails")
