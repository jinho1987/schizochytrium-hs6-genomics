"""
Render static top-docked-pose PNGs (+ 360-degree orbiting GIFs) for the full
18-enzyme lysis panel: receptor cartoon (colored by pLDDT, same palette as
render_lysis_structures.py) + top Vina pose ligand as sticks (orange carbons),
white ray-traced background.

Two kinds of figure:
  1. PRIMARY_JOBS -- every one of the 18 enzymes' top pose vs. the shared
     galactan-probe ligand, always with BOTH a static PNG and a 36-frame GIF
     (this used to be limited to 2 "headline" enzymes; now every enzyme gets
     the same treatment per the site's full-coverage requirement).
  2. CANONICAL_JOBS -- an additional static-only pose vs. the enzyme's own
     canonical substrate, for the enzymes where that pose is scientifically
     interesting to show as a real, on-target mechanism: beta-agarase A
     (already rendered in an earlier run, kept as-is) plus the 3 new
     "new mechanistic class" enzymes added in the wet-lab-mirroring expansion
     (pectinase, beta-galactosidase, the B. subtilis neutral protease).

Output naming:
  - {receptor}_docking_top_pose.png / .gif   -> galactan-probe pose (all 18)
  - {receptor}_canonical_docking_top_pose.png -> canonical-substrate pose
    (agarase A's pre-existing canonical figure keeps its ORIGINAL filename,
    {label}_docking_top_pose.png, untouched -- see AGARASE_A_LEGACY note below)

GIF frames now go to their own per-enzyme subdirectory under
results/lysis_enzymes/docking/_gif_frames/<label>/ instead of one shared flat
directory. This was a deliberate fix: running two render passes concurrently
previously raced on the same flat _gif_frames/ dir and one process's cleanup
os.remove() could hit FileNotFoundError because the other process had already
deleted a same-named frame file. Per-enzyme subdirectories make this robust
even if the script is ever invoked concurrently again -- though the
recommended practice is still to only run one instance at a time.

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

# AGARASE_A_LEGACY: betaagaraseA_Zg's original job rendered its CANONICAL pose
# (vs agarase_substrate / neoagarohexaose) at the bare filename
# betaagaraseA_Zg_docking_top_pose.png, and section 02 of lysis-enzymes.html
# references that exact file with a caption describing it as the canonical
# pose. That file/caption is left untouched. Its NEW galactan-probe pose
# (required now that all 18 enzymes get one) is rendered below under the
# distinct out_label "betaagaraseA_Zg_galactan" so it cannot collide with or
# overwrite the legacy canonical figure.

# receptor_label -> (ligand_key, make_gif, out_label_override_or_None)
PRIMARY_JOBS = [
    ("subtilisin_carlsberg",  "galactan_probe", True, None),
    ("subtilisin_BPN",        "galactan_probe", True, None),
    ("subtilisin_E",          "galactan_probe", True, None),   # already had a GIF; re-run is a no-op (files exist)
    ("subtilisin_pumilus",    "galactan_probe", True, None),
    ("lysozyme_C_chicken",    "galactan_probe", True, None),
    ("endolysin_T4",          "galactan_probe", True, None),
    ("lysozyme_M1_strepto",   "galactan_probe", True, None),
    ("endoglucanase_I_Tr",    "galactan_probe", True, None),
    ("endoglucanase_B_An",    "galactan_probe", True, None),   # already had a GIF; re-run is a no-op (files exist)
    ("xylanase2_Tr",          "galactan_probe", True, None),
    ("xylanaseB_An",          "galactan_probe", True, None),
    ("betaagaraseA_Zg",       "galactan_probe", True, "betaagaraseA_Zg_galactan"),  # see AGARASE_A_LEGACY above
    ("betaagaraseB_Zg",       "galactan_probe", True, None),
    # 5 new enzymes (wet-lab Experiment 4 mirror) -- galactan-probe pose + GIF
    ("pectinase_An",             "galactan_probe", True, None),
    ("betagalactosidase_Ao",     "galactan_probe", True, None),
    ("subtilisin_savinase_Bsp",  "galactan_probe", True, None),
    ("neutralprotease_Bs",       "galactan_probe", True, None),
    ("xylanaseC_Bs",              "galactan_probe", True, None),
]

# Additional static-only pose vs the enzyme's OWN canonical substrate, for
# enzymes where that's a scientifically interesting on-target mechanism to
# show (mirrors the pre-existing beta-agarase A treatment).
CANONICAL_JOBS = [
    ("pectinase_An",         "pectinase_substrate",     "pectinase_An_canonical"),
    ("betagalactosidase_Ao", "galactosidase_substrate", "betagalactosidase_Ao_canonical"),
    ("neutralprotease_Bs",   "protease_substrate",      "neutralprotease_Bs_canonical"),
]


def setup_scene(receptor_label, ligand_key):
    top_pose_pdb = f"{DOCK_DIR}/{receptor_label}_{ligand_key}_top_pose.pdb"
    receptor_pdb = f"{STRUCT_DIR}/{receptor_label}.pdb"
    if not (os.path.exists(top_pose_pdb) and os.path.exists(receptor_pdb)):
        print(f"SKIP {receptor_label}/{ligand_key}: missing {top_pose_pdb if not os.path.exists(top_pose_pdb) else receptor_pdb}")
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


def render_static(receptor_label, ligand_key, out_label=None):
    out_label = out_label or receptor_label
    out_png = f"{OUT_DIR}/{out_label}_docking_top_pose.png"
    if os.path.exists(out_png):
        print(f"SKIP static {out_label}: already rendered")
        return
    if not setup_scene(receptor_label, ligand_key):
        return
    cmd.ray(1400, 1000)
    cmd.png(out_png, dpi=150)
    print(f"saved {out_png}")


def render_gif(receptor_label, ligand_key, out_label=None, n_frames=36):
    out_label = out_label or receptor_label
    out_gif = f"{OUT_DIR}/{out_label}_docking_top_pose.gif"
    out_png = f"{OUT_DIR}/{out_label}_docking_top_pose.png"
    if os.path.exists(out_gif) and os.path.exists(out_png):
        print(f"SKIP gif {out_label}: already rendered")
        return
    if not setup_scene(receptor_label, ligand_key):
        return
    frame_dir = f"{TMP_DIR}/{out_label}"
    os.makedirs(frame_dir, exist_ok=True)
    frame_paths = []
    for i in range(n_frames):
        cmd.turn("y", 360.0 / n_frames)
        cmd.ray(700, 500)
        fp = f"{frame_dir}/{out_label}_{i:03d}.png"
        cmd.png(fp, dpi=100)
        frame_paths.append(fp)
    frames = [Image.open(p).convert("RGB") for p in frame_paths]
    frames[0].save(
        out_gif, save_all=True, append_images=frames[1:],
        duration=90, loop=0, optimize=True,
    )
    print(f"saved {out_gif} ({len(frames)} frames)")
    for p in frame_paths:
        try:
            os.remove(p)
        except FileNotFoundError:
            pass  # another concurrent run may have already cleaned this frame
    try:
        os.rmdir(frame_dir)
    except OSError:
        pass  # not empty / already removed by a concurrent run -- harmless
    # also save a static PNG of the starting (pre-rotation) orientation
    setup_scene(receptor_label, ligand_key)
    cmd.ray(1400, 1000)
    cmd.png(out_png, dpi=150)
    print(f"saved {out_png}")


for receptor_label, ligand_key, make_gif, out_label in PRIMARY_JOBS:
    label = out_label or receptor_label
    print(f"\n=== {receptor_label} vs {ligand_key} -> {label} (gif={make_gif}) ===")
    if make_gif:
        render_gif(receptor_label, ligand_key, out_label)
    else:
        render_static(receptor_label, ligand_key, out_label)

for receptor_label, ligand_key, out_label in CANONICAL_JOBS:
    print(f"\n=== {receptor_label} vs {ligand_key} (canonical, static-only) -> {out_label} ===")
    render_static(receptor_label, ligand_key, out_label)

print("\nDONE rendering docking pose figures")
