"""
Score and rank every detected environmental colony by how closely its
color/texture profile matches the Schizochytrium reference strains.

Reference signal comes from two sources, pooled:
  - reference_manifest.csv: true single-colony crops (currently just 1,
    from strain7 -- isolation was poor across the reference plates, per
    the user's note that lab strains weren't necessarily under optimal
    conditions for this).
  - reference_lawn_features.csv: whole-growth-region features for all 5
    strains, used as a fallback signature where single colonies aren't
    available. This is the majority of the reference signal here.

This is a similarity/anomaly score, not a trained classifier -- with this
little reference data (5 plates, mostly whole-region rather than per-
colony), a real classifier would be overfit noise. A per-feature z-score
distance against the pooled reference distribution is a more honest match
to how much signal actually exists.
"""
import csv
import os
import sys

import cv2
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import segment_schizo_colonies as sc

ROOT = sc.ROOT
ENV_MANIFEST = os.path.join(ROOT, "environmental_manifest.csv")
REF_MANIFEST = os.path.join(ROOT, "reference_manifest.csv")
REF_LAWN = os.path.join(ROOT, "reference_lawn_features.csv")
OUT_CSV = os.path.join(ROOT, "environmental_colony_ranking.csv")

FEATURES = ["mean_hue", "mean_sat", "mean_val", "std_hue", "std_sat", "std_val", "texture_std"]


def extract_features_from_crop(path):
    img = sc.imread_unicode(path)
    if img is None:
        return None
    hsv = cv2.cvtColor(img, cv2.COLOR_BGR2HSV)
    gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    mean_hsv, std_hsv = cv2.meanStdDev(hsv)
    texture_std = cv2.meanStdDev(gray)[1].flatten()[0]
    h, w = img.shape[:2]
    equiv_diam = (h + w) / 2  # crop is already tight to the colony; use crop size as a diameter proxy
    return dict(
        mean_hue=mean_hsv[0, 0], mean_sat=mean_hsv[1, 0], mean_val=mean_hsv[2, 0],
        std_hue=std_hsv[0, 0], std_sat=std_hsv[1, 0], std_val=std_hsv[2, 0],
        texture_std=texture_std, equiv_diameter=equiv_diam,
    )


def build_reference_distribution():
    """Pool per-colony reference features (n=1) and per-strain lawn features (n=5) into one reference sample set."""
    samples = []

    if os.path.exists(REF_MANIFEST):
        ref_rows = list(csv.DictReader(open(REF_MANIFEST, encoding="utf-8")))
        for r in ref_rows:
            feat = extract_features_from_crop(os.path.join(ROOT, r["crop_path"]))
            if feat:
                samples.append(feat)

    if os.path.exists(REF_LAWN):
        lawn_rows = list(csv.DictReader(open(REF_LAWN, encoding="utf-8")))
        for r in lawn_rows:
            samples.append({f: float(r[f]) for f in FEATURES})

    return samples


def main():
    ref_samples = build_reference_distribution()
    print(f"Reference samples pooled: {len(ref_samples)} "
          f"(single colonies + per-strain lawn regions)")
    if len(ref_samples) < 2:
        print("Not enough reference data to build a distribution. Aborting.")
        return

    ref_mean = {f: np.mean([s[f] for s in ref_samples]) for f in FEATURES}
    ref_std = {f: max(np.std([s[f] for s in ref_samples]), 1e-6) for f in FEATURES}
    print("Reference mean/std per feature:")
    for f in FEATURES:
        print(f"  {f:12s} mean={ref_mean[f]:7.2f} std={ref_std[f]:7.2f}")

    env_rows = list(csv.DictReader(open(ENV_MANIFEST, encoding="utf-8")))
    print(f"\nScoring {len(env_rows)} environmental colony crops...")

    scored = []
    for r in env_rows:
        feat = extract_features_from_crop(os.path.join(ROOT, r["crop_path"]))
        if feat is None:
            continue
        z_sq_sum = sum(((feat[f] - ref_mean[f]) / ref_std[f]) ** 2 for f in FEATURES)
        distance = np.sqrt(z_sq_sum / len(FEATURES))  # RMS z-score distance
        scored.append(dict(
            plate_id=r["plate_id"], crop_path=r["crop_path"],
            schizo_likeness_distance=round(distance, 3),
            **{f: round(feat[f], 1) for f in FEATURES},
        ))

    scored.sort(key=lambda d: d["schizo_likeness_distance"])
    fieldnames = ["plate_id", "crop_path", "schizo_likeness_distance"] + FEATURES
    with open(OUT_CSV, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(scored)

    print(f"\nWrote ranking to {OUT_CSV}")
    print("\nTop 15 most Schizochytrium-like candidates:")
    for row in scored[:15]:
        print(f"  plate {row['plate_id']:>3} dist={row['schizo_likeness_distance']:.2f}  {row['crop_path']}")


if __name__ == "__main__":
    main()
