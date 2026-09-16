"""
Fallback reference signature for strains where single colonies are too
scarce/absent (thick lawn/streak growth, per the user's note that lab
strains may not be under optimal isolation conditions).

Instead of individual colonies, this computes the same style of
color/texture descriptors over the WHOLE growth mass (streak + lawn +
colonies, whatever is present) as one pooled region per strain. This can't
support per-colony classification, but gives a "what does this strain's
growth generally look like" signature that can still contribute to the
Schizochytrium-likeness reference profile.
"""
import csv
import os
import sys

import cv2
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import segment_schizo_colonies as sc

ROOT = sc.ROOT
REF_DIR = os.path.join(ROOT, "reference")
OUT_CSV = os.path.join(ROOT, "reference_lawn_features.csv")


def extract_lawn_features(path):
    img = sc.imread_unicode(path)
    hsv = cv2.cvtColor(img, cv2.COLOR_BGR2HSV)
    dish_roi, _ = sc.find_dish_roi(img)

    v = hsv[:, :, 2].astype(np.float32)
    blur = cv2.GaussianBlur(v, (0, 0), sigmaX=sc.LOCAL_CONTRAST_BLUR_SIGMA)
    local_contrast = v - blur
    growth_mask = np.where(local_contrast > sc.LOCAL_CONTRAST_CUTOFF, 255, 0).astype(np.uint8)
    growth_mask = cv2.bitwise_and(growth_mask, growth_mask, mask=dish_roi)
    growth_mask = cv2.morphologyEx(growth_mask, cv2.MORPH_OPEN, np.ones((5, 5), np.uint8))
    growth_mask = cv2.morphologyEx(growth_mask, cv2.MORPH_CLOSE, np.ones((9, 9), np.uint8))

    growth_area = int((growth_mask > 0).sum())
    dish_area = int((dish_roi > 0).sum())
    if growth_area < 500:
        return None

    mean_hsv, std_hsv = cv2.meanStdDev(hsv, mask=growth_mask)
    gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    texture_std = cv2.meanStdDev(gray, mask=growth_mask)[1].flatten()[0]

    return dict(
        growth_fraction=round(growth_area / dish_area, 4),
        mean_hue=round(mean_hsv[0, 0], 1), mean_sat=round(mean_hsv[1, 0], 1), mean_val=round(mean_hsv[2, 0], 1),
        std_hue=round(std_hsv[0, 0], 1), std_sat=round(std_hsv[1, 0], 1), std_val=round(std_hsv[2, 0], 1),
        texture_std=round(float(texture_std), 1),
    )


def main():
    rows = []
    for strain in sorted(os.listdir(REF_DIR)):
        sdir = os.path.join(REF_DIR, strain)
        files = sorted(os.listdir(sdir))
        lid_off = files[1]
        feat = extract_lawn_features(os.path.join(sdir, lid_off))
        if feat is None:
            print(f"{strain}: no usable growth region found")
            continue
        feat["strain"] = strain
        feat["source_image"] = os.path.join("reference", strain, lid_off)
        rows.append(feat)
        print(f"{strain}: growth_fraction={feat['growth_fraction']:.3f} "
              f"mean_hue={feat['mean_hue']} mean_sat={feat['mean_sat']} mean_val={feat['mean_val']} "
              f"texture_std={feat['texture_std']}")

    fieldnames = ["strain", "source_image", "growth_fraction",
                  "mean_hue", "mean_sat", "mean_val", "std_hue", "std_sat", "std_val", "texture_std"]
    with open(OUT_CSV, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)
    print(f"\nWrote {len(rows)} rows to {OUT_CSV}")


if __name__ == "__main__":
    main()
