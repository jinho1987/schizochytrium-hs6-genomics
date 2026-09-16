"""
Sanity check on the Schizochytrium-likeness ranking (rank_environmental_colonies.py):
does the scorer actually rank a KNOWN true positive favorably?

Leave-one-out: score the one genuine isolated reference colony (strain7) using
a reference distribution built WITHOUT it (the other 4 strains' whole-growth-
region lawn features only), then see where that known true positive would
land among the 555 environmental candidates.

Result (as of this run): rank 552 of 556 (99.3rd percentile) -- the scorer
would place a known real Schizochytrium colony almost dead last. Traced to
pooling whole-lawn/streak-region features with a single-colony feature as if
they were the same kind of measurement, and to having only one single-colony
reference sample in existence (not enough to build a same-type reference
distribution). See colony-screening.html Section 03 for the write-up.
"""
import csv
import os
import sys

import cv2
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import segment_schizo_colonies as sc

ROOT = sc.ROOT
FEATURES = ["mean_hue", "mean_sat", "mean_val", "std_hue", "std_sat", "std_val", "texture_std"]
TRUE_POSITIVE_CROP = os.path.join(ROOT, "cropped_colonies_reference", "strain7", "strain7__colony001.jpg")


def extract(path):
    img = sc.imread_unicode(path)
    hsv = cv2.cvtColor(img, cv2.COLOR_BGR2HSV)
    gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    mean_hsv, std_hsv = cv2.meanStdDev(hsv)
    texture_std = cv2.meanStdDev(gray)[1].flatten()[0]
    return dict(
        mean_hue=mean_hsv[0, 0], mean_sat=mean_hsv[1, 0], mean_val=mean_hsv[2, 0],
        std_hue=std_hsv[0, 0], std_sat=std_hsv[1, 0], std_val=std_hsv[2, 0],
        texture_std=texture_std,
    )


def main():
    lawn_rows = list(csv.DictReader(open(os.path.join(ROOT, "reference_lawn_features.csv"), encoding="utf-8")))
    other_strains = [r for r in lawn_rows if r["strain"] != "strain7"]
    ref_mean = {f: np.mean([float(r[f]) for r in other_strains]) for f in FEATURES}
    ref_std = {f: max(np.std([float(r[f]) for r in other_strains]), 1e-6) for f in FEATURES}

    feat = extract(TRUE_POSITIVE_CROP)
    dist = np.sqrt(np.mean([((feat[f] - ref_mean[f]) / ref_std[f]) ** 2 for f in FEATURES]))
    print("Known true-positive (strain7 colony) leave-one-out distance:", round(dist, 3))

    env_rows = list(csv.DictReader(open(os.path.join(ROOT, "environmental_colony_ranking_cleaned.csv"), encoding="utf-8")))
    env_dists = sorted(float(r["schizo_likeness_distance"]) for r in env_rows)
    rank = sum(1 for d in env_dists if d < dist)
    pct = 100 * rank / len(env_dists)
    print(f"Would rank #{rank + 1} of {len(env_dists) + 1} ({pct:.1f}th percentile) among environmental candidates")
    print("(lower distance = more Schizochytrium-like; lower percentile = better)")


if __name__ == "__main__":
    main()
