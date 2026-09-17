"""
Real validation of the "Schizochytrium colony" reference profile, now that
we have 297 same-type (single-colony) crops of a confirmed strain (#7)
instead of just one. Yesterday's leave-one-out check failed (n=1, and mixed
colony-type with lawn-type features) -- this repeats that logic properly:

For each of the ~7 plates contributing colonies, build a reference from the
OTHER plates' colonies, then check whether the held-out plate's own colonies
score closer to that reference than a random sample of environmental
colonies does. If real #7 colonies consistently score closer to a
#7-built reference than environmental colonies do, that's evidence the
feature vector carries real, self-consistent colony-identity signal (a
necessary condition for the whole ranking approach to mean anything) --
not proof any specific environmental colony IS Schizochytrium, but proof
the measurement isn't noise.
"""
import csv
import os
import random
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import segment_schizo_colonies as sc

ROOT = sc.ROOT
GFP7_MANIFEST = os.path.join(ROOT, "gfp7_manifest.csv")
ENV_RANKING = os.path.join(ROOT, "environmental_colony_ranking_cleaned.csv")
FEATURES = ["mean_hue", "mean_sat", "mean_val", "std_hue", "std_sat", "std_val", "texture_std"]


def extract(path):
    import cv2
    img = sc.imread_unicode(path)
    hsv = cv2.cvtColor(img, cv2.COLOR_BGR2HSV)
    gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    mean_hsv, std_hsv = cv2.meanStdDev(hsv)
    texture_std = cv2.meanStdDev(gray)[1].flatten()[0]
    return dict(mean_hue=mean_hsv[0, 0], mean_sat=mean_hsv[1, 0], mean_val=mean_hsv[2, 0],
                std_hue=std_hsv[0, 0], std_sat=std_hsv[1, 0], std_val=std_hsv[2, 0],
                texture_std=texture_std)


def dist_to_ref(feat, ref_mean, ref_std):
    return np.sqrt(np.mean([((feat[f] - ref_mean[f]) / ref_std[f]) ** 2 for f in FEATURES]))


def main():
    rows = list(csv.DictReader(open(GFP7_MANIFEST, encoding="utf-8")))
    by_plate = {}
    for r in rows:
        plate = os.path.basename(r["source_image"]).split(".")[0]
        feat = extract(os.path.join(ROOT, r["crop_path"]))
        by_plate.setdefault(plate, []).append(feat)

    plates = {p: v for p, v in by_plate.items() if len(v) >= 5}
    print(f"Plates with >=5 colonies: {list(plates.keys())} (n={[len(v) for v in plates.values()]})")

    env_rows = list(csv.DictReader(open(ENV_RANKING, encoding="utf-8")))
    random.seed(0)
    env_sample = random.sample(env_rows, min(200, len(env_rows)))
    env_feats = [{f: float(r[f]) for f in FEATURES} for r in env_sample]

    print(f"\n{'held-out plate':16s} {'n':>4s} {'median dist (own)':>18s} {'median dist (env)':>18s} {'separation?'}")
    for held_out in plates:
        train_feats = [f for p, feats in plates.items() if p != held_out for f in feats]
        ref_mean = {f: np.mean([x[f] for x in train_feats]) for f in FEATURES}
        ref_std = {f: max(np.std([x[f] for x in train_feats]), 1e-6) for f in FEATURES}

        own_dists = [dist_to_ref(f, ref_mean, ref_std) for f in plates[held_out]]
        env_dists = [dist_to_ref(f, ref_mean, ref_std) for f in env_feats]

        med_own = np.median(own_dists)
        med_env = np.median(env_dists)
        sep = "YES (own closer)" if med_own < med_env else "no (env closer or tied)"
        print(f"{held_out:16s} {len(own_dists):4d} {med_own:18.3f} {med_env:18.3f}   {sep}")


if __name__ == "__main__":
    main()
