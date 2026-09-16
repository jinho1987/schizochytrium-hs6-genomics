"""
Draw the top-ranked Schizochytrium-like candidates directly on their
original full plate photo, labeled with rank + distance, so a candidate
found in the ranking can be located back on the physical plate.
"""
import csv
import os
import sys
from collections import defaultdict

import cv2

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import segment_schizo_colonies as sc

ROOT = sc.ROOT
RANKING_CSV = os.path.join(ROOT, "environmental_colony_ranking_cleaned.csv")
ENV_MANIFEST = os.path.join(ROOT, "environmental_manifest.csv")
OUT_DIR = os.path.join(ROOT, "candidate_collages", "annotated_plates")
TOP_N = 60


def main():
    ranking_rows = list(csv.DictReader(open(RANKING_CSV, encoding="utf-8")))[:TOP_N]
    manifest_rows = list(csv.DictReader(open(ENV_MANIFEST, encoding="utf-8")))
    bbox_by_crop = {r["crop_path"]: r for r in manifest_rows}

    by_plate = defaultdict(list)
    for rank, r in enumerate(ranking_rows, start=1):
        bbox_row = bbox_by_crop.get(r["crop_path"])
        if bbox_row is None:
            continue
        by_plate[r["plate_id"]].append((rank, r, bbox_row))

    os.makedirs(OUT_DIR, exist_ok=True)
    for plate_id, items in sorted(by_plate.items(), key=lambda kv: -len(kv[1])):
        source_image = os.path.join(ROOT, items[0][2]["source_image"])
        img = sc.imread_unicode(source_image)
        if img is None:
            print(f"could not read {source_image}")
            continue

        for rank, r, bbox in items:
            x, y, w, h = int(bbox["x"]), int(bbox["y"]), int(bbox["w"]), int(bbox["h"])
            pad = 30
            cv2.rectangle(img, (x - pad, y - pad), (x + w + pad, y + h + pad), (0, 255, 0), 8)
            label = f"#{rank} d={float(r['schizo_likeness_distance']):.2f}"
            (tw, th), _ = cv2.getTextSize(label, cv2.FONT_HERSHEY_SIMPLEX, 1.6, 4)
            ly = max(0, y - pad - 15)
            cv2.rectangle(img, (x - pad, ly - th - 10), (x - pad + tw + 10, ly + 5), (0, 255, 0), -1)
            cv2.putText(img, label, (x - pad + 5, ly), cv2.FONT_HERSHEY_SIMPLEX, 1.6, (0, 0, 0), 4, cv2.LINE_AA)

        small = cv2.resize(img, (img.shape[1] // 3, img.shape[0] // 3))
        out_path = os.path.join(OUT_DIR, f"plate{plate_id}_annotated.jpg")
        sc.imwrite_unicode(out_path, small)
        print(f"plate {plate_id}: {len(items)} candidates marked -> {out_path}")


if __name__ == "__main__":
    main()
