"""
Segment the new GFP-transformant (strain #7) batch. Same dish-detection and
shape/isolation-filter pipeline as the original steel-bench dataset, but a
different colony-detection mask: this batch was shot under different
lighting where colonies are DARKER than agar (opposite of the steel-bench
relationship) with a yellow-green hue vs. the agar's blue-cyan hue. Local
contrast (the steel-bench approach) found essentially no usable signal here
even at generous cutoffs; direct HSV thresholding calibrated from sampled
pixel values worked cleanly. See segment_schizo_colonies.compute_colony_mask_hsv_dark.
"""
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import segment_schizo_colonies as sc

sc.MASK_FN = sc.compute_colony_mask_hsv_dark
sc.MIN_AREA = 100  # colonies in this batch skew smaller than the steel-bench dataset

ROOT = sc.ROOT
SRC_DIR = os.path.join(ROOT, "new_batch_202609", "main")
CROP_DIR = os.path.join(ROOT, "cropped_colonies_gfp7")
QC_DIR = os.path.join(ROOT, "qc_overlays_gfp7")
MANIFEST_PATH = os.path.join(ROOT, "gfp7_manifest.csv")


def main():
    files = sorted(os.listdir(SRC_DIR))
    manifest_rows = []
    for f in files:
        label = os.path.splitext(f)[0]
        path = os.path.join(SRC_DIR, f)
        qc_path = os.path.join(QC_DIR, f"{label}_overlay.jpg")
        records = sc.process_image(path, CROP_DIR, qc_path, label_prefix=label)
        manifest_rows.extend(records)
        print(f"  {label}: kept {len(records)} colonies")

    with open(MANIFEST_PATH, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=[
            "crop_path", "source_image", "x", "y", "w", "h",
            "area", "circularity", "solidity",
        ])
        writer.writeheader()
        writer.writerows(manifest_rows)

    print(f"\nTotal colonies: {len(manifest_rows)} from {len(files)} plates")
    print(f"Manifest: {MANIFEST_PATH}")


if __name__ == "__main__":
    main()
