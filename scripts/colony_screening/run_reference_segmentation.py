import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import segment_schizo_colonies as sc

ROOT = sc.ROOT
REF_DIR = os.path.join(ROOT, "reference")
CROP_DIR = os.path.join(ROOT, "cropped_colonies_reference")
QC_DIR = os.path.join(ROOT, "qc_overlays_reference")
MANIFEST_PATH = os.path.join(ROOT, "reference_manifest.csv")


def main():
    strains = sorted(os.listdir(REF_DIR))
    print("Reference strains:", strains)
    manifest_rows = []
    for strain in strains:
        sdir = os.path.join(REF_DIR, strain)
        files = sorted(os.listdir(sdir))
        # process the lid-off photo only (2nd file, index 1) -- less glare
        lid_off = files[1]
        path = os.path.join(sdir, lid_off)
        out_crop_dir = os.path.join(CROP_DIR, strain)
        qc_path = os.path.join(QC_DIR, f"{strain}_overlay.jpg")
        records = sc.process_image(path, out_crop_dir, qc_path, label_prefix=strain)
        for r in records:
            r["strain"] = strain
        manifest_rows.extend(records)
        print(f"  {strain} ({lid_off}): kept {len(records)} colonies")

    with open(MANIFEST_PATH, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "strain", "crop_path", "source_image", "x", "y", "w", "h",
            "area", "circularity", "solidity",
        ])
        writer.writeheader()
        writer.writerows(manifest_rows)
    print(f"\nTotal reference colonies: {len(manifest_rows)}")
    print(f"Manifest: {MANIFEST_PATH}")


if __name__ == "__main__":
    main()
