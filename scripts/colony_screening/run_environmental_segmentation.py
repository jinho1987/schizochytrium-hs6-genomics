import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import segment_schizo_colonies as sc

ROOT = sc.ROOT
ENV_DIR = os.path.join(ROOT, "environmental")
PLATES_CSV = os.path.join(ROOT, "environmental_plates.csv")
CROP_DIR = os.path.join(ROOT, "cropped_colonies_environmental")
QC_DIR = os.path.join(ROOT, "qc_overlays_environmental")
MANIFEST_PATH = os.path.join(ROOT, "environmental_manifest.csv")


def main():
    plates = list(csv.DictReader(open(PLATES_CSV, encoding="utf-8")))
    manifest_rows = []
    total = 0
    for p in plates:
        plate_id = int(p["plate_id"])
        lid_off_path = os.path.join(ROOT, p["lid_off"])
        label = f"plate{plate_id:03d}"
        qc_path = os.path.join(QC_DIR, f"{label}_overlay.jpg")
        records = sc.process_image(lid_off_path, CROP_DIR, qc_path, label_prefix=label)
        for r in records:
            r["plate_id"] = plate_id
        manifest_rows.extend(records)
        total += len(records)
        print(f"  plate {plate_id:3d}: kept {len(records)} colonies")

    with open(MANIFEST_PATH, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "plate_id", "crop_path", "source_image", "x", "y", "w", "h",
            "area", "circularity", "solidity",
        ])
        writer.writeheader()
        writer.writerows(manifest_rows)

    print(f"\nTotal environmental colonies: {total} from {len(plates)} plates")
    print(f"Manifest: {MANIFEST_PATH}")


if __name__ == "__main__":
    main()
