"""
Determine which of the 120 environmental photos belong to the same physical
plate (lid-on + lid-off pair, or repeated shots) vs different plates, using
actual image similarity rather than trusting capture-time gaps alone (which
turned out to merge multiple distinct plates into one apparent "burst").
"""
import os
import re

import cv2
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ENV_DIR = os.path.join(ROOT, "environmental")


def imread_unicode(path):
    return cv2.imdecode(np.fromfile(path, dtype=np.uint8), cv2.IMREAD_COLOR)


def thumbnail_gray(path, size=96):
    img = imread_unicode(path)
    img = cv2.resize(img, (size, size))
    return cv2.cvtColor(img, cv2.COLOR_BGR2GRAY).astype(np.float32)


def main():
    files = sorted(os.listdir(ENV_DIR))
    thumbs = {f: thumbnail_gray(os.path.join(ENV_DIR, f)) for f in files}

    print(f"{'file_a':30s} {'file_b':30s} {'gap_s':>6s} {'similarity':>10s}")
    sims = []
    for i in range(len(files) - 1):
        a, b = files[i], files[i + 1]
        ta = int(re.search(r"_(\d{6})\.jpg$", a).group(1))
        tb = int(re.search(r"_(\d{6})\.jpg$", b).group(1))
        def to_sec(t):
            h, m, s = t // 10000, (t // 100) % 100, t % 100
            return h * 3600 + m * 60 + s
        gap = to_sec(tb) - to_sec(ta)

        diff = np.abs(thumbs[a] - thumbs[b])
        similarity = 1.0 - (diff.mean() / 255.0)
        sims.append(similarity)
        marker = "  <-- SAME PLATE?" if similarity > 0.90 else ""
        print(f"{a:30s} {b:30s} {gap:6d} {similarity:10.3f}{marker}")

    sims = np.array(sims)
    print(f"\nsimilarity distribution: min={sims.min():.3f} p25={np.percentile(sims,25):.3f} "
          f"median={np.median(sims):.3f} p75={np.percentile(sims,75):.3f} max={sims.max():.3f}")


if __name__ == "__main__":
    main()
