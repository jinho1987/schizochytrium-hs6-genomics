"""
Detect well-isolated single Schizochytrium-type colonies on plate photos
(reference lab strains, and environmental sample plates) and crop them out.

Adapted from the Chlorella pipeline (../../colony_pics/scripts/segment_colonies.py)
for a different organism and imaging setup:
  - Background is a reflective, scratched stainless-steel bench (highly
    variable brightness/color), NOT a uniform pale backdrop -- so the dish
    is found via Hough circle detection on grayscale, not color thresholding.
  - Colonies are cream/tan and BRIGHT (high value, moderate-high saturation,
    warm hue ~0-35), vs. the duller, lower-saturation agar background --
    the opposite brightness relationship from Chlorella's dark-on-light
    green colonies.
  - Images are much higher resolution (3000x4000) with a physically larger
    dish fraction of the frame.

Same shape-based philosophy as before: accept only blobs that are round,
solid, appropriately sized, and -- critically -- whose crop window contains
no other detected blob (rejects streak fragments, handwriting, and
neighboring colonies all in one check).
"""
import csv
import os
from dataclasses import dataclass

import cv2
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # schizo_pics/

# ---- colony detection ---------------------------------------------------
# Plain HSV thresholding doesn't work: the agar has a warm cast under this
# lighting, overlapping colony color across differently-exposed shots.
# Global Otsu on raw brightness within the dish doesn't work either: a
# handheld phone shot has real illumination gradients across the dish
# (one side brighter than the other) independent of colonies, and a global
# cutoff picks that up as false positives.
# Local contrast (difference-of-Gaussians) handles both problems: it's
# invariant to the photo's overall exposure AND to smooth illumination
# gradients across the dish, as long as the blur radius is well above
# typical colony size (otherwise large colonies fragment into edge-only
# rings instead of filled blobs).
LOCAL_CONTRAST_BLUR_SIGMA = 250
LOCAL_CONTRAST_CUTOFF = 20  # fixed; background/agar noise sits near 0 at this sigma

# ---- dish detection ---------------------------------------------------
DISH_ERODE_PX = 130  # the plastic dish rim itself catches a bright highlight band;
                      # needs a wide erosion at this resolution to fully exclude it

# ---- shape filters (same logic as the Chlorella pipeline) -----------------
MIN_AREA = 150
MAX_AREA = 90000
MIN_CIRCULARITY = 0.40  # local-contrast mask edges are rougher than a clean color mask, lowered vs. the Chlorella pipeline
MIN_SOLIDITY = 0.78  # real colony edges have small bumps/dents; lowered vs. the Chlorella pipeline
MAX_ASPECT = 1.8
MIN_DIAM = 18
CROP_MARGIN_PX = 25
MIN_CROP_SIZE = 90
CROP_CONFLICT_PAD = 4


def imread_unicode(path):
    data = np.fromfile(path, dtype=np.uint8)
    return cv2.imdecode(data, cv2.IMREAD_COLOR)


def imwrite_unicode(path, img):
    ext = os.path.splitext(path)[1]
    ok, buf = cv2.imencode(ext, img)
    if ok:
        buf.tofile(path)
    return ok


@dataclass
class Candidate:
    x: int
    y: int
    w: int
    h: int
    cx: float
    cy: float
    area: float
    circularity: float
    solidity: float
    accepted: bool = False
    reason: str = ""
    crop_box: tuple = None


def find_dish_roi(img):
    h, w = img.shape[:2]
    gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    gray_blur = cv2.medianBlur(gray, 9)
    circles = cv2.HoughCircles(
        gray_blur, cv2.HOUGH_GRADIENT, dp=1.5, minDist=w,
        param1=80, param2=60, minRadius=int(w * 0.25), maxRadius=int(w * 0.48),
    )
    roi = np.zeros((h, w), np.uint8)
    if circles is None:
        # fallback: assume a centered circle filling most of the frame
        cx, cy, r = w // 2, int(h * 0.55), int(w * 0.40)
    else:
        cx, cy, r = circles[0][0]
        cx, cy, r = int(cx), int(cy), int(r)
    cv2.circle(roi, (cx, cy), max(0, r - DISH_ERODE_PX), 255, thickness=cv2.FILLED)
    return roi, (cx, cy, r)


def compute_colony_mask_local_contrast(img, hsv, dish_roi):
    """Default mask: local contrast (difference-of-Gaussians) on brightness.
    Works when colonies are BRIGHTER than agar and exposure/illumination
    varies shot-to-shot (the original steel-bench dataset)."""
    v = hsv[:, :, 2].astype(np.float32)
    blur = cv2.GaussianBlur(v, (0, 0), sigmaX=LOCAL_CONTRAST_BLUR_SIGMA)
    local_contrast = v - blur
    mask = np.where(local_contrast > LOCAL_CONTRAST_CUTOFF, 255, 0).astype(np.uint8)
    mask = cv2.bitwise_and(mask, mask, mask=dish_roi)
    mask = cv2.morphologyEx(mask, cv2.MORPH_OPEN, np.ones((5, 5), np.uint8))
    mask = cv2.morphologyEx(mask, cv2.MORPH_CLOSE, np.ones((5, 5), np.uint8))
    return mask


# module-level hook so a batch-specific driver can swap in a different mask
# function (e.g. compute_colony_mask_hsv_dark below) without duplicating the
# shape-filter / crop-isolation / output logic in process_image.
MASK_FN = compute_colony_mask_local_contrast


def compute_colony_mask_hsv_dark(img, hsv, dish_roi):
    """Alternate mask for the white-background GFP-transformant batch: under
    that lighting, colonies are DARKER than agar (opposite of the steel-bench
    dataset) with a yellow-green hue (~20-65) vs. the agar's blue-cyan hue
    (~95-105). Local contrast found essentially no signal on this batch even
    at generous percentile cutoffs -- direct HSV thresholding, calibrated by
    sampling real colony vs. agar pixels, worked far better."""
    mask = cv2.inRange(hsv, (20, 0, 40), (65, 255, 155))
    mask = cv2.bitwise_and(mask, mask, mask=dish_roi)
    mask = cv2.morphologyEx(mask, cv2.MORPH_OPEN, np.ones((3, 3), np.uint8))
    mask = cv2.morphologyEx(mask, cv2.MORPH_CLOSE, np.ones((3, 3), np.uint8))
    return mask


def rects_overlap_with_pad(a, b, pad):
    ax1, ay1, aw, ah = a
    ax2, ay2 = ax1 + aw, ay1 + ah
    bx1, by1, bw, bh = b
    bx2, by2 = bx1 + bw, by1 + bh
    ax1 -= pad; ay1 -= pad; ax2 += pad; ay2 += pad
    return not (ax2 < bx1 or bx2 < ax1 or ay2 < by1 or by2 < ay1)


def process_image(path, out_crop_dir, out_qc_path, label_prefix):
    img = imread_unicode(path)
    if img is None:
        print(f"  !! could not read {path}")
        return []
    hsv = cv2.cvtColor(img, cv2.COLOR_BGR2HSV)
    h, w = img.shape[:2]

    dish_roi, (dcx, dcy, dr) = find_dish_roi(img)

    colony_mask = MASK_FN(img, hsv, dish_roi)

    contours, _ = cv2.findContours(colony_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

    candidates = []
    all_bboxes = []
    for c in contours:
        area = cv2.contourArea(c)
        if area < 15:
            continue
        x, y, bw, bh = cv2.boundingRect(c)
        if area >= 100:  # ignore dust-speck-scale contours for crop-conflict purposes; lowered from 300
                          # so small satellite colonies inside a crop window are caught, not just dust
            all_bboxes.append((x, y, bw, bh))
        perim = cv2.arcLength(c, True)
        circularity = 4 * np.pi * area / (perim * perim) if perim > 0 else 0
        hull = cv2.convexHull(c)
        hull_area = cv2.contourArea(hull)
        solidity = area / hull_area if hull_area > 0 else 0
        M = cv2.moments(c)
        if M["m00"] == 0:
            continue
        cx, cy = M["m10"] / M["m00"], M["m01"] / M["m00"]
        aspect = max(bw, bh) / max(1, min(bw, bh))

        cand = Candidate(x, y, bw, bh, cx, cy, area, circularity, solidity)

        if area < MIN_AREA:
            cand.reason = "too small"
        elif min(bw, bh) < MIN_DIAM:
            cand.reason = "too thin/small"
        elif area > MAX_AREA:
            cand.reason = "too large / merged"
        elif circularity < MIN_CIRCULARITY:
            cand.reason = f"not round (circ={circularity:.2f})"
        elif solidity < MIN_SOLIDITY:
            cand.reason = f"irregular shape (sol={solidity:.2f})"
        elif aspect > MAX_ASPECT:
            cand.reason = f"elongated (aspect={aspect:.2f})"
        elif x <= 1 or y <= 1 or x + bw >= w - 1 or y + bh >= h - 1:
            cand.reason = "touches image edge"
        else:
            cand.accepted = True
        candidates.append(cand)

    for cand in candidates:
        if not cand.accepted:
            continue
        side = max(cand.w, cand.h) + 2 * CROP_MARGIN_PX
        side = max(side, MIN_CROP_SIZE)
        x0 = int(cand.cx - side / 2); y0 = int(cand.cy - side / 2)
        x1 = int(cand.cx + side / 2); y1 = int(cand.cy + side / 2)
        x0, y0 = max(0, x0), max(0, y0)
        x1, y1 = min(w, x1), min(h, y1)
        cand.crop_box = (x0, y0, x1, y1)
        crop_rect = (x0, y0, x1 - x0, y1 - y0)

        my_box = (cand.x, cand.y, cand.w, cand.h)
        for other_box in all_bboxes:
            if other_box == my_box:
                continue
            if rects_overlap_with_pad(crop_rect, other_box, -CROP_CONFLICT_PAD):
                cand.accepted = False
                cand.reason = "other object inside crop frame"
                break

    overlay = img.copy()
    cv2.circle(overlay, (dcx, dcy), dr, (255, 128, 0), 6)
    os.makedirs(out_crop_dir, exist_ok=True)
    accepted_records = []
    colony_idx = 0

    for cand in candidates:
        color = (0, 200, 0) if cand.accepted else (0, 0, 220)
        cv2.rectangle(overlay, (cand.x, cand.y), (cand.x + cand.w, cand.y + cand.h), color, 4)
        if not cand.accepted and cand.reason:
            cv2.putText(overlay, cand.reason, (cand.x, max(0, cand.y - 8)),
                        cv2.FONT_HERSHEY_SIMPLEX, 1.0, color, 2, cv2.LINE_AA)

        if cand.accepted:
            x0, y0, x1, y1 = cand.crop_box
            crop = img[y0:y1, x0:x1]
            colony_idx += 1
            crop_name = f"{label_prefix}__colony{colony_idx:03d}.jpg"
            crop_path = os.path.join(out_crop_dir, crop_name)
            imwrite_unicode(crop_path, crop)
            accepted_records.append(dict(
                crop_path=os.path.relpath(crop_path, ROOT),
                source_image=os.path.relpath(path, ROOT),
                x=cand.x, y=cand.y, w=cand.w, h=cand.h,
                area=round(cand.area, 1), circularity=round(cand.circularity, 3),
                solidity=round(cand.solidity, 3),
            ))

    cv2.putText(overlay, f"kept {colony_idx} / {len(candidates)} candidates",
                (40, 100), cv2.FONT_HERSHEY_SIMPLEX, 2.5, (0, 140, 255), 4, cv2.LINE_AA)
    os.makedirs(os.path.dirname(out_qc_path), exist_ok=True)
    overlay_small = cv2.resize(overlay, (overlay.shape[1] // 3, overlay.shape[0] // 3))
    imwrite_unicode(out_qc_path, overlay_small)

    return accepted_records
