"""Feature-based TEM segmentation, calibrated particle measurements, and plots."""
from __future__ import annotations
import argparse
import json
import heapq
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from scipy.spatial import ConvexHull
from PIL import Image
from ppt_size_distribution import load_data, plot_figures


def otsu_threshold(values):
    """Between-class variance maximization on a 256-bin intensity histogram."""
    counts, edges = np.histogram(values, bins=256)
    centers = (edges[:-1]+edges[1:])/2
    weight = np.cumsum(counts)
    moment = np.cumsum(counts*centers)
    denom = weight*(weight[-1]-weight)
    score = np.full(256, -np.inf)
    ok = denom > 0
    score[ok] = (moment[-1]*weight[ok]-moment[ok]*weight[-1])**2 / denom[ok]
    return float(centers[np.argmax(score)])


def watershed(distance, markers, mask):
    """Flood negative distance using a priority queue and four-connected pixels.

    Queue order breaks equal-height ties by arrival time. Original marker pixels
    remain fixed. This separates distance peaks without assuming circular shapes.
    """
    labels = markers.copy()
    queued = labels > 0
    heap = []
    counter = 0
    for y, x in np.argwhere(queued):
        heapq.heappush(heap, (-distance[y, x], counter, int(y), int(x), int(labels[y, x])))
        counter += 1
    height, width = mask.shape
    while heap:
        level, _, y, x, label = heapq.heappop(heap)
        labels[y, x] = label
        for ny, nx in ((y-1,x),(y+1,x),(y,x-1),(y,x+1)):
            if 0 <= ny < height and 0 <= nx < width and mask[ny,nx] and not queued[ny,nx]:
                queued[ny,nx] = True
                heapq.heappush(heap, (max(level, -distance[ny,nx]), counter, ny, nx, label))
                counter += 1
    return labels


def detect(image, nm_per_pixel, *, polarity="bright", sigma_px=1.0,
           min_area_nm2=100.0, seed_distance_px=2.0, min_solidity=0.0,
           max_aspect_ratio=None, exclude_border=True, roi=None):
    """Return labelled candidates and audited morphology measurements.

    Phase identity is supplied by the user through TEM imaging context, not
    inferred chemically. ROI is a boolean valid-specimen mask.
    """
    values = (nm_per_pixel, sigma_px, min_area_nm2, seed_distance_px, min_solidity)
    if not np.isfinite(values).all() or nm_per_pixel <= 0 or sigma_px < 0 or min_area_nm2 <= 0 or seed_distance_px <= 0:
        raise ValueError("Calibration, area, and seed distance must be positive; smoothing must be nonnegative")
    if not 0 <= min_solidity <= 1 or polarity not in ("bright", "dark"):
        raise ValueError("Invalid solidity or polarity")
    if max_aspect_ratio is not None and (not np.isfinite(max_aspect_ratio) or max_aspect_ratio < 1):
        raise ValueError("Maximum aspect ratio must be finite and >= 1")
    image = np.asarray(image)
    if image.ndim == 3 and image.shape[-1] in (3, 4):
        image = np.dot(image[..., :3], [0.2126, 0.7152, 0.0722])
    if image.ndim != 2:
        raise ValueError("Use a single grayscale or RGB micrograph, not an image stack")
    gray = image.astype(float)
    if not np.isfinite(gray).all():
        raise ValueError("Image contains nonfinite pixels")
    valid = np.ones(gray.shape, bool) if roi is None else np.asarray(roi, dtype=bool)
    if valid.shape != gray.shape or not valid.any():
        raise ValueError("ROI mask must match the image and contain valid specimen pixels")
    smoothed = ndi.gaussian_filter(gray, sigma=sigma_px) if sigma_px else gray
    if np.ptp(smoothed[valid]) == 0:
        raise ValueError("No image contrast inside the ROI")
    threshold = otsu_threshold(smoothed[valid])
    foreground = (smoothed > threshold if polarity == "bright" else smoothed < threshold) & valid
    distance = ndi.distance_transform_edt(np.pad(foreground, 1))[1:-1, 1:-1]
    # Local distance maxima seed a marker-controlled watershed. The user controls
    # peak suppression distance; it is not an automatic estimate of particle size.
    radius = max(1, int(np.ceil(seed_distance_px)))
    peaks = (distance == ndi.maximum_filter(distance, size=2*radius+1)) & foreground
    markers, _ = ndi.label(peaks)
    components, count = ndi.label(foreground)
    for component in range(1, count + 1):
        region = components == component
        if not np.any(markers[region]):
            index = np.argmax(np.where(region, distance, -1))
            markers.flat[index] = markers.max() + 1
    labels = watershed(distance, markers, foreground)
    boundary = np.zeros(gray.shape, bool)
    boundary[[0, -1], :] = True
    boundary[:, [0, -1]] = True
    boundary |= ndi.binary_dilation(~valid)
    border_ids = set(np.unique(labels[boundary]))
    rows = []
    for particle_id, box in enumerate(ndi.find_objects(labels), 1):
        if box is None:
            continue
        region = labels[box] == particle_id
        yy, xx = np.nonzero(region)
        area_px = len(xx)
        area = area_px * nm_per_pixel**2
        # Moment axes include the variance 1/12 of each square pixel.
        coords = np.column_stack((xx, yy)).astype(float)
        centered = coords - coords.mean(axis=0)
        eig = np.linalg.eigvalsh(centered.T @ centered / area_px + np.eye(2)/12)
        minor, major = 4*np.sqrt(eig) * nm_per_pixel
        aspect = major / minor
        padded = np.pad(region.astype(int), 1)
        perimeter = (np.abs(np.diff(padded, axis=0)).sum() +
                     np.abs(np.diff(padded, axis=1)).sum()) * nm_per_pixel
        corners = np.concatenate([coords + offset for offset in
                                  [(-.5,-.5),(-.5,.5),(.5,-.5),(.5,.5)]])
        solidity = min(1.0, area_px / ConvexHull(corners).volume)
        reasons = []
        if area < min_area_nm2:
            reasons.append("below_min_area")
        if exclude_border and particle_id in border_ids:
            reasons.append("touches_image_or_roi_border")
        if solidity < min_solidity:
            reasons.append("below_min_solidity")
        if max_aspect_ratio is not None and aspect > max_aspect_ratio:
            reasons.append("above_max_aspect_ratio")
        rows.append(dict(particle_id=particle_id, area_nm2=area,
                         diameter_nm=2*np.sqrt(area/np.pi), perimeter_nm=perimeter,
                         major_axis_nm=major, minor_axis_nm=minor, aspect_ratio=aspect,
                         solidity=solidity, circularity=4*np.pi*area/perimeter**2,
                         centroid_x_px=xx.mean()+box[1].start,
                         centroid_y_px=yy.mean()+box[0].start,
                         accepted=not reasons, rejection_reason=";".join(reasons)))
    if not rows:
        raise ValueError("No candidate particles detected; inspect contrast, polarity, and ROI")
    return gray, labels, pd.DataFrame(rows), float(threshold)


def run_manifest(manifest: Path, output: Path, make_plots=True):
    config = json.loads(manifest.read_text())
    images = config.get("images", [])
    if not images:
        raise ValueError("Manifest needs a nonempty images list")
    output.mkdir(parents=True, exist_ok=True)
    all_rows, runs = [], []
    for index, item in enumerate(images, 1):
        path = (manifest.parent / item["path"]).resolve()
        time = item["aging_time_min"]
        if time not in (15, 60, 165):
            raise ValueError("Aging time must be 15, 60, or 165 minutes")
        image = np.asarray(Image.open(path))
        roi = None
        if item.get("roi_mask"):
            roi = np.asarray(Image.open(manifest.parent / item["roi_mask"]))
            if roi.ndim != 2:
                raise ValueError("ROI mask must be grayscale: nonzero = valid specimen")
            roi = roi > 0
        settings = {**config.get("segmentation", {}), **item.get("segmentation", {})}
        gray, labels, table, threshold = detect(image, item["nm_per_pixel"], roi=roi, **settings)
        table.insert(0, "image_id", f"image_{index:03d}")
        table.insert(1, "source_image", item["path"])
        table.insert(2, "aging_time_min", time)
        all_rows.append(table)
        stem = f"image_{index:03d}"
        Image.fromarray(labels.astype(np.int32)).save(output / f"{stem}_labels.tif")
        accepted_ids = table.loc[table.accepted, "particle_id"].to_numpy()
        accepted_mask = np.isin(labels, accepted_ids) & (labels > 0)
        Image.fromarray((accepted_mask*255).astype(np.uint8)).save(output / f"{stem}_mask.png")
        base = np.repeat(gray[..., None], 3, axis=2)
        lo, hi = base.min(), base.max()
        base = (base-lo)/(hi-lo) if hi > lo else base
        edges = (ndi.maximum_filter(labels, size=3) != ndi.minimum_filter(labels, size=3)) & (labels > 0)
        base[edges & accepted_mask] = (0, 1, 0.6)
        base[edges & ~accepted_mask] = (1, 0.25, 0.25)
        Image.fromarray(np.round(base*255).astype(np.uint8)).save(output / f"{stem}_overlay.png")
        runs.append(dict(image_id=stem, source_image=item["path"], aging_time_min=time,
                         nm_per_pixel=item["nm_per_pixel"], roi_mask=item.get("roi_mask"),
                         settings=settings, otsu_threshold=threshold,
                         candidates=len(table), accepted=int(table.accepted.sum())))
    audit = pd.concat(all_rows, ignore_index=True)
    audit.to_csv(output / "all_candidates.csv", index=False)
    accepted = audit[audit.accepted]
    accepted.to_csv(output / "particle_measurements.csv", index=False)
    accepted[["aging_time_min", "diameter_nm"]].to_csv(output / "diameters.csv", index=False)
    (output / "run_metadata.json").write_text(json.dumps({"configuration": config, "runs": runs}, indent=2))
    if make_plots:
        # Validate minimum sample counts and variability before fitting curves.
        frame = load_data(output / "diameters.csv")
        plot_figures(frame, output / "figures", label="Automatic segmentation · unreviewed")
    return audit


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path("results"))
    parser.add_argument("--segment-only", action="store_true", help="Inspect one image without requiring three conditions")
    args = parser.parse_args()
    try:
        audit = run_manifest(args.manifest, args.output, not args.segment_only)
        print(f"Detected {len(audit)} candidates; accepted {audit.accepted.sum()}. Review overlays in {args.output}")
    except (ValueError, OSError, KeyError, TypeError) as exc:
        parser.exit(2, f"Error: {exc}\n")


if __name__ == "__main__":
    main()
