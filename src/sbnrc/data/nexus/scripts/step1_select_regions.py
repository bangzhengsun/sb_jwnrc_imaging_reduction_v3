#!/usr/bin/env python3
"""Automatically select sparse, commonly covered blank-sky HRMS regions."""

import argparse
import json
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from matplotlib import pyplot as plt
from scipy.ndimage import binary_dilation, gaussian_filter

import tomllib
from pathlib import Path
_CONFIG = tomllib.loads(Path(__file__).with_name('nexus.toml').read_text())
FILTERS = _CONFIG['filters']['sw'] + _CONFIG['filters']['lw']
DETECTION_FILTERS = _CONFIG['hrms']['detection_filters']


def product(collection, filt, suffix):
    matches = sorted(collection.glob(f"JWST_NRC_*_{filt}_60mas_{suffix}.fits"))
    if len(matches) != 1:
        raise RuntimeError(f"Expected one {filt} {suffix} file, found {len(matches)}")
    return matches[0]


def robust_scale(data, valid):
    sample = data[::8, ::8][valid[::8, ::8]]
    _, median, std = sigma_clipped_stats(sample, sigma=3, maxiters=8)
    return float(median), float(std)


def nonoverlapping_ranked(candidates, count, min_distance):
    selected = []
    for candidate in sorted(candidates, key=lambda item: item["score"]):
        y, x = candidate["center"]
        if all(np.hypot(y-q["center"][0], x-q["center"][1]) >= min_distance
               for q in selected):
            selected.append(candidate)
            if len(selected) == count:
                break
    return selected


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--collection", type=Path, default=Path("collection/60mas"))
    parser.add_argument("--box-size", type=int, default=301)
    parser.add_argument("--num-candidates", type=int, default=20)
    parser.add_argument("--num-select", type=int, default=10)
    parser.add_argument("--stride", type=int, default=75)
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    parser.add_argument(
        "--select-indices",
        help="Optional comma-separated indices from the ranked-20 plot; records a visual veto",
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    reference = fits.getdata(product(args.collection, "F444W", "sci"), memmap=True)
    shape = reference.shape
    common_valid = np.ones(shape, dtype=bool)
    for filt in FILTERS:
        wht = fits.getdata(product(args.collection, filt, "wht"), memmap=True)
        common_valid &= np.isfinite(wht) & (wht > 0)

    source_mask = np.zeros(shape, dtype=bool)
    normalized = {}
    for filt in DETECTION_FILTERS:
        sci = np.asarray(fits.getdata(product(args.collection, filt, "sci"), memmap=True))
        err = np.asarray(fits.getdata(product(args.collection, filt, "rms"), memmap=True))
        valid = common_valid & np.isfinite(sci) & np.isfinite(err) & (err > 0)
        median, fallback_std = robust_scale(sci, valid)
        finite_noise = err[valid & (err > 0)]
        typical_err = np.nanmedian(finite_noise)
        floor = min(typical_err, fallback_std) if fallback_std > 0 else typical_err
        denom = np.where(valid, np.maximum(err, floor), np.nan)
        sn = (sci - median) / denom
        sn_smooth = gaussian_filter(np.nan_to_num(sn, nan=0.0), sigma=1.2)
        source_mask |= binary_dilation((sn > 2.5) | (sn_smooth > 1.8), iterations=5)
        normalized[filt] = sn

    half = args.box_size // 2
    candidates = []
    for cy in range(half, shape[0]-half, args.stride):
        for cx in range(half, shape[1]-half, args.stride):
            ys = slice(cy-half, cy-half+args.box_size)
            xs = slice(cx-half, cx-half+args.box_size)
            valid = common_valid[ys, xs]
            valid_fraction = float(valid.mean())
            if valid_fraction < 0.995:
                continue
            occupancy = float(source_mask[ys, xs].mean())
            scatters, tails = [], []
            for filt in DETECTION_FILTERS:
                values = normalized[filt][ys, xs][valid]
                med = np.nanmedian(values)
                scatters.append(1.4826*np.nanmedian(np.abs(values-med)))
                tails.append(np.nanpercentile(values, 99))
            scatter = float(np.nanmedian(scatters))
            tail = float(np.nanmedian(tails))
            score = 12.0*occupancy + 0.35*scatter + 0.025*max(tail, 0)
            candidates.append({"center": [cy, cx], "corner": [cy-half, cx-half],
                               "score": score, "occupancy": occupancy,
                               "scatter": scatter, "tail99": tail,
                               "valid_fraction": valid_fraction})

    ranked20 = nonoverlapping_ranked(candidates, args.num_candidates, args.box_size*1.15)
    if args.select_indices:
        indices = [int(value) for value in args.select_indices.split(",")]
        if len(indices) != args.num_select or len(set(indices)) != args.num_select:
            raise ValueError(f"--select-indices must contain {args.num_select} unique values")
        if min(indices) < 0 or max(indices) >= len(ranked20):
            raise ValueError("--select-indices contains an out-of-range value")
        selected = [ranked20[index] for index in indices]
    else:
        indices = list(range(args.num_select))
        selected = ranked20[:args.num_select]
    if len(selected) != args.num_select:
        raise RuntimeError(f"Only found {len(selected)} suitable regions")

    def plot_regions(regions, filename, columns):
        rows = int(np.ceil(len(regions)/columns))
        fig, axes = plt.subplots(rows, columns, figsize=(4*columns, 3.5*rows))
        axes = np.atleast_1d(axes).ravel()
        finite = reference[np.isfinite(reference)]
        vmin, vmax = np.nanpercentile(finite, [2, 98])
        for index, (axis, region) in enumerate(zip(axes, regions)):
            y, x = region["corner"]
            axis.imshow(reference[y:y+args.box_size, x:x+args.box_size],
                        origin="lower", cmap="gray_r", vmin=vmin, vmax=vmax)
            axis.set_title(f"{index}: score={region['score']:.3f}\n"
                           f"src={region['occupancy']:.3f}, MAD={region['scatter']:.2f}")
            axis.set_xticks([]); axis.set_yticks([])
        for axis in axes[len(regions):]: axis.axis("off")
        fig.tight_layout(); fig.savefig(args.output_dir/filename, dpi=150); plt.close(fig)

    plot_regions(ranked20, "candidate_regions.png", 5)
    plot_regions(selected, "selected_regions_preview.png", 5)
    payload = {"method": "multi-band sparsity ranking (F200W+F444W), common 2-filter coverage, visual artifact veto",
               "box_size_60mas": args.box_size, "detection_filters": DETECTION_FILTERS,
               "selected_candidate_indices": indices,
               "corners": [item["corner"] for item in selected],
               "selected_metrics": selected, "candidate_metrics": ranked20}
    (args.output_dir/"selected_regions.json").write_text(json.dumps(payload, indent=2))
    print("Selected corners:", payload["corners"])
    print(f"Wrote diagnostics and selected_regions.json to {args.output_dir}")


if __name__ == "__main__":
    main()
