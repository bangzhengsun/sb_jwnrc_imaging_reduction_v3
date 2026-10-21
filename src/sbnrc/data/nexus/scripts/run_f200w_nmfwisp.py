#!/usr/bin/env python3
"""Fit/subtract NMFwisp from F200W flat-field products and make QA panels."""

from pathlib import Path
import csv
import os
import shutil
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.stats import SigmaClip
from photutils.segmentation import detect_sources, detect_threshold
from photutils.utils import circular_footprint
from nmfwisp import fit_wisp
import nmfwisp

FILTER = sys.argv[1].upper() if len(sys.argv) > 1 else "F200W"
BASE = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(os.environ.get("BASE", "/mnt/Data_Analysis_6/Data_Reduction/JWST_NRC_SMACS0723_v2026Aug"))
WORK = BASE / FILTER / "stage2_medstep_wisprm"
ORIGINALS = WORK / "original_flatfieldstep"
MODELS = WORK / "nmfwisp_models"
DIAGNOSTICS = BASE / FILTER / "nmfwisp_diagnostics"
ALIASES = BASE / FILTER / "nmfwisp_template_aliases"
AFFECTED = {"NRCA3", "NRCA4", "NRCB3", "NRCB4"}


def make_template_aliases() -> None:
    package_templates = Path(nmfwisp.__file__).parent / "templates"
    for source in package_templates.glob("**/*_wisp.fits.gz"):
        rel = source.relative_to(package_templates)
        target_dir = ALIASES / rel.parent
        target_dir.mkdir(parents=True, exist_ok=True)
        # NMFwisp uppercases filter_name before lookup, but v1.1.4 files are lowercase.
        bits = source.name.split("_")
        bits[-2] = bits[-2].upper()
        target = target_dir / "_".join(bits)
        if not target.exists():
            target.symlink_to(source)


def template_available(detector: str) -> bool:
    det = detector.lower()
    stem = f"{det}_{FILTER}_wisp"
    return any(
        (ALIASES / directory / f"{stem}{extension}").exists()
        for directory in (det, f"{det}_org")
        for extension in (".fits", ".fits.gz")
    )


def source_and_bad_mask(data: np.ndarray, dq: np.ndarray, err: np.ndarray) -> np.ndarray:
    invalid = (~np.isfinite(data)) | (~np.isfinite(err)) | (err <= 0) | (dq != 0)
    clean = data.copy()
    clean[invalid] = np.nan
    threshold = detect_threshold(clean, nsigma=5.5, sigma_clip=SigmaClip(sigma=3, maxiters=10))
    segment = detect_sources(clean, threshold, npixels=10)
    if segment is None:
        sources = np.zeros(data.shape, dtype=bool)
    else:
        sources = segment.make_source_mask(footprint=circular_footprint(radius=10))
    return invalid | sources


def save_diagnostic(path: Path, before: np.ndarray, after: np.ndarray, detector: str) -> None:
    valid = np.isfinite(before)
    lo, hi = np.nanpercentile(before[valid], [5, 95])
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.8), constrained_layout=True)
    for ax, image, title in zip(axes, [before, after], ["Data", "Data - NMF WISP"]):
        shown = np.nan_to_num(image, nan=lo, posinf=hi, neginf=lo)
        artist = ax.imshow(shown, origin="lower", cmap="gray", vmin=lo, vmax=hi)
        ax.set_title(title)
        ax.set_xlabel("x [pixel]")
        ax.set_ylabel("y [pixel]")
        fig.colorbar(artist, ax=ax, fraction=0.046, pad=0.04, label="DN/s")
    fig.suptitle(f"{path.name}  ({detector}, {FILTER})", fontsize=10)
    fig.savefig(DIAGNOSTICS / path.with_suffix(".png").name, dpi=150)
    plt.close(fig)


def main() -> None:
    make_template_aliases()
    for directory in (ORIGINALS, MODELS, DIAGNOSTICS):
        directory.mkdir(parents=True, exist_ok=True)
    files = sorted(WORK.glob("*_flatfieldstep.fits"))
    if not files:
        raise RuntimeError(f"No flatfieldstep inputs found under {WORK}")

    rows = []
    for index, path in enumerate(files, 1):
        backup = ORIGINALS / path.name
        if not backup.exists():
            shutil.copy2(path, backup)
        # Always derive the test product from the preserved original, making reruns idempotent.
        shutil.copy2(backup, path)
        with fits.open(backup, memmap=False) as hdul:
            before = hdul["SCI"].data.astype(np.float64)
            err = hdul["ERR"].data.astype(np.float64)
            dq = hdul["DQ"].data
            detector = hdul[0].header["DETECTOR"].upper()
            mask = source_and_bad_mask(before, dq, err)

        if detector in AFFECTED and template_available(detector):
            wisp, wisp_err = fit_wisp(
                before, err, mask,
                wisp_path=str(ALIASES), detector_name=detector.lower(),
                filter_name=FILTER, correct_1f=True,
            )
        else:
            wisp = np.zeros_like(before)
            wisp_err = np.zeros_like(before)
        after = before - wisp

        with fits.open(path, mode="update", memmap=False) as hdul:
            hdul["SCI"].data[...] = after.astype(hdul["SCI"].data.dtype)
            hdul[0].header["HISTORY"] = f"NMFwisp 1.1.4 {FILTER} model subtracted; correct_1f=True"
            hdul.flush()
        fits.HDUList([
            fits.PrimaryHDU(header=fits.Header({"DETECTOR": detector, "FILTER": FILTER})),
            fits.ImageHDU(wisp.astype(np.float32), name="WISP"),
            fits.ImageHDU(wisp_err.astype(np.float32), name="ERR"),
            fits.ImageHDU(mask.astype(np.uint8), name="MASK"),
        ]).writeto(MODELS / path.name.replace("_flatfieldstep.fits", "_nmfwisp_model.fits"), overwrite=True)
        save_diagnostic(path, before, after, detector)
        rows.append({
            "filename": path.name, "detector": detector,
            "masked_fraction": float(mask.mean()),
            "wisp_median": float(np.nanmedian(wisp)),
            "wisp_p95": float(np.nanpercentile(wisp, 95)),
            "wisp_max": float(np.nanmax(wisp)),
        })
        print(f"[{index:02d}/{len(files)}] {detector}: {path.name}", flush=True)

    with (DIAGNOSTICS / "nmfwisp_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print(f"Completed {len(files)} files; diagnostics: {DIAGNOSTICS}")


if __name__ == "__main__":
    main()
