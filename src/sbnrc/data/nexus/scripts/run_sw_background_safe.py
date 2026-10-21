#!/usr/bin/env python3
"""Idempotent 2-D background subtraction for SW calibrated images."""

from pathlib import Path
import os
import shutil
import sys

import numpy as np
from astropy import stats
from astropy.io import fits
from photutils.background import Background2D, MedianBackground
from photutils.segmentation import detect_sources, detect_threshold
from photutils.utils import circular_footprint

FILTER = sys.argv[1].upper()
BASE = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(os.environ.get("BASE", "/mnt/Data_Analysis_6/Data_Reduction/JWST_NRC_SMACS0723_v2026Aug"))
CAL = BASE / FILTER / "stage2_cal"
BACKUP = CAL / "Orig_cal"


def background(data):
    image = data.copy()
    image[image == 0] = np.nan
    clip = stats.SigmaClip(sigma=3.0, maxiters=10)
    threshold = detect_threshold(image, nsigma=3.0, sigma_clip=clip)
    segment = detect_sources(image, threshold, npixels=10)
    mask = None if segment is None else segment.make_source_mask(footprint=circular_footprint(radius=10))
    return Background2D(
        image, (128, 128), filter_size=(5, 5), mask=mask,
        sigma_clip=clip, bkg_estimator=MedianBackground(), exclude_percentile=80.0,
    ).background


def main():
    BACKUP.mkdir(parents=True, exist_ok=True)
    files = sorted(CAL.glob("*_cal.fits"))
    if not files:
        raise RuntimeError(f"No calibrated files under {CAL}")
    for number, path in enumerate(files, 1):
        original = BACKUP / path.name
        if original.exists():
            print(f"[{number:03d}/{len(files)}] ALREADY SUBTRACTED: {path.name}", flush=True)
            continue
        shutil.copy2(path, original)
        with fits.open(original, memmap=False) as hdul:
            data = hdul["SCI"].data
            bkg = background(data)
        with fits.open(path, mode="update", memmap=False) as hdul:
            hdul["SCI"].data[...] = data - bkg
            hdul[0].header["HISTORY"] = "2-D Background2D model subtracted"
            hdul.flush()
        print(f"[{number:03d}/{len(files)}] SUBTRACTED: {path.name}", flush=True)


if __name__ == "__main__":
    main()
