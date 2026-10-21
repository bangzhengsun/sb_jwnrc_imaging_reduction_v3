"""
find_mosaic_center.py
---------------------
Given a list of astrometrically-registered JWST NIRCam stage-2 calibrated
images (_cal.fits), find the optimal central coordinate (RA, Dec) to use
as the output frame center for stage-3 drizzling (sky_method in
calwebb_image3 / drizzlepac).

Strategy
--------
Three methods are implemented, in increasing sophistication:

  1. footprint_mean   – Mean RA/Dec of all image corner vertices (fast,
                        works even when images only partially overlap).
  2. overlap_centroid – Centroid of the TRUE sky intersection polygon
                        (requires shapely; best for heterogeneous pointings).
  3. weighted_mean    – Exposure-time-weighted mean of image centers
                        (simplest; useful when all frames cover similar area).

The script prints all three and recommends the most appropriate one.
It also writes a DS9 region file so you can visually inspect the result.

Usage
-----
  python find_mosaic_center.py file1_cal.fits file2_cal.fits ...
  python find_mosaic_center.py --filelist my_asn.txt   # one path per line
  python find_mosaic_center.py --glob "path/to/*_cal.fits"

Dependencies
------------
  astropy, numpy, (optional) shapely
"""

import argparse
import glob
import sys
import warnings
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning

warnings.filterwarnings("ignore", category=FITSFixedWarning)

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def load_wcs_and_shape(fpath):
    """Return (wcs, naxis2, naxis1, exptime) from a _cal.fits file."""
    with fits.open(fpath) as hdul:
        # Stage-2 products: science in extension 1 ('SCI') or primary
        hdr = None
        for ext in hdul:
            if ext.data is not None and ext.data.ndim == 2:
                hdr = ext.header
                break
        if hdr is None:
            raise ValueError(f"No 2-D science extension found in {fpath}")

        wcs = WCS(hdr, naxis=2)
        ny, nx = hdr["NAXIS2"], hdr["NAXIS1"]
        exptime = hdr.get("XPOSURE", hdr.get("EFFEXPTM", hdr.get("EXPTIME", 1.0)))
    return wcs, ny, nx, float(exptime)


def image_corners_sky(wcs, ny, nx):
    """
    Return sky coordinates (RA, Dec in degrees) of the 4 corners
    and the centre of a 2-D image described by *wcs*.
    Pixel convention: 0-indexed, corner = (0,0) etc.
    """
    corners_pix = np.array([
        [0,    0   ],
        [nx-1, 0   ],
        [nx-1, ny-1],
        [0,    ny-1],
    ], dtype=float)
    center_pix = np.array([[nx / 2.0, ny / 2.0]])

    sky_corners = wcs.all_pix2world(corners_pix, 0)   # shape (4,2)
    sky_center  = wcs.all_pix2world(center_pix,  0)   # shape (1,2)
    return sky_corners, sky_center[0]


def mean_ra_dec(ra_arr, dec_arr):
    """
    Compute the mean sky position robustly (handles RA wrap at 0/360).
    Uses unit-vector averaging.
    """
    ra_rad  = np.deg2rad(ra_arr)
    dec_rad = np.deg2rad(dec_arr)
    x = np.cos(dec_rad) * np.cos(ra_rad)
    y = np.cos(dec_rad) * np.sin(ra_rad)
    z = np.sin(dec_rad)
    xm, ym, zm = x.mean(), y.mean(), z.mean()
    ra_mean  = np.rad2deg(np.arctan2(ym, xm)) % 360.0
    dec_mean = np.rad2deg(np.arctan2(zm, np.hypot(xm, ym)))
    return ra_mean, dec_mean


# ---------------------------------------------------------------------------
# Method 1 – footprint mean
# ---------------------------------------------------------------------------

def method_footprint_mean(images):
    """Average RA/Dec of all corner vertices across all images."""
    all_ra, all_dec = [], []
    for fpath, wcs, ny, nx, _ in images:
        corners, _ = image_corners_sky(wcs, ny, nx)
        all_ra.extend(corners[:, 0])
        all_dec.extend(corners[:, 1])
    return mean_ra_dec(np.array(all_ra), np.array(all_dec))


# ---------------------------------------------------------------------------
# Method 2 – overlap centroid (requires shapely)
# ---------------------------------------------------------------------------

def method_overlap_centroid(images):
    """
    Compute the centroid of the polygon that is the intersection of all
    image footprints projected onto a tangent plane centred on the rough
    footprint mean.
    """
    try:
        from shapely.geometry import Polygon
        from shapely.ops import unary_union
    except ImportError:
        return None, None, "shapely not installed; skipping overlap_centroid"

    # Rough centre for gnomonic projection
    all_ra, all_dec = [], []
    for _, wcs, ny, nx, _ in images:
        c, _ = image_corners_sky(wcs, ny, nx)
        all_ra.extend(c[:, 0]); all_dec.extend(c[:, 1])
    ra0, dec0 = mean_ra_dec(np.array(all_ra), np.array(all_dec))

    def gnomonic_xy(ra, dec, ra0=ra0, dec0=dec0):
        """Project (ra, dec) onto tangent plane at (ra0, dec0). Returns (x,y) deg."""
        ra_rad  = np.deg2rad(ra);  dec_rad  = np.deg2rad(dec)
        ra0_rad = np.deg2rad(ra0); dec0_rad = np.deg2rad(dec0)
        cosc = (np.sin(dec0_rad)*np.sin(dec_rad)
                + np.cos(dec0_rad)*np.cos(dec_rad)*np.cos(ra_rad - ra0_rad))
        x = np.cos(dec_rad)*np.sin(ra_rad - ra0_rad) / cosc
        y = (np.cos(dec0_rad)*np.sin(dec_rad)
             - np.sin(dec0_rad)*np.cos(dec_rad)*np.cos(ra_rad - ra0_rad)) / cosc
        return np.rad2deg(x), np.rad2deg(y)

    polys = []
    for _, wcs, ny, nx, _ in images:
        corners, _ = image_corners_sky(wcs, ny, nx)
        xs, ys = gnomonic_xy(corners[:, 0], corners[:, 1])
        polys.append(Polygon(zip(xs, ys)))

    intersection = polys[0]
    for p in polys[1:]:
        intersection = intersection.intersection(p)
        if intersection.is_empty:
            return None, None, "Images have no common overlap region"

    cx_deg, cy_deg = intersection.centroid.x, intersection.centroid.y

    # Back-project from tangent plane to sky
    ra0_rad  = np.deg2rad(ra0);  dec0_rad = np.deg2rad(dec0)
    x_rad    = np.deg2rad(cx_deg); y_rad   = np.deg2rad(cy_deg)
    rho      = np.hypot(x_rad, y_rad)
    c        = np.arctan(rho)
    dec_out  = np.rad2deg(np.arcsin(
        np.cos(c)*np.sin(dec0_rad) + y_rad*np.sin(c)*np.cos(dec0_rad)/rho
    )) if rho > 0 else dec0
    ra_out   = ra0 + np.rad2deg(np.arctan2(
        x_rad*np.sin(c),
        rho*np.cos(dec0_rad)*np.cos(c) - y_rad*np.sin(dec0_rad)*np.sin(c)
    )) if rho > 0 else ra0
    ra_out %= 360.0
    return ra_out, dec_out, "OK"


# ---------------------------------------------------------------------------
# Method 3 – exposure-time-weighted mean of image centres
# ---------------------------------------------------------------------------

def method_weighted_mean(images):
    """Exposure-time-weighted mean of individual image centres."""
    ra_list, dec_list, wt_list = [], [], []
    for _, wcs, ny, nx, exptime in images:
        _, cen = image_corners_sky(wcs, ny, nx)
        ra_list.append(cen[0]); dec_list.append(cen[1]); wt_list.append(exptime)
    wt = np.array(wt_list)
    ra_arr  = np.array(ra_list)
    dec_arr = np.array(dec_list)
    # Weighted unit-vector mean
    ra_rad  = np.deg2rad(ra_arr);  dec_rad = np.deg2rad(dec_arr)
    x = np.average(np.cos(dec_rad)*np.cos(ra_rad), weights=wt)
    y = np.average(np.cos(dec_rad)*np.sin(ra_rad), weights=wt)
    z = np.average(np.sin(dec_rad),                weights=wt)
    ra_m  = np.rad2deg(np.arctan2(y, x)) % 360.0
    dec_m = np.rad2deg(np.arctan2(z, np.hypot(x, y)))
    return ra_m, dec_m


# ---------------------------------------------------------------------------
# DS9 region file
# ---------------------------------------------------------------------------

def write_ds9_regions(images, results, outfile="mosaic_center.reg"):
    lines = [
        "# Region file format: DS9 version 4.1",
        "global color=green dashlist=8 3 width=1 font='helvetica 10 normal roman'",
        "fk5",
    ]
    # Image footprints
    for fpath, wcs, ny, nx, _ in images:
        corners, cen = image_corners_sky(wcs, ny, nx)
        verts = " ".join(f"{r:.6f},{d:.6f}" for r, d in corners)
        lines.append(f"polygon({verts}) # color=cyan text={{{Path(fpath).name}}}")
        lines.append(f"point({cen[0]:.6f},{cen[1]:.6f}) # point=x color=cyan")

    # Candidate centres
    colors = {"footprint_mean": "green", "overlap_centroid": "magenta", "weighted_mean": "yellow"}
    for method, (ra, dec) in results.items():
        if ra is None:
            continue
        c = colors.get(method, "white")
        lines.append(f"point({ra:.6f},{dec:.6f}) # point=circle 15 color={c} text={{{method}}}")

    with open(outfile, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"\nDS9 region file written → {outfile}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("files", nargs="*", help="_cal.fits file paths")
    p.add_argument("--filelist", help="Text file with one _cal.fits path per line")
    p.add_argument("--glob",    help="Glob pattern, e.g. 'data/*_cal.fits'")
    p.add_argument("--outregion", default="mosaic_center.reg",
                   help="Output DS9 region file (default: mosaic_center.reg)")
    return p.parse_args()


def collect_files(args):
    paths = list(args.files)
    if args.filelist:
        with open(args.filelist) as f:
            paths += [l.strip() for l in f if l.strip() and not l.startswith("#")]
    if args.glob:
        paths += glob.glob(args.glob)
    paths = sorted(set(paths))
    if not paths:
        sys.exit("ERROR: No input files provided. Use positional args, --filelist, or --glob.")
    return paths


def main():
    args = parse_args()
    filepaths = collect_files(args)

    print(f"\nLoading WCS from {len(filepaths)} image(s)...")
    images = []  # list of (fpath, wcs, ny, nx, exptime)
    for fp in filepaths:
        try:
            wcs, ny, nx, exptime = load_wcs_and_shape(fp)
            images.append((fp, wcs, ny, nx, exptime))
            print(f"  OK  {fp}  ({nx}×{ny} pix, exptime={exptime:.1f}s)")
        except Exception as e:
            print(f"  SKIP {fp}: {e}")

    if not images:
        sys.exit("ERROR: No valid images could be loaded.")

    # -----------------------------------------------------------------------
    print("\n" + "="*60)
    print("  Computing centre coordinates")
    print("="*60)

    results = {}

    # Method 1
    ra1, dec1 = method_footprint_mean(images)
    results["footprint_mean"] = (ra1, dec1)
    print(f"\n[1] Footprint mean (all corners averaged)")
    print(f"    RA  = {ra1:.6f} deg  ({_deg2hms(ra1)})")
    print(f"    Dec = {dec1:+.6f} deg  ({_deg2dms(dec1)})")

    # Method 2
    ra2, dec2, msg2 = method_overlap_centroid(images)
    results["overlap_centroid"] = (ra2, dec2)
    print(f"\n[2] Overlap centroid (intersection polygon centroid)")
    if ra2 is not None:
        print(f"    RA  = {ra2:.6f} deg  ({_deg2hms(ra2)})")
        print(f"    Dec = {dec2:+.6f} deg  ({_deg2dms(dec2)})")
    else:
        print(f"    Not available: {msg2}")

    # Method 3
    ra3, dec3 = method_weighted_mean(images)
    results["weighted_mean"] = (ra3, dec3)
    print(f"\n[3] Exposure-time-weighted mean of image centres")
    print(f"    RA  = {ra3:.6f} deg  ({_deg2hms(ra3)})")
    print(f"    Dec = {dec3:+.6f} deg  ({_deg2dms(dec3)})")

    # -----------------------------------------------------------------------
    print("\n" + "="*60)
    print("  RECOMMENDATION")
    print("="*60)
    if ra2 is not None:
        best_ra, best_dec = ra2, dec2
        best_method = "overlap_centroid"
        reason = "Centroid of the true sky overlap polygon — minimises wasted output pixels."
    else:
        best_ra, best_dec = ra1, dec1
        best_method = "footprint_mean"
        reason = "Overlap centroid unavailable (install shapely). Using footprint mean."

    print(f"\n  Method  : {best_method}")
    print(f"  RA      : {best_ra:.6f} deg  ({_deg2hms(best_ra)})")
    print(f"  Dec     : {best_dec:+.6f} deg  ({_deg2dms(best_dec)})")
    print(f"  Reason  : {reason}")

    print(f"\n  Use in calwebb_image3 association JSON or pipeline call:")
    print(f"    ra_center  = {best_ra:.6f}")
    print(f"    dec_center = {best_dec:.6f}")

    # DS9 region
    write_ds9_regions(images, results, outfile=args.outregion)


# ---------------------------------------------------------------------------
# Formatting helpers
# ---------------------------------------------------------------------------

def _deg2hms(deg):
    h_total = deg / 15.0
    h  = int(h_total)
    m  = int((h_total - h) * 60)
    s  = ((h_total - h) * 60 - m) * 60
    return f"{h:02d}h{m:02d}m{s:06.3f}s"


def _deg2dms(deg):
    sign = "+" if deg >= 0 else "-"
    deg  = abs(deg)
    d    = int(deg)
    m    = int((deg - d) * 60)
    s    = ((deg - d) * 60 - m) * 60
    return f"{sign}{d:02d}d{m:02d}m{s:05.2f}s"


if __name__ == "__main__":
    main()