"""
step2_calc_hrms.py
------------------
Compute hrms maps for all filters and both resolutions (30mas / 60mas) using
the sky regions selected by step1_select_regions.py.

Usage
-----
    python step2_calc_hrms.py

Prerequisites
-------------
    selected_regions.json   (produced by step1_select_regions.py)

Output
------
    ./collection/{res}/JWST_NRC_*_{FILT}_{res}_hrms.fits
    for every combination of filter × resolution.
"""

import json
import glob
import sys
import warnings
import numpy as np
from numpy.fft import fftn, ifftn, ifftshift
from astropy.io import fits
from astropy.stats import biweight_location, biweight_midvariance, sigma_clipped_stats
from scipy.ndimage import median_filter, binary_dilation
from scipy import ndimage as ndimg
from copy import deepcopy

warnings.filterwarnings("ignore")


# ─────────────────────────────────────────────────────────────────────────────
# All helper functions (verbatim from the notebook)
# ─────────────────────────────────────────────────────────────────────────────

def measure_autocorrelation(image_data, aper_radius=5, bg_annulus_radii=(5, 7)):
    fft_data = fftn(image_data)
    autocorr_image = ifftn(fft_data * np.conjugate(fft_data)).real
    autocorr_image = ifftshift(autocorr_image)
    im_shape = autocorr_image.shape
    center = np.array(im_shape) / 2 - 0.5

    x, y = np.meshgrid(range(im_shape[1]), range(im_shape[0]))
    sq_dist = (x - center[1]) ** 2 + (y - center[0]) ** 2

    sky_mask = sq_dist > min(bg_annulus_radii) ** 2
    sky_mask &= sq_dist < max(bg_annulus_radii) ** 2
    sky_fix_mask = ~sky_mask
    sky_fix_mask &= sq_dist > (min(bg_annulus_radii) - 1) ** 2
    sky_fix_mask &= sq_dist < (max(bg_annulus_radii) + 1) ** 2

    sky_fix_area = (np.pi * (max(bg_annulus_radii) ** 2 -
                             min(bg_annulus_radii) ** 2) - np.sum(sky_mask))
    fix_pixels_weight = sky_fix_area / np.sum(sky_fix_mask)
    sky_wts = 1.0 * sky_mask + fix_pixels_weight * sky_fix_mask

    flux_mask = sq_dist < aper_radius ** 2
    flux_fix_mask = ~flux_mask & (sq_dist < (aper_radius + 1) ** 2)
    flux_fix_area = np.pi * aper_radius ** 2 - np.sum(flux_mask)
    fix_pixels_weight = flux_fix_area / np.sum(flux_fix_mask)
    flux_wts = 1.0 * flux_mask + fix_pixels_weight * flux_fix_mask

    peak_val = np.max(autocorr_image[flux_mask])
    bg_val = np.average(autocorr_image, weights=sky_wts)
    total_corr = np.sum((autocorr_image - bg_val) * flux_wts)
    corr_rms = np.sqrt((peak_val - bg_val) / autocorr_image.size)
    corr_fac = np.sqrt(total_corr / (peak_val - bg_val))
    autocorr_image -= bg_val

    return corr_rms, corr_fac


def _binned_median_filter(image_data, bin_size=4, filt_size=25):
    img_shape = image_data.shape
    if img_shape[0] % bin_size != 0 or img_shape[1] % bin_size != 0:
        raise ValueError('Size of image_data must be a multiple of bin_size in both dimensions')
    if bin_size > 1:
        locs_y = np.arange(0, img_shape[0], bin_size)
        locs_x = np.arange(0, img_shape[1], bin_size)
        binned = np.add.reduceat(np.add.reduceat(image_data, locs_y, axis=0),
                                 locs_x, axis=1) / bin_size ** 2
    else:
        binned = image_data
    binned = ndimg.median_filter(binned, size=filt_size, mode='nearest')
    return np.repeat(np.repeat(binned, bin_size, axis=0), bin_size, axis=1)


def _pixel_stats(pixels, clip_sig=3, n_clip=10, min_pix=50):
    clip_iter = 0
    mode, scale = 0, 1
    mask = np.ones(pixels.shape, dtype=bool)
    while True:
        mode = biweight_location(pixels[mask])
        scale = biweight_midvariance(pixels[mask])
        mask &= np.abs(pixels - mode) < clip_sig * scale
        clip_iter += 1
        if np.sum(mask) < min_pix or clip_iter >= n_clip:
            break
    return mode, scale


def make_mask(img_data, subsample=1, filtsize=32, nsmooth=3,
              threshtype="sigma", nsigthresh=3, constthresh=32,
              ngrow=1, verbose=True):
    if subsample > 1:
        img_data = img_data[::subsample, ::subsample]
    if filtsize > 0:
        sky_image = median_filter(img_data, size=filtsize)
        img_data -= sky_image
    if subsample > 1:
        img_data = np.kron(img_data, np.ones((subsample, subsample)))
    if nsmooth > 0:
        img_data = median_filter(img_data, size=nsmooth)
    img_data_nan = deepcopy(img_data)
    img_data_nan[img_data_nan == 0.0] = np.nan
    mean, median, std = sigma_clipped_stats(img_data_nan, sigma=3, maxiters=10)
    thresh = median + (nsigthresh * std if threshtype == "sigma" else constthresh)
    mask = (img_data > thresh).astype(int)
    if ngrow > 0:
        mask = binary_dilation(mask, iterations=ngrow)
    mask[np.isnan(img_data_nan)] = True
    return mask


def calc_rms(image_data, weight_data=None, sky_bin=4, sky_size=25,
             autocorr_aper=5, autocorr_bg_annulus_radii=(5, 7),
             mask_smooth_size=3, mask_sigma=1.5, mask_grow_size=3,
             use_rms_filt=False, rms_filt_size=7, log_file=None):
    working_data = image_data.copy()
    bp_mask = ~np.isfinite(image_data)

    if weight_data is not None:
        bp_mask |= (weight_data <= 0) | ~np.isfinite(weight_data)
        weight_data[bp_mask] = 0.0
        weight_avg = np.mean(weight_data)
        working_data *= np.sqrt(weight_data / weight_avg)

    working_data[bp_mask] = 0.0
    working_data[bp_mask] = 0.0

    if mask_sigma > 0:
        bp_mask |= make_mask(working_data, verbose=False)

    frac_flagged = bp_mask.sum() / image_data.size
    working_data[bp_mask] = 0.0

    autocorr_rms, corr_factor = measure_autocorrelation(
        working_data, aper_radius=autocorr_aper,
        bg_annulus_radii=autocorr_bg_annulus_radii)
    autocorr_rms /= np.sqrt(1.0 - frac_flagged)
    measured_rms = autocorr_rms

    weight_nominal = 1.0 / (measured_rms * corr_factor) ** 2

    rms_info = {'frac_flagged': frac_flagged,
                'autocorr_rms': autocorr_rms,
                'corr_factor': corr_factor,
                'weight_nominal': weight_nominal}
    if weight_data is not None:
        weight_scale = weight_nominal / weight_avg
        rms_info['weight_avg'] = weight_avg
        rms_info['weight_scale'] = weight_scale

    return rms_info


# ─────────────────────────────────────────────────────────────────────────────
# Configuration
# ─────────────────────────────────────────────────────────────────────────────

import tomllib
from pathlib import Path
_CONFIG = tomllib.loads(Path(__file__).with_name('nexus.toml').read_text())
FILTERS = sys.argv[1:] or (_CONFIG['filters']['sw'] + _CONFIG['filters']['lw'])
assert set(FILTERS) <= set(_CONFIG['filters']['sw'] + _CONFIG['filters']['lw'])

# box_size per resolution (pixels)
BOX_SIZES = {"30mas": 601, "60mas": 301}

RESOLUTIONS = ['60mas']


# ─────────────────────────────────────────────────────────────────────────────
# Load selected region corners from JSON file
# ─────────────────────────────────────────────────────────────────────────────

def load_corners(filt):
    per_filter = Path("automation/hrms/per_filter/selected_regions_per_filter.json")
    if per_filter.exists():
        return json.loads(per_filter.read_text())["filters"][filt]["corners"]
    json_file = "automation/hrms/selected_regions.json"
    try:
        with open(json_file) as f:
            data = json.load(f)
    except FileNotFoundError:
        sys.exit(
            f"ERROR: {json_file} not found.\n"
            f"  Run:  python step1_select_regions.py"
        )
    return data["corners"]   # list of [y0, x0]


def corners_to_slices(corners, box_size, scale=1.0):
    """Reconstruct slice tuples from saved corner coordinates.

    Parameters
    ----------
    corners  : list of [y0, x0] derived from the 60mas image
    box_size : side length of the square region in pixels
    scale    : coordinate scaling factor (1.0 for 60mas, 2.0 for 30mas)
    """
    return [
        (slice(int(y0 * scale), int(y0 * scale) + box_size),
         slice(int(x0 * scale), int(x0 * scale) + box_size))
        for y0, x0 in corners
    ]


# ─────────────────────────────────────────────────────────────────────────────
# Main loop
# ─────────────────────────────────────────────────────────────────────────────

def main():
    # Load corners once

    for res in RESOLUTIONS:
        box_size = BOX_SIZES[res]
        scale = 2.0 if res == "30mas" else 1.0

        print(f"\n{'='*60}")
        print(f"  Resolution: {res}   box_size: {box_size}")
        print(f"{'='*60}")

        for filt in FILTERS:
            corners = load_corners(filt)
            selected_regions = corners_to_slices(corners, box_size, scale=scale)
            # Locate science image
            pattern = f"./collection/{res}/JWST_NRC_*_{filt}_{res}_sci.fits"
            matches = sorted(glob.glob(pattern))
            if not matches:
                print(f"  WARNING: no file found for {filt} {res} – skipping")
                continue

            master_sci = matches[0]
            master_wht = master_sci.replace("sci.fits", "wht.fits")

            sci_image = fits.getdata(master_sci)
            wht_image = fits.getdata(master_wht)

            bpmask = wht_image <= 0.0
            bpmask |= ~np.isfinite(wht_image)
            bpmask |= ~np.isfinite(fits.getdata(master_sci))

            rmss = []
            for slice_y, slice_x in selected_regions:
                rmsdict = calc_rms(
                    image_data=sci_image[slice_y, slice_x],
                    weight_data=wht_image[slice_y, slice_x],
                    sky_bin=1,
                )
                rmss.append(rmsdict["weight_scale"])

            w2i_factor = biweight_location(rmss, ignore_nan=True)
            print(f"  {filt}  {res}  autoreg  {np.round(w2i_factor, 4)}")

            wht_image = fits.getdata(master_wht)   # reload (was modified in loop)
            bpmask = fits.getdata(master_wht) <= 0.0
            bpmask |= ~np.isfinite(fits.getdata(master_wht))
            bpmask |= ~np.isfinite(fits.getdata(master_sci))

            wht_image *= w2i_factor
            wht_image[bpmask] = 0
            wht_image = np.sqrt(wht_image)
            irms_image = np.divide(1.0, wht_image)
            irms_image[bpmask] = 1.0e6

            out_file = master_sci.replace("sci.fits", "hrms.fits")
            hdu = fits.PrimaryHDU(irms_image, header=fits.getheader(master_sci))
            hdu.writeto(out_file, overwrite=True)

    print("\n✓ All hrms maps written.")


if __name__ == "__main__":
    main()
