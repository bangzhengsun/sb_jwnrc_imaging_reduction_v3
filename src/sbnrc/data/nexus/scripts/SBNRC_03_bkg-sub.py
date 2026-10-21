######################## CHOOSE WHETHER TO SUBTRACT BACKGROUND HERE #####################
bkg_subtract = 'Y'; 
astrometry = 'N'; 
##### DEFINE REFERENCE CATALOG BELOW #####
refcat = '../../../IDF_F356W_refcat.txt'; 
#########################################################################################

print('Initiating... importing packages'); 
from jwst.pipeline import calwebb_image3; 

from jwst.associations import asn_from_list; 
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base; 
from jwst.tweakreg import TweakRegStep; 

from astropy.io import ascii, fits; 
from astropy.visualization import ImageNormalize, ZScaleInterval; 
from astropy import stats; 
from astropy.wcs import WCS; 
from astropy.convolution import convolve, Gaussian2DKernel; 
from astropy.modeling.models import Gaussian2D; 
from astropy.table import Table; 

from matplotlib import pyplot as plt; 
import numpy as np; 
import glob; 
from tqdm import tqdm; 
import os, sys, shutil; 
from IPython.display import clear_output; 

from photutils.background import Background2D, MedianBackground; 
from photutils.segmentation import detect_threshold, detect_sources; 
from photutils.utils import circular_footprint; 
from scipy.ndimage import gaussian_filter; 

import jhat; 

filt = str(sys.argv[1]); 
os.chdir(f'./{filt}'); 

def cal_bkg_estimate(jwcal, sigma_thresh=5.0, min_area=10, clip_sigma=3.0, box_size=128, exc_percent=20.0): 
    f1 = fits.open(jwcal); 
    sci = f1[1].data; 
    sci[sci==0.0] = np.nan; # in order for background estimator to automatically mask those values

    sigma_clip = stats.SigmaClip(sigma=clip_sigma, maxiters=10); 
    footprint = circular_footprint(radius=10); 
    bkg_estimator = MedianBackground(); 

    thresh = detect_threshold(sci, nsigma=sigma_thresh, sigma_clip=sigma_clip); 
    segm = detect_sources(sci, thresh, npixels=min_area); 
    maskk = segm.make_source_mask(footprint=footprint); 

    bkg = Background2D(sci, (box_size,box_size), filter_size=(5,5), mask=maskk, 
                       sigma_clip=sigma_clip, bkg_estimator=bkg_estimator,
                       exclude_percentile=exc_percent); 
    med = bkg.background_median; 
    background = bkg.background; 
    f1.close(); 
    return med, background; 

if bkg_subtract == 'Y': 
	cal_dir = './stage2_cal/'; 
	cal_list = sorted(glob.glob(f'{cal_dir}/*cal.fits')); 
	os.makedirs(f'{cal_dir}/Orig_cal', exist_ok=True); 

	for cal in tqdm(cal_list): 
	    shutil.copy(cal, f'{cal_dir}/Orig_cal/'); 
	    med, bkg = cal_bkg_estimate(cal, sigma_thresh=3.0, min_area=10, 
		                        clip_sigma=3.0, box_size=128, 
		                        exc_percent=80.0); 
	    with fits.open(cal) as f1: 
		    data = f1[1].data; hdr = f1[1].header; 
	    new_data = np.subtract(data, bkg); 
	    fits.update(cal, new_data, header=hdr, ext=1); 

if astrometry == 'Y': 
    output_dir = './stage2_cal_jhat_aligned'; 
    os.makedirs(output_dir, exist_ok=True); 
    align_imgs = sorted(glob.glob('./stage2_cal/*cal.fits')); 
    wcs_align = jhat.st_wcs_align(); 
    for align_image in tqdm(align_imgs): 
        wcs_align.run_all(align_image, 
                        telescope='jwst', 
                        outsubdir=output_dir, 
                        overwrite=True, 
                        d2d_max=1.0, 
                        refcatname=refcat, 
                        refcat_racol='ra', 
                        refcat_deccol='dec', 
                        refcat_magcol='mag', 
                        refcat_magerrcol='magerr', 
                        histocut_order='dxdy', 
                        SNR_min=5, 
                        objmag_lim=(14,26), 
                        #sharpness_lim=(0.3,0.9), 
                        #roundness1_lim=(-0.7,0.7),
                        ); 

os.chdir('../'); 
