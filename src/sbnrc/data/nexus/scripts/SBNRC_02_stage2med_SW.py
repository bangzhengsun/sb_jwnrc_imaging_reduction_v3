import jwst; 
print(f'Current JWST package version: {jwst.__version__}'); 

from jwst.pipeline import calwebb_image2; 

from jwst.background import BackgroundStep;
from jwst.assign_wcs import AssignWcsStep;
from jwst.flatfield import FlatFieldStep;
from jwst.photom import PhotomStep;
from jwst.resample import ResampleStep;
from jwst import datamodels;

from astropy.io import fits, ascii; 
from astropy.visualization import ImageNormalize, ZScaleInterval; 
from astropy import stats; 
from astropy.wcs import WCS; 

from photutils.background import Background2D, MedianBackground; 
from photutils.segmentation import detect_threshold, detect_sources; 
from photutils.utils import circular_footprint; 

from jwst.pipeline import calwebb_image3; 

from jwst.associations import asn_from_list; 
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base; 

from astropy.io import ascii, fits; 
from astropy.visualization import ImageNormalize, ZScaleInterval; 

from matplotlib import pyplot as plt; 
import numpy as np; 
import glob; 
import shutil; 
import os, sys; 
from tqdm import tqdm; 
from IPython.display import clear_output; 
import time; 

time.sleep(5); 

def cal_bkg_estimate(jwcal, sigma_thresh=5.0, min_area=5, clip_sigma=3.0, box_size=128, exc_percent=20.0): 
    f1 = fits.open(jwcal); 
    sci = f1[1].data; 
    sci[sci==0.0] = np.nan; # in order for background estimator to automatically mask those values

    sigma_clip = stats.SigmaClip(sigma=clip_sigma, maxiters=10); 
    footprint = circular_footprint(radius=25); 
    bkg_estimator = MedianBackground(); 

    thresh = detect_threshold(sci, nsigma=sigma_thresh, sigma_clip=sigma_clip); 
    segm = detect_sources(sci, thresh, npixels=min_area); 
    maskk = segm.make_source_mask(footprint=footprint); 

    bkg = Background2D(sci, (box_size,box_size), filter_size=(3,3), mask=maskk, 
                       sigma_clip=sigma_clip, bkg_estimator=bkg_estimator,
                       exclude_percentile=exc_percent); 
    med = bkg.background_median; 
    background = bkg.background; 
    f1.close(); 
    return med, background; 

def source_mask(data, sigma_thresh=5.5, min_area=10.0, clip_sigma=3.0, radius=10): 
    sigma_clip = stats.SigmaClip(sigma=clip_sigma, maxiters=10); 
    footprint = circular_footprint(radius=radius); 
    bkg_estimator = MedianBackground(); 
    
    thresh = detect_threshold(data, nsigma=sigma_thresh, sigma_clip=sigma_clip); 
    segm = detect_sources(data, thresh, npixels=min_area); 
    maskk = segm.make_source_mask(footprint=footprint); 
    return maskk; 

filt = sys.argv[1]; 

s2fails_wcs = []; 
stage1cal_dir = f'./{filt}/stage1_cal/'; 
s2rfiles = [os.path.basename(s) for s in sorted(glob.glob(f'./{filt}/stage1_cal/*rate.fits'))]; 
output_dir = f'./{filt}/stage2_medstep_ff'; 
os.makedirs(output_dir, exist_ok=True); 

# do assign_wcs +flatfield step
for f in s2rfiles: 
    try: 
        assign_wcs_step = AssignWcsStep(); 
        assign_wcs_step.output_dir = output_dir; 
        assign_wcs_step.save_results = True; 
        assign_wcs_step.run(f'{stage1cal_dir}/{f}'); 
        
        flatfield_step = FlatFieldStep(); 
        flatfield_step.output_dir = output_dir; 
        flatfield_step.save_results = True; 
        wcs_fn = f.replace('rate.fits', 'assignwcsstep.fits'); 
        flatfield_step.run(f'{output_dir}/{wcs_fn}'); 
    except: 
        s2fails_wcs.append(f); 
    if len(s2fails_wcs): 
        with open(f'{output_dir}/ffstep_fails.txt', 'w+') as f: 
            for fail in s2fails_wcs: 
                f.write(fail); 

rate_dir = f'./{filt}/stage1_cal/'; 
s2ffiles = sorted(glob.glob(f'./{filt}/stage2_medstep_ff/*flatfieldstep.fits')); 
output_dir = f'./{filt}/stage2_medstep_bkgsubtract'; 
os.makedirs(output_dir, exist_ok=True); 
recal_fails = []; 
for f in s2ffiles: 
    try: 
        med, bkg = cal_bkg_estimate(f, sigma_thresh=3.0, min_area=5.0, clip_sigma=3.0, box_size=128, exc_percent=80.0); 
        aa = f.index('jw'); bb = f.index('flatfieldstep.fits'); 
        rate_fn = f[aa:bb]+'rate.fits'; 
        shutil.copy(f'{rate_dir}/{rate_fn}', output_dir); 
        with fits.open(f'{rate_dir}/{rate_fn}') as f1: 
            data = f1[1].data; hdr = f1[1].header; 
        newdata = np.subtract(data, med); 
        fits.update(f'{output_dir}/{rate_fn}', newdata, header=hdr, ext=1); 
    except: 
        recal_fails.append(f); 
    if len(recal_fails): 
        with open(f'{output_dir}/recal_fails.txt', 'w+') as f: 
            for fail in recal_fails: 
                f.write(fail); 


