import jwst; 
print(f'Current JWST package version: {jwst.__version__}'); 

from jwst.pipeline import calwebb_image2; 

from jwst.background import BackgroundStep;
from jwst.assign_wcs import AssignWcsStep;
from jwst.flatfield import FlatFieldStep;
from jwst.photom import PhotomStep;
from jwst.resample import ResampleStep;
from jwst import datamodels;

from jwst.pipeline import calwebb_image3; 

from jwst.associations import asn_from_list; 
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base; 

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
import time; 
import os, sys, shutil; 
from IPython.display import clear_output; 

RD = os.environ['Research_Data']; 
sys.path.insert(1, f'{RD}/python_modules'); 
from image1overf import sub1fimaging; 

from photutils.background import Background2D, MedianBackground, BiweightLocationBackground; 
from photutils.segmentation import detect_threshold, detect_sources; 
from photutils.utils import circular_footprint; 
from scipy.ndimage import gaussian_filter; 

#import jhat; 

time.sleep(5); 
filt = str(sys.argv[1]); 
os.chdir(f'./{filt}'); 

print(f'Correction 1/f noise...'); 
time.sleep(5); 
flist = sorted(glob.glob('./stage2_cal_aligned/*tweakregstep.fits')); 
for cal2file in tqdm(flist): 
    cal21overffile = cal2file.replace('_tweakregstep.fits','_cal_tweakreg_1overf.fits'); 
    #print ('Running 1/f correction on {} to produce {}'.format(cal2file,cal21overffile)); 
    with fits.open(cal2file) as cal2hdulist:
        if cal2hdulist['PRIMARY'].header['SUBARRAY']=='FULL' or cal2hdulist['PRIMARY'].header['SUBARRAY']=='SUB256':
            sigma_bgmask=3.0; 
            sigma_1fmask=2.0; 
            splitamps=False;   #Set to True only in a sparse field so each amplifier will be fit separately. 
            correcteddata = sub1fimaging(cal2hdulist,sigma_bgmask,sigma_1fmask,splitamps); 
            if cal2hdulist['PRIMARY'].header['SUBARRAY']=='FULL':
                cal2hdulist['SCI'].data[4:2044,4:2044] = correcteddata;   
            elif cal2hdulist['PRIMARY'].header['SUBARRAY']=='SUB256':
                cal2hdulist['SCI'].data[:252,:252] = correcteddata; 
            cal2hdulist.writeto(cal21overffile, overwrite=True); 
    clear_output(wait=True); 
print('Job complete! Images ready for final drizzling!'); 
os.chdir('../'); 
