# run stage2 for NGDEEP F277W and F356W
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

filt = sys.argv[1]; 

input_dir = f'./{filt}/stage1_cal'; 
output_dir = f'./{filt}/stage2_cal'; 
os.makedirs(output_dir, exist_ok=True); 

file_list = sorted(glob.glob(f'{input_dir}/*rate.fits')); 
fails = []; 
for f in tqdm(file_list): 
    try: 
        image2 = calwebb_image2.Image2Pipeline(); 
        image2.output_dir = output_dir; 
        image2.save_results = True; 
        image2.resample.skip = True; 
        image2.run(f); 
    except Exception:
        import traceback
        traceback.print_exc()
        fails.append(f)
        raise 
    clear_output(wait=True); 
if fails: 
    print(f'List of failed calibrations: {fails}'); 
    ascii.write(fails, f'{output_dir}/stage2_fails.txt', format='basic', overwrite=True); 
