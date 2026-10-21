##### IMPORTING PACKAGES #####
print('Initializing...'); 
import time; 
time.sleep(3); 
from astropy.wcs import WCS;
from astropy.wcs.utils import skycoord_to_pixel
from astropy.table import Table; 

from collections import defaultdict; 
import astropy.units as u; 
from astropy.coordinates import SkyCoord; 
from astropy.stats import biweight_location, sigma_clip; 
from astropy.io import ascii, fits; 
import math; 
import matplotlib.pyplot as plt; 
import glob; 
from tqdm import tqdm; 

from photutils.background import Background2D, MedianBackground; 
from photutils.segmentation import detect_threshold, detect_sources; 
from photutils.detection import DAOStarFinder, find_peaks; 
from photutils.utils import circular_footprint; 
from scipy.ndimage import gaussian_filter; 
from astropy import stats; 

import os, sys; 
###################################

filt = sys.argv[1]; 

print('Running Sex...?????'); 
# Running SExtractor for each exposure for Astrometric Alignment
# We supply custom catalogs for tweakreg for better alignment results 
output_dir = f'./{filt}/stage2_cal/'; 
#os.makedirs(output_dir, exist_ok=True); 
align_imgs = sorted(glob.glob(f'./{filt}/stage2_cal/*cal.fits')); 

trim_edge = 10; 
detect_thresh = 4.0; 

for fn in align_imgs[:]: 
    with fits.open(fn) as f1: 
        sci = f1[1].data; shdr = f1[1].header; 
        rms = f1[2].data; rhdr = f1[2].header; 
        shdu = fits.PrimaryHDU(sci, header=shdr); 
        rhdu = fits.PrimaryHDU(rms, header=rhdr); 
        shdu.writeto(fn.replace('cal.fits', 'cal_sci.fits'), overwrite=True); 
        rhdu.writeto(fn.replace('cal.fits', 'cal_rms.fits'), overwrite=True); 

    sex_command = f'sex -c sextractor.conf {fn.replace("cal.fits", "cal_sci.fits")} -DETECT_THRESH {detect_thresh} -ANALYSIS_THRESH {detect_thresh} -CATALOG_NAME {fn.replace("_cal.fits", ".sex.txt")} -WEIGHT_IMAGE {fn.replace("cal.fits", "cal_rms.fits")} -MAG_ZEROPOINT 21.0'; 
    os.system(sex_command); 

    print('Processing photometry table...'); 
    tab = ascii.read(fn.replace("_cal.fits", ".sex.txt")); 
    selects = []; 
    for i, row in enumerate(tab): 
        if trim_edge <= row['X_IMAGE'] <= sci.shape[0] and trim_edge <= row['Y_IMAGE'] <= sci.shape[0] and row['MAGERR_ISO'] <= 0.2 and row['MAG_ISO'] != 99.0 and row['ISOAREA_IMAGE'] >= 10: 
            selects.append(i); 
    catalog = tab[selects]; 
    catxy = Table([catalog['NUMBER'], catalog['XWIN_IMAGE'], 
                   catalog['YWIN_IMAGE'], catalog['FLUX_AUTO']],  
                  names=['id', 'xcentroid', 'ycentroid', 'flux']); 
    
    catxy['xcentroid']  -=  1.0
    catxy['ycentroid']  -=  1.0

    catxy.write(output_dir+os.path.basename(fn).replace("_cal.fits", ".phot.ecsv"), format='ascii.ecsv', overwrite=True); 
    os.system(f"rm {fn.replace('cal.fits', 'cal_sci.fits')}"); 
    os.system(f"rm {fn.replace('cal.fits', 'cal_rms.fits')}"); 
    #clear_output(wait=True)


print('Done! Ready for astrometry alignment!'); 


