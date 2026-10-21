##### IMPORTING PACKAGES #####
print('Initializing...'); 
import time; 
time.sleep(3); 
import numpy as np; 
import os, sys, shutil; 
from astropy.io import fits; 

from jwst.pipeline import calwebb_detector1; 
from jwst.pipeline import calwebb_image2; 
from jwst.pipeline import calwebb_image3; 
from jwst.background import BackgroundStep;
from jwst.assign_wcs import AssignWcsStep;
from jwst.flatfield import FlatFieldStep;
from jwst.photom import PhotomStep;
from jwst.resample import ResampleStep;
from jwst import datamodels;
from jwst.associations import asn_from_list; 
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base; 
from jwst.tweakreg import TweakRegStep; 

import logging
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

import glob; 
from tqdm import tqdm; 
#from astroscrappy import detect_cosmics; 

from astropy.wcs import WCS; 
from astropy.table import Table; 

from collections import defaultdict; 
import astropy.units as u; 
from astropy.coordinates import SkyCoord; 
import math; 
import matplotlib.pyplot as plt; 

from photutils.background import Background2D, MedianBackground; 
from photutils.segmentation import detect_threshold, detect_sources; 
from photutils.utils import circular_footprint; 
from scipy.ndimage import gaussian_filter; 
from astropy import stats; 

from IPython.display import clear_output; 

import jhat; 
###################################

filt = sys.argv[1]; 

print('Running astrometric alignment...'); 
# Astrometric Alignment
basepath = '/mnt/SBHPCData09/Data_Reduction/JWST_NRC_NEXUS_v2026Oct'; 

pwd = os.getcwd(); 
input_dir = f'{basepath}/{filt}/stage2_cal/'; 
output_dir = f'{basepath}/{filt}/stage2_cal_aligned/'; 
os.makedirs(output_dir, exist_ok=True); 
import tomllib
with open(os.path.join(basepath, 'nexus.toml'), 'rb') as config_handle:
    abs_refcat = tomllib.load(config_handle)['astrometry']['reference_catalog'] 

os.chdir(input_dir); 
def run_tweakregstep(f): 
    tweakregstep = TweakRegStep(); 
    tweakregstep.log.setLevel(logging.INFO)
    tweakregstep.output_dir = output_dir; 
    tweakregstep.save_results = True; 
    tweakregstep.save_catalogs = False; 
    tweakregstep.use_custom_catalogs = True; 
    #tweakregstep.catfile = f'custom_phot_alignment.ecsv'; 
    ### Source detection parameters, useless if custom catalos are provided 
    tweakregstep.kernel_fwhm = 2.0; 
    tweakregstep.snr_threshold = 5.0; 
    tweakregstep.enforce_user_order = True; 
    tweakregstep.expand_refcat = False; 
    ### relative alignment parameters, note that tolerance should be 
    ### smaller than tolerance to avoid duplicate matches 
    tweakregstep.minobj = 10; 
    tweakregstep.use2dhist = False; # no need to use 2D histogram 
    tweakregstep.fitgeometry = 'shift'; 
    tweakregstep.separation = 0.2; 
    tweakregstep.tolerance = 0.1; 
    tweakregstep.searchrad = 3.0; 
    tweakregstep.nclip = 1; 
    tweakregstep.sigma = 0.7; 
    ### absolute alignment 
    tweakregstep.abs_minobj = int(os.environ.get('ASTROM_ABS_MINOBJ', '10')); 
    tweakregstep.abs_searchrad = 3.0; 
    tweakregstep.abs_nclip = 1; 
    tweakregstep.abs_sigma = 0.7; 
    tweakregstep.abs_separation = 0.2; 
    tweakregstep.abs_tolerance = 0.1;   
    tweakregstep.abs_refcat = abs_refcat; 
    #tweakregstep.abs_minobj = 10; 
    tweakregstep.abs_use2dhist = True; 
    tweakregstep.run(f); 
    
def execute_tweakreg_on_file(img): 
    asn_name = img.replace('cal.fits', 'tweak_asn')
    asn = asn_from_list.asn_from_list(list([img]), rule=DMS_Level3_Base, product_name=asn_name); 
    asn_file = '.'.join([asn_name, 'json']); 
    for m, member in enumerate(asn['products'][0]['members']): 
        cat = asn['products'][0]['members'][m]['expname'].replace('_cal.fits', '.phot.ecsv'); 
        asn['products'][0]['members'][m]['tweakreg_catalog'] = cat; 
    with open(asn_file, 'w') as outfile: 
        name, serialized = asn.dump(format='json'); 
        outfile.write(serialized); 
    run_tweakregstep(asn_file); 
    clear_output(wait=True); 

#align_imgs = sorted(glob.glob('jw01837002007*cal.fits')) + sorted(glob.glob('jw01837002009*cal.fits')); 
#align_imgs = sorted(glob.glob('jw01837004002*cal.fits'))
align_imgs = sorted(glob.glob('*cal.fits'))
print(f"{len(align_imgs)} images to be aligned..."); 
time.sleep(3); 
#n_procs = 16; 

failures = []
for img in tqdm(align_imgs[:]): 
    try:
        execute_tweakreg_on_file(img); 
    except Exception as exc:
        failures.append((img, repr(exc)))
        logging.exception("Astrometry failed for %s", img)

failure_path = os.path.join(basepath, 'automation', filt, 'astrometry_failures.txt')
os.makedirs(os.path.dirname(failure_path), exist_ok=True)
with open(failure_path, 'w') as handle:
    for img, exc in failures:
        handle.write(f'{img}\t{exc}\n')
if failures:
    raise RuntimeError(f'{len(failures)} astrometry failures; see {failure_path}')

os.chdir(pwd); 
