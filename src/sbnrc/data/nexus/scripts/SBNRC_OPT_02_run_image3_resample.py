import tomllib
from pathlib import Path
_config = tomllib.loads(Path('/mnt/SBHPCData09/Data_Reduction/JWST_NRC_NEXUS_v2026Oct/nexus.toml').read_text())
if not _config['mosaic'].get('geometry_finalized', False):
    raise RuntimeError('Derive and validate the shared WCS from all aligned LW exposures first.')
import jwst; 
from jwst.pipeline import calwebb_image3;
from jwst.associations import asn_from_list; 
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base; 
from jwst.resample import ResampleStep; 

import logging
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

from astropy.io import fits; 

import numpy as np; 
import glob; 
from tqdm import tqdm; 
import os, sys, shutil; 
import time; 
from datetime import datetime; 

filt = sys.argv[1].upper()

os.chdir(f"./{filt}")

#wcsref = '/mnt/Data_Storage2/Data_Reduction_Archive/JWST_EGS/others/hlsp_candels_hst_wfc3_egs-tot-60mas_f160w_v1.0_drz.fits'; 
#hdr = fits.getheader(wcsref); 
#crp1 = hdr['CRPIX1']; crp2 = hdr['CRPIX2']; 
#crv1 = hdr['CRVAL1']; crv2 = hdr['CRVAL2']; 
#nax1 = hdr['NAXIS1']; nax2 = hdr['NAXIS2']; 

start_time = datetime.now() 
print(f"Drizzling started at {start_time.strftime('%Y-%m-%d %H:%M:%S')}")

resample = ResampleStep(); 

resample.save_results = True; 
resample.kernel = 'square'; 
resample.pixfrac = 1.0; 
resample.rotation = 0.0; 
resample.pixel_scale = 0.03; 
resample.crpix = _config['mosaic']['crpix_step_30']; 
resample.crval = [_config['mosaic']['center_ra_deg'], _config['mosaic']['center_dec_deg']]; 
resample.output_shape = _config['mosaic']['shape_30']; 
resample.in_memory = True; 

flist = sorted(glob.glob('stage3_cal/*crf.fits')); 
print(f'Total number of files to be drizzled: {len(flist)}'); 
if not flist:
    raise RuntimeError(f'No CRF inputs found for {filt}')
time.sleep(5); 
resample.run(flist); 

end_time = datetime.now()
print(f"Drizzling started at {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
print(f"Drizzling completed at {end_time.strftime('%Y-%m-%d %H:%M:%S')}")

duration_sec = (end_time - start_time).total_seconds(); 
hours, rem = divmod(duration_sec, 3600)
minutes, seconds = divmod(rem, 60)
print(f"Total duration is {int(hours):02d}h {int(minutes):02d}m {seconds:05.2f}s")
