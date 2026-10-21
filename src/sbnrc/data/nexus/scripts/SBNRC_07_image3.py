import tomllib
from pathlib import Path
_config = tomllib.loads(Path('/mnt/SBHPCData09/Data_Reduction/JWST_NRC_NEXUS_v2026Oct/nexus.toml').read_text())
if not _config['mosaic'].get('geometry_finalized', False):
    raise RuntimeError('Derive and validate the shared WCS from all aligned LW exposures first.')
import jwst; 
from jwst.pipeline import calwebb_image3;
from jwst.associations import asn_from_list; 
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base; 

import logging
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

from astropy.io import fits; 

import numpy as np; 
import glob; 
from tqdm import tqdm; 
import os, sys, shutil; 
import time; 
from datetime import datetime; 

base_path = '/mnt/SBHPCData09/Data_Reduction/JWST_NRC_NEXUS_v2026Oct'; 
#abs_refcat = '/mnt/Data_Analysis_1/Data_Reduction/JWST_NRC_ACO370/sex/aco370_f356w_refcat.ecsv'; 

#filt = 'F277W'; 
#wcsref = '/mnt/Data_Analysis_1/Data_Reduction/JWST_NRC_NEXUS_Ep1/sex/JWST_NEXUS_Ep1_F444W_60mas_sci.fits'; 

filt = str(sys.argv[1]); 
file_format = str(sys.argv[2]); 
#module = str(sys.argv[3]); 

os.chdir(f'{base_path}/{filt}/stage2_cal_aligned'); 

file_list = [os.path.basename(s) for s in sorted(glob.glob(f'jw*{file_format}.fits'))]; 
#file_list = []
#for part in ['002003', '002004', '002005', '002007', '002009', '034', '035']: 
#    file_list.extend(sorted(glob.glob(f'jw01837{part}*{file_format}.fits')))
print(f'Number of images to be drizzled together: {len(file_list)}'); 

asn = asn_from_list.asn_from_list(file_list, 
                                  rule=DMS_Level3_Base, 
                                  product_name = f'JWST_NRC_NEXUS_v2026Oct_{filt}_60mas'); 
with open('asn_60mas.json', 'w+') as outfile: 
    name, serialized = asn.dump(format='json'); 
    outfile.write(serialized); 


#hdr = fits.getheader(wcsref); 
#crp1 = hdr['CRPIX1']; crp2 = hdr['CRPIX2']; 
#crv1 = hdr['CRVAL1']; crv2 = hdr['CRVAL2']; 
#nax1 = hdr['NAXIS1']; nax2 = hdr['NAXIS2']; 

output_dir = f'{base_path}/{filt}/stage3_cal'; 
os.makedirs(output_dir, exist_ok=True); 
#os.chdir(f'{base_path}/{filt}/stage2_cal'); 

# make a copy if not done yet. 
#if not os.path.exists('./Orig_aligned/'): 
#    print('No backup of calibrated images, making a backup now...'); 
#    os.makedirs('./Orig_aligned', exist_ok=True); 
#    os.system(f'cp *{file_format}.fits Orig_aligned'); 

time.sleep(5)

start_time = datetime.now() 
print(f"calwebb_image3 started at {start_time.strftime('%Y-%m-%d %H:%M:%S')}")

image3 = calwebb_image3.Image3Pipeline(); 
image3.log.setLevel(logging.INFO)
image3.in_memory = _config['execution'].get('image3_in_memory_sw', True) if filt in _config['filters']['sw'] else True

image3.output_dir = output_dir; 
image3.save_results = True; 

image3.tweakreg.skip = True; 
#image3.tweakreg.save_catalogs = False; 
#image3.tweakreg.use_custom_catalogs = True; 
#image3.tweakreg.catfile = 'custom_phot_alignment.ecsv'; 
#image3.tweakreg.catalog_format = 'ecsv'; 
#image3.tweakreg.kernel_fwhm = 2.0; 
#image3.tweakreg.snr_threshold = 10.0; 
#image3.tweakreg.enforce_user_order = True; 
#image3.tweakreg.expand_refcat = False; 
#image3.tweakreg.minobj = 10; 
#image3.tweakreg.use2dhist = False; 
#image3.tweakreg.abs_refcat = abs_refcat; 
#image3.tweakreg.abs_minobj = 10; 
#image3.tweakreg.abs_use2dhist = True; 
#image3.tweakreg.save_results = True; 


image3.skymatch.skip = True; 

image3.outlier_detection.in_memory = image3.in_memory; 

image3.resample.kernel = 'square'; 
image3.resample.pixfrac = 1.0; 
image3.resample.rotation = 0.0; 
image3.resample.pixel_scale = 0.06; 
image3.resample.crpix = _config['mosaic']['crpix_step_60']; 
image3.resample.crval = [_config['mosaic']['center_ra_deg'], _config['mosaic']['center_dec_deg']]; 
image3.resample.output_shape = _config['mosaic']['shape_60']; 
image3.resample.in_memory = image3.in_memory; 

image3.source_catalog.skip = True; 

asn_file = 'asn_60mas.json'; 
image3.run(asn_file); 

end_time = datetime.now()
print(f"calwebb_image3 started at {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
print(f"calwebb_image3 completed at {end_time.strftime('%Y-%m-%d %H:%M:%S')}")

duration_sec = (end_time - start_time).total_seconds(); 
hours, rem = divmod(duration_sec, 3600)
minutes, seconds = divmod(rem, 60)
print(f"Total duration is {int(hours):02d}h {int(minutes):02d}m {seconds:05.2f}s")
