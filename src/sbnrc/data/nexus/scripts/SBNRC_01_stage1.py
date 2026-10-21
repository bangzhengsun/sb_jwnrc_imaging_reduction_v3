import os, sys; 
################# DEFINE PARAMETERS HERE ####################
base_path = '.'; 
filt = str(sys.argv[1]); 
uncal_path = f'{base_path}/{filt}/stage1'; 
output_dir = f'{base_path}/{filt}/stage1_cal'; 

#jump_rejection_threshold = 3.0; 
#jump_expand_large_events = True; 
#jump_expand_factor = 1.5; 
#jump_max_cores = 'all'; 
#snowball_sat_required = True; 
#snowmall_min_sat_area = 10.0; 
cpu_nap = False; 
cpu_nap_time = 10; # seconds

################## DO NOT EDIT ANYTHING BELOW THIS LINE ##################

import jwst; 

import time; 
print(f'Installed JWST package version: {jwst.__version__}'); 
time.sleep(5); 

try: 
    print(f'Current CRDS_PATH: {os.environ["CRDS_PATH"]}'); 
    print(f'Current CRDS_SERVER_URL: {os.environ["CRDS_SERVER_URL"]}'); 
except: 
    raise Exception(f'At least one of CRDS_PATH and CRDS_SERVER_URL environment variables is missing.'); 

import glob; 
from tqdm import tqdm; 

import numpy as np; 

from jwst.pipeline import calwebb_detector1; 

if not os.path.exists(uncal_path): 
    raise Exception('Directory containing uncal files does not exist.'); 

os.makedirs(output_dir, exist_ok=True); 

file_list = sorted(glob.glob(f'{uncal_path}/*uncal.fits')); 
if len(file_list) == 0: 
    raise Exception(f'No uncal files in {uncal_path}'); 

print(f'Total number of files to be calibrated: {len(file_list)}'); 
time.sleep(5); 

fail_list = []; 

for i, f in enumerate(file_list): 
    detector1 = calwebb_detector1.Detector1Pipeline(); 

    detector1.output_dir = output_dir; 
    detector1.save_results = True; 
    detector1.linearity.save_results = False; 
    detector1.jump.rejection_threshold = 3.0; 
    detector1.jump.expand_large_events = True; 
    detector1.jump.expand_factor = 1.5; 
    #detector1.jump.maximum_cores = 'quarter'; 
    detector1.jump.sat_required_snowball = True; 
    detector1.jump.min_sat_area = 10.0; 
    detector1.jump.max_jump_to_flag_neighbors = 1.0; 
    detector1.jump.min_jump_to_flag_neighbors = 100000.0; 
    try: 
        run_output = detector1.run(f); 
    except Exception:
        import traceback
        traceback.print_exc()
        fail_list.append(f)
        with open(f'{output_dir}/calwebb_detector1_failed_files.txt', 'a') as failed:
            failed.write(f + '\n')
        raise 

    if cpu_nap and (i+1) % 10 == 0.0: 
        time.sleep(cpu_rest_time); 
        print(f'Resting for {cpu_nap_time} seconds to prevent CPU from dying.'); 

print(f'calwebb_detector1 complete. Number of files failed during this step: {len(fail_list)}'); 
if len(fail_list): 
    with open(f'{output_dir}/calwebb_detector1_failed_files.txt', 'w+') as f: 
        f.write('\n'.join(fail_list) + '\n'); 
    print('A file has been saved in the output directory showing all files failed in this step.'); 
