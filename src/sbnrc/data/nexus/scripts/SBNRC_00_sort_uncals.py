import os, sys; 
import shutil; 
import glob; 
from astropy.io import fits; 

uncal_dir = './uncals'; 
print('Sorting all uncals...'); 

flist = sorted(glob.glob(f'{uncal_dir}/*uncal.fits')); 
for f in flist: 
    hdr = fits.getheader(f, ext=0); 
    filt_name = hdr['FILTER']; 
    os.makedirs(f'./{filt_name}/stage1', exist_ok=True); 
    shutil.move(f, f'./{filt_name}/stage1'); 

flist = glob.glob(f'{uncal_dir}/*uncal.fits'); 
if len(flist) == 0: 
    os.system(f'rm -r {uncal_dir}'); 

print('Job done!'); 
