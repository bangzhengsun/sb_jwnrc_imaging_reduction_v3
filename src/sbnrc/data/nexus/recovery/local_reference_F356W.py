from pathlib import Path
import os,json,logging,datetime,tomllib,re
import numpy as np
from astropy.io import fits,ascii
from astropy.table import Table
from astropy.wcs import WCS
from scipy.ndimage import binary_erosion
from jwst.resample import ResampleStep
from jwst.tweakreg import TweakRegStep
from jwst.associations.asn_from_list import asn_from_list
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base
import subprocess
B=Path.cwd();cfg=tomllib.loads((B/'nexus.toml').read_text());A=B/'automation/F356W/local_reference';A.mkdir(exist_ok=True)
logging.basicConfig(level=logging.INFO,format='%(levelname)s: %(message)s')
(A/'pid').write_text(str(os.getpid()));(A/'status').write_text('RUNNING\n')
def event(stage,**kw):
 with (B/'automation/reduction_log.jsonl').open('a') as f:f.write(json.dumps({'timestamp':datetime.datetime.now(datetime.timezone.utc).isoformat(),'filter':'F356W','stage':stage,**kw})+'\n')
qc={r['prefix']:r for r in json.loads((B/'automation/F356W/paired_retry/qc.json').read_text())}
s=(B/'SBNRC_05_astrometry.py').read_text();exec(s[s.index('def run_tweakregstep(f):'):s.index('def execute_tweakreg_on_file')])
for block in cfg['local_reference_f356w']:
 visit=block['visit'];d=A/visit;d.mkdir(exist_ok=True);paths=[]
 for prefix in block['seed_groups']:
  assert qc[prefix]['pass'] and qc[prefix]['final_matches']>=16
  paths += [str(B/'F356W/stage2_cal_aligned_paired_retry'/f'{prefix}_{m}_tweakregstep.fits') for m in ['nrcalong','nrcblong']]
 assert len(paths)==len(set(paths)) and all(Path(p).exists() for p in paths)
 asn=asn_from_list(paths,rule=DMS_Level3_Base,product_name='local_'+visit);ap=d/'drizzle_asn.json';ap.write_text(asn.dump(format='json')[1]);out=d/f'local_{visit}_i2d.fits';assert not out.exists()
 event('local_reference_drizzle',state='RUNNING',visit=visit,inputs=paths,pixel_scale=.06,kernel='square',pixfrac=1.,note='User-authorized temporary astrometric reference; directly resample aligned seed groups, no target exposure included.')
 step=ResampleStep();step.output_dir=str(d);step.output_file='local_'+visit;step.save_results=True;step.pixel_scale=.06;step.rotation=0.;step.pixfrac=1.;step.kernel='square';step.in_memory=True
 result=step.run(str(ap));result.save(str(out));result.close()
 with fits.open(out,memmap=True) as h:
  h.verify('exception');sci=h['SCI'].data.copy();err=h['ERR'].data.copy();weight=h['WHT'].data;hdr=h['SCI'].header.copy();w=WCS(hdr).celestial;valid=np.isfinite(sci)&np.isfinite(err)&(err>0)&np.isfinite(weight)&(weight>0)
  assert sci.shape[0]<10000 and sci.shape[1]<10000
  assert len(h['HDRTAB'].data)==len(paths)
  bounds=[]
  for p in paths:
   ih=fits.getheader(p,'SCI');iw=WCS(ih).celestial;nx=ih['NAXIS1'];ny=ih['NAXIS2'];t=np.linspace(0,1,50);x=np.r_[t*(nx-1),np.full(50,nx-1),t[::-1]*(nx-1),np.zeros(50)];y=np.r_[np.zeros(50),t*(ny-1),np.full(50,ny-1),t[::-1]*(ny-1)];ra,de=iw.all_pix2world(x,y,0);px,py=w.all_world2pix(ra,de,0);assert px.min()>-2 and py.min()>-2 and px.max()<sci.shape[1]+1 and py.max()<sci.shape[0]+1;bounds.append([float(px.min()),float(px.max()),float(py.min()),float(py.max())])
 good=binary_erosion(valid,iterations=10);sci[~valid]=0;err[~valid]=1e10
 sf=d/'sci.fits';ef=d/'rms.fits';fits.writeto(sf,sci,hdr,overwrite=True);fits.writeto(ef,err,hdr,overwrite=True)
 cat=d/'sources.sex.txt';cmd=['sex','-c','sextractor.conf',str(sf),'-DETECT_THRESH','4','-ANALYSIS_THRESH','4','-CATALOG_NAME',str(cat),'-WEIGHT_IMAGE',str(ef),'-MAG_ZEROPOINT','21','-CHECKIMAGE_TYPE','NONE']
 with (d/'sex.log').open('w') as log:subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
 tab=ascii.read(cat);x=np.asarray(tab['XWIN_IMAGE'])-1;y=np.asarray(tab['YWIN_IMAGE'])-1;ok=np.isfinite(x)&np.isfinite(y)&(x>=0)&(y>=0)&(x<sci.shape[1])&(y<sci.shape[0])
 ix=np.clip(np.rint(np.nan_to_num(x)).astype(int),0,sci.shape[1]-1);iy=np.clip(np.rint(np.nan_to_num(y)).astype(int),0,sci.shape[0]-1)
 ok &= good[iy,ix]&(tab['MAGERR_ISO']<=.2)&(tab['MAG_ISO']!=99)&(tab['ISOAREA_IMAGE']>=10)
 ra,de=w.all_pix2world(x[ok],y[ok],0);ref=Table([ra,de,np.asarray(tab['MAG_AUTO'])[ok]],names=['RA','DEC','MAG']);assert len(ref)>=10 and np.isfinite(ra).all() and np.isfinite(de).all();ref.meta['seed_groups']=block['seed_groups'];ref.meta['note']='Local relative reference anchored by paired fits to nexus_f444w_refcat; no retry targets included';rp=d/'reference.ecsv';ref.write(rp,overwrite=True)
 event('local_reference_catalog',state='COMPLETE',visit=visit,sources=len(ref),shape=list(sci.shape),footprint_bounds=bounds,reference=str(rp));print('REFERENCE_COMPLETE',visit,len(ref),sci.shape,flush=True)
 abs_refcat=str(rp);output_dir=str(d/'aligned');Path(output_dir).mkdir(exist_ok=True)
 for prefix in block['target_groups']:
  assert prefix not in block['seed_groups'];members=[str(B/'F356W/stage2_cal'/f'{prefix}_{m}_cal.fits') for m in ['nrcalong','nrcblong']];asn=asn_from_list(members,rule=DMS_Level3_Base,product_name=prefix+'_local')
  for member in asn['products'][0]['members']:member['group_id']=prefix;member['tweakreg_catalog']=member['expname'].replace('_cal.fits','.phot.ecsv')
  ap=d/(prefix+'.json');ap.write_text(asn.dump(format='json')[1]);os.environ['ASTROM_ABS_MINOBJ']='5';print('LOCAL_PAIR_BEGIN',prefix,flush=True);run_tweakregstep(str(ap));print('LOCAL_PAIR_END',prefix,flush=True)
(A/'status').write_text('PROCESS_COMPLETE numerical QC pending\n')
