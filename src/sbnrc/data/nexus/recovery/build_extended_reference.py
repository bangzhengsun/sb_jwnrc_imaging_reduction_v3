from pathlib import Path
import os,json,datetime,subprocess
import numpy as np
from astropy.io import fits,ascii
from astropy.table import Table,vstack
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
b=Path.cwd();a=b/'automation/extended_reference';a.mkdir(exist_ok=True);(a/'status').write_text('RUNNING photometry\n');(a/'pid').write_text(str(os.getpid()))
prefix=b/'collection/60mas/JWST_NRC_NEXUS_v2026Oct_F356W_60mas';sci=Path(str(prefix)+'_sci.fits');rms=Path(str(prefix)+'_rms.fits');cat=a/'F356W.sex.txt'
with fits.open(sci,memmap=True) as h:
 hdr=h[0].header.copy();scale=abs(np.linalg.det(WCS(hdr).celestial.pixel_scale_matrix))*(np.pi/180)**2
 zp= -2.5*np.log10(1e6*scale/3631.)
cmd=['sex','-c','sextractor.conf',str(sci),'-DETECT_THRESH','4','-ANALYSIS_THRESH','4','-CATALOG_NAME',str(cat),'-WEIGHT_IMAGE',str(rms),'-MAG_ZEROPOINT',str(zp),'-CHECKIMAGE_TYPE','NONE']
with (a/'sex.log').open('w') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
t=ascii.read(cat);x=np.asarray(t['XWIN_IMAGE'])-1;y=np.asarray(t['YWIN_IMAGE'])-1;ok=np.isfinite(x)&np.isfinite(y)&(x>=10)&(y>=10)&(x<29990)&(y<29990)&(t['MAGERR_ISO']<=.2)&(t['MAG_ISO']!=99)&(t['ISOAREA_IMAGE']>=10)&np.isfinite(t['MAG_AUTO'])&(t['MAG_AUTO']!=99)
x=x[ok];y=y[ok];mag=np.asarray(t['MAG_AUTO'])[ok];ix=np.rint(x).astype(int);iy=np.rint(y).astype(int)
with fits.open(str(prefix)+'_wht.fits',memmap=True) as h:wt=h[0].data[iy,ix];good=np.isfinite(wt)&(wt>0)
with fits.open(b/'collection/60mas/JWST_NRC_NEXUS_v2026Oct_F444W_60mas_wht.fits',memmap=True) as h:wt=h[0].data[iy,ix];extra=~(np.isfinite(wt)&(wt>0))
keep=good&extra;x=x[keep];y=y[keep];mag=mag[keep];ra,de=WCS(hdr).celestial.all_pix2world(x,y,0)
old=Table.read(b/'nexus_f444w_refcat.ecsv');new=Table([ra,de,mag],names=['RA','DEC','MAG']);coords=SkyCoord(ra*u.deg,de*u.deg);orig=SkyCoord(np.asarray(old['RA'])*u.deg,np.asarray(old['DEC'])*u.deg);idx,sep,_=coords.match_to_catalog_sky(orig);new=new[sep.arcsec>.2]
# Keep original columns and values intact; reference consumer requires RA/DEC/MAG.
assert set(old.colnames)=={'RA','DEC','MAG'},old.colnames
old['REF_FILTER']=np.full(len(old),'F444W');new['REF_FILTER']=np.full(len(new),'F356W');combined=vstack([old,new],metadata_conflicts='silent');combined.meta['extension_note']='F356W sources outside F444W positive-weight footprint; existing coordinates and magnitudes preserved; 0.2 arcsec deduplication.'
out=b/'nexus_f444w_f356w_extended_refcat.ecsv';assert not out.exists();combined.write(out);new.write(a/'added_F356W_sources.ecsv',overwrite=True)
report={'original_sources':len(old),'added_sources':len(new),'total_sources':len(combined),'zeropoint_AB':zp,'dedup_radius_arcsec':.2,'catalog':str(out)};(a/'report.json').write_text(json.dumps(report,indent=2));(a/'status').write_text('COMPLETE pending spatial QC\n')
with (b/'automation/reduction_log.jsonl').open('a') as f:f.write(json.dumps({'timestamp':datetime.datetime.now(datetime.timezone.utc).isoformat(),'stage':'extended_reference','state':'PHOTOMETRY_COMPLETE',**report})+'\n')
print(json.dumps(report),flush=True)
