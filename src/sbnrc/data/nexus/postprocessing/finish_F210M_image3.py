from pathlib import Path
import json,datetime,os
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from run_sci_rms_wht_extraction import extract
b=Path.cwd();a=b/'automation/F210M';expected={p.name.removesuffix('_uncal.fits') for p in (b/'F210M/stage1').glob('*_uncal.fits')};assert len(expected)==384
p=b/'F210M/stage3_cal/JWST_NRC_NEXUS_v2026Oct_F210M_60mas_i2d.fits'
crfs=list(p.parent.glob('*_crf.fits'));assert {x.name.split('_cal_tweakreg')[0] for x in crfs}==expected
for x in crfs:
 with fits.open(x,memmap=True,do_not_scale_image_data=True) as h:
  h.verify('exception');assert h[0].header['S_OUTLIR']=='COMPLETE'
  for e in h:
   i=e.fileinfo();assert i['datLoc']+i['datSpan']<=x.stat().st_size
with fits.open(p,memmap=True) as h:
 h.verify('exception');assert h['SCI'].shape==(30000,30000);assert h[0].header['S_RESAMP']=='COMPLETE'
 for e in h:
  i=e.fileinfo();assert i['datLoc']+i['datSpan']<=p.stat().st_size
 assert {str(n).split('_cal_tweakreg')[0] for n in h['HDRTAB'].data['FILENAME']}==expected
 w=WCS(h['SCI'].header).celestial;ref=WCS(fits.getheader(b/'collection/60mas/JWST_NRC_NEXUS_v2026Oct_F444W_60mas_sci.fits')).celestial
 assert np.allclose(w.wcs.crpix,ref.wcs.crpix,rtol=0,atol=1e-7);assert np.allclose(w.wcs.crval,ref.wcs.crval,rtol=0,atol=1e-10);assert np.allclose(w.pixel_scale_matrix,ref.pixel_scale_matrix,rtol=0,atol=1e-14)
 valid_count=0
 for y in range(0,30000,300):
  wt=h['WHT'].data[y:y+300];v=np.isfinite(wt)&(wt>0);s=h['SCI'].data[y:y+300];r=h['ERR'].data[y:y+300];assert np.isfinite(s[v]).all() and np.isfinite(r[v]).all() and (r[v]>0).all();valid_count+=int(v.sum())
report={'state':'PASS','exposures':384,'shape':[30000,30000],'WCS':'PASS matches F444W','valid_pixels':valid_count};(a/'image3_qc.json').write_text(json.dumps(report,indent=2));(a/'image3.status').write_text('COMPLETE expected=384 actual=384 WCS PASS\n')
with (b/'automation/reduction_log.jsonl').open('a') as f:f.write(json.dumps({'timestamp':datetime.datetime.now(datetime.timezone.utc).isoformat(),'filter':'F210M','stage':'image3_qc',**report})+'\n')
dest=b/'collection/60mas'/p.name
if not dest.exists():os.link(p,dest)
assert os.path.samefile(p,dest)
(a/'extraction.status').write_text('RUNNING\n');extract(dest,dest.parent);(a/'extraction.status').write_text('COMPLETE SCI RMS WHT; HRMS pending\n')
