from pathlib import Path
import json,datetime,subprocess,sys,zipfile
import numpy as np
from astropy.io import fits
b=Path.cwd();a=b/'automation/F210M';col=b/'collection/60mas';prefix='JWST_NRC_NEXUS_v2026Oct_F210M_60mas';(a/'delivery.status').write_text('RUNNING HRMS validation and ZIP\n')
with fits.open(col/(prefix+'_hrms.fits'),memmap=True) as h, fits.open(col/(prefix+'_wht.fits'),memmap=True) as w:
 assert h[0].data.shape==w[0].data.shape==(30000,30000)
 for k in ['CRVAL1','CRVAL2','CRPIX1','CRPIX2','CTYPE1','CTYPE2']:assert h[0].header[k]==w[0].header[k]
 for y in range(0,30000,300):
  wt=w[0].data[y:y+300];v=np.isfinite(wt)&(wt>0);r=h[0].data[y:y+300];assert np.isfinite(r[v]).all() and (r[v]>0).all()
paths=[col/(prefix+'_'+s+'.fits') for s in ['sci','hrms','wht']];subprocess.run([sys.executable,'bulk_zipfiles.py',*[str(p) for p in paths]],check=True)
for p in paths:
 with zipfile.ZipFile(str(p)+'.zip') as z:assert z.namelist()==[p.name] and z.getinfo(p.name).file_size==p.stat().st_size and z.testzip() is None
(a/'delivery.status').write_text('COMPLETE 60mas SCI RMS WHT HRMS; SCI HRMS WHT ZIP CRC PASS\n')
with (b/'automation/reduction_log.jsonl').open('a') as f:f.write(json.dumps({'timestamp':datetime.datetime.now(datetime.timezone.utc).isoformat(),'filter':'F210M','stage':'delivery','state':'COMPLETE','archives':3,'HRMS':'finite positive where weight positive; WCS PASS','ZIP':'CRC PASS'})+'\n')
