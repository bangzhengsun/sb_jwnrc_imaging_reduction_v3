from pathlib import Path
import json,os,sys,logging,tomllib,subprocess,re,datetime,shutil,fcntl
from jwst.tweakreg import TweakRegStep
from jwst.associations.asn_from_list import asn_from_list
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base
from jwst.datamodels import ModelLibrary
from astropy.io import fits
B=Path.cwd();C=tomllib.loads((B/'nexus.toml').read_text());A=B/'automation/F090W/joint8_retry'; O=B/'F090W/stage2_cal_aligned_joint8_retry'
def event(state,**kw):
 with (B/'automation/reduction_log.jsonl').open('a') as f:f.write(json.dumps(dict(timestamp=datetime.datetime.now(datetime.timezone.utc).isoformat(),filter='F090W',stage='joint8_astrometry',state=state,**kw))+'\n')
if len(sys.argv)>1:
 logging.basicConfig(level=logging.INFO,format='%(levelname)s: %(message)s')
 output_dir=str(O);abs_refcat=C['astrometry']['reference_catalog'];os.environ['ASTROM_ABS_MINOBJ']=str(C['astrometry']['abs_minobj'])
 s=(B/'SBNRC_05_astrometry.py').read_text();exec(s[s.index('def run_tweakregstep(f):'):s.index('def execute_tweakreg_on_file')]);run_tweakregstep(sys.argv[1]);sys.exit()
A.mkdir(exist_ok=True);O.mkdir(exist_ok=True)
lock=(A/'lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
assert not (A/'status').exists(),'Retry already exists; inspect before restarting'
q=json.loads((A.parent/'astrometry_qc.json').read_text());prefixes=sorted({r['file'].rsplit('_',1)[0] for r in q['failures']})
manifest=[]
for prefix in prefixes:
 files=[B/'F090W/stage2_cal'/f'{prefix}_nrc{m}{d}_cal.fits' for m in 'ab' for d in range(1,5)]
 assert all(p.exists() for p in files)
 asn=asn_from_list([str(p) for p in files],rule=DMS_Level3_Base,product_name=prefix+'_joint8')
 for m in asn['products'][0]['members']:
  m['group_id']=prefix;m['tweakreg_catalog']=m['expname'].replace('_cal.fits','.phot.ecsv');assert Path(m['tweakreg_catalog']).exists()
 p=A/(prefix+'.json');p.write_text(asn.dump(format='json')[1]);assert len(ModelLibrary(str(p)).group_names)==1
 manifest.append(dict(group=prefix,association=str(p),members=[str(p) for p in files]))
(A/'manifest.json').write_text(json.dumps(manifest,indent=2));(A/'pid').write_text(str(os.getpid()));(A/'status').write_text('RUNNING\n');event('RUNNING',groups=len(manifest),frames=8*len(manifest),minimum_matches=10,reference_catalog=C['astrometry']['reference_catalog'])
rows=[]
try:
 for r in manifest:
  log=A/(r['group']+'.log')
  with log.open('w') as f:ret=subprocess.run([sys.executable,'-u',__file__,r['association']],stdout=f,stderr=subprocess.STDOUT)
  s=log.read_text();n=re.findall(r'Final solution based on (\d+) objects',s);rms=re.findall(r'FIT RMSE:\s*([\d.eE+-]+)',s)
  good=ret.returncode==0 and len(n)==len(rms)==1 and int(n[0])>=10 and float(rms[0])<=C['astrometry']['max_rmse_arcsec']
  outs=[O/(Path(p).name.replace('_cal.fits','_tweakregstep.fits')) for p in r['members']]
  if good:
   for p in outs:
    with fits.open(p,memmap=True,do_not_scale_image_data=True) as h:
     h.verify('exception');assert h[0].header['S_TWKREG']=='COMPLETE'
     for e in h:
      i=e.fileinfo();assert i['datLoc']+i['datSpan']<=p.stat().st_size
  rows.append(dict(group=r['group'],matches=[int(x) for x in n],rmse_arcsec=[float(x) for x in rms],passed=good))
  (A/'qc.json').write_text(json.dumps(rows,indent=2));print(rows[-1],flush=True)
 if not all(r['passed'] for r in rows):
  (A/'status').write_text('FAILED_QC original aligned outputs unchanged\n');event('FAILED_QC',results=rows);sys.exit(2)
 backup=B/'F090W/astrometry_before_joint8';backup.mkdir(exist_ok=True)
 # Preflight all replacement/backup paths before any mutation.
 for r in manifest:
  for p in r['members']:
   name=Path(p).name.replace('_cal.fits','_tweakregstep.fits');assert (B/'F090W/stage2_cal_aligned'/name).exists();assert not (backup/name).exists()
 for r in manifest:
  for p in r['members']:
   name=Path(p).name.replace('_cal.fits','_tweakregstep.fits');dest=B/'F090W/stage2_cal_aligned'/name
   shutil.copy2(dest,backup/name);shutil.copy2(O/name,dest)
 (A/'status').write_text('COMPLETE adopted all eight detectors for each passing group\n')
 (A.parent/'astrometry_qc.status').write_text(f'PASS_WITH_JOINT8_RETRY {len(manifest)} groups adopted; original QC JSON preserved\n')
 event('COMPLETE',groups=len(manifest),replaced_frames=8*len(manifest),results=rows,backup=str(backup))
except Exception as e:
 (A/'status').write_text('FAILED '+repr(e)+'\n');event('FAILED',error=repr(e));raise
