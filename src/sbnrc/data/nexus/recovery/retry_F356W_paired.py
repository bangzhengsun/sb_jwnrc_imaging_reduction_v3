from pathlib import Path
import os,json,datetime,logging,re,sys
from jwst.tweakreg import TweakRegStep
from jwst.associations.asn_from_list import asn_from_list
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base
from jwst.datamodels import ModelLibrary
B=Path.cwd(); A=B/'automation/F356W/paired_retry';A.mkdir(exist_ok=True)
output_dir=str(B/'F356W/stage2_cal_aligned_paired_retry');Path(output_dir).mkdir(exist_ok=True)
abs_refcat=str(B/'nexus_f444w_refcat.ecsv')
logging.basicConfig(level=logging.INFO,format='%(levelname)s: %(message)s')
s=(B/'SBNRC_05_astrometry.py').read_text();func=s[s.index('def run_tweakregstep(f):'):s.index('def execute_tweakreg_on_file')];exec(func)
failed=json.loads((B/'automation/F356W/astrometry_qc_failures.json').read_text());prefixes=sorted({r['exposure'].rsplit('_',1)[0] for r in failed})
manifest=[]
for prefix in prefixes:
 files=[B/'F356W/stage2_cal'/f'{prefix}_{m}_cal.fits' for m in ['nrcalong','nrcblong']]
 assert all(p.exists() for p in files),files
 minimum=5 if prefix[7:13] in ['001033','001016'] else 10
 asn=asn_from_list([str(p) for p in files],rule=DMS_Level3_Base,product_name=prefix+'_paired')
 for member in asn['products'][0]['members']:
  member['group_id']=prefix;member['tweakreg_catalog']=member['expname'].replace('_cal.fits','.phot.ecsv');assert Path(member['tweakreg_catalog']).exists()
 path=A/(prefix+'.json');path.write_text(asn.dump(format='json')[1]);lib=ModelLibrary(str(path));assert len(lib.group_names)==1
 manifest.append({'prefix':prefix,'minimum':minimum,'association':str(path),'members':[str(p) for p in files]})
(A/'manifest.json').write_text(json.dumps(manifest,indent=2));(A/'status').write_text('RUNNING\n');(A/'pid').write_text(str(os.getpid()))
for row in manifest:
 print('PAIR_BEGIN '+json.dumps(row),flush=True);os.environ['ASTROM_ABS_MINOBJ']=str(row['minimum']);run_tweakregstep(row['association']);print('PAIR_END '+row['prefix'],flush=True)
(A/'status').write_text('PROCESS_COMPLETE pending numerical and output QC\n')
