from pathlib import Path
import json,logging,os,datetime
from jwst.tweakreg import TweakRegStep
from jwst.associations.asn_from_list import asn_from_list
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base
from jwst.datamodels import ModelLibrary
b=Path.cwd();logging.basicConfig(level=logging.INFO,format='%(levelname)s: %(message)s');s=(b/'SBNRC_05_astrometry.py').read_text();exec(s[s.index('def run_tweakregstep(f):'):s.index('def execute_tweakreg_on_file')]);abs_refcat=str(b/'nexus_f444w_refcat.ecsv');os.environ['ASTROM_ABS_MINOBJ']='10'
for filt in ['F277W','F360M']:
 a=b/'automation'/filt/'paired_retry';a.mkdir(exist_ok=True);assert not (a/'status').exists();(a/'status').write_text('RUNNING\n');(a/'pid').write_text(str(os.getpid()));output_dir=str(b/filt/'stage2_cal_aligned_paired_retry');Path(output_dir).mkdir(exist_ok=True)
 failures=json.loads((a.parent/'astrometry_qc_failures.json').read_text());names=failures['failed_exposures'] if isinstance(failures,dict) else [r['exposure'] for r in failures];prefixes=sorted({n.rsplit('_',1)[0] for n in names})
 for prefix in prefixes:
  members=[str(b/filt/'stage2_cal'/f'{prefix}_{m}_cal.fits') for m in ['nrcalong','nrcblong']];assert all(Path(p).exists() for p in members);asn=asn_from_list(members,rule=DMS_Level3_Base,product_name=prefix+'_paired')
  for m in asn['products'][0]['members']:m['group_id']=prefix;m['tweakreg_catalog']=m['expname'].replace('_cal.fits','.phot.ecsv');assert Path(m['tweakreg_catalog']).exists()
  p=a/(prefix+'.json');p.write_text(asn.dump(format='json')[1]);assert len(ModelLibrary(str(p)).group_names)==1
  print('PAIR_BEGIN',filt,prefix,flush=True);run_tweakregstep(str(p));print('PAIR_END',filt,prefix,flush=True)
 (a/'status').write_text('PROCESS_COMPLETE pending QC\n')
