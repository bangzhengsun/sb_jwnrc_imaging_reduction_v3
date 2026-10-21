import datetime, fcntl, json, os, re, sys, tomllib
from pathlib import Path
from astropy.io import fits
base=Path(os.environ['NEXUS_BASE']);band=sys.argv[1];os.chdir(base)
cfg=tomllib.loads((base/'nexus.toml').read_text());auto=base/'automation'/band
lock=(auto/'image3_driver.lock').open('a')
fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
status=auto/'image3.status'
def event(state,**details):
    row=dict(timestamp=datetime.datetime.now(datetime.timezone.utc).isoformat(),filter=band,stage='image3_disk_retry',state=state,**details)
    with (base/'automation/reduction_log.jsonl').open('a') as stream:stream.write(json.dumps(row)+'\n')
try:
    assert (auto/'astrometry_qc.status').read_text().startswith('PASS')
    assert (auto/'oneoverf.status').read_text().startswith('COMPLETE')
    expected={p.name.removesuffix('_uncal.fits') for p in (base/band/'stage1').glob('*_uncal.fits')}
    suffix='_cal_tweakreg_1overf.fits'
    inputs=sorted((base/band/'stage2_cal_aligned').glob('*'+suffix))
    assert len(expected)==cfg['expected_exposures'][band]
    assert {p.name.removesuffix(suffix) for p in inputs}==expected
    output=base/band/'stage3_cal'
    if output.exists() and any(output.glob('*.fits')):
        raise RuntimeError(f'Existing FITS outputs require review: {output}')
    print(f'Validating {len(inputs)} inputs for {band}',flush=True)
    for path in inputs:
        with fits.open(path,memmap=True,do_not_scale_image_data=True) as hdul:
            hdul.verify('exception')
            assert hdul[0].header.get('S_TWKREG')=='COMPLETE'
            for hdu in hdul:
                info=hdu.fileinfo()
                assert info['datLoc']+info['datSpan']<=path.stat().st_size
    script=base/'SBNRC_07_image3.py';code=script.read_text()
    code,count=re.subn(r'(?m)^image3\.in_memory\s*=.*$','image3.in_memory = False',code)
    assert count==1,'Pipeline script changed; inspect memory settings'
    assert 'image3.outlier_detection.in_memory = image3.in_memory' in code
    assert 'image3.resample.in_memory = image3.in_memory' in code
    (auto/'image3.pid').write_text(str(os.getpid())+'\n')
    status.write_text(f'RUNNING disk-backed expected_inputs={len(inputs)}\n')
    event('RUNNING',inputs=len(inputs),in_memory=False,scratch=os.environ['TMPDIR'])
    sys.argv=[str(script),band,'cal_tweakreg_1overf']
    exec(compile(code,str(script),'exec'),{'__name__':'__main__','__file__':str(script)})
    status.write_text('PROCESS_COMPLETE pending output/WCS QC\n');event('PROCESS_COMPLETE')
except Exception as exc:
    status.write_text(f'FAILED disk-backed retry: {exc!r}\n');event('FAILED',error=repr(exc));raise
