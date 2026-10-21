#!/usr/bin/env python
"""Resample retained NEXUS CRFs only; never perform outlier rejection."""
import datetime, fcntl, json, logging, os, sys, tomllib
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from jwst.resample import ResampleStep
from jwst.associations.asn_from_list import asn_from_list
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base

B = Path(__file__).resolve().parent
C = tomllib.loads((B / 'nexus.toml').read_text())
assert B == Path(C['field']['base_path'])
os.chdir(B)
filters = C['filters']['lw'] + C['filters']['sw']
assert set(filters) == set(C['expected_exposures'])
for band in filters:
    assert (B / 'automation' / band / 'delivery.status').read_text().startswith('COMPLETE'), f'{band}: 60mas delivery incomplete'
A = B / 'automation/resample_30mas'
A.mkdir(exist_ok=True)
lock = (A / 'driver.lock').open('a')
fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
(A / 'driver.pid').write_text(str(os.getpid()))
logging.basicConfig(level=logging.INFO)
m = C['mosaic']
shape = [2 * n for n in m['shape_60']]
scale = m['pixel_scale_60_arcsec'] / 2
crpix = [n / 2 - 1 for n in shape]
assert shape == m['shape_30']
assert np.allclose(crpix, m['crpix_step_30'], rtol=0, atol=1e-7)
assert np.isclose(scale, m['pixel_scale_30_arcsec'], rtol=0, atol=1e-12)

def event(band, state, **extra):
    with (B / 'automation/reduction_log.jsonl').open('a') as h:
        h.write(json.dumps(dict(timestamp=datetime.datetime.now(datetime.timezone.utc).isoformat(), filter=band, stage='resample_only_30mas', state=state, **extra)) + '\n')

def validate(path, band, expected):
    with fits.open(path, memmap=True) as h:
        h.verify('exception')
        assert h['SCI'].shape == tuple(shape[::-1])
        assert h[0].header['S_RESAMP'] == 'COMPLETE'
        for ext in h:
            info = ext.fileinfo()
            assert info['datLoc'] + info['datSpan'] <= path.stat().st_size
        names = [str(n).split('_cal_tweakreg')[0] for n in h['HDRTAB'].data['FILENAME']]
        assert len(names) == len(expected) and set(names) == expected
        w = WCS(h['SCI'].header).celestial
        ref = WCS(fits.getheader(B / 'collection/60mas' / f"{C['field']['product_prefix']}_{band}_60mas_sci.fits")).celestial
        assert np.allclose(w.wcs.crval, ref.wcs.crval, rtol=0, atol=1e-10)
        assert np.allclose(w.wcs.crpix, 2 * ref.wcs.crpix, rtol=0, atol=1e-7)
        assert np.allclose(w.pixel_scale_matrix, ref.pixel_scale_matrix / 2, rtol=0, atol=1e-14)
        valid = 0
        for y in range(0, shape[1], 100):
            wt = h['WHT'].data[y:y+100]
            mask = np.isfinite(wt) & (wt > 0)
            sci, err = h['SCI'].data[y:y+100], h['ERR'].data[y:y+100]
            assert np.isfinite(sci[mask]).all() and np.isfinite(err[mask]).all() and (err[mask] > 0).all()
            valid += int(mask.sum())
        assert valid > 0

out = B / 'collection/30mas'
out.mkdir(exist_ok=True)
for band in filters:
    status = A / f'{band}.status'
    try:
        expected = {p.name.removesuffix('_uncal.fits') for p in (B / band / 'stage1').glob('*_uncal.fits')}
        assert len(expected) == C['expected_exposures'][band]
        name = f"{C['field']['product_prefix']}_{band}_30mas"
        target = out / f'{name}_i2d.fits'
        if target.exists():
            validate(target, band, expected)
            status.write_text('COMPLETE validated existing mosaic\n')
            event(band, 'SKIPPED', reason='Existing mosaic passed full QC')
            continue
        files = sorted((B / band / 'stage3_cal').glob('*_crf.fits'))
        assert len(files) == len(expected)
        assert {p.name.split('_cal_tweakreg')[0] for p in files} == expected
        for p in files:
            with fits.open(p, memmap=True) as h:
                h.verify('exception')
                assert h[0].header['S_OUTLIR'] == 'COMPLETE'
                for ext in h:
                    info = ext.fileinfo()
                    assert info['datLoc'] + info['datSpan'] <= p.stat().st_size
        asn = asn_from_list([str(p) for p in files], rule=DMS_Level3_Base, product_name=name)
        ap = A / f'{band}_asn.json'
        ap.write_text(asn.dump(format='json')[1])
        step = ResampleStep()
        step.save_results = True
        step.output_dir = str(out)
        step.output_file = name
        step.suffix = 'i2d'
        step.in_memory = True
        step.kernel, step.pixfrac = m['kernel'], m['pixfrac']
        step.rotation = m['rotation_deg']
        step.pixel_scale, step.crpix = scale, crpix
        step.crval = [m['center_ra_deg'], m['center_dec_deg']]
        step.output_shape = shape
        status.write_text('RUNNING disk-backed ResampleStep only\n')
        event(band, 'RUNNING', exposures=len(files), shape=shape, pixel_scale=scale, crpix_step=crpix, in_memory=False)
        print(f'Starting {band}: {len(files)} CRFs -> {shape}, {scale} arcsec/pixel', flush=True)
        step.run(str(ap))
        validate(target, band, expected)
        status.write_text('COMPLETE mosaic and WCS QC PASS\n')
        event(band, 'COMPLETE', output=str(target))
    except Exception as exc:
        status.write_text(f'FAILED {exc!r}\n')
        event(band, 'FAILED', error=repr(exc))
        raise
print('All nine 30mas i2d mosaics complete and validated.', flush=True)
