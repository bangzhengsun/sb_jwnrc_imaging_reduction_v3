"""Inspect and export the versioned workflow without executing science at import."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import tomllib

DATA = Path(__file__).with_name('data')

def grid(config, resolution):
    m = config['mosaic']
    factor = {'60mas': 1, '30mas': 2}[resolution]
    shape = [factor * n for n in m['shape_60']]
    if len(shape) != 2 or any(n <= 0 or int(n) != n for n in shape):
        raise ValueError('shape_60 must contain two positive integer dimensions')
    return dict(shape=shape, pixel_scale_arcsec=m['pixel_scale_60_arcsec']/factor,
                crpix_step=[n/2-1 for n in shape],
                crval=[m['center_ra_deg'],m['center_dec_deg']], rotation=m['rotation_deg'])

def verify():
    manifest = json.loads((DATA/'manifest.json').read_text())
    for item in manifest['files']:
        path = DATA / item['path'].split('src/sbnrc/data/', 1)[1]
        if hashlib.sha256(path.read_bytes()).hexdigest() != item['sha256']:
            raise ValueError(f'Snapshot hash mismatch: {path}')
    return len(manifest['files'])

def main(argv=None):
    p=argparse.ArgumentParser(description='SBNRC JWST 3.0.0 workflow package and reduction record')
    sub=p.add_subparsers(dest='command', required=True)
    sub.add_parser('verify', help='Verify all bundled source and QC snapshot hashes')
    e=sub.add_parser('export', help='Export the full workflow record into a new directory')
    e.add_argument('destination',type=Path)
    g=sub.add_parser('grid',help='Print the configured WCS grid without running resampling')
    g.add_argument('config',type=Path);g.add_argument('--resolution',choices=['60mas','30mas'],default='60mas')
    a=p.parse_args(argv)
    if a.command=='verify':
        print(f'{verify()} bundled files verified')
    elif a.command=='export':
        verify()
        if a.destination.exists():
            p.error('Destination already exists; choose a new directory to preserve existing work')
        shutil.copytree(DATA,a.destination)
        print(f'Exported workflow record to {a.destination.resolve()}')
        print('Read README.md in the repository before adapting archived field scripts.')
    elif a.command=='grid':
        print(json.dumps(grid(tomllib.loads(a.config.read_text()),a.resolution),indent=2))

if __name__=='__main__':
    main()
