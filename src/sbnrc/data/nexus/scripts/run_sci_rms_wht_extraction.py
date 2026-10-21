#!/usr/bin/env python3
"""Extract SCI, ERR (RMS), and WHT images from JWST i2d mosaics."""

import argparse
from pathlib import Path

from astropy.io import fits


def extract(i2d: Path, output_dir: Path) -> None:
    with fits.open(i2d, memmap=True) as hdul:
        # JWST i2d ERR/WHT extension headers may omit WCS cards. Start from
        # SCI so every standalone map is sky-aware, then overlay the source
        # extension metadata (notably EXTNAME and BUNIT).
        science_header = hdul["SCI"].header.copy()
        for extension, suffix in (("SCI", "sci"), ("ERR", "rms"), ("WHT", "wht")):
            source = hdul[extension]
            header = science_header.copy()
            header.update(source.header)
            output = output_dir / i2d.name.replace("i2d.fits", f"{suffix}.fits")
            fits.PrimaryHDU(data=source.data, header=header).writeto(
                output, overwrite=True
            )
            print(f"WROTE {output}", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("directories", nargs="*", default=["."], type=Path)
    args = parser.parse_args()
    total = 0
    for directory in args.directories:
        directory.mkdir(parents=True, exist_ok=True)
        inputs = sorted(directory.glob("*i2d.fits"))
        if not inputs:
            raise RuntimeError(f"No i2d products found in {directory}")
        for i2d in inputs:
            extract(i2d, directory)
            total += 1
    print(f"COMPLETE: extracted {total} mosaics", flush=True)


if __name__ == "__main__":
    main()
