#!/usr/bin/env python3
"""Prepare median-background-subtracted F200W flat-field products for NMFwisp."""

from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import os
import sys

FILTER = sys.argv[1].upper() if len(sys.argv) > 1 else "F200W"
BASE = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(os.environ.get("BASE", "/mnt/Data_Analysis_6/Data_Reduction/JWST_NRC_SMACS0723_v2026Aug"))
INPUT = BASE / FILTER / "stage2_medstep_bkgsubtract"
OUTPUT = BASE / FILTER / "stage2_medstep_wisprm"


def process_one(rate_path_str: str) -> tuple[str, str]:
    from jwst.assign_wcs import AssignWcsStep
    from jwst.flatfield import FlatFieldStep

    rate_path = Path(rate_path_str)
    assign_name = rate_path.name.replace("_rate.fits", "_assignwcsstep.fits")
    flat_name = rate_path.name.replace("_rate.fits", "_flatfieldstep.fits")
    flat_path = OUTPUT / flat_name
    if flat_path.exists():
        return rate_path.name, "EXISTS"

    assign = AssignWcsStep()
    assign.output_dir = str(OUTPUT)
    assign.save_results = True
    assign.suffix = "assignwcsstep"
    assign.run(str(rate_path))

    flat = FlatFieldStep()
    flat.output_dir = str(OUTPUT)
    flat.save_results = True
    flat.suffix = "flatfieldstep"
    flat.run(str(OUTPUT / assign_name))
    return rate_path.name, "CREATED"


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    inputs = sorted(INPUT.glob("*_rate.fits"))
    if not inputs:
        raise RuntimeError(f"No rate inputs found under {INPUT}")
    print(f"Preparing {len(inputs)} {FILTER} inputs with 1 worker", flush=True)
    failures = []
    with ProcessPoolExecutor(max_workers=1) as pool:
        futures = {pool.submit(process_one, str(path)): path for path in inputs}
        for index, future in enumerate(as_completed(futures), 1):
            path = futures[future]
            try:
                name, state = future.result()
                print(f"[{index:02d}/{len(inputs)}] {state}: {name}", flush=True)
            except Exception as exc:
                failures.append((path.name, repr(exc)))
                print(f"[{index:02d}/{len(inputs)}] FAILED: {path.name}: {exc!r}", flush=True)
    if failures:
        raise RuntimeError(f"{len(failures)} preparation failures: {failures}")


if __name__ == "__main__":
    main()
