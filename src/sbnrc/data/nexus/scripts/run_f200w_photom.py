#!/usr/bin/env python3
"""Run JWST PhotomStep on the approved F200W NMFwisp flat-field products."""

from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import os
import sys

FILTER = sys.argv[1].upper() if len(sys.argv) > 1 else "F200W"
BASE = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(os.environ.get("BASE", "/mnt/Data_Analysis_6/Data_Reduction/JWST_NRC_SMACS0723_v2026Aug"))
INPUT = BASE / FILTER / "stage2_medstep_wisprm"
OUTPUT = BASE / FILTER / "stage2_cal"


def photom_one(path_string: str) -> tuple[str, str]:
    from jwst.photom import PhotomStep

    path = Path(path_string)
    output = OUTPUT / path.name.replace("_flatfieldstep.fits", "_cal.fits")
    if output.exists():
        return output.name, "EXISTS"
    step = PhotomStep()
    step.output_dir = str(OUTPUT)
    step.save_results = True
    step.suffix = "cal"
    step.run(str(path))
    return output.name, "CREATED"


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    inputs = sorted(INPUT.glob("*_flatfieldstep.fits"))
    expected = len(inputs)
    if expected == 0:
        raise RuntimeError(f"No inputs found under {INPUT}")
    failures = []
    with ProcessPoolExecutor(max_workers=1) as pool:
        jobs = {pool.submit(photom_one, str(path)): path for path in inputs}
        for number, job in enumerate(as_completed(jobs), 1):
            path = jobs[job]
            try:
                name, state = job.result()
                print(f"[{number:03d}/{expected}] {state}: {name}", flush=True)
            except Exception as exc:
                failures.append((path.name, repr(exc)))
                print(f"[{number:03d}/{expected}] FAILED: {path.name}: {exc!r}", flush=True)
    if failures:
        raise RuntimeError(f"Photom failures: {failures}")
    outputs = list(OUTPUT.glob("*_cal.fits"))
    if len(outputs) != expected:
        raise RuntimeError(f"Expected {expected} cal outputs, found {len(outputs)}")


if __name__ == "__main__":
    main()
