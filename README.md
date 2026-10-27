Disclaimer: Generative AI tool ChatGPT (GPT-6 Astra) was used to prepare this pipeline for GitHub repository. The pipeline itself was mostly human-written and parts involving GAI were all validated. 

# SBNRC: JWST NIRCam imaging reduction, pipeline 3.0.0

Pipeline scripts for reducing NIRCam imaging from uncalibrated exposures through
astrometric alignment, mosaics, SCI/RMS/WHT extraction, empirical HRMS maps,
and verified delivery archives. This version preserves the implementation used
for the NEXUS full-depth reduction with **JWST 3.0.0**, **CRDS jwst_1596.pmap**,
and **NMFwisp 1.1.4**.

## Installation

```bash
git clone https://github.com/jerryfinn6/sb_jwnrc_imaging_reduction_v3.git
cd sb_jwnrc_imaging_reduction_v3
conda env create -f environment.yml
conda activate jwst_v3
pip install --no-deps .
sbnrc verify
```

Alternatively, create a Python 3.12 environment and run `pip install '.[reduction]'`.
Install SExtractor separately and make sure its executable is available as `sex`.
`requirements-observed.txt` records the working environment; it is an inventory,
not a portable dependency lock. Science data, astrometric catalogs and CRDS
reference files must be supplied separately.

## Prepare a working directory

Export the pipeline into a new directory:

```bash
sbnrc export /path/to/pipeline-copy
mkdir -p /path/to/field/uncals
cp /path/to/pipeline-copy/nexus/scripts/* /path/to/field/
cd /path/to/field
export BASE="$PWD"
export Research_Data=/path/to/pipeline-copy/support
export CRDS_SERVER_URL=https://jwst-crds.stsci.edu
export CRDS_CONTEXT=jwst_1596.pmap
export CRDS_PATH=/path/to/crds_cache
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
```

**Before executing a science stage, adapt the scripts and `nexus.toml` to your
field.** This is the tested NEXUS implementation, not a fully portable one-command
runner. Some scripts contain absolute NEXUS paths or older fallback paths.
Editing only the TOML is insufficient. Review paths and product names with:

```bash
rg -n '/mnt/|/home/|NEXUS|SMACS0723' *.py *.sh nexus.toml
```

Set the base directory, filter list, exposure counts, reference catalog,
output prefix, WCS, and HRMS settings explicitly. Keep the configuration filename
`nexus.toml` unless you also change its references in the scripts. Adapt the shell
environment and launcher paths before using them. Do not inherit NEXUS-specific
fit exceptions or HRMS box coordinates for a different field.

Place only imaging `*_uncal.fits` in `uncals/`. Exclude `EXP_TYPE=NRC_WFSS` and
other non-imaging data; confirm FITS headers, unique basenames and expected counts.
The sorter uses the header FILTER, so inspect pupil-selected medium-band data
and confirm the resulting directories match the effective imaging filters.

## Run the calibration stages

Run commands from the field directory in the activated environment. First obtain
references and sort inputs:

```bash
crds sync --fetch-references --dataset-files uncals/*_uncal.fits
python SBNRC_00_sort_uncals.py
```

For each long-wavelength filter, for example F444W:

```bash
python SBNRC_01_stage1.py F444W
python SBNRC_02_stage2_LW.py F444W
python SBNRC_03_bkg-sub.py F444W
```

For each short-wavelength filter, for example F200W:

```bash
python SBNRC_01_stage1.py F200W
python SBNRC_02_stage2med_SW.py F200W
python prepare_f200w_nmfwisp_inputs.py F200W "$BASE"
python run_f200w_nmfwisp.py F200W "$BASE"
# Inspect NMF outputs/diagnostics and confirm complete, valid corrections.
python run_f200w_photom.py F200W "$BASE"
python run_sw_background_safe.py F200W "$BASE"
```

Despite their historical filenames, the NMF preparation and photometry helpers
accept the filter argument. Preserve original backups and verify output membership
before advancing. Do not repeat background subtraction over already modified data.

## Astrometry and 1/f correction

For each band:

```bash
python SBNRC_04_RunSextractor.py F444W
bash SBNRC_05_do_astrometry.sh F444W 1
# Review match counts and fitted RMS for every exposure before continuing.
python SBNRC_06_1overf.py F444W
```

Keep `sextractor.conf`, `outparam.param`, `gauss_2.0_5x5.conv` and `default.nnw` available
in the working directory. Stage 06 uses the bundled `image1overf.py` through
`Research_Data/python_modules`. The final 1/f correction is still required after
NMFwisp. Grouped-detector and local-reference recovery recipes are supplied under
`nexus/recovery/`; adapt their selection and acceptance criteria to your field.

## Build mosaics

Finalize a shared WCS from the aligned footprint or adopt a supplied reference
WCS. Set `geometry_finalized` only after checking coverage. For a centered grid,
JWST step CRPIX is `shape/2 - 1`. Inspect the configured grid with:

```bash
sbnrc grid nexus.toml --resolution 60mas
sbnrc grid nexus.toml --resolution 30mas
```

For the full 60 mas Image3 stage:

```bash
python SBNRC_07_image3.py F444W cal_tweakreg_1overf
```

Review the Image3 memory configuration before launching. The archived TOML retains
an earlier in-memory preference. The final NEXUS SW runs used disk-backed
processing instead: see `nexus/execution/run_band.py` and
`run_remaining_on_disk.sh`. Their filter selections, locks, paths, and QC
prerequisites are specific to that run. Do not launch overlapping jobs or assume
that a quiet log means a stalled process.

For optional 30 mas mosaics, use retained CRFs and **ResampleStep only**:

```bash
bash run_nexus_30mas.sh
```

Adapt the launcher's environment paths first. It requires every configured 60 mas
band's `delivery.status` to be COMPLETE, validates the CRFs, runs sequentially on
disk, and checks each output WCS. It doubles both dimensions and halves pixel
scale while preserving center and rotation. It does not repeat outlier detection.
The NEXUS example is 30000 square at 60 mas and 60000 square at 30 mas;
these dimensions are not defaults for other fields.

## Extract maps, calculate HRMS and make ZIP files

Validate the i2d FITS structure, CRF/HDRTAB membership, shared WCS, and finite
SCI/ERR values on positive-weight coverage. Place validated i2d mosaics in
`collection/60mas/` (and optionally `collection/30mas/`), then run:

```bash
python run_sci_rms_wht_extraction.py collection/60mas
python step1_select_regions.py --collection collection/60mas
# Inspect candidate and selected-region previews before calculating HRMS.
python step2_calc_hrms.py
python bulk_zipfiles.py collection/60mas/*_sci.fits collection/60mas/*_hrms.fits collection/60mas/*_wht.fits
```

For 30 mas extraction, pass `collection/30mas` to the extractor. Enable the desired
resolutions in the configuration before HRMS calculation. HRMS region selection
must avoid artifacts and use valid coverage. If common boxes are insufficient,
use `select_hrms_regions_per_filter.py` and record the per-band corners; the HRMS
calculator supports per-filter selections. Do not copy the archived NEXUS box
coordinates into a new field. Check each HRMS map for finite positive values on
valid coverage and verify each ZIP's member name, size, and CRC.

The `nexus/postprocessing/` scripts show the full validation and delivery checks
used for NEXUS; adapt their band names, expected counts and paths before reuse.
A successful process exit alone does not establish complete scientific products.

## Repository contents and validation

- `src/sbnrc/data/nexus/scripts/`: calibration, astrometry, mosaics, HRMS and delivery.
- `src/sbnrc/data/nexus/execution/`: disk-backed execution recipe.
- `src/sbnrc/data/nexus/recovery/`: field-specific astrometric recovery examples.
- `src/sbnrc/data/nexus/postprocessing/`: output validation and delivery recipes.
- `src/sbnrc/data/nexus/qc/`: NEXUS completion evidence and region selections.
- `src/sbnrc/data/support/`: local 1/f helper used in this reduction.
- `src/sbnrc/data/manifest.json`: source locations and SHA-256 provenance.

All nine NEXUS 60 mas bands passed mosaic and delivery QC (6624 exposures and
27 CRC-checked archives). The 30 mas launcher was prepared separately; that is
not evidence of completed 30 mas products. No science FITS, source catalogs,
calibration cache, or credentials are included.

```bash
python -m unittest discover -s tests -v
python -m pip wheel --no-deps . -w dist
```

Package tests cover grid conventions, source integrity, and safe export behavior;
they do not replace end-to-end science validation. See `docs/VALIDATION.md`.

## Attribution

SBNRC workflow: Bangzheng Sun. External dependencies retain their authorship and
licenses. The supplied local `image1overf.py` has no explicit license header;
this repository does not assign a new license to third-party material.
