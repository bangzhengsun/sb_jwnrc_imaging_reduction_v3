# Release validation

Version 0.1.0 was checked locally with:

- Python parsing for every bundled Python source.
- Three unit tests covering non-square 30 mas WCS conventions, SHA-256 integrity
  of 80 preserved source/QC files, and refusal to overwrite an existing export.
- Wheel build using the working jwst_v3 environment without resolving new packages.
- Installation of the wheel into an isolated target directory, followed by
  `python -m sbnrc.cli verify` (80 files passed).
- Inspection of wheel contents for the pipeline scripts and field configuration.
- Credential-pattern scan of the publication file set (no matches).

The NEXUS field was reduced separately using JWST 3.0.0 and jwst_1596.pmap.
All nine 60 mas bands passed the retained image3 and delivery checks. No multi-day
science stage was rerun merely to package the record. The environment YAML and
optional dependencies were not tested via a fresh network install; the observed
working versions are separately captured in requirements-observed.txt.

The optional 30 mas launcher is a prepared recipe, not an already validated
30 mas science dataset. The CLI is an inspection/export interface, not an
unconditional automatic execution engine. Do not confuse packaging tests with
scientific end-to-end validation on a new field.
