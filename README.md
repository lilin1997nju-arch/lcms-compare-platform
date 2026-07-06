# LC-MS Compare Platform

This repository contains the LC-MS data reading and comparison tools developed from the Empower chromatogram workflow.

The current focus is a Peak-first LC-MS comparison workflow:

- Convert vendor/raw data to mzML.
- Build a local SQLite API database from scan-level LC-MS data.
- Align chromatograms and compare selected 2- or 3-sample datasets.
- Detect TIC peaks, compare TIC-peak-level summed spectra, extract XICs, and align LCMSFeatures with local RT correction.
- Review low-similarity feature groups through tables, spectra, XIC views, and feature-level RT-m/z heatmaps.

## Repository Layout

- `lcms_feature_mvp/`: Python core algorithms, parsers, local API servers, and tests.
- `lcms_realdata_platform/`: user-facing PowerShell wrappers for converting, building, and serving real-data LC-MS comparison projects.
- `lcms_realdata_platform/outputs*/`: lightweight HTML frontends are kept as examples; SQLite databases and generated matrices are ignored.
- Data directories are kept with `.gitkeep` placeholders only. Raw files, mzML files, and generated databases should stay local.

## Run Peak-first Compare

Place mzML files under `lcms_realdata_platform/data/mzML/`, then run from the repository root:

```powershell
.\lcms_realdata_platform\build_peak_first_compare.ps1 -Group PTM -OutputDir .\outputs_peak_first_ptm
.\lcms_realdata_platform\build_peak_first_compare.ps1 -Group SVA -OutputDir .\outputs_peak_first_sva
.\lcms_realdata_platform\serve_peak_first_compare.ps1 -OutputDir .\outputs_peak_first_ptm -Port 8768
```

For custom 2- or 3-sample comparison, pass `-IncludeSample sample_id` repeatedly to `build_peak_first_compare.ps1`.

## Convert RAW To mzML

Thermo RAW conversion requires ThermoRawFileParser or an equivalent local converter:

```powershell
.\lcms_realdata_platform\convert_raw.ps1
```

Generated `.raw`, `.mzML`, `.sqlite`, `.csv`, and `.json` result files are intentionally excluded from Git.
