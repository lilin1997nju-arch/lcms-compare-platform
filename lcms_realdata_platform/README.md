# LC-MS Real Data Platform

This workspace contains the local LC-MS website built from the five user RAW files.

## Layout

- `data/mzML/`: converted mzML files generated from Thermo RAW data.
- `outputs/lcms_workbench.sqlite`: local API database for scans, chromatograms, spectra, heatmaps, and saved features.
- `outputs/lcms_xcalibur_workbench.html`: frontend only. Open it through the local server, not by double-clicking.
- `outputs_peak_first/lcms_peak_first_compare.sqlite`: separate Peak-first Compare V2 local API database.
- `outputs_peak_first/lcms_peak_first_compare.html`: separate Peak-first Compare V2 frontend.
- `convert_raw.ps1`: convert RAW files to mzML with ThermoRawFileParser.
- `build_workbench.ps1`: rebuild the SQLite database and frontend from mzML.
- `serve.ps1`: start the local website/API server.
- `build_peak_first_compare.ps1`: rebuild the Peak-first Compare V2 database and frontend from mzML.
- `serve_peak_first_compare.ps1`: start the Peak-first Compare V2 website/API server.

## Run

```powershell
.\lcms_realdata_platform\serve.ps1
```

Then open the printed localhost URL.

The current server can host several global RT x m/z comparison datasets on the
same port. Pass repeated `-Comparison "id|label|output_dir"` entries and switch
them from the page header.

## Run Peak-first Compare V2

This is a separate copy of the LC-MS compare workflow. It keeps the original
global RT x m/z heatmap workbench intact while adding the TIC-peak-first
screening flow.

```powershell
.\lcms_realdata_platform\build_peak_first_compare.ps1
.\lcms_realdata_platform\serve_peak_first_compare.ps1 -Port 8768
```

Open `http://127.0.0.1:8768/`.

The default V2 build analyzes up to 60 confirmed TIC peaks, sorted by retention
time in the table. Increase or reduce the scope with `-TopNPeaks` after checking
the first result.

For the current samples, keep the same comparison method on one port and switch
datasets from the page header:

```powershell
.\lcms_realdata_platform\build_peak_first_compare.ps1 -Group PTM -OutputDir .\outputs_peak_first_ptm
.\lcms_realdata_platform\build_peak_first_compare.ps1 -Group SVA -OutputDir .\outputs_peak_first_sva
.\lcms_realdata_platform\build_peak_first_compare.ps1 -InputDir .\data\mzML_new_ptm_20240911 -OutputDir .\outputs_peak_first_newptm_20240911

.\lcms_realdata_platform\serve_peak_first_compare.ps1 -OutputDir .\outputs_peak_first_ptm -Port 8768 -Comparison @(
  "ptm_old|PTM_old|$PWD\lcms_realdata_platform\outputs_peak_first_ptm",
  "sva_old|SVA_old|$PWD\lcms_realdata_platform\outputs_peak_first_sva",
  "ptm_20240911|PTM_20240911_boyoupei_vs_mailishu|$PWD\lcms_realdata_platform\outputs_peak_first_newptm_20240911"
)
```

Use `-IncludeSample sample_id` repeatedly for a custom 2- or 3-sample compare.

The 20240911 RAW files were converted into `data/mzML_new_ptm_20240911/` and
also built for the global RT x m/z workbench in `outputs_newptm_20240911/`.

## Rebuild From RAW

```powershell
.\lcms_realdata_platform\convert_raw.ps1
.\lcms_realdata_platform\build_workbench.ps1
.\lcms_realdata_platform\serve.ps1
```

The frontend calls `/api/*` endpoints from the local server. Opening the HTML file directly will show `Failed to fetch`.
