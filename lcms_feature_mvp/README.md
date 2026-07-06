# LC-MS RT-m/z Feature Comparison MVP

This directory is an isolated LC-MS prototype. It does not modify the existing Empower chromatogram app.

## Data

Downloaded source files:

- `data/raw/zenodo_5005513/SCX-HPLC-MS_Intact_MabThera_Deamidated.raw`
- `data/raw/zenodo_5005513/SCX-HPLC-MS_Intact_MabThera_Glycated.raw`
- `data/raw/zenodo_5005513/SCX-HPLC-MS_Intact_MabThera_Untreated_1.raw`
- `data/raw/zenodo_5005513/SCX-HPLC-MS_Intact_MabThera_Untreated_2.raw`
- `data/raw/zenodo_5005513/SCX-HPLC-MS_Intact_MabThera_Untreated_3.raw`
- `data/raw/zenodo_5005513/SCX-HPLC-MS_Intact_Reditux_Deamidated.raw`
- `data/raw/zenodo_5005513/SCX-HPLC-MS_Intact_Reditux_Glycated.raw`
- `data/raw/zenodo_5005513/SCX-HPLC-MS_Intact_Reditux_Untreated_1.raw`
- `data/raw/zenodo_5005513/SCX-HPLC-MS_Intact_Reditux_Untreated_2.raw`
- `data/raw/zenodo_5005513/SCX-HPLC-MS_Intact_Reditux_Untreated_3.raw`

Source dataset: https://zenodo.org/records/5005513

ThermoRawFileParser has been installed locally under:

```text
../.local-tools/ThermoRawFileParser/current/ThermoRawFileParser.exe
```

Converted mzML files are stored under:

```text
data/converted/mzML
```

The parser can now read `.mzML` directly. If the input is still `.raw`, the MVP keeps the original RAW as a source artifact and falls back to deterministic mock centroid scans marked `mock_from_vendor_raw`.

## Convert RAW to mzML

Download several additional lightweight SCX-HPLC-MS intact RAW files from the same Zenodo record. The script supports resume and verifies file size plus MD5 checksum:

```powershell
powershell.exe -ExecutionPolicy Bypass -File .\lcms_feature_mvp\download_additional_rawdata.ps1
```

```powershell
powershell.exe -ExecutionPolicy Bypass -File .\lcms_feature_mvp\convert_raw_to_mzml.ps1
```

## Run

```powershell
python .\lcms_feature_mvp\run_lcms_mvp.py
```

## Xcalibur-style workbench

Generate the new browser-based review page requested for LC-MS raw-data inspection:

```powershell
python .\lcms_feature_mvp\run_xcalibur_workbench.py --input-dir .\lcms_feature_mvp\data\converted\mzML --output-dir .\lcms_feature_mvp\outputs_workbench --rt-min 0 --rt-max 60 --mz-min 2500 --mz-max 8000
```

The generator writes the heavy workbench payload to SQLite instead of embedding it in the HTML. Start the local server:

```powershell
python .\lcms_feature_mvp\serve_xcalibur_workbench.py --output-dir .\lcms_feature_mvp\outputs_workbench
```

Open the printed localhost URL, usually:

```text
http://127.0.0.1:8765/
```

Implemented in this workbench:

- Top pane chromatogram browser, similar to Xcalibur.
- Bottom pane MS spectrum browser.
- Single-RT full scan spectrum by default.
- Summed spectrum over a selected RT range.
- Apex spectrum from the selected RT range.
- Average spectrum over the selected RT range.
- Box zoom for chromatogram and spectrum canvases.
- Multi-sample selection.
- Main-peak RT alignment before comparison. The reference sample first defines the main peak from the BPC trace by smoothing, baseline estimation, peak-boundary search, and area ranking; other samples then match the corresponding peak within `--alignment-match-window-min` before calculating `reference_rt - sample_rt`. By default the first 2 min are ignored to avoid injection/front artifacts.
- RT x m/z binned similarity/difference heatmap after alignment, with selectable `max`, `tic`, `none`, and `median` intensity normalization.
- Pairwise heatmaps for reference-vs-sample comparisons plus an all-samples CV heatmap. The cohort heatmap uses `similarity = 1 / (1 + CV)` at each aligned RT-m/z bin, so one view can quickly expose differences across all selected samples.
- Box zoom, reset, and undo for the RT x m/z heatmap.
- Local heatmap windows can be requested from the SQLite-backed API with `/api/heatmap-window?key=...&method=...&rt_min=...&rt_max=...&mz_min=...&mz_max=...`. Heatmap box zoom now calls this endpoint and draws the returned local RT-m/z window; reset or comparison changes return to the full cached heatmap. The first implementation slices cached matrices and keeps intensity grids plus low-similarity points aligned, providing the API shape needed for later zoom-triggered high-resolution recalculation.
- Top difference region table with rank, aligned RT, m/z, difference score, similarity score, max raw intensity, sample presence, normalized mean intensity, fold change, and difference type.
- Low-similarity bin table rows can be clicked to jump to the corresponding XIC and spectrum region.
- The selected RT-m/z intensity table now shows per-sample nearest-scan `scan_id`, raw RT, aligned RT, RT delta, nearest-scan raw intensity, heatmap-bin raw intensity, normalized intensity, relative percentages, and present/missing status. The selected intensity table can be exported as CSV for traceable drill-down review.
- Selected RT-m/z regions can be saved as LCMSFeature regions into the local SQLite `saved_lcms_features` table and exported as CSV for later feature matrix/database integration. Saved regions now preserve per-sample raw and normalized intensities plus source scan trace metadata for both pairwise and cohort heatmaps, and they are restored when the workbench is reopened.
- Top difference regions from the current heatmap can be batch-saved as LCMSFeature regions. A feature matrix table is built from saved regions with one feature per row and per-sample raw intensity, normalized intensity, and present/missing status columns. The matrix can be exported as `lcms_feature_matrix.csv`.
- The saved-feature matrix is also available directly from the local SQLite-backed API as `/api/feature-matrix` and `/api/feature-matrix.csv`, including per-sample raw intensity, normalized intensity, present/missing status, and source scan trace columns.
- The workbench generator stores scan/heatmap payloads in `lcms_workbench.sqlite`; the HTML is now a thin frontend. It first loads `/api/bootstrap`, then fetches selected sample spectra through `/api/spectra?sample=...`, heatmaps through `/api/heatmap?key=...&method=...`, saved LCMSFeature regions through `/api/features`, and DB-backed feature matrices through `/api/feature-matrix.csv`. The legacy `/api/payload` endpoint is kept for compatibility and is assembled from the split local artifacts only when requested. This avoids packaging all LC-MS data into the HTML as sample count grows.
- The workbench generator also writes automatic backend artifacts from the default-normalization Top difference regions: `lcms_auto_difference_regions.json`, `lcms_auto_difference_regions.csv`, and `lcms_auto_feature_matrix.csv`. These files are meant as the first bridge toward an impurity/feature database.

Run against the real converted mzML files:

```powershell
python .\lcms_feature_mvp\run_lcms_mvp.py --input-dir .\lcms_feature_mvp\data\converted\mzML --output-dir .\lcms_feature_mvp\outputs_real_mzml --mz-min 2500 --mz-max 8000 --rt-start 0 --rt-end 60 --intensity-threshold 1500 --min-scan-count 5 --top-n-mz 30 --min-peak-height 1500 --min-peak-area 10 --min-snr 2
```

Outputs:

- `outputs/lcms_mvp_report.html`
- `outputs/lcms_mvp_results.json`
- `outputs/feature_matrix.csv`
- `outputs/candidates_<sample>.csv`
- `outputs/mock_centroid_scans.csv`
- `outputs_real_mzml/lcms_mvp_report.html`
- `outputs_real_mzml/feature_matrix.csv`
- `outputs_workbench/lcms_xcalibur_workbench.html`
- `outputs_workbench/lcms_workbench.sqlite`
- `outputs_workbench/lcms_auto_difference_regions.json`
- `outputs_workbench/lcms_auto_difference_regions.csv`
- `outputs_workbench/lcms_auto_feature_matrix.csv`

## Implemented MVP

- LC-MS raw file metadata structure.
- Scan-level `RT -> mz_array + intensity_array` model.
- TIC/BPC calculation from scans.
- RT-window candidate m/z screening.
- XIC extraction.
- XIC peak detection and trapezoid integration.
- Cross-sample feature matching by m/z and RT tolerance.
- Feature matrix and difference type classification.
- Standalone HTML report with TIC, feature table, XIC overlay, and matrix.

## Next converter step

To parse real Thermo `.raw` scans, add one of these before the current parser:

- ProteoWizard `msconvert` to produce `.mzML`.
- ThermoRawFileParser to produce `.mzML` or `.mgf`.
- A validated vendor SDK conversion service under the company Empower/Waters/Thermo environment.

After conversion, the existing pipeline can consume centroid CSV immediately, and an mzML reader can be added behind `core/lcms_parser.py`.
