# LC-MS RT-m/z Feature Comparison MVP

This directory is an isolated LC-MS prototype. It does not modify the existing Empower chromatogram app.

## Vedolizumab peptide-map LC-MS/MS

The current project extends the Peak-first MS1 comparison with a standard-library MS/MS search:

```text
Peak-first differential Feature
-> MS2 precursor/RT link
-> verified FASTA tryptic peptide
-> b/y fragment match
-> common PTM candidate
-> target-decoy q-value
-> JSON/CSV/HTML evidence report
```

Peak-first MS1 feature areas are normalized per run to the median total-TIC area before difference typing. Quantitative presence is determined from the normalized area ratio, while S/N-based detection is retained as a separate confidence field. With the defaults, `common_feature` is below 2-fold, `moderate_difference` is 2-fold to below 4-fold, `area_changed` is at least 4-fold, and only an at-most-2% weak side is `presence_absence`. If comparable areas are present but one side misses the S/N threshold, the fold-based type is retained and ranking is multiplied by 0.35. The report's fold column is the true symmetric normalized-area ratio rather than a ratio of log-transformed values. These settings are configurable with `--feature-presence-relative-area-fraction`, `--feature-common-fold-change-threshold`, `--feature-strong-fold-change-threshold`, and `--feature-partial-detection-rank-weight`.

The RT-m/z feature heatmap defaults to a selected reference/test pair. Red means the test sample is higher and blue means it is lower; color depth uses the capped signed `ln(test/reference)`. Square area is independent of fold change and represents the higher normalized abundance in the selected sample pair. Abundance is transformed with `ln(1 + area)` and scaled between the 5th and 95th percentiles, preventing a few extreme signals from making all other Features invisible. Reference and test can be selected independently (or swapped), so the same cohort result supports pairwise views when more than two samples are loaded. The optional fold-magnitude mode retains the unsigned `ln(fold) / ln(max fold)` color palette while keeping abundance-driven square sizes. Both the local and global Feature tables show the pairwise type, signed direction, and sample-colored abundance-order dots; clicking any sortable header toggles ascending/descending order.

`core/lcms_parser.py` keeps MS1-only behavior by default and can now read selected MS levels plus precursor m/z, charge, isolation window, activation method, and collision energy. `core/lcms_msms.py` implements tryptic digestion, configurable carbamidomethyl Cys, one variable PTM per peptide, target-decoy competition, fragment scoring, and Peak-first Feature linking. The linker now maps direct precursors, charge-state envelopes, and `1.003355 / z` isotope envelopes separately, so isotope/charge peaks are not mislabeled as independent modifications. An MS1-first evidence table reports every significant Feature as direct identification, isotope/charge inference, selected-but-unidentified precursor, co-isolated unresolved MS2, or no acquired MS2. Use `--no-carbamidomethyl-cys` when the confirmed sample preparation did not alkylate Cys.

The antibody-focused candidate catalog covers Met/Trp oxidation and dioxidation, Asn/Gln deamidation, Asn succinimide, N-terminal pyroglutamate, Lys/N-terminal glycation, Fc-sequon G0F/G1F/G2F glycopeptides, and heavy-chain C-terminal Lys clipping. Isolation-window-guided re-search is permitted only for unresolved significant MS1 Features and is capped at tentative confidence because co-isolation can produce chimeric spectra. The generated `lcms_ms1_ms2_feature_evidence.csv` is the traceable Feature-first export.

Unidentified MS1 component groups can also be searched jointly. Candidate peptides are first restricted by the component neutral mass (20 ppm in the current workflow), then b/y evidence is combined across selected precursor scans assigned to different charge-state or isotope members. Isolation-window-only scans are not used. A joint result requires target-decoy q <= 1%, score >= 45, at least five unique matched ions, at least 15% fragment coverage, and support from at least two scans and two member Features. It is reported as `tentative_component_consensus_identification`, never as a direct identification. Features already supported at B/C confidence are skipped; D-level single-spectrum links remain eligible so independent group members can rescue an otherwise insufficient spectrum.

Components still unexplained by the strict neutral-mass search enter a separate open sequence-tag search. Its sequence catalog includes fully tryptic peptides with up to four missed cleavages plus explicit one-direction N- or C-terminal truncations of up to 20 residues. A component mass offset from -250 to +2500 Da is allowed, and regular/offset b/y ions are combined across genuinely selected charge-state or isotope-member spectra. Acceptance requires q <= 1%, score >= 45, at least five ions, at least 6% coverage, a sequence tag of at least two residues, two scans, two component members, and a score margin of at least three over competing sequences/sites. A peptide shorter than eight residues with a large unexplained mass offset cannot be promoted unless the offset matches a known modification; an explicit truncation with near-zero residual mass is also allowed. Accepted results are reported separately as `backbone_sequence_supported` or `truncation_sequence_supported`. Structurally useful but non-unique evidence remains `sequence_region_candidate` at D level and is displayed as a candidate rather than an identification.

For every peptide-form family containing at least one significant MS1 Feature, the MS/MS workflow now recalls all B/C-level identified members, including members classified as `common_feature`, and pairs modified or terminally processed forms with their unmodified counterpart. It extracts the full-run MS1 XIC from the Peak-first SQLite spectra, sums M/M+1/M+2 isotope signals, merges supported charge states into one peptide-form result, and reports normalized form composition plus pairwise percentage-point differences. This family trigger detects patterns such as differential retained heavy-chain C-terminal Lys with a stable clipped counterpart. Consistent direction across at least two comparable charge states is recorded as supporting evidence; charge states and isotope peaks are not ranked as independent differential components. Isobaric site alternatives such as `Oxidation@1 / Oxidation@12` remain grouped when the current fragments cannot localize the site. Traceable exports are `lcms_modification_pairs.csv` and `lcms_modification_level_quantitation.csv`. B-level complete families enter formal relative quantitation, C-level families remain tentative, and D-level or mass-only hypotheses are excluded. Reported proportions are semi-quantitative relative MS-response composition, not response-factor-corrected absolute occupancy.

Run after supplying a user-verified heavy/light-chain FASTA:

```powershell
python .\lcms_feature_mvp\run_msms_compare.py `
  --mzml .\data\mzML\20260407_QL2519_20260316_T5-3G6-ProA_Trypsin_PTM.mzML `
  --mzml .\data\mzML\20260407_Vedolizumab_12756432_Trypsin_PTM.mzML `
  --fasta .\data\sequences\vedolizumab_verified.fasta `
  --feature-sqlite .\outputs\vedolizumab_peak_first\lcms_peak_first_compare.sqlite `
  --output-dir .\outputs\vedolizumab_peak_first
```

The report never assigns confirmation level A from MS/MS alone. A verified product sequence is required; the program does not download or silently substitute a public antibody sequence.

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

## Current boundary

Thermo RAW conversion and MS1/MS2 mzML parsing are operational. Real sequence-level identification still requires the verified product FASTA and confirmed sample-preparation settings. The search intentionally limits candidates to one variable PTM per peptide; expand this only for unexplained features because combinatorial PTM search rapidly enlarges the search space. A Feature covered only by an isolation window is not considered identified, and MS1-only mass differences are never promoted to a named PTM without qualifying fragment evidence.

## Sequence and 3D structure mapping

The Peak-first report now includes a final “差异组分的序列与三级结构定位” module:

- All differential Features with a main or candidate peptide sequence are deduplicated into peptide regions and rendered on the complete HC and LC FASTA sequences. Red means the selected test sample is higher than the reference, blue means lower, purple means overlapping evidence has conflicting directions, and semi-transparent color means tentative evidence.
- Clicking a colored residue selects the linked Feature, scrolls the Feature table to its group, refreshes the XIC/MS2 detail, and focuses the selected region in 3D.
- Selecting a Feature reads its linked MS2 peptide candidate and highlights the corresponding 1-based residue interval in the product FASTA chain.
- A local PDB or mmCIF file can be uploaded directly in the browser.
- A PDB ID can be fetched through the local server from RCSB and cached under `data/structures/pdb_cache`.
- Structure mapping first assigns the full product HC/LC sequence to the best matching PDB antibody chains and then performs local peptide mapping only inside those chains. This avoids short peptides being assigned to an antigen or receptor chain in antibody-complex structures.
- All mappable differential regions are colored in the structure using the same red/blue/purple direction scheme. The selected Feature is orange, and a localized `@position` modification is additionally shown as an orange sphere/stick site.
- Vedolizumab currently auto-loads PDB `3V4P` (parent ACT-1 Fab bound to α4β7), and denosumab auto-loads PDB `5I1C` (a close human Fab variable-region template). RCSB has no experimental entry with a ≥99% exact match to either drug's HC/LC variable regions, so both defaults are explicitly labelled as homologous templates rather than exact drug structures.
- Structure mapping never upgrades the original MS2 identification confidence. Missing coordinates, construct truncation, sequence variants, and ambiguous antibody chains can all prevent a reliable structural match.

The bundled viewer is served from `lcms_feature_mvp/ui/vendor/3Dmol-min.js`, so uploaded structures still work without a CDN. Online PDB retrieval requires outbound HTTPS access to `files.rcsb.org`.
