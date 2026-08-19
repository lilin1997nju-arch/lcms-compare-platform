# UI Notes

The current MVP writes a standalone report at `../outputs/lcms_mvp_report.html`.

The report already contains the four requested regions:

- TIC/BPC overview.
- Candidate/feature list.
- XIC overlay viewer.
- Feature matrix.

When integrating with the existing chromatogram web app later, split the report into components equivalent to:

- `LCMSOverview`
- `LCMSFeatureTable`
- `LCMSXICViewer`
- `LCMSFeatureMatrix`
- `LCMSParameterPanel`

