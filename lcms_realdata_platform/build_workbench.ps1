param(
    [string]$InputDir = ".\data\mzML",
    [string]$OutputDir = ".\outputs",
    [double]$RtMin = 0,
    [double]$RtMax = 120,
    [double]$MzMin = 250,
    [double]$MzMax = 2000,
    [int]$RtBinCount = 220,
    [int]$MzBinCount = 180
)

$ErrorActionPreference = "Stop"
$Root = Split-Path -Parent $MyInvocation.MyCommand.Path
$Workspace = Split-Path -Parent $Root
$Runner = Join-Path $Workspace "lcms_feature_mvp\run_xcalibur_workbench.py"
$ResolvedInput = Join-Path $Root $InputDir
$ResolvedOutput = Join-Path $Root $OutputDir

python $Runner `
    --input-dir $ResolvedInput `
    --output-dir $ResolvedOutput `
    --project-id user_raw_5samples `
    --rt-min $RtMin `
    --rt-max $RtMax `
    --mz-min $MzMin `
    --mz-max $MzMax `
    --alignment-rt-start 0 `
    --alignment-signal bpc `
    --alignment-match-window-min 3.0 `
    --max-peaks-per-scan 220 `
    --spectrum-min-intensity 0 `
    --rt-bin-count $RtBinCount `
    --mz-bin-count $MzBinCount `
    --heatmap-normalization-methods tic_log1p `
    --default-heatmap-normalization tic_log1p `
    --auto-feature-top-n 30
