param(
    [string]$InputDir = ".\data\mzML",
    [string]$OutputDir = ".\outputs_peak_first",
    [ValidateSet("All", "PTM", "SVA")]
    [string]$Group = "All",
    [string[]]$SampleContains = @(),
    [string[]]$IncludeSample = @(),
    [string]$ReferenceSample = "",
    [int]$TopNPeaks = 60,
    [int]$TopNMz = 40,
    [int]$TopNChangedMz = 15,
    [double]$MzToleranceDa = 0.5
)

$ErrorActionPreference = "Stop"
$Root = Split-Path -Parent $MyInvocation.MyCommand.Path
$Workspace = Split-Path -Parent $Root
$Runner = Join-Path $Workspace "lcms_feature_mvp\run_peak_first_compare.py"
$ResolvedInput = Join-Path $Root $InputDir
$ResolvedOutput = Join-Path $Root $OutputDir

$argsList = @(
    $Runner,
    "--input-dir", $ResolvedInput,
    "--output-dir", $ResolvedOutput,
    "--project-id", "user_raw_5samples_peak_first",
    "--top-n-peaks", $TopNPeaks,
    "--top-n-mz", $TopNMz,
    "--top-n-changed-mz", $TopNChangedMz,
    "--max-spectrum-points-per-scan", 200,
    "--mz-tolerance-mode", "da",
    "--mz-tolerance-da", $MzToleranceDa,
    "--max-peaks-per-scan", 500,
    "--spectrum-min-intensity", 0
)

if ($ReferenceSample -ne "") {
    $argsList += @("--reference-sample", $ReferenceSample)
}

if ($Group -eq "PTM") {
    $argsList += @("--sample-contains", "PTM")
} elseif ($Group -eq "SVA") {
    $argsList += @("--sample-contains", "SVA")
}

foreach ($token in $SampleContains) {
    if ($token -ne "") {
        $argsList += @("--sample-contains", $token)
    }
}

foreach ($sample in $IncludeSample) {
    if ($sample -ne "") {
        $argsList += @("--include-sample", $sample)
    }
}

python @argsList
