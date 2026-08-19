param(
    [string]$InputDir = "data/raw/zenodo_5005513",
    [string]$OutputDir = "data/converted/mzML"
)

$ErrorActionPreference = "Stop"
$Root = Split-Path -Parent $MyInvocation.MyCommand.Path
$Workspace = Split-Path -Parent $Root
$Parser = Join-Path $Workspace ".local-tools/ThermoRawFileParser/current/ThermoRawFileParser.exe"

if (-not (Test-Path $Parser)) {
    throw "ThermoRawFileParser not found: $Parser"
}

$ResolvedInput = Join-Path $Root $InputDir
$ResolvedOutput = Join-Path $Root $OutputDir
New-Item -ItemType Directory -Force -Path $ResolvedOutput | Out-Null

Get-ChildItem -Path $ResolvedInput -Filter "*.raw" | ForEach-Object {
    $targetMzml = Join-Path $ResolvedOutput ($_.BaseName + ".mzML")
    if (Test-Path -LiteralPath $targetMzml) {
        Write-Host "Skipping existing $($_.BaseName).mzML"
        return
    }
    Write-Host "Converting $($_.Name)"
    & $Parser -i $_.FullName -o $ResolvedOutput -f 2 -m 0 -l 3
}

Write-Host "mzML output: $ResolvedOutput"
