param(
    [string]$InputDir = "..\lcms_feature_mvp\质谱数据",
    [string]$OutputDir = ".\data\mzML"
)

$ErrorActionPreference = "Stop"
$Root = Split-Path -Parent $MyInvocation.MyCommand.Path
$Workspace = Split-Path -Parent $Root
$Parser = Join-Path $Workspace ".local-tools\ThermoRawFileParser\current\ThermoRawFileParser.exe"

if (-not (Test-Path -LiteralPath $Parser)) {
    throw "ThermoRawFileParser not found: $Parser"
}

$ResolvedInput = Resolve-Path -LiteralPath (Join-Path $Root $InputDir)
$ResolvedOutput = Join-Path $Root $OutputDir
New-Item -ItemType Directory -Force -Path $ResolvedOutput | Out-Null

Get-ChildItem -LiteralPath $ResolvedInput -Filter "*.raw" | ForEach-Object {
    $TargetMzml = Join-Path $ResolvedOutput ($_.BaseName + ".mzML")
    if (Test-Path -LiteralPath $TargetMzml) {
        Write-Host "Skipping existing $($_.BaseName).mzML"
        return
    }
    Write-Host "Converting $($_.Name)"
    & $Parser -i $_.FullName -o $ResolvedOutput -f 2 -m 0 -l 3
    if ($LASTEXITCODE -ne 0) {
        throw "ThermoRawFileParser failed for $($_.FullName)"
    }
}

Write-Host "mzML output: $ResolvedOutput"
