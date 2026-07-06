param(
    [int]$Port = 8767,
    [string]$OutputDir = ".\outputs",
    [string[]]$Comparison = @()
)

$ErrorActionPreference = "Stop"
$Root = Split-Path -Parent $MyInvocation.MyCommand.Path
$Workspace = Split-Path -Parent $Root
$Server = Join-Path $Workspace "lcms_feature_mvp\serve_xcalibur_workbench.py"
$ResolvedOutput = Join-Path $Root $OutputDir

$argsList = @($Server, "--output-dir", $ResolvedOutput, "--port", $Port)
foreach ($item in $Comparison) {
    if ($item -ne "") {
        $argsList += @("--comparison", $item)
    }
}

python @argsList
