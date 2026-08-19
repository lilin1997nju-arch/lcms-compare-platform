param(
    [int]$Port = 8770,
    [switch]$OpenBrowser,
    [switch]$Visible
)

$ErrorActionPreference = "Stop"
$InstallRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$Workspace = Join-Path $InstallRoot "app"
$PlatformRoot = Join-Path $Workspace "lcms_department_platform"
$Server = Join-Path $PlatformRoot "server.py"
$Python = Join-Path $InstallRoot "runtime\python.exe"
$Parser = Join-Path $PlatformRoot "tools\ThermoRawFileParser\ThermoRawFileParser.exe"

foreach ($required in @($Server, $Python, $Parser)) {
    if (-not (Test-Path -LiteralPath $required -PathType Leaf)) {
        throw "Required platform file is missing: $required"
    }
}

$logDir = Join-Path $PlatformRoot "logs"
$stateDir = Join-Path $PlatformRoot "state"
New-Item -ItemType Directory -Force -Path $logDir, $stateDir | Out-Null

$arguments = @(
    $Server,
    "--host", "127.0.0.1",
    "--port", $Port,
    "--root", $PlatformRoot,
    "--parser-path", $Parser,
    "--python", $Python
)

if ($Visible) {
    & $Python @arguments
    exit $LASTEXITCODE
}

$proc = Start-Process -FilePath $Python `
    -ArgumentList $arguments `
    -WorkingDirectory $Workspace `
    -WindowStyle Hidden `
    -PassThru `
    -RedirectStandardOutput (Join-Path $logDir "lcms_department_platform_$Port.out.log") `
    -RedirectStandardError (Join-Path $logDir "lcms_department_platform_$Port.err.log")

Set-Content -LiteralPath (Join-Path $stateDir "lcms_department_platform_$Port.pid") -Value $proc.Id -Encoding ascii
Start-Sleep -Seconds 2
$proc.Refresh()
if ($proc.HasExited) {
    $errorLog = Join-Path $logDir "lcms_department_platform_$Port.err.log"
    $details = if (Test-Path -LiteralPath $errorLog) { Get-Content -LiteralPath $errorLog -Raw } else { "no error log was written" }
    throw "LC-MS Department Platform failed to start. $details"
}

Write-Host "LC-MS Department Platform started at http://127.0.0.1:$Port/"
if ($OpenBrowser) { Start-Process "http://127.0.0.1:$Port/" }
