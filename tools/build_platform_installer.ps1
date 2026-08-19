param(
    [string]$OutputDir = ""
)

$ErrorActionPreference = "Stop"
$ProgressPreference = "SilentlyContinue"
$Workspace = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
if (-not $OutputDir) {
    $OutputDir = Join-Path $Workspace "release"
}
$OutputDir = [System.IO.Path]::GetFullPath($OutputDir)
$BundleRoot = Join-Path $OutputDir "LCMS_Department_Platform"
$PayloadZip = Join-Path $OutputDir "LCMS-Payload.zip"
$IExpressStage = Join-Path $OutputDir "iexpress_stage"
$Installer = Join-Path $OutputDir "LCMS_Department_Platform_Installer.exe"

foreach ($path in @($BundleRoot, $IExpressStage)) {
    if (Test-Path -LiteralPath $path) {
        Remove-Item -LiteralPath $path -Recurse -Force
    }
}
foreach ($path in @($PayloadZip, $Installer)) {
    if (Test-Path -LiteralPath $path) {
        Remove-Item -LiteralPath $path -Force
    }
}
New-Item -ItemType Directory -Force -Path $OutputDir, $BundleRoot | Out-Null

function Copy-RequiredFile([string]$Source, [string]$Destination) {
    if (-not (Test-Path -LiteralPath $Source -PathType Leaf)) {
        throw "Required file is missing: $Source"
    }
    $parent = Split-Path -Parent $Destination
    New-Item -ItemType Directory -Force -Path $parent | Out-Null
    Copy-Item -LiteralPath $Source -Destination $Destination -Force
}

function Copy-RequiredTree([string]$Source, [string]$Destination) {
    if (-not (Test-Path -LiteralPath $Source -PathType Container)) {
        throw "Required directory is missing: $Source"
    }
    New-Item -ItemType Directory -Force -Path $Destination | Out-Null
    Copy-Item -Path (Join-Path $Source "*") -Destination $Destination -Recurse -Force
}

$appRoot = Join-Path $BundleRoot "app"
$platformRoot = Join-Path $appRoot "lcms_department_platform"
$featureRoot = Join-Path $appRoot "lcms_feature_mvp"
$parserRoot = Join-Path $platformRoot "tools\ThermoRawFileParser"
New-Item -ItemType Directory -Force -Path $platformRoot, $featureRoot | Out-Null

# Platform service and the MCP adapter.
Copy-RequiredFile (Join-Path $Workspace "lcms_department_platform\server.py") (Join-Path $platformRoot "server.py")
Copy-RequiredFile (Join-Path $Workspace "lcms_department_platform\mcp_server.py") (Join-Path $platformRoot "mcp_server.py")

# Only modules invoked by the current department task pipeline are included.
foreach ($name in @("run_peak_first_compare.py", "run_msms_compare.py", "serve_peak_first_compare.py")) {
    Copy-RequiredFile (Join-Path $Workspace "lcms_feature_mvp\$name") (Join-Path $featureRoot $name)
}
Copy-RequiredTree (Join-Path $Workspace "lcms_feature_mvp\core") (Join-Path $featureRoot "core")
Copy-RequiredFile (Join-Path $Workspace "lcms_feature_mvp\ui\vendor\3Dmol-min.js") (Join-Path $featureRoot "ui\vendor\3Dmol-min.js")

# ThermoRawFileParser is a self-contained .NET deployment; keep all of its
# assemblies because the executable resolves them from this directory.
Copy-RequiredTree (Join-Path $Workspace "lcms_department_platform\tools\ThermoRawFileParser") $parserRoot

# Empty runtime data directories are created by the service when first used.
New-Item -ItemType Directory -Force -Path (Join-Path $platformRoot "jobs"), (Join-Path $platformRoot "logs"), (Join-Path $platformRoot "state") | Out-Null

$runtimeRoot = Join-Path $BundleRoot "runtime"
New-Item -ItemType Directory -Force -Path $runtimeRoot, (Join-Path $runtimeRoot "DLLs"), (Join-Path $runtimeRoot "Lib") | Out-Null

# Prefer the already-built portable runtime when rebuilding in a restricted
# workspace. Otherwise, source it from an installed Python distribution.
$runtimeCandidates = @(
    (Join-Path $Workspace "release\LCMS_Department_Platform\runtime"),
    (Join-Path $Workspace "release\installer_test_stage\runtime"),
    (Join-Path $Workspace "release\installer_test_install\runtime")
)
$existingRuntime = $runtimeCandidates |
    Where-Object { Test-Path -LiteralPath (Join-Path $_ "python.exe") -PathType Leaf } |
    Select-Object -First 1
if ($existingRuntime) {
    Copy-RequiredTree $existingRuntime $runtimeRoot
} else {
    $pythonCandidates = @()
    try { $pythonCandidates += (Get-Command python.exe -ErrorAction Stop).Source } catch { }
    $pythonCandidates += @(
        (Join-Path $env:LOCALAPPDATA "Programs\Python\Python312\python.exe"),
        (Join-Path $env:LOCALAPPDATA "Programs\Python\Python311\python.exe"),
        (Join-Path $env:LOCALAPPDATA "Programs\Python\Python310\python.exe")
    )
    $pythonExe = $pythonCandidates | Where-Object { $_ -and (Test-Path -LiteralPath $_ -PathType Leaf) } | Select-Object -First 1
    if (-not $pythonExe) {
        throw "Python 3.10+ was not found for building the self-contained runtime."
    }
    $pythonHome = Split-Path $pythonExe -Parent

    foreach ($source in Get-ChildItem -LiteralPath $pythonHome -File | Where-Object {
        $_.Name -match "^(python\.exe|pythonw\.exe|python\d*\.dll|vcruntime.*\.dll|LICENSE\.txt)$"
    }) {
        Copy-Item -LiteralPath $source.FullName -Destination (Join-Path $runtimeRoot $source.Name) -Force
    }
    Copy-RequiredTree (Join-Path $pythonHome "DLLs") (Join-Path $runtimeRoot "DLLs")

    # Python standard-library-only runtime. Test and cache directories are not
    # needed by this application and are excluded to reduce the package size.
    & robocopy (Join-Path $pythonHome "Lib") (Join-Path $runtimeRoot "Lib") /E /XF *.pyc /XD __pycache__ test tests idlelib tkinter | Out-Null
    if ($LASTEXITCODE -ge 8) {
        throw "Failed to copy Python standard library (robocopy exit code $LASTEXITCODE)."
    }
}

foreach ($name in @("Start-LCMS-Platform.cmd", "Start-LCMS-Platform.ps1", "Stop-LCMS-Platform.cmd", "README.txt")) {
    Copy-RequiredFile (Join-Path $PSScriptRoot "..\packaging\$name") (Join-Path $BundleRoot $name)
}

# Remove interpreter caches that may exist in the development workspace.
Get-ChildItem -LiteralPath $BundleRoot -Recurse -File -Filter "*.pyc" -ErrorAction SilentlyContinue | Remove-Item -Force
Get-ChildItem -LiteralPath $BundleRoot -Recurse -Directory -Filter "__pycache__" -ErrorAction SilentlyContinue |
    Sort-Object FullName -Descending |
    Remove-Item -Recurse -Force

$runtimeCheck = & (Join-Path $runtimeRoot "python.exe") -c "import csv, json, sqlite3, http.server, pathlib; print('portable-python-ok')" 2>&1
if ($LASTEXITCODE -ne 0 -or $runtimeCheck -notmatch "portable-python-ok") {
    throw "Bundled Python runtime failed validation: $runtimeCheck"
}

# The validation import may recreate bytecode caches; remove them once more
# immediately before archiving the final payload.
Get-ChildItem -LiteralPath $BundleRoot -Recurse -File -Filter "*.pyc" -ErrorAction SilentlyContinue | Remove-Item -Force
Get-ChildItem -LiteralPath $BundleRoot -Recurse -Directory -Filter "__pycache__" -ErrorAction SilentlyContinue |
    Sort-Object FullName -Descending |
    Remove-Item -Recurse -Force

if (Test-Path -LiteralPath $PayloadZip) {
    Remove-Item -LiteralPath $PayloadZip -Force
}
Compress-Archive -Path (Join-Path $BundleRoot "*") -DestinationPath $PayloadZip -CompressionLevel Optimal

# IExpress receives only three files. The payload remains a ZIP so the
# installer preserves all subdirectories without flattening them.
New-Item -ItemType Directory -Force -Path $IExpressStage | Out-Null
foreach ($name in @("Install-LCMSPlatform.ps1", "Install-LCMSPlatform.cmd")) {
    Copy-RequiredFile (Join-Path $PSScriptRoot "..\packaging\$name") (Join-Path $IExpressStage $name)
}
Copy-RequiredFile $PayloadZip (Join-Path $IExpressStage "LCMS-Payload.zip")

$targetName = $Installer
$stagePath = $IExpressStage
$sed = @"
[Version]
Class=IEXPRESS
SEDVersion=3
[Options]
PackagePurpose=InstallApp
ShowInstallProgramWindow=1
HideExtractAnimation=1
UseLongFileName=1
InsideCompressed=0
CAB_FixedSize=0
CAB_ResvCodeSigning=0
RebootMode=N
InstallPrompt=%InstallPrompt%
DisplayLicense=
FinishMessage=%FinishMessage%
TargetName=$targetName
FriendlyName=%FriendlyName%
AppLaunched=%AppLaunched%
PostInstallCmd=<None>
AdminQuietInstCmd=
UserQuietInstCmd=
SourceFiles=SourceFiles
[Strings]
InstallPrompt=
FinishMessage=LC-MS Department Platform installation completed.
FriendlyName=LC-MS Department Platform
AppLaunched=cmd.exe /c Install-LCMSPlatform.cmd
FILE0="Install-LCMSPlatform.ps1"
FILE1="Install-LCMSPlatform.cmd"
FILE2="LCMS-Payload.zip"
[SourceFiles]
SourceFiles0=$stagePath
[SourceFiles0]
%FILE0%=
%FILE1%=
%FILE2%=
"@
$sedPath = Join-Path $IExpressStage "LCMS_Department_Platform.sed"
Set-Content -LiteralPath $sedPath -Value $sed -Encoding ascii

$iexpressProcess = Start-Process -FilePath "iexpress.exe" -ArgumentList @("/N", $sedPath) -PassThru -Wait
for ($attempt = 0; $attempt -lt 60 -and -not (Test-Path -LiteralPath $Installer -PathType Leaf); $attempt++) {
    Start-Sleep -Milliseconds 500
}
if ($iexpressProcess.ExitCode -ne 0 -or -not (Test-Path -LiteralPath $Installer -PathType Leaf)) {
    throw "IExpress failed to create the installer. Exit code: $($iexpressProcess.ExitCode)"
}

$bundleBytes = (Get-ChildItem $BundleRoot -Recurse -File | Measure-Object Length -Sum).Sum
$installerBytes = (Get-Item $Installer).Length
Write-Host "Clean bundle: $BundleRoot"
Write-Host "Installer:    $Installer"
Write-Host ("Bundle size:  {0:N1} MB" -f ($bundleBytes / 1MB))
Write-Host ("Installer:    {0:N1} MB" -f ($installerBytes / 1MB))
