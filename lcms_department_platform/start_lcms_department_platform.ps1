param(
    [int]$Port = 8770,
    [string]$HostName = "127.0.0.1",
    [string]$ParserPath = "",
    [string]$PythonPath = "",
    [string]$RootPath = "",
    [switch]$OpenBrowser,
    [switch]$Visible
)

$ErrorActionPreference = "Stop"
$Root = Split-Path -Parent $MyInvocation.MyCommand.Path
$Workspace = Split-Path -Parent $Root
if ($RootPath -ne "") {
    $Root = (Resolve-Path -LiteralPath $RootPath).Path
}
$Server = Join-Path (Split-Path -Parent $MyInvocation.MyCommand.Path) "server.py"
if ($ParserPath -eq "") {
    $parserCandidates = @(
        (Join-Path $Root "tools\ThermoRawFileParser\ThermoRawFileParser.exe"),
        (Join-Path $Root "ThermoRawFileParser.exe"),
        (Join-Path $Workspace ".local-tools\ThermoRawFileParser\current\ThermoRawFileParser.exe")
    )
    $ParserPath = $parserCandidates | Where-Object { Test-Path -LiteralPath $_ } | Select-Object -First 1
    if (-not $ParserPath) { $ParserPath = $parserCandidates[0] }
}
$pythonCandidates = @()
if ($PythonPath -ne "") { $pythonCandidates += $PythonPath }
$pythonCandidates += @(
    (Join-Path $Workspace ".venv\Scripts\python.exe"),
    (Join-Path $Workspace "venv\Scripts\python.exe")
)
$codexPythonRoots = Get-ChildItem -Path "C:\Users" -Directory -ErrorAction SilentlyContinue | ForEach-Object {
    Join-Path $_.FullName ".cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe"
}
$pythonCandidates += $codexPythonRoots
$pythonCommand = $null
try { $pythonCommand = (Get-Command python.exe -ErrorAction Stop).Source } catch { }
if ($pythonCommand) { $pythonCandidates += $pythonCommand }
try { $pythonCandidates += (Get-Command py.exe -ErrorAction Stop).Source } catch { }

function Test-PythonCandidate([string]$Candidate) {
    if (-not $Candidate -or -not (Test-Path -LiteralPath $Candidate -PathType Leaf)) { return $false }
    try {
        $versionText = (& $Candidate --version 2>&1 | Out-String).Trim()
        if ($LASTEXITCODE -ne 0 -or $versionText -notmatch '^Python 3\.(\d+)') { return $false }
        return ([int]$Matches[1] -ge 10)
    } catch {
        return $false
    }
}

$PythonExe = $null
foreach ($candidate in $pythonCandidates) {
    if (Test-PythonCandidate $candidate) {
        $PythonExe = (Resolve-Path -LiteralPath $candidate).Path
        break
    }
}
if (-not $PythonExe) { throw "Python 3.10+ was not found. Install Python or pass -PythonPath 'C:\path\python.exe'." }
$LogDir = Join-Path $Root "logs"
$StateDir = Join-Path $Root "state"
New-Item -ItemType Directory -Force -Path $LogDir | Out-Null
New-Item -ItemType Directory -Force -Path $StateDir | Out-Null

$argsList = @(
    $Server,
    "--host", $HostName,
    "--port", $Port,
    "--root", $Root,
    "--parser-path", $ParserPath
)

if ($Visible) {
    & $PythonExe @argsList
} else {
    $proc = Start-Process -FilePath $PythonExe `
        -ArgumentList $argsList `
        -WorkingDirectory $Workspace `
        -WindowStyle Hidden `
        -PassThru `
        -RedirectStandardOutput (Join-Path $LogDir "lcms_department_platform_$Port.out.log") `
        -RedirectStandardError (Join-Path $LogDir "lcms_department_platform_$Port.err.log")
    Set-Content -LiteralPath (Join-Path $StateDir "lcms_department_platform_$Port.pid") -Value $proc.Id -Encoding ascii
    Start-Sleep -Seconds 2
    $proc.Refresh()
    if ($proc.HasExited) {
        $errorLog = Join-Path $LogDir "lcms_department_platform_$Port.err.log"
        $details = if (Test-Path -LiteralPath $errorLog) { Get-Content -LiteralPath $errorLog -Raw } else { "no error log was written" }
        throw "LC-MS department platform failed to start with $PythonExe. $details"
    }
    Write-Host "LC-MS department platform started: http://$HostName`:$Port/"
    Write-Host "PID: $($proc.Id)"
    if ($OpenBrowser) { Start-Process "http://$HostName`:$Port/" }
}
