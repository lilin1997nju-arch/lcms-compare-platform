param(
    [string]$PythonPath = ""
)

$ErrorActionPreference = "Stop"
$Root = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
if ($PythonPath -eq "") {
    $candidates = @(
        (Join-Path $Root ".venv\Scripts\python.exe"),
        (Join-Path $Root "venv\Scripts\python.exe")
    )
    try { $candidates += (Get-Command python.exe -ErrorAction Stop).Source } catch { }
    try { $candidates += (Get-Command py.exe -ErrorAction Stop).Source } catch { }
    $PythonPath = $candidates | Where-Object { $_ -and (Test-Path -LiteralPath $_) } | Select-Object -First 1
}
if (-not $PythonPath) { throw "Python 3.10+ was not found. Pass -PythonPath explicitly." }

$required = @(
    "lcms_feature_mvp\core\lcms_peak_first.py",
    "lcms_feature_mvp\core\lcms_parser.py",
    "lcms_feature_mvp\run_peak_first_compare.py",
    "lcms_feature_mvp\serve_peak_first_compare.py",
    "lcms_department_platform\server.py",
    "lcms_department_platform\start_lcms_department_platform.ps1",
    "start_department_platform.cmd",
    "requirements.txt"
)
foreach ($relativePath in $required) {
    if (-not (Test-Path -LiteralPath (Join-Path $Root $relativePath))) {
        throw "Required file is missing: $relativePath"
    }
}

Push-Location $Root
try {
    & $PythonPath -m compileall -q lcms_feature_mvp lcms_department_platform examples
    if ($LASTEXITCODE -ne 0) { throw "Python compilation failed." }
    if (Test-Path -LiteralPath (Join-Path $Root "lcms_feature_mvp\tests")) {
        & $PythonPath -m unittest discover -s lcms_feature_mvp\tests -p "test_*.py"
        if ($LASTEXITCODE -ne 0) { throw "LC-MS tests failed." }
    }
} finally {
    Pop-Location
}
Write-Host "LC-MS migration package verification passed."
