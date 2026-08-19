@echo off
setlocal
cd /d "%~dp0"
powershell.exe -NoProfile -ExecutionPolicy Bypass -Command ^
  "$root=Join-Path (Get-Location) 'lcms_department_platform\state'; $pidFile=Join-Path $root 'lcms_department_platform_8770.pid'; if(Test-Path -LiteralPath $pidFile){$id=[int](Get-Content -LiteralPath $pidFile); $p=Get-Process -Id $id -ErrorAction SilentlyContinue; if($p){Stop-Process -Id $id -Force; Write-Host ('Stopped PID '+$id)}; Remove-Item -LiteralPath $pidFile -Force -ErrorAction SilentlyContinue} else {Write-Host 'No recorded department platform process.'}"
endlocal
