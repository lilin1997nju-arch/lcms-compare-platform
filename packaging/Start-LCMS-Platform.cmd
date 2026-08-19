@echo off
setlocal
cd /d "%~dp0"
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "%~dp0Start-LCMS-Platform.ps1" -OpenBrowser
if errorlevel 1 (
  echo.
  echo LC-MS Department Platform failed to start.
  pause
)
endlocal
