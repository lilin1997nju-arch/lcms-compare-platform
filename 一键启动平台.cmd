@echo off
setlocal
cd /d "%~dp0"
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "%~dp0lcms_department_platform\start_lcms_department_platform.ps1" -OpenBrowser
if errorlevel 1 (
  echo.
  echo LC-MS department platform failed to start.
  pause
)
endlocal
