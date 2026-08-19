@echo off
setlocal
cd /d "%~dp0"
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "%~dp0Install-LCMSPlatform.ps1"
if errorlevel 1 (
  echo.
  echo LC-MS Department Platform installation failed.
  pause
  exit /b 1
)
endlocal
