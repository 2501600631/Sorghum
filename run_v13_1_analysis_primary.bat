@echo off
setlocal
cd /d "%~dp0"
python run_v13_1_analysis_primary.py
if errorlevel 1 (
  echo.
  echo V13.2.0-98 primary preflight FAILED.
  exit /b 1
)
echo.
echo V13.2.0-98 primary preflight completed. NOT FOR SUBMISSION.
endlocal
