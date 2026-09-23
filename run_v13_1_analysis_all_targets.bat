@echo off
setlocal
cd /d "%~dp0"
python run_v13_1_analysis.py
if errorlevel 1 (
  echo.
  echo V13.2.0-98 formal all-target analysis FAILED.
  exit /b 1
)
echo.
echo V13.2.0-98 formal all-target analysis completed successfully.
endlocal
