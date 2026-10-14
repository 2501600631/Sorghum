@echo off
setlocal
cd /d "%~dp0"
python run_v13_1_validate.py --verify-checksums
if errorlevel 1 (
  echo.
  echo V13.2.0-98 validation FAILED.
  exit /b 1
)
echo.
echo V13.2.0-98 validation completed successfully.
endlocal
