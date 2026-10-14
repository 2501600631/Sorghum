@echo off
setlocal
cd /d "%~dp0"
python run_v13_1_submission_readiness.py
if errorlevel 1 (
  echo.
  echo V13.2.0-98 submission-readiness check FAILED. Complete only genuine provenance/reference-method records and a valid formal four-target run; do not invent metadata.
  exit /b 1
)
echo.
echo V13.2.0-98 submission-readiness checks passed.
endlocal
