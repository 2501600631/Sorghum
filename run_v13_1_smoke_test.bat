@echo off
setlocal
cd /d "%~dp0"
python run_v13_1_smoke_test.py
if errorlevel 1 (
  echo.
  echo V13.2.0-98 smoke test FAILED.
  exit /b 1
)
echo.
echo V13.2.0-98 smoke test completed successfully.
endlocal
