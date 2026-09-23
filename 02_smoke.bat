@echo off
setlocal
cd /d "%~dp0"
python run_v13_1_smoke_test.py
if errorlevel 1 (
 echo FAILED. Please fix the error before the next step.
 pause
 exit /b 1
)
echo PASSED.
pause
