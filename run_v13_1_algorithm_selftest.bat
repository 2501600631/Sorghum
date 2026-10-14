@echo off
setlocal
cd /d "%~dp0"
python run_v13_1_algorithm_selftest.py
if errorlevel 1 exit /b 1
echo V13.2.0-98 algorithm-family self-test completed successfully.
endlocal
