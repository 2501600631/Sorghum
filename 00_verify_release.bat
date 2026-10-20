@echo off
setlocal
cd /d "%~dp0"
python verify_release.py
if errorlevel 1 exit /b 1
