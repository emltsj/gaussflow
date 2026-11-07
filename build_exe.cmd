@echo off
rem Build the single-file Windows GUI executable (dist\GaussFlow.exe).
rem Requires: python + pip install pyinstaller customtkinter rdkit cclib openpyxl
cd /d "%~dp0"
python assets\make_icon.py
python -m PyInstaller --noconfirm --clean --onefile --windowed ^
  --name GaussFlow ^
  --icon assets\icon.ico ^
  --collect-data rdkit ^
  --collect-data customtkinter ^
  GaussFlow.py
echo.
echo Done. Executable: %~dp0dist\GaussFlow.exe
pause
