@echo off
echo ============================================
echo   Instalacja pakietow Python dla AI Analysis
echo ============================================
echo.

echo Sprawdzanie instalacji Python...
python --version
if errorlevel 1 (
    echo BLAD: Python nie jest zainstalowany lub nie ma go w PATH.
    pause
    exit /b 1
)

echo.
echo Aktualizacja pip...
python -m pip install --upgrade pip

echo.
echo Instalacja bibliotek (moze to potrwac kilka minut)...
echo   - torch, transformers (AI models)
echo   - pandas, numpy (Data manipulation)
echo   - scikit-learn, umap-learn (Dimensionality reduction)
echo   - matplotlib, seaborn (Visualization)
echo   - requests (API calls)

python -m pip install torch transformers pandas numpy scikit-learn umap-learn matplotlib seaborn requests pydantic plotly

echo.
echo ============================================
echo   Instalacja pakietow Python zakonczona!
echo ============================================
pause
