@echo off
echo ============================================
echo   Instalacja pakietow R dla proteomiki
echo ============================================
echo.
echo WAZNE: Zamknij RStudio przed uruchomieniem!
echo.
pause

echo.
echo === Instaluje scales (1.4.0+) ===
Rscript -e "install.packages('scales', repos='https://cloud.r-project.org', dependencies=TRUE)"

echo.
echo === Instaluje ggrepel ===
Rscript -e "install.packages('ggrepel', repos='https://cloud.r-project.org', dependencies=TRUE)"

echo.
echo === Instaluje patchwork ===
Rscript -e "install.packages('patchwork', repos='https://cloud.r-project.org', dependencies=TRUE)"

echo.
echo === Instaluje BiocManager ===
Rscript -e "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager', repos='https://cloud.r-project.org')"

echo.
echo === Instaluje vsn (Bioconductor) ===
Rscript -e "BiocManager::install('vsn', ask=FALSE, update=FALSE)"

echo.
echo === Instaluje imputeLCMD (Bioconductor) ===
Rscript -e "BiocManager::install('imputeLCMD', ask=FALSE, update=FALSE)"

echo.
echo === Instaluje limma (Bioconductor) ===
Rscript -e "BiocManager::install('limma', ask=FALSE, update=FALSE)"

echo.
echo ============================================
echo   Instalacja zakonczona!
echo ============================================
echo.
echo Mozesz teraz otworzyc RStudio i uruchomic:
echo   source("run_analysis.R")
echo.
pause
