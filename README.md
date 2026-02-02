# E. coli Proteomics Analysis - UPS2 vs UPS1 Benchmark

Analiza danych proteomicznych z eksperymentu "Dynamic Range Benchmark" (Cox et al. 2014, PRIDE PXD000279).

## Cel projektu

1. **Analiza różnicowa UPS2 vs UPS1** w tle białek E. coli
2. **Empiryczny FDR** - ocena fałszywych trafień na białkach E. coli (powinny być niezmienne)
3. **Dokładność fold-change** - porównanie obserwowanych vs oczekiwanych zmian dla białek UPS
4. **Porównanie wariantów pipeline'u** - median vs VSN normalizacja, brak vs QRILC imputacja

## Struktura projektu

```
ecoli_proteomics/
├── data/                          # Dane wejściowe (pobrane automatycznie)
│   └── proteinGroups.txt
├── scripts/
│   ├── 01_download_data.R         # Pobieranie danych z PRIDE
│   ├── 02_preprocessing.R         # Filtrowanie, normalizacja, imputacja
│   ├── 03_differential_analysis.R # Analiza różnicowa z limma
│   ├── 04_benchmark.R             # Obliczenie FDR i dokładności FC
│   └── 05_report.Rmd              # Generowanie raportu
├── results/
│   ├── figures/                   # Wykresy (volcano, PCA, heatmap)
│   └── tables/                    # Tabele CSV z wynikami
├── run_analysis.R                 # Główny skrypt uruchamiający
├── report.html                    # Raport końcowy
└── README.md
```

## Wymagania

### R pakiety (CRAN)
- tidyverse
- pheatmap
- ggrepel
- RColorBrewer
- knitr, kableExtra, rmarkdown
- patchwork

### R pakiety (Bioconductor)
- limma
- vsn
- imputeLCMD

## Uruchomienie analizy

### Opcja 1: Cała analiza naraz

```r
# W R/RStudio, ustaw working directory na folder projektu
setwd("c:/Users/konrad_guest/Documents/GitHub/ecoli_proteomics")

# Uruchom główny skrypt
source("run_analysis.R")
```

### Opcja 2: Krok po kroku

```r
setwd("c:/Users/konrad_guest/Documents/GitHub/ecoli_proteomics")

# 1. Pobierz dane
source("scripts/01_download_data.R")

# 2. Preprocessing
source("scripts/02_preprocessing.R")

# 3. Analiza różnicowa
source("scripts/03_differential_analysis.R")

# 4. Benchmark
source("scripts/04_benchmark.R")

# 5. Generuj raport
rmarkdown::render("scripts/05_report.Rmd", output_file = "../report.html")
```

## Opis danych

- **UPS1**: 48 białek ludzkich w równym stężeniu (10,000 fmol każde)
- **UPS2**: Te same białka w 5 grupach o różnych stężeniach (dynamic range)
- **Tło E. coli**: ~4000 białek bakteryjnych

| Grupa | Stężenie UPS2 | Expected log2FC |
|-------|---------------|-----------------|
| A (5) | 50,000 fmol   | +2.32           |
| B (10)| 5,000 fmol    | -1.00           |
| C (8) | 500 fmol      | -4.32           |
| D (15)| 50 fmol       | -7.64           |
| E (10)| 5 fmol        | -10.97          |

## Warianty pipeline'u

| Wariant | Normalizacja | Imputacja |
|---------|--------------|-----------|
| 1       | Median centering | Brak (analiza z NA) |
| 2       | VSN | QRILC |

## Wyniki

Po uruchomieniu analizy znajdziesz:

1. **report.html** - pełny raport z wykresami i wnioskami
2. **results/figures/** - wykresy PNG:
   - `volcano_*.png` - volcano plots
   - `pca_*.png` - PCA
   - `heatmap_*.png` - heatmapy top białek
   - `expected_vs_observed_*.png` - dokładność FC
3. **results/tables/** - tabele CSV:
   - `differential_results_*.csv` - wyniki analizy różnicowej
   - `benchmark_comparison.csv` - porównanie wariantów

## Referencje

- Cox J, et al. (2014) Accurate proteome-wide label-free quantification. *Mol Cell Proteomics* 13(9):2513-26. DOI: [10.1074/mcp.M113.031591](https://doi.org/10.1074/mcp.M113.031591)
- PRIDE: [PXD000279](https://www.ebi.ac.uk/pride/archive/projects/PXD000279)
- UPS Standards: [Sigma-Aldrich](https://www.sigmaaldrich.com/technical-documents/articles/biology/ups1-and-ups2-proteomic.html)
