# ==============================================================================
# run_analysis.R - Główny skrypt uruchamiający całą analizę
# ==============================================================================
#
# PRZED URUCHOMIENIEM:
# 1. Zamknij wszystkie sesje R/RStudio
# 2. Otwórz nową sesję R
# 3. Uruchom: source("scripts/00_install_packages.R")
# 4. Zrestartuj R (Ctrl+Shift+F10)
# 5. Dopiero teraz uruchom ten skrypt
#
# ==============================================================================

.libPaths("C:/Users/konrad_guest/AppData/Local/Programs/R/R-4.4.3/library")
if (file.exists(".RData")) file.remove(".RData")


message("========================================")
message("  E. coli Proteomics Analysis Pipeline  ")
message("       UPS2 vs UPS1 Benchmark           ")
message("========================================")
message("")

# Ustaw working directory jeśli potrzebne
# setwd("c:/Users/konrad_guest/Documents/GitHub/ecoli_proteomics")

# ==============================================================================
# Sprawdź czy pakiety zainstalowane (NIE instaluj tutaj!)
# ==============================================================================

message("=== Sprawdzanie pakietów ===")

required_packages <- c(
  "tidyverse", "pheatmap", "ggrepel", "RColorBrewer", "scales",
  "knitr", "kableExtra", "rmarkdown", "patchwork",
  "limma", "vsn", "imputeLCMD"
)

missing <- c()
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    missing <- c(missing, pkg)
  }
}

if (length(missing) > 0) {
  stop(
    "\n\nBrakujące pakiety: ", paste(missing, collapse = ", "), 
    "\n\nUruchom najpierw:",
    "\n  source('scripts/00_install_packages.R')",
    "\n\nPotem zrestartuj R i spróbuj ponownie."
  )
}

message("Wszystkie pakiety dostępne.\n")

# Utwórz katalogi wynikowe
dir.create("data", showWarnings = FALSE)
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# KROK 1: Pobieranie danych
# ==============================================================================

message("\n========================================")
message("KROK 1: Pobieranie danych z PRIDE")
message("========================================")

source("scripts/01_download_data.R")

# ==============================================================================
# KROK 2: Preprocessing
# ==============================================================================

message("\n========================================")
message("KROK 2: Preprocessing danych")
message("========================================")

source("scripts/02_preprocessing.R")

# ==============================================================================
# KROK 3: Analiza różnicowa
# ==============================================================================

message("\n========================================")
message("KROK 3: Analiza różnicowa UPS2 vs UPS1")
message("========================================")

source("scripts/03_differential_analysis.R")

# ==============================================================================
# KROK 4: Benchmark
# ==============================================================================

message("\n========================================")
message("KROK 4: Benchmark na E. coli i UPS")
message("========================================")

source("scripts/04_benchmark.R")

# ==============================================================================
# KROK 5: Generowanie raportu
# ==============================================================================

message("\n========================================")
message("KROK 5: Generowanie raportu HTML")
message("========================================")

rmarkdown::render(
  input = "scripts/05_report.Rmd",
  output_file = "../report.html",
  output_format = "html_document",
  quiet = FALSE
)

message("\n========================================")
message("  ANALIZA ZAKOŃCZONA POMYŚLNIE!        ")
message("========================================")
message("")
message("Wyniki zapisane w:")
message("  - results/figures/ (wykresy)")
message("  - results/tables/ (tabele CSV)")
message("  - report.html (raport końcowy)")
message("")