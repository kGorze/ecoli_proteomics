# ==============================================================================
# 00_install_packages.R - Instalacja wszystkich wymaganych pakietów
# ==============================================================================
# 
# PROBLEM: RStudio automatycznie restartuje R przy install.packages()
#          co powoduje utratę zmiennych.
#
# ROZWIĄZANIE: Instaluj każdy pakiet osobno, sprawdzając czy jest zainstalowany.
#              Uruchom ten skrypt WIELOKROTNIE aż wszystko się zainstaluje.
#
# ==============================================================================

message("==============================================")
message("  Instalacja pakietów dla analizy proteomiki  ")
message("==============================================\n")

# Wymuś użycie CRAN bez pytań
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Wyłącz kompilację ze źródeł (szybsza instalacja)
options(install.packages.compile.from.source = "never")

# ==============================================================================
# Definicja pakietów (musi być na początku żeby przetrwać restart)
# ==============================================================================

cran_packages <- c(
  "BiocManager",
  "tidyverse",
  "pheatmap", 
  "ggrepel", 
  "RColorBrewer", 
  "scales",
  "knitr", 
  "kableExtra", 
  "rmarkdown", 
  "patchwork",
  "pROC"  # dla ROC/AUC analysis
)

bioc_packages <- c("limma", "vsn", "imputeLCMD", "sva")  # sva dla wykrywania batch effects

# ==============================================================================
# Funkcja do sprawdzania i instalacji CRAN
# ==============================================================================

install_cran_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(">>> Instaluję ", pkg, "...")
    install.packages(pkg, dependencies = TRUE, ask = FALSE)
    return(TRUE)  # Zainstalowano
  } else {
    message("    ", pkg, " - OK (już zainstalowany)")
    return(FALSE)  # Już był
  }
}

# ==============================================================================
# Funkcja do sprawdzania i instalacji Bioconductor
# ==============================================================================

install_bioc_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(">>> Instaluję ", pkg, " (Bioconductor)...")
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
    return(TRUE)
  } else {
    message("    ", pkg, " - OK (już zainstalowany)")
    return(FALSE)
  }
}

# ==============================================================================
# GŁÓWNA LOGIKA - instaluj pierwszy brakujący pakiet
# ==============================================================================

message("=== Sprawdzam pakiety CRAN ===\n")

installed_something <- FALSE

# Najpierw sprawdź CRAN
for (pkg in cran_packages) {
  tryCatch({
    if (install_cran_if_missing(pkg)) {
      installed_something <- TRUE
      # Po instalacji przerwij - RStudio może zrestartować
      message("\n>>> Zainstalowano ", pkg)
      message(">>> Uruchom skrypt ponownie żeby kontynuować!\n")
      break
    }
  }, error = function(e) {
    message("!!! Błąd przy ", pkg, ": ", e$message)
  })
}

# Jeśli CRAN OK, sprawdź Bioconductor
if (!installed_something) {
  message("\n=== Sprawdzam pakiety Bioconductor ===\n")
  
  for (pkg in bioc_packages) {
    tryCatch({
      if (install_bioc_if_missing(pkg)) {
        installed_something <- TRUE
        message("\n>>> Zainstalowano ", pkg)
        message(">>> Uruchom skrypt ponownie żeby kontynuować!\n")
        break
      }
    }, error = function(e) {
      message("!!! Błąd przy ", pkg, ": ", e$message)
    })
  }
}

# ==============================================================================
# Podsumowanie
# ==============================================================================

if (!installed_something) {
  message("\n==============================================")
  message("  WSZYSTKIE PAKIETY SĄ ZAINSTALOWANE!")
  message("==============================================")
  message("")
  message("Możesz teraz uruchomić analizę:")
  message("  source('run_analysis.R')")
  message("")
} else {
  message("==============================================")
  message("  Niektóre pakiety wymagają jeszcze instalacji")
  message("==============================================")
  message("")
  message("Uruchom ten skrypt ponownie:")
  message("  source('scripts/00_install_packages.R')")
  message("")
}
