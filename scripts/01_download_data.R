# ==============================================================================
# 01_download_data.R - Pobieranie danych z PRIDE PXD000279
# ==============================================================================
# 
# UWAGA: W PXD000279 są DWA różne benchmarki:
#   1. proteomebenchmark.zip - HeLa + E. coli (E. coli zmienia się 3:1, human stałe)
#   2. dynamicrangebenchmark.zip - UPS1/UPS2 w tle E. coli (E. coli stałe, UPS zmienne)
#
# Ten skrypt pobiera DYNAMIC RANGE BENCHMARK (UPS1/UPS2)
# ==============================================================================

# Ścieżki
data_dir <- "data"
zip_file <- file.path(data_dir, "dynamicrangebenchmark.zip")

# URL do Dynamic Range Benchmark
ftp_urls <- c(
  "https://ftp.pride.ebi.ac.uk/pride/data/archive/2014/09/PXD000279/dynamicrangebenchmark.zip",
  "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2014/09/PXD000279/dynamicrangebenchmark.zip"
)

# Sprawdź czy katalog data istnieje
if (!dir.exists(data_dir)) {
  dir.create(data_dir, recursive = TRUE)
}

# Pobierz dane jeśli nie istnieją lub są uszkodzone
if (!file.exists(zip_file) || file.size(zip_file) < 1000000) {  # < 1MB = uszkodzony
  
  # Usuń stare pliki proteomebenchmark jeśli istnieją
  old_files <- c(
    file.path(data_dir, "proteomebenchmark.zip"),
    file.path(data_dir, "proteinGroups.txt"),
    file.path(data_dir, "experimentalDesign.txt"),
    file.path(data_dir, "peptides.txt"),
    file.path(data_dir, "modificationSpecificPeptides.txt"),
    file.path(data_dir, "parameters.txt")
  )
  for (f in old_files) {
    if (file.exists(f)) {
      message("Usuwam stary plik: ", basename(f))
      file.remove(f)
    }
  }
  
  if (file.exists(zip_file)) {
    message("Usuwam uszkodzony plik...")
    file.remove(zip_file)
  }
  
  message("Pobieranie danych Dynamic Range Benchmark z PRIDE FTP...")
  message("(UPS1/UPS2 spike-in w tle E. coli)")
  message("UWAGA: Może potrwać kilka minut.")
  
  # Zwiększ timeout
  old_timeout <- getOption("timeout")
  options(timeout = 600)  # 10 minut
  
  download_success <- FALSE
  
  for (url in ftp_urls) {
    message("\nPróbuję: ", url)
    
    tryCatch({
      download.file(
        url = url,
        destfile = zip_file,
        mode = "wb",
        method = "auto",
        quiet = FALSE
      )
      
      # Sprawdź czy pobrany plik jest wystarczająco duży
      if (file.exists(zip_file) && file.size(zip_file) > 1000000) {
        message("\nPobrano pomyślnie: ", zip_file)
        download_success <- TRUE
        break
      } else {
        message("Plik za mały, próbuję inny URL...")
        if (file.exists(zip_file)) file.remove(zip_file)
      }
      
    }, error = function(e) {
      message("Błąd: ", e$message)
    })
  }
  
  # Przywróć stary timeout
  options(timeout = old_timeout)
  
  if (!download_success) {
    stop(
      "\n\nNie udało się pobrać danych automatycznie.",
      "\n\nPobierz ręcznie:",
      "\n1. Otwórz w przeglądarce: https://www.ebi.ac.uk/pride/archive/projects/PXD000279",
      "\n2. Pobierz plik: dynamicrangebenchmark.zip",
      "\n3. Umieść go w folderze: ", normalizePath(data_dir),
      "\n4. Uruchom skrypt ponownie"
    )
  }
  
} else {
  message("Plik już istnieje: ", zip_file)
  message("Rozmiar: ", round(file.size(zip_file) / 1024 / 1024, 1), " MB")
}

# Rozpakuj archiwum
message("\nRozpakowanie archiwum...")
unzip(zip_file, exdir = data_dir, overwrite = TRUE)

# Znajdź plik proteinGroups.txt
protein_files <- list.files(
  data_dir, 
  pattern = "proteinGroups.txt", 
  recursive = TRUE, 
  full.names = TRUE
)

if (length(protein_files) == 0) {
  stop("Nie znaleziono pliku proteinGroups.txt!")
}

message("\nZnalezione pliki proteinGroups.txt:")
for (f in protein_files) {
  message("  - ", f, " (", round(file.size(f) / 1024 / 1024, 1), " MB)")
}

# Sprawdź experimental design
exp_design_files <- list.files(
  data_dir,
  pattern = "experimentalDesign|summary",
  recursive = TRUE,
  full.names = TRUE
)

if (length(exp_design_files) > 0) {
  message("\nPliki z designem eksperymentu:")
  for (f in exp_design_files) {
    message("  - ", f)
  }
}

message("\n=== Pobieranie Dynamic Range Benchmark zakończone ===")
