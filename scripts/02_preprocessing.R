# ==============================================================================
# 02_preprocessing.R - Preprocessing danych proteomicznych
# ==============================================================================

# Załaduj pakiety (muszą być już zainstalowane przez 00_install_packages.R)
# Używamy suppressPackageStartupMessages aby uniknąć konfliktów przy przeładowaniu
suppressPackageStartupMessages({
  library(tidyverse)
  library(vsn)      # vsn MUSI być przed limma (vsn importuje limma)
  library(limma)
  library(imputeLCMD)
})

# ==============================================================================
# 1. WCZYTANIE DANYCH
# ==============================================================================

message("=== Wczytywanie danych ===")

# Znajdź plik proteinGroups.txt
data_dir <- "data"
protein_files <- list.files(data_dir, pattern = "proteinGroups.txt", 
                             recursive = TRUE, full.names = TRUE)

if (length(protein_files) == 0) {
  stop("Brak pliku proteinGroups.txt! Uruchom najpierw 01_download_data.R")
}

# Wybierz plik dla UPS/Dynamic Range Benchmark
# Sprawdź który plik zawiera dane UPS
for (f in protein_files) {
  header <- read_tsv(f, n_max = 0, show_col_types = FALSE)
  cols <- names(header)
  if (any(grepl("UPS|ups", cols))) {
    data_file <- f
    message("Znaleziono dane UPS w: ", f)
    break
  }
}

if (!exists("data_file")) {
  # Jeśli nie znaleziono UPS w nazwie kolumn, użyj pierwszego pliku
  data_file <- protein_files[1]
  message("Używam pliku: ", data_file)
}

# Wczytaj dane
raw_data <- read_tsv(data_file, show_col_types = FALSE)
message("Wczytano ", nrow(raw_data), " białek")
message("Liczba kolumn: ", ncol(raw_data))

# ==============================================================================
# 2. FILTROWANIE
# ==============================================================================

message("\n=== Filtrowanie danych ===")

# Sprawdź dostępne kolumny filtrowania
filter_cols <- c("Reverse", "Potential contaminant", "Only identified by site",
                 "Contaminant", "REV__", "CON__")
available_filters <- intersect(filter_cols, names(raw_data))
message("Dostępne kolumny filtrowania: ", paste(available_filters, collapse = ", "))

# Filtrowanie
filtered_data <- raw_data

# Usuń reverse hits
if ("Reverse" %in% names(filtered_data)) {
  n_before <- nrow(filtered_data)
  filtered_data <- filtered_data %>% filter(is.na(Reverse) | Reverse != "+")
  message("Usunięto ", n_before - nrow(filtered_data), " reverse hits")
}

# Usuń contaminanty
if ("Potential contaminant" %in% names(filtered_data)) {
  n_before <- nrow(filtered_data)
  filtered_data <- filtered_data %>% filter(is.na(`Potential contaminant`) | `Potential contaminant` != "+")
  message("Usunięto ", n_before - nrow(filtered_data), " potencjalnych kontaminantów")
} else if ("Contaminant" %in% names(filtered_data)) {
  n_before <- nrow(filtered_data)
  filtered_data <- filtered_data %>% filter(is.na(Contaminant) | Contaminant != "+")
  message("Usunięto ", n_before - nrow(filtered_data), " kontaminantów")
}

# Usuń "only identified by site"
if ("Only identified by site" %in% names(filtered_data)) {
  n_before <- nrow(filtered_data)
  filtered_data <- filtered_data %>% filter(is.na(`Only identified by site`) | `Only identified by site` != "+")
  message("Usunięto ", n_before - nrow(filtered_data), " białek 'only identified by site'")
}

message("Pozostało ", nrow(filtered_data), " białek po filtrowaniu")

# ==============================================================================
# 3. EKSTRAKCJA INTENSYWNOŚCI LFQ
# ==============================================================================

message("\n=== Ekstrakcja intensywności LFQ ===")

# Znajdź kolumny z intensywnościami LFQ
lfq_cols <- names(filtered_data)[grepl("^LFQ intensity", names(filtered_data))]
message("Znaleziono ", length(lfq_cols), " kolumn LFQ intensity")

if (length(lfq_cols) == 0) {
  # Spróbuj znaleźć kolumny Intensity
  lfq_cols <- names(filtered_data)[grepl("^Intensity ", names(filtered_data))]
  message("Używam kolumn Intensity: ", length(lfq_cols))
}

# Wyświetl nazwy kolumn
message("Kolumny: ", paste(lfq_cols, collapse = ", "))

# Identyfikuj grupy UPS1 i UPS2
# UWAGA: W tym datasecie kolumny nazywają się:
#   "LFQ intensity L1", "L2", "L3" -> Light = UPS1 (niższe stężenie spike-in)
#   "LFQ intensity H1", "H2", "H3" -> Heavy = UPS2 (wyższe stężenie spike-in)
ups1_cols <- lfq_cols[grepl(" L[0-9]", lfq_cols)]  # Light = UPS1
ups2_cols <- lfq_cols[grepl(" H[0-9]", lfq_cols)]  # Heavy = UPS2

# Jeśli nadal puste, spróbuj alternatywnych wzorców
if (length(ups1_cols) == 0) {
  ups1_cols <- lfq_cols[grepl("UPS1|ups1", lfq_cols, ignore.case = TRUE)]
}
if (length(ups2_cols) == 0) {
  ups2_cols <- lfq_cols[grepl("UPS2|ups2", lfq_cols, ignore.case = TRUE)]
}

message("\nKolumny UPS1/Light (", length(ups1_cols), "): ", paste(ups1_cols, collapse = ", "))
message("Kolumny UPS2/Heavy (", length(ups2_cols), "): ", paste(ups2_cols, collapse = ", "))

# Przygotuj macierz intensywności
intensity_matrix <- filtered_data %>%
  select(all_of(c(ups1_cols, ups2_cols))) %>%
  as.matrix()

# Zamień 0 na NA
intensity_matrix[intensity_matrix == 0] <- NA

# Dodaj nazwy wierszy (identyfikatory białek)
if ("Protein IDs" %in% names(filtered_data)) {
  rownames(intensity_matrix) <- filtered_data$`Protein IDs`
} else if ("Majority protein IDs" %in% names(filtered_data)) {
  rownames(intensity_matrix) <- filtered_data$`Majority protein IDs`
}

# Przygotuj metadane białek
protein_info <- filtered_data %>%
  select(any_of(c("Protein IDs", "Majority protein IDs", "Protein names", 
                  "Gene names", "Fasta headers")))

# ==============================================================================
# 4. KLASYFIKACJA BIAŁEK (E. coli vs UPS/Human)
# ==============================================================================

message("\n=== Klasyfikacja białek ===")

# Użyj Fasta headers do klasyfikacji - tam są informacje o gatunku
# Format: >sp|P12345|NAME_SPECIES Description OS=Organism Name
fasta_headers <- filtered_data$`Fasta headers`

# Klasyfikacja na podstawie Fasta headers (lepsze niż same Protein IDs)
protein_species <- case_when(
  grepl("ECOLI|Escherichia coli", fasta_headers, ignore.case = TRUE) ~ "ECOLI",
  grepl("HUMAN|Homo sapiens|UPS", fasta_headers, ignore.case = TRUE) ~ "UPS",
  TRUE ~ "Unknown"
)

message("Klasyfikacja białek:")
message("  E. coli: ", sum(protein_species == "ECOLI"))
message("  UPS (Human): ", sum(protein_species == "UPS"))
message("  Unknown: ", sum(protein_species == "Unknown"))

# Zapisz klasyfikację
protein_info$Species <- protein_species

# ==============================================================================
# 5. TRANSFORMACJA LOG2
# ==============================================================================

message("\n=== Transformacja log2 ===")

log2_matrix <- log2(intensity_matrix)
message("Zakres wartości log2: ", round(min(log2_matrix, na.rm = TRUE), 2), 
        " - ", round(max(log2_matrix, na.rm = TRUE), 2))

# ==============================================================================
# 6. FILTROWANIE "MIN 2 Z 3 POWTÓRZEŃ"
# ==============================================================================

message("\n=== Filtrowanie min 2/3 powtórzeń (lub 3 z 4) ===")

n_ups1 <- length(ups1_cols)
n_ups2 <- length(ups2_cols)

# Funkcja sprawdzająca czy białko ma wystarczająco danych
valid_in_group <- function(row, group_cols, min_valid = 2) {
  sum(!is.na(row[group_cols])) >= min_valid
}

# Oblicz minimalną liczbę (2/3 replik)
min_ups1 <- ceiling(n_ups1 * 2/3)
min_ups2 <- ceiling(n_ups2 * 2/3)

message("Wymagane: min ", min_ups1, " z ", n_ups1, " dla UPS1")
message("Wymagane: min ", min_ups2, " z ", n_ups2, " dla UPS2")

# Sprawdź dla każdego białka
valid_ups1 <- apply(log2_matrix[, ups1_cols, drop = FALSE], 1, 
                    function(x) sum(!is.na(x)) >= min_ups1)
valid_ups2 <- apply(log2_matrix[, ups2_cols, drop = FALSE], 1, 
                    function(x) sum(!is.na(x)) >= min_ups2)

# Białko musi mieć dane w OBIE grupy dla porównania
valid_both <- valid_ups1 & valid_ups2

message("Białek z wystarczającymi danymi w UPS1: ", sum(valid_ups1))
message("Białek z wystarczającymi danymi w UPS2: ", sum(valid_ups2))
message("Białek z danymi w obu grupach: ", sum(valid_both))

# Filtruj
log2_filtered <- log2_matrix[valid_both, ]
protein_info_filtered <- protein_info[valid_both, ]
protein_species_filtered <- protein_species[valid_both]

message("\nPo filtrowaniu pozostało ", nrow(log2_filtered), " białek")
message("  E. coli: ", sum(protein_species_filtered == "ECOLI"))
message("  UPS: ", sum(protein_species_filtered == "UPS"))

# ==============================================================================
# 7. ANALIZA BRAKUJĄCYCH WARTOŚCI
# ==============================================================================

message("\n=== Analiza brakujących wartości ===")

missing_pct <- mean(is.na(log2_filtered)) * 100
message("Procent brakujących wartości: ", round(missing_pct, 1), "%")

# Zapisz wizualizację missingness
missing_per_sample <- colMeans(is.na(log2_filtered)) * 100
missing_df <- data.frame(
  Sample = names(missing_per_sample),
  Missing_pct = missing_per_sample,
  Group = ifelse(grepl("UPS1", names(missing_per_sample)), "UPS1", "UPS2")
)

# Generowanie Missingness Heatmap - PORÓWNANIE PRZED I PO FILTROWANIU
message("Generowanie Missingness Heatmap (przed i po filtrowaniu)...")

# Utwórz katalog na wykresy jeśli nie istnieje
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

# --- PRZED FILTROWANIEM ---
missing_matrix_before <- is.na(log2_matrix) * 1

# Wybierz losową próbkę białek (te same indeksy dla obu wykresów)
set.seed(42)
n_proteins_before <- nrow(missing_matrix_before)
if (n_proteins_before > 100) {
  sample_idx_before <- sample(1:n_proteins_before, 100)
  missing_before_plot <- missing_matrix_before[sample_idx_before, ]
} else {
  missing_before_plot <- missing_matrix_before
}

# --- PO FILTROWANIU ---
missing_matrix_after <- is.na(log2_filtered) * 1

n_proteins_after <- nrow(missing_matrix_after)
if (n_proteins_after > 100) {
  sample_idx_after <- sample(1:n_proteins_after, 100)
  missing_after_plot <- missing_matrix_after[sample_idx_after, ]
} else {
  missing_after_plot <- missing_matrix_after
}

# Przygotuj annotacje kolumn
col_annotation <- data.frame(
  Group = factor(ifelse(grepl("UPS1", colnames(missing_before_plot)), "UPS1", "UPS2"), 
                 levels = c("UPS1", "UPS2")),
  row.names = colnames(missing_before_plot)
)
ann_colors <- list(Group = c("UPS1" = "#E41A1C", "UPS2" = "#377EB8"))  # czerwony/niebieski - spójne z innymi wykresami

# Statystyki do tytułów
pct_missing_before <- round(mean(is.na(log2_matrix)) * 100, 1)
pct_missing_after <- round(mean(is.na(log2_filtered)) * 100, 1)

# Generuj PORÓWNANIE - dwa wykresy obok siebie
png("results/figures/missingness_heatmap.png", width = 1600, height = 800, res = 100)
par(mfrow = c(1, 2))

# Heatmapa 1: PRZED filtrowaniem
p1 <- pheatmap::pheatmap(
  missing_before_plot,
  color = c("white", "#2c3e50"),  # białe = dane obecne, ciemnoszary = NA
  main = paste0("PRZED filtrowaniem\n(", nrow(log2_matrix), " białek, ", pct_missing_before, "% NA)"),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  fontsize_col = 9,
  labels_col = gsub("LFQ intensity ", "", colnames(missing_before_plot)),
  annotation_col = col_annotation,
  annotation_colors = ann_colors,
  legend = FALSE,
  silent = TRUE
)

# Heatmapa 2: PO filtrowaniu
p2 <- pheatmap::pheatmap(
  missing_after_plot,
  color = c("white", "#2c3e50"),  # białe = dane obecne, ciemnoszary = NA
  main = paste0("PO filtrowaniu (min 3/4 replik)\n(", nrow(log2_filtered), " białek, ", pct_missing_after, "% NA)"),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  fontsize_col = 9,
  labels_col = gsub("LFQ intensity ", "", colnames(missing_after_plot)),
  annotation_col = col_annotation,
  annotation_colors = ann_colors,
  legend = TRUE,
  silent = TRUE
)

# Połącz wykresy obok siebie
gridExtra::grid.arrange(p1$gtable, p2$gtable, ncol = 2, 
                         top = grid::textGrob("Porównanie rozkładu brakujących wartości", 
                                               gp = grid::gpar(fontsize = 14, fontface = "bold")))
dev.off()
message("Zapisano: results/figures/missingness_heatmap.png")
message("  Przed filtrowaniem: ", nrow(log2_matrix), " białek, ", pct_missing_before, "% NA")
message("  Po filtrowaniu: ", nrow(log2_filtered), " białek, ", pct_missing_after, "% NA")

# ==============================================================================
# 8. WARIANT 1: MEDIAN NORMALIZATION (BEZ IMPUTACJI)
# ==============================================================================

message("\n=== WARIANT 1: Median normalization (bez imputacji) ===")

# Median normalization na surowych danych
sample_medians <- apply(log2_filtered, 2, median, na.rm = TRUE)
global_median <- median(sample_medians)
norm_factors <- global_median - sample_medians

log2_v1 <- sweep(log2_filtered, 2, norm_factors, "+")

message("Normalizacja zakończona")
message("Mediany przed: ", paste(round(sample_medians, 2), collapse = ", "))
message("Mediany po: ", paste(round(apply(log2_v1, 2, median, na.rm = TRUE), 2), collapse = ", "))

# ==============================================================================
# 9. WARIANT 2: VSN NORMALIZATION (bez imputacji - na kompletnych przypadkach)
# ==============================================================================

message("\n=== WARIANT 2: VSN normalization ===")

# VSN (Variance Stabilizing Normalization) to transformacja glog (generalized log)
# która stabilizuje wariancję i już zawiera transformację logarytmiczną.
# NIE należy robić dodatkowego log2 po VSN!
#
# WAŻNE: W tym benchmarku mamy tylko 1% brakujących wartości.
# Imputacja typu QRILC może bardziej zaszkodzić niż pomóc, bo "ściąga" wartości
# do dolnej granicy i kompresuje różnice. Dlatego tutaj:
# - Używamy VSN na danych z wypełnionymi NA medianą wiersza (minimalna interwencja)
# - NIE używamy agresywnej imputacji QRILC

set.seed(123)  # dla reproducibility

# Przygotuj dane - wypełnij NA medianą wiersza (minimalna interwencja)
log2_for_vsn <- log2_filtered
for (i in 1:nrow(log2_for_vsn)) {
  row_median <- median(log2_for_vsn[i, ], na.rm = TRUE)
  log2_for_vsn[i, is.na(log2_for_vsn[i, ])] <- row_median
}

tryCatch({
  # VSN działa na danych w skali liniowej (nie log)
  # i zwraca dane w skali glog (generalized log) - już zlogarytmowane!
  linear_data <- 2^log2_for_vsn
  log2_v2 <- normalizeVSN(linear_data)
  # NIE robimy log2() - VSN już dał nam dane w skali logarytmicznej (glog)
  
  message("VSN zakończone pomyślnie")
  message("  Uwaga: VSN zwraca dane w skali glog (generalized log)")
  message("  Zakres: ", round(min(log2_v2), 2), " - ", round(max(log2_v2), 2))
}, error = function(e) {
  message("Błąd w VSN: ", e$message)
  message("Używam alternatywnej metody: quantile normalization")
  
  # Alternatywa: quantile normalization (bez imputacji)
  log2_v2 <<- normalizeBetweenArrays(log2_for_vsn, method = "quantile")
})

# ==============================================================================
# 10. ZAPISZ WYNIKI PREPROCESSINGU
# ==============================================================================

message("\n=== Zapisywanie wyników ===")

# Utwórz listę z wynikami
preprocessing_results <- list(
  # Dane surowe
  raw_matrix = intensity_matrix,
  log2_raw = log2_matrix,
  
  # Dane po filtrowaniu
  log2_filtered = log2_filtered,
  
  # Warianty
  log2_v1_median = log2_v1,  # Median norm, bez imputacji
  log2_v2_vsn = log2_v2,  # VSN
  
  # Metadane
  protein_info = protein_info_filtered,
  protein_species = protein_species_filtered,
  
  # Grupy próbek
  ups1_cols = ups1_cols,
  ups2_cols = ups2_cols,
  
  # Statystyki
  missing_stats = missing_df
)

# Zapisz do pliku RDS
saveRDS(preprocessing_results, file = "results/preprocessing_results.rds")
message("Wyniki zapisane do: results/preprocessing_results.rds")

# Zapisz także jako CSV
write_csv(
  bind_cols(protein_info_filtered, as.data.frame(log2_v1)),
  "results/tables/log2_median_normalized.csv"
)

write_csv(
  bind_cols(protein_info_filtered, as.data.frame(log2_v2)),
  "results/tables/log2_vsn.csv"
)

message("\n=== Preprocessing zakończony ===")
message("Wariant 1 (Median): ", nrow(log2_v1), " białek")
message("Wariant 2 (VSN+QRILC): ", nrow(log2_v2), " białek")
