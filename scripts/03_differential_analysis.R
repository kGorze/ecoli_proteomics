# ==============================================================================
# 03_differential_analysis.R - Analiza różnicowa UPS2 vs UPS1
# ==============================================================================

library(tidyverse)
library(limma)
library(pheatmap)
library(ggrepel)
library(RColorBrewer)

# ==============================================================================
# 1. WCZYTANIE DANYCH Z PREPROCESSINGU
# ==============================================================================

message("=== Wczytywanie danych z preprocessingu ===")

if (!file.exists("results/preprocessing_results.rds")) {
  stop("Brak wyników preprocessingu! Uruchom najpierw 02_preprocessing.R")
}

results <- readRDS("results/preprocessing_results.rds")

# Wyciągnij dane
log2_v1 <- results$log2_v1_median
log2_v2 <- results$log2_v2_vsn
protein_info <- results$protein_info
protein_species <- results$protein_species
ups1_cols <- results$ups1_cols
ups2_cols <- results$ups2_cols

message("Wczytano ", nrow(log2_v1), " białek")
message("Grupy: UPS1 (", length(ups1_cols), " próbek), UPS2 (", length(ups2_cols), " próbek)")

# ==============================================================================
# 2. DIAGNOSTYKA BATCH EFFECTS I UKRYTYCH EFEKTÓW
# ==============================================================================

message("\n=== DIAGNOSTYKA: Wykrywanie batch effects i ukrytych efektów ===")

# ---------------------------------------------------------------------------
# 2.1. Analiza struktury próbek z nazw
# ---------------------------------------------------------------------------

message("\n--- 2.1. Analiza nazw próbek pod kątem potencjalnych batch ---")

all_samples <- c(ups1_cols, ups2_cols)
sample_info <- data.frame(
  Sample = all_samples,
  Group = c(rep("UPS1", length(ups1_cols)), rep("UPS2", length(ups2_cols))),
  stringsAsFactors = FALSE
)

# Spróbuj wyciągnąć potencjalne batch z nazw (np. numer repliki)
sample_info$Replicate <- gsub(".*([0-9]+)$", "\\1", sample_info$Sample)
sample_info$Prefix <- gsub("LFQ intensity ", "", sample_info$Sample)
sample_info$Prefix <- gsub("[0-9]+$", "", sample_info$Prefix)

message("Struktura próbek:")
print(sample_info)

# ---------------------------------------------------------------------------
# 2.2. Hierarchiczne klastrowanie próbek
# ---------------------------------------------------------------------------

message("\n--- 2.2. Hierarchiczne klastrowanie próbek ---")

# Użyj log2_v2 (kompletne dane po VSN) do klastrowania
expr_complete <- log2_v2[complete.cases(log2_v2), all_samples]

# Oblicz macierz korelacji
cor_matrix <- cor(expr_complete, method = "pearson")

message("Macierz korelacji między próbkami:")
print(round(cor_matrix, 3))

# Hierarchiczne klastrowanie
hc <- hclust(as.dist(1 - cor_matrix), method = "complete")

# Zapisz dendrogram z lepszą wizualizacją używając dendextend
png("results/figures/diagnostic_dendrogram.png", width = 900, height = 700, res = 150)

# Sprawdź czy dendextend jest dostępny
if (requireNamespace("dendextend", quietly = TRUE)) {
  library(dendextend)
  
  # Stwórz dendrogram
  dend <- as.dendrogram(hc)
  
  # Przypisz kolory do etykiet (UPS1 = czerwony, UPS2 = niebieski)
  sample_labels <- labels(dend)
  label_colors <- ifelse(grepl("UPS1", sample_labels), "#E41A1C", "#377EB8")
  
  # Koloruj etykiety i gałęzie
  dend <- dend %>%
    set("labels_col", label_colors) %>%
    set("labels_cex", 0.8) %>%
    set("branches_k_color", k = 2, value = c("#E41A1C", "#377EB8"))
  
  # Rysuj
  par(mar = c(8, 5, 4, 2))
  plot(dend, 
       main = "Hierarchiczne klastrowanie próbek\n(oparte na korelacji Pearsona)",
       ylab = "Dystans (1 - r)",
       sub = "")
  
  # Dodaj legendę
  legend("topright", 
         legend = c("UPS1", "UPS2"), 
         col = c("#E41A1C", "#377EB8"), 
         pch = 15, 
         cex = 0.9,
         title = "Grupa")
  

  
} else {
  # Fallback do podstawowego wykresu
  par(mar = c(8, 5, 4, 2))
  plot(hc, main = "Hierarchiczne klastrowanie próbek\n(oparte na korelacji Pearsona)",
       xlab = "", ylab = "Dystans (1 - r)",
       sub = "Czy próbki klastrują się według grupy (UPS1 vs UPS2)?")
  rect.hclust(hc, k = 2, border = c("#E41A1C", "#377EB8"))
}

dev.off()
message("Zapisano: results/figures/diagnostic_dendrogram.png")

# Sprawdź czy klastrowanie jest zgodne z grupami
clusters <- cutree(hc, k = 2)
cluster_vs_group <- table(sample_info$Group, clusters[sample_info$Sample])
message("\nKlastrowanie vs grupy biologiczne:")
print(cluster_vs_group)

# Ocena: czy klastry odpowiadają grupom?
cluster_purity <- max(
  sum(diag(cluster_vs_group)) / sum(cluster_vs_group),
  sum(diag(cluster_vs_group[, 2:1])) / sum(cluster_vs_group)
)
if (cluster_purity < 0.8) {
  message("⚠️  UWAGA: Próbki nie klastrują się czysto według grup!")
  message("   Może to wskazywać na batch effect lub inny czynnik zakłócający.")
} else {
  message("✓ Próbki klastrują się zgodnie z grupami biologicznymi.")
}

# ---------------------------------------------------------------------------
# 2.3. Korelacja wewnątrz- vs między-grupowa
# ---------------------------------------------------------------------------

message("\n--- 2.3. Korelacja wewnątrz- vs między-grupowa ---")

# Korelacje wewnątrz grup
ups1_cors <- cor_matrix[ups1_cols, ups1_cols]
ups2_cors <- cor_matrix[ups2_cols, ups2_cols]
within_cors <- c(ups1_cors[upper.tri(ups1_cors)], ups2_cors[upper.tri(ups2_cors)])

# Korelacje między grupami
between_cors <- as.vector(cor_matrix[ups1_cols, ups2_cols])

message("Średnia korelacja WEWNĄTRZ grup: ", round(mean(within_cors), 4))
message("Średnia korelacja MIĘDZY grupami: ", round(mean(between_cors), 4))
message("Różnica (wewnątrz - między): ", round(mean(within_cors) - mean(between_cors), 4))

if (mean(within_cors) - mean(between_cors) < 0.01) {
  message("⚠️  UWAGA: Korelacje wewnątrz i między grupami są podobne!")
  message("   To może wskazywać na słaby efekt biologiczny lub dominujący batch effect.")
}

# ---------------------------------------------------------------------------
# 2.4. Wykrywanie potencjalnych powtórzeń technicznych
# ---------------------------------------------------------------------------

message("\n--- 2.4. Wykrywanie potencjalnych powtórzeń technicznych ---")

# Bardzo wysokie korelacje mogą wskazywać na powtórzenia techniczne
very_high_cor_threshold <- 0.98
high_cor_pairs <- which(cor_matrix > very_high_cor_threshold & cor_matrix < 1, arr.ind = TRUE)

if (nrow(high_cor_pairs) > 0) {
  # Usuń duplikaty (górny trójkąt)
  high_cor_pairs <- high_cor_pairs[high_cor_pairs[,1] < high_cor_pairs[,2], , drop = FALSE]
  
  if (nrow(high_cor_pairs) > 0) {
    message("⚠️  Znaleziono pary próbek o bardzo wysokiej korelacji (>", very_high_cor_threshold, "):")
    for (i in 1:nrow(high_cor_pairs)) {
      s1 <- rownames(cor_matrix)[high_cor_pairs[i, 1]]
      s2 <- colnames(cor_matrix)[high_cor_pairs[i, 2]]
      cor_val <- cor_matrix[high_cor_pairs[i, 1], high_cor_pairs[i, 2]]
      message("   ", s1, " <-> ", s2, ": ", round(cor_val, 4))
    }
    message("   Mogą to być powtórzenia techniczne - rozważ duplicateCorrelation()")
  }
} else {
  message("✓ Nie znaleziono par o podejrzanie wysokiej korelacji.")
}

# ---------------------------------------------------------------------------
# 2.5. PCA - wizualna diagnostyka batch effects
# ---------------------------------------------------------------------------

message("\n--- 2.5. PCA - wizualna diagnostyka ---")

pca_result <- prcomp(t(expr_complete), scale. = TRUE)
var_explained <- summary(pca_result)$importance[2, ] * 100

pca_data <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  PC3 = pca_result$x[, 3],
  Sample = rownames(pca_result$x),
  stringsAsFactors = FALSE
)
pca_data <- merge(pca_data, sample_info, by = "Sample")

# Test: czy PC1 koreluje z grupą biologiczną?
pc1_group_cor <- cor(pca_data$PC1, as.numeric(factor(pca_data$Group)))
message("Korelacja PC1 z grupą biologiczną: ", round(abs(pc1_group_cor), 3))

if (abs(pc1_group_cor) < 0.5) {
  message("⚠️  PC1 słabo koreluje z grupami biologicznymi!")
  message("   To może wskazywać na batch effect dominujący nad efektem biologicznym.")
  message("   Sprawdź, czy PC1 lub PC2 koreluje z innymi czynnikami (dzień, operator, itp.)")
} else {
  message("✓ PC1 silnie koreluje z grupami biologicznymi (oczekiwane).")
}

# Zapisz PCA z adnotacjami
library(ggrepel)
pca_diagnostic <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point(size = 4) +
  geom_text_repel(size = 3) +
  scale_color_manual(values = c("UPS1" = "#E41A1C", "UPS2" = "#377EB8")) +
  labs(
    title = "DIAGNOSTYKA: PCA próbek",
    subtitle = paste0("PC1 wyjaśnia ", round(var_explained[1], 1), "% wariancji | ",
                      "Korelacja PC1-Group: ", round(abs(pc1_group_cor), 2)),
    x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "%)")
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("results/figures/diagnostic_pca.png", pca_diagnostic, width = 10, height = 8, dpi = 150)
message("Zapisano: results/figures/diagnostic_pca.png")

# ---------------------------------------------------------------------------
# 2.6. Wykrywanie ukrytych zmiennych z SVA (Surrogate Variable Analysis)
# ---------------------------------------------------------------------------

message("\n--- 2.6. SVA - Surrogate Variable Analysis ---")

# SVA wykrywa ukryte źródła zmienności (batch effects) z samych danych
# Wymaga pakietu sva

if (requireNamespace("sva", quietly = TRUE)) {
  library(sva)
  
  # Przygotuj dane dla SVA (musi być macierzą bez NA)
  mod <- model.matrix(~ Group, data = sample_info)
  mod0 <- model.matrix(~ 1, data = sample_info)  # null model
  
  # SVA
  tryCatch({
    svobj <- sva(as.matrix(expr_complete), mod, mod0, method = "irw")
    
    n_sv <- svobj$n.sv
    message("SVA znalazło ", n_sv, " surrogate variable(s)")
    
    if (n_sv > 0) {
      message("⚠️  Wykryto ", n_sv, " ukrytych źródeł zmienności (potential batch effects)!")
      message("   Rozważ dodanie surrogate variables do modelu:")
      message("   design <- cbind(design, svobj$sv)")
      
      # Korelacje SV z replikatami
      for (i in 1:n_sv) {
        sv_rep_cor <- cor(svobj$sv[, i], as.numeric(factor(sample_info$Replicate)))
        message("   SV", i, " korelacja z numerem repliki: ", round(sv_rep_cor, 3))
      }
      
      # Zapisz SV dla późniejszego użycia
      sva_results <- list(
        n.sv = n_sv,
        sv = svobj$sv,
        pprob.b = svobj$pprob.b
      )
      saveRDS(sva_results, "results/sva_results.rds")
      message("   Zapisano surrogate variables do: results/sva_results.rds")
    } else {
      message("✓ SVA nie wykryło żadnych ukrytych źródeł zmienności.")
      sva_results <- NULL
    }
    
  }, error = function(e) {
    message("SVA błąd: ", e$message)
    message("   Kontynuuję bez SVA.")
    sva_results <<- NULL
  })
  
} else {
  message("Pakiet 'sva' nie jest zainstalowany.")
  message("Aby włączyć SVA, zainstaluj: BiocManager::install('sva')")
  sva_results <- NULL
}

# ---------------------------------------------------------------------------
# 2.7. Podsumowanie diagnostyki
# ---------------------------------------------------------------------------

message("\n=== PODSUMOWANIE DIAGNOSTYKI ===")
message("1. Klastrowanie próbek: ", ifelse(cluster_purity >= 0.8, "✓ OK", "⚠️ UWAGA"))
message("2. Korelacja wewnątrz/między grupami: różnica = ", 
        round(mean(within_cors) - mean(between_cors), 4))
message("3. Korelacja PC1-Grupa: ", round(abs(pc1_group_cor), 3), 
        ifelse(abs(pc1_group_cor) >= 0.5, " ✓", " ⚠️"))
message("4. Potencjalne powtórzenia techniczne: ", 
        ifelse(nrow(high_cor_pairs) > 0, "TAK", "NIE"))
message("5. SVA ukryte zmienne: ", 
        ifelse(exists("sva_results") && !is.null(sva_results), 
               paste0(sva_results$n.sv, " wykrytych"), "0 lub brak analizy"))

message("\nWykresy diagnostyczne zapisane w results/figures/diagnostic_*.png")
message("=" |> rep(60) |> paste(collapse = ""))

# ==============================================================================
# 3. FUNKCJA DO ANALIZY RÓŻNICOWEJ Z LIMMA
# ==============================================================================

run_limma_analysis <- function(log2_matrix, ups1_cols, ups2_cols, 
                                protein_species, variant_name) {
  
  message("\n=== Analiza różnicowa: ", variant_name, " ===")
  
  # ===========================================================================
  # METODA STATYSTYCZNA:
  # ===========================================================================
    # 1. Model liniowy: limma (Linear Models for Microarray Data)
    # 2. Moderacja wariancji: empirical Bayes (eBayes)
    # 3. Korekta wielokrotnych testów: Benjamini-Hochberg (BH) FDR
    # 4. Kryterium istotności: adj.P.Val < 0.05 AND |log2FC| > 1
  #
  # Design matrix: ~0 + group (bez intercept, grupowe średnie)
  # Contrast: UPS2 - UPS1 (log2FC = log2(UPS2/UPS1))
  # ===========================================================================
  
  # Design matrix
  all_cols <- c(ups1_cols, ups2_cols)
  group <- factor(c(rep("UPS1", length(ups1_cols)), 
                    rep("UPS2", length(ups2_cols))))
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  
  # Contrast: UPS2 - UPS1 (czyli log2FC = log2(UPS2/UPS1))
  contrast_matrix <- makeContrasts(
    UPS2_vs_UPS1 = UPS2 - UPS1,
    levels = design
  )
  
  # Dopasuj macierz
  expr_matrix <- log2_matrix[, all_cols]
  
  # Fit model (limma automatycznie obsługuje NA)
  fit <- lmFit(expr_matrix, design)
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)  # empirical Bayes moderation
  
  # Wyniki z korektą BH FDR
  results_table <- topTable(fit2, coef = "UPS2_vs_UPS1", 
                             number = Inf, sort.by = "none",
                             adjust.method = "BH")  # Benjamini-Hochberg
  
  # Dodaj informacje
  results_table$Protein_ID <- rownames(results_table)
  results_table$Species <- protein_species
  
  # Kryterium istotności: adj.P.Val < 0.05 AND |logFC| > 1
  results_table$Significant <- results_table$adj.P.Val < 0.05 & 
                               abs(results_table$logFC) > 1
  
  message("Test statystyczny: limma + empirical Bayes")
  message("Korekta wielokrotnych testów: Benjamini-Hochberg (BH)")
  message("Kryterium istotności: adj.P.Val < 0.05 AND |logFC| > 1")
  message("Białek istotnych: ", sum(results_table$Significant, na.rm = TRUE))
  message("  - E. coli: ", sum(results_table$Significant & 
                                results_table$Species == "ECOLI", na.rm = TRUE))
  message("  - UPS: ", sum(results_table$Significant & 
                            results_table$Species == "UPS", na.rm = TRUE))
  
  return(results_table)
}

# ==============================================================================
# 3. ANALIZA DLA OBU WARIANTÓW
# ==============================================================================

# Wariant 1: Median normalization (bez imputacji)
results_v1 <- run_limma_analysis(log2_v1, ups1_cols, ups2_cols, 
                                  protein_species, "Median (bez imputacji)")

# Wariant 2: VSN
results_v2 <- run_limma_analysis(log2_v2, ups1_cols, ups2_cols, 
                                  protein_species, "VSN")

# ==============================================================================
# 4. VOLCANO PLOT
# ==============================================================================

message("\n=== Generowanie Volcano Plots ===")

create_volcano_plot <- function(results_table, title) {
  
  # Przygotuj dane
  plot_data <- results_table %>%
    mutate(
      neg_log10_pval = -log10(P.Value),
      Color = case_when(
        Species == "UPS" & Significant ~ "UPS (significant)",
        Species == "UPS" & !Significant ~ "UPS (not significant)",
        Species == "ECOLI" & Significant ~ "E. coli (FALSE POSITIVE)",
        Species == "ECOLI" & !Significant ~ "E. coli (true negative)",
        TRUE ~ "Other"
      )
    )
  
  # Definiuj kolory
  color_palette <- c(
    "UPS (significant)" = "#E41A1C",
    "UPS (not significant)" = "#FB9A99",
    "E. coli (FALSE POSITIVE)" = "#377EB8",
    "E. coli (true negative)" = "#A6CEE3",
    "Other" = "grey50"
  )
  
  # Stwórz wykres
  p <- ggplot(plot_data, aes(x = logFC, y = neg_log10_pval, color = Color)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
    scale_color_manual(values = color_palette) +
    labs(
      title = title,
      subtitle = paste0("UPS2 vs UPS1 | ", 
                        sum(plot_data$Significant, na.rm = TRUE), " significant"),
      x = "log2 Fold Change (UPS2/UPS1)",
      y = "-log10(p-value)",
      color = "Group"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    ) +
    guides(color = guide_legend(nrow = 2))
  
  return(p)
}

# Volcano dla wariantu 1
volcano_v1 <- create_volcano_plot(results_v1, "Volcano Plot - Median Normalization")
ggsave("results/figures/volcano_median.png", volcano_v1, width = 10, height = 8, dpi = 150)

# Volcano dla wariantu 2
volcano_v2 <- create_volcano_plot(results_v2, "Volcano Plot - VSN")
ggsave("results/figures/volcano_vsn.png", volcano_v2, width = 10, height = 8, dpi = 150)

message("Zapisano: results/figures/volcano_median.png")
message("Zapisano: results/figures/volcano_vsn.png")

# ==============================================================================
# 5. PCA PLOT
# ==============================================================================

message("\n=== Generowanie PCA Plots ===")

create_pca_plot <- function(log2_matrix, ups1_cols, ups2_cols, title) {
  
  # Przygotuj dane - usuń białka z brakującymi wartościami
  complete_rows <- complete.cases(log2_matrix)
  pca_matrix <- log2_matrix[complete_rows, c(ups1_cols, ups2_cols)]
  
  # PCA
  pca_result <- prcomp(t(pca_matrix), scale. = TRUE)
  
  # Wyciągnij % wyjaśnionej wariancji
  var_explained <- summary(pca_result)$importance[2, 1:2] * 100
  
  # Dane do wykresu
  pca_data <- data.frame(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    Sample = rownames(pca_result$x),
    Group = c(rep("UPS1", length(ups1_cols)), rep("UPS2", length(ups2_cols)))
  )
  
  # Wykres
  p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
    geom_point(size = 4) +
    geom_text_repel(size = 3, max.overlaps = 20) +
    scale_color_manual(values = c("UPS1" = "#E41A1C", "UPS2" = "#377EB8")) +
    labs(
      title = title,
      x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
      y = paste0("PC2 (", round(var_explained[2], 1), "%)")
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold")
    )
  
  return(p)
}

# PCA dla wariantu 2 (lepsze dla pełnych danych po imputacji)
pca_v2 <- create_pca_plot(log2_v2, ups1_cols, ups2_cols, 
                           "PCA - VSN Normalized Data")
ggsave("results/figures/pca_vsn.png", pca_v2, width = 10, height = 8, dpi = 150)

# Spróbuj PCA dla wariantu 1 (może mieć problemy z NA)
tryCatch({
  pca_v1 <- create_pca_plot(log2_v1, ups1_cols, ups2_cols, 
                             "PCA - Median Normalized Data")
  ggsave("results/figures/pca_median.png", pca_v1, width = 10, height = 8, dpi = 150)
  message("Zapisano: results/figures/pca_median.png")
}, error = function(e) {
  message("PCA dla wariantu 1 niemożliwe z powodu brakujących wartości")
})

message("Zapisano: results/figures/pca_vsn.png")

# ==============================================================================
# 6. HEATMAP TOP BIAŁEK
# ==============================================================================

message("\n=== Generowanie Heatmap ===")

create_heatmap <- function(log2_matrix, results_table, ups1_cols, ups2_cols, 
                           title, filename, top_n = 50) {
  
  # Wybierz top białek według p-value
  top_proteins <- results_table %>%
    arrange(P.Value) %>%
    head(top_n) %>%
    pull(Protein_ID)
  
  # Przygotuj dane
  heatmap_data <- log2_matrix[top_proteins, c(ups1_cols, ups2_cols)]
  
  # Dla NA zamień na medianę wiersza
  for (i in 1:nrow(heatmap_data)) {
    row_median <- median(heatmap_data[i, ], na.rm = TRUE)
    heatmap_data[i, is.na(heatmap_data[i, ])] <- row_median
  }
  
  # Skaluj wiersze
  heatmap_scaled <- t(scale(t(heatmap_data)))
  
  # Annotacje
  annotation_col <- data.frame(
    Group = c(rep("UPS1", length(ups1_cols)), rep("UPS2", length(ups2_cols))),
    row.names = c(ups1_cols, ups2_cols)
  )
  
  annotation_row <- data.frame(
    Species = results_table$Species[match(top_proteins, results_table$Protein_ID)],
    row.names = top_proteins
  )
  
  # Kolory
  annotation_colors <- list(
    Group = c(UPS1 = "#E41A1C", UPS2 = "#377EB8"),
    Species = c(ECOLI = "#A6CEE3", UPS = "#FB9A99", Unknown = "grey50")
  )
  
  # Zapisz heatmapę
  png(filename, width = 1200, height = 1000, res = 150)
  pheatmap(
    heatmap_scaled,
    main = title,
    annotation_col = annotation_col,
    annotation_row = annotation_row,
    annotation_colors = annotation_colors,
    cluster_cols = FALSE,
    show_rownames = FALSE,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    border_color = NA
  )
  dev.off()
}

# Heatmap dla wariantu 1 (Median)
create_heatmap(log2_v1, results_v1, ups1_cols, ups2_cols,
               "Top 50 Differential Proteins - Median",
               "results/figures/heatmap_median.png")
message("Zapisano: results/figures/heatmap_median.png")

# Heatmap dla wariantu 2 (VSN)
create_heatmap(log2_v2, results_v2, ups1_cols, ups2_cols,
               "Top 50 Differential Proteins - VSN",
               "results/figures/heatmap_vsn.png")
message("Zapisano: results/figures/heatmap_vsn.png")

# --- PORÓWNANIE SIDE-BY-SIDE ---
message("Tworzenie porównania heatmap Median vs VSN...")

# Użyj magick do połączenia obrazów jeśli dostępny
tryCatch({
  library(magick)
  
  img1 <- image_read("results/figures/heatmap_median.png")
  img2 <- image_read("results/figures/heatmap_vsn.png")
  
  # Połącz obok siebie
  combined <- image_append(c(img1, img2), stack = FALSE)
  
  # Dodaj tytuł
  combined <- image_annotate(combined, "Porównanie Heatmap: Median vs VSN", 
                              size = 40, gravity = "north", 
                              color = "black", weight = 700)
  
  image_write(combined, "results/figures/comparison_heatmap.png")
  message("Zapisano: results/figures/comparison_heatmap.png")
  
}, error = function(e) {
  message("Pakiet magick niedostępny - zachowano osobne heatmapy")
  message("Instalacja: install.packages('magick')")
})

# ==============================================================================
# 7. ZAPISZ WYNIKI
# ==============================================================================

message("\n=== Zapisywanie wyników analizy różnicowej ===")

# Zapisz tabele wyników
write_csv(results_v1, "results/tables/differential_results_median.csv")
write_csv(results_v2, "results/tables/differential_results_vsn.csv")

# Zapisz do RDS dla dalszych analiz
differential_results <- list(
  results_v1 = results_v1,
  results_v2 = results_v2,
  log2_v1 = log2_v1,
  log2_v2 = log2_v2,
  protein_info = protein_info,
  protein_species = protein_species,
  ups1_cols = ups1_cols,
  ups2_cols = ups2_cols
)

saveRDS(differential_results, "results/differential_results.rds")

message("Zapisano: results/tables/differential_results_median.csv")
message("Zapisano: results/tables/differential_results_vsn.csv")
message("Zapisano: results/differential_results.rds")

message("\n=== Analiza różnicowa zakończona ===")
