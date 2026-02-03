# ==============================================================================
# 04_benchmark.R - Ocena jakości pipeline'u na E. coli i UPS
# ==============================================================================

library(tidyverse)
library(ggrepel)

# Sprawdź czy pROC jest dostępny (dla ROC/AUC analysis)
if (requireNamespace("pROC", quietly = TRUE)) {
  library(pROC)
  has_pROC <- TRUE
} else {
  message("Pakiet pROC nie jest zainstalowany - pomijam analizę ROC/AUC")
  has_pROC <- FALSE
}

# ==============================================================================
# 1. WCZYTANIE WYNIKÓW
# ==============================================================================

message("=== Wczytywanie wyników analizy różnicowej ===")

if (!file.exists("results/differential_results.rds")) {
  stop("Brak wyników! Uruchom najpierw 03_differential_analysis.R")
}

diff_results <- readRDS("results/differential_results.rds")

results_v1 <- diff_results$results_v1
results_v2 <- diff_results$results_v2

# ==============================================================================
# 2. OCZEKIWANE FOLD-CHANGE DLA UPS
# ==============================================================================

message("\n=== Definiowanie oczekiwanych fold-change dla UPS ===")

# UPS Dynamic Range Standard (Sigma-Aldrich):
# UPS1: wszystkie białka w stężeniu 5,000 fmol
# UPS2: 6 grup stężeń (48 białek, po 8 w każdej grupie)
# 
# log2FC = log2(UPS2_fmol / UPS1_fmol)
#
# Grupa A: 50,000 fmol -> log2(50000/5000) = log2(10) = 3.32
# Grupa B: 5,000 fmol  -> log2(5000/5000)  = log2(1)  = 0.00
# Grupa C: 500 fmol    -> log2(500/5000)   = log2(0.1) = -3.32
# Grupa D: 50 fmol     -> log2(50/5000)    = log2(0.01) = -5.64
# Grupa E: 5 fmol      -> log2(5/5000)     = log2(0.001) = -9.97
# Grupa F: 0.5 fmol    -> log2(0.5/5000)   = log2(0.0001) = -13.29

# Tabela oczekiwanych FC (Sigma-Aldrich UPS1/UPS2 Product Information)
ups_expected_fc <- tribble(
  ~Protein, ~Group, ~UPS2_fmol, ~UPS1_fmol, ~Expected_log2FC,
  
  # Grupa A: 50,000 fmol (8 białek)
  "P00915", "A", 50000, 5000, log2(50000/5000),  # Carbonic anhydrase 1
  "P00918", "A", 50000, 5000, log2(50000/5000),  # Carbonic anhydrase 2
  "P01031", "A", 50000, 5000, log2(50000/5000),  # Complement C5
  "P69905", "A", 50000, 5000, log2(50000/5000),  # Hemoglobin alpha
  "P68871", "A", 50000, 5000, log2(50000/5000),  # Hemoglobin beta
  "P41159", "A", 50000, 5000, log2(50000/5000),  # Leptin
  "P02768", "A", 50000, 5000, log2(50000/5000),  # Serum Albumin
  "P62988", "A", 50000, 5000, log2(50000/5000),  # Ubiquitin
  
  # Grupa B: 5,000 fmol (8 białek) - brak zmiany!
  "P04040", "B", 5000, 5000, log2(5000/5000),    # Catalase
  "P00167", "B", 5000, 5000, log2(5000/5000),    # Cytochrome b5
  "P01133", "B", 5000, 5000, log2(5000/5000),    # EGF
  "P02144", "B", 5000, 5000, log2(5000/5000),    # Myoglobin
  "P15559", "B", 5000, 5000, log2(5000/5000),    # NQO1
  "P62937", "B", 5000, 5000, log2(5000/5000),    # Cyclophilin A
  "Q06830", "B", 5000, 5000, log2(5000/5000),    # Peroxiredoxin 1
  "P63165", "B", 5000, 5000, log2(5000/5000),    # SUMO-1
  
  # Grupa C: 500 fmol (8 białek)
  "P00709", "C", 500, 5000, log2(500/5000),      # Alpha-lactalbumin
  "P06732", "C", 500, 5000, log2(500/5000),      # Creatine kinase M
  "P12081", "C", 500, 5000, log2(500/5000),      # Histidyl-tRNA synthetase
  "P61626", "C", 500, 5000, log2(500/5000),      # Lysozyme C
  "Q15843", "C", 500, 5000, log2(500/5000),      # Neddylin
  "P02753", "C", 500, 5000, log2(500/5000),      # Retinol-binding protein
  "P16083", "C", 500, 5000, log2(500/5000),      # NQO2
  "P63279", "C", 500, 5000, log2(500/5000),      # UbcH9
  
  # Grupa D: 50 fmol (8 białek)
  "P01008", "D", 50, 5000, log2(50/5000),        # Antithrombin-III
  "P61769", "D", 50, 5000, log2(50/5000),        # Beta-2-microglobulin
  "P55957", "D", 50, 5000, log2(50/5000),        # BID
  "O76070", "D", 50, 5000, log2(50/5000),        # Gamma-synuclein
  "P08263", "D", 50, 5000, log2(50/5000),        # GST A1
  "P01344", "D", 50, 5000, log2(50/5000),        # IGF-II
  "P01127", "D", 50, 5000, log2(50/5000),        # PDGF B
  "P10599", "D", 50, 5000, log2(50/5000),        # Thioredoxin
  
  # Grupa E: 5 fmol (8 białek)
  "P99999", "E", 5, 5000, log2(5/5000),          # Cytochrome c
  "P06396", "E", 5, 5000, log2(5/5000),          # Gelsolin
  "P09211", "E", 5, 5000, log2(5/5000),          # GST P
  "P01112", "E", 5, 5000, log2(5/5000),          # HRas
  "P01579", "E", 5, 5000, log2(5/5000),          # IFN-gamma
  "P02787", "E", 5, 5000, log2(5/5000),          # Serotransferrin
  "O00762", "E", 5, 5000, log2(5/5000),          # UbcH10
  "P51965", "E", 5, 5000, log2(5/5000),          # UbcH6
  
  # Grupa F: 0.5 fmol (8 białek) - najniższe stężenie
  "P08758", "F", 0.5, 5000, log2(0.5/5000),      # Annexin A5
  "P02741", "F", 0.5, 5000, log2(0.5/5000),      # CRP
  "P05413", "F", 0.5, 5000, log2(0.5/5000),      # FABP
  "P10145", "F", 0.5, 5000, log2(0.5/5000),      # IL-8
  "P02788", "F", 0.5, 5000, log2(0.5/5000),      # Lactotransferrin
  "P10636", "F", 0.5, 5000, log2(0.5/5000),      # Tau protein
  "P00441", "F", 0.5, 5000, log2(0.5/5000),      # SOD1
  "P01375", "F", 0.5, 5000, log2(0.5/5000)       # TNF-alpha
)

message("Zdefiniowano ", nrow(ups_expected_fc), " białek UPS z oczekiwanymi FC")
message("Grupy stężeń (log2FC):")
message("  A: 50,000 fmol -> log2FC = ", round(log2(10), 2))
message("  B: 5,000 fmol  -> log2FC = ", round(log2(1), 2), " (brak zmiany!)")
message("  C: 500 fmol    -> log2FC = ", round(log2(0.1), 2))
message("  D: 50 fmol     -> log2FC = ", round(log2(0.01), 2))
message("  E: 5 fmol      -> log2FC = ", round(log2(0.001), 2))
message("  F: 0.5 fmol    -> log2FC = ", round(log2(0.0001), 2))

# ==============================================================================
# 3. FUNKCJA BENCHMARK
# ==============================================================================

run_benchmark <- function(results_table, ups_expected_fc, variant_name) {
  
  message("\n=== BENCHMARK: ", variant_name, " ===")
  
  # -------------------------
  # A) EMPIRYCZNY FDR NA E. COLI
  # -------------------------
  # DEFINICJA:
  # W tym benchmarku E. coli jest tłem 1:1 (stałe między UPS1 a UPS2).
  # Każde białko E. coli oznaczone jako "istotne" jest false positive.
  # 
  # Empiryczny FDR = (liczba E. coli z adj.P.Val < 0.05 AND |logFC| > 1) / (całkowita liczba E. coli)
  #
  # To NIE jest standardowy FDR (odsetek FP wśród wszystkich pozytywów),
  # ale odsetek E. coli błędnie uznanych za różnicowe.
  # -------------------------
  
  ecoli_results <- results_table %>%
    filter(Species == "ECOLI")
  
  n_ecoli_total <- nrow(ecoli_results)
  n_ecoli_significant <- sum(ecoli_results$Significant, na.rm = TRUE)
  
  # Empiryczny FDR = % E. coli białek oznaczonych jako istotne
  # (E. coli powinny być niezmienne między UPS1 a UPS2)
  empirical_fdr <- n_ecoli_significant / n_ecoli_total * 100
  
  message("\n--- Empiryczny FDR na E. coli ---")
  message("DEFINICJA: Odsetek białek E. coli błędnie uznanych za różnicowe")
  message("  (kryterium: adj.P.Val < 0.05 AND |logFC| > 1)")
  message("Całkowita liczba białek E. coli: ", n_ecoli_total)
  message("E. coli oznaczone jako istotne (FP): ", n_ecoli_significant)
  message("Empiryczny FDR: ", round(empirical_fdr, 2), "%")
  
  # Również policz przy różnych progach adj.P.Val
  fdr_thresholds <- c(0.01, 0.05, 0.1, 0.2)
  fdr_at_thresholds <- sapply(fdr_thresholds, function(thr) {
    sum(ecoli_results$adj.P.Val < thr, na.rm = TRUE) / n_ecoli_total * 100
  })
  
  fdr_summary <- data.frame(
    Threshold = fdr_thresholds,
    N_FalsePositive = sapply(fdr_thresholds, function(thr) {
      sum(ecoli_results$adj.P.Val < thr, na.rm = TRUE)
    }),
    Empirical_FDR_pct = fdr_at_thresholds
  )
  
  message("\nEmZpiryczny FDR przy różnych progach:")
  print(fdr_summary)
  
  # -------------------------
  # B) DOKŁADNOŚĆ FOLD-CHANGE DLA UPS
  # -------------------------
  
  message("\n--- Dokładność Fold-Change dla UPS ---")
  
  # Dopasuj białka UPS z wyników do oczekiwanych
  ups_results <- results_table %>%
    filter(Species == "UPS")
  
  # Wyciągnij identyfikatory UniProt z protein IDs
  # Typowy format: "sp|P02768|ALBU_HUMAN" lub "P02768"
  extract_uniprot <- function(protein_id) {
    # Szukaj wzorca P/Q/O + 5 cyfr/liter
    matches <- str_extract_all(protein_id, "[PQO][0-9][A-Z0-9]{4}")[[1]]
    if (length(matches) > 0) {
      return(matches[1])
    }
    return(NA_character_)
  }
  
  ups_results$UniProt <- sapply(ups_results$Protein_ID, extract_uniprot)
  
  # Dopasuj do oczekiwanych FC
  ups_comparison <- ups_results %>%
    inner_join(ups_expected_fc, by = c("UniProt" = "Protein")) %>%
    select(Protein_ID, UniProt, Group, logFC, Expected_log2FC, P.Value, adj.P.Val)
  
  # RAPORT O DOSTĘPNOŚCI UPS
  message("\n--- Dostępność białek UPS ---")
  message("Białka UPS w danych: ", nrow(ups_results))
  message("Dopasowano do oczekiwanych FC: ", nrow(ups_comparison), " z 48 teoretycznych")
  
  # Sprawdź które grupy są obecne
  expected_groups <- table(ups_expected_fc$Group)
  if (nrow(ups_comparison) > 0) {
    observed_groups <- table(ups_comparison$Group)
    message("\nGrupy UPS (oczekiwane -> obecne):")
    for (g in names(expected_groups)) {
      n_exp <- expected_groups[g]
      n_obs <- ifelse(g %in% names(observed_groups), observed_groups[g], 0)
      status <- ifelse(n_obs > 0, "✓", "✗ BRAK")
      message("  Grupa ", g, " (", ups_expected_fc$UPS2_fmol[ups_expected_fc$Group == g][1], 
              " fmol): ", n_obs, "/", n_exp, " ", status)
    }
  }
  
  if (nrow(ups_comparison) > 0) {
    # Oblicz metryki
    correlation <- cor(ups_comparison$logFC, ups_comparison$Expected_log2FC, 
                       use = "complete.obs")
    rmse <- sqrt(mean((ups_comparison$logFC - ups_comparison$Expected_log2FC)^2, 
                       na.rm = TRUE))
    mae <- mean(abs(ups_comparison$logFC - ups_comparison$Expected_log2FC), 
                 na.rm = TRUE)
    
    message("\n--- Metryki dokładności fold-change ---")
    message("UWAGA: Niskie wartości mogą wynikać z:")
    message("  - Filtrowaniu (grupy D i E o niskim stężeniu często nie są wykrywane)")
    message("  - Kompresji sygnału przez normalizację")
    message("  - Błędów w mapowaniu identyfikatorów UniProt")
    message("")
    message("  Korelacja Pearsona (expected vs observed): ", round(correlation, 3))
    message("  RMSE (log2 scale): ", round(rmse, 3))
    message("  MAE (log2 scale): ", round(mae, 3))
    
    # Metryki per grupa
    group_metrics <- ups_comparison %>%
      group_by(Group) %>%
      summarise(
        N = n(),
        Expected_FC = first(Expected_log2FC),
        Mean_Observed_FC = mean(logFC, na.rm = TRUE),
        SD_Observed_FC = sd(logFC, na.rm = TRUE),
        Bias = mean(logFC - Expected_log2FC, na.rm = TRUE),
        .groups = "drop"
      )
    
    message("\nMetryki per grupa (KLUCZOWE dla interpretacji):")
    print(group_metrics)
    
  } else {
    correlation <- NA
    rmse <- NA
    mae <- NA
    ups_comparison <- NULL
    message("BŁĄD: Nie udało się dopasować białek UPS!")
    message("  Sprawdź czy identyfikatory UniProt są poprawnie wyodrębniane.")
  }
  
  # Zwróć wyniki
  return(list(
    variant = variant_name,
    n_ecoli = n_ecoli_total,
    n_ecoli_fp = n_ecoli_significant,
    empirical_fdr = empirical_fdr,
    fdr_summary = fdr_summary,
    n_ups_matched = nrow(ups_comparison),
    correlation = correlation,
    rmse = rmse,
    mae = mae,
    ups_comparison = ups_comparison,
    ecoli_results = ecoli_results
  ))
}

# ==============================================================================
# 4. BENCHMARK DLA OBU WARIANTÓW
# ==============================================================================

benchmark_v1 <- run_benchmark(results_v1, ups_expected_fc, "Median (bez imputacji)")
benchmark_v2 <- run_benchmark(results_v2, ups_expected_fc, "VSN")

# ==============================================================================
# 5. PORÓWNANIE WARIANTÓW
# ==============================================================================

message("\n\n========================================")
message("=== PODSUMOWANIE PORÓWNANIA WARIANTÓW ===")
message("========================================\n")

comparison_summary <- data.frame(
  Wariant = c("Median (bez imputacji)", "VSN"),
  N_Ecoli = c(benchmark_v1$n_ecoli, benchmark_v2$n_ecoli),
  N_Ecoli_FP = c(benchmark_v1$n_ecoli_fp, benchmark_v2$n_ecoli_fp),
  Empirical_FDR_pct = c(benchmark_v1$empirical_fdr, benchmark_v2$empirical_fdr),
  N_UPS_Matched = c(benchmark_v1$n_ups_matched, benchmark_v2$n_ups_matched),
  Correlation = c(benchmark_v1$correlation, benchmark_v2$correlation),
  RMSE = c(benchmark_v1$rmse, benchmark_v2$rmse)
)

print(comparison_summary)

# Który wariant lepszy?
message("\n--- REKOMENDACJA ---")

# Logika wniosków
fdr_diff <- abs(benchmark_v1$empirical_fdr - benchmark_v2$empirical_fdr)

if (fdr_diff < 0.1) {
  # Oba FDR są praktycznie takie same
  message("Oba warianty mają podobny empiryczny FDR (~", 
          round(benchmark_v1$empirical_fdr, 1), "%).")
  message("W tym przypadku rozstrzygające powinny być inne kryteria.")
} else if (benchmark_v1$empirical_fdr < benchmark_v2$empirical_fdr) {
  message("Wariant 1 (Median) ma niższy empiryczny FDR: ", 
          round(benchmark_v1$empirical_fdr, 2), "% vs ", 
          round(benchmark_v2$empirical_fdr, 2), "%")
} else {
  message("Wariant 2 (VSN) ma niższy empiryczny FDR: ", 
          round(benchmark_v2$empirical_fdr, 2), "% vs ", 
          round(benchmark_v1$empirical_fdr, 2), "%")
}

if (!is.na(benchmark_v1$rmse) && !is.na(benchmark_v2$rmse)) {
  if (benchmark_v1$rmse < benchmark_v2$rmse) {
    message("Wariant 1 (Median) ma lepszą dokładność FC (niższe RMSE): ", 
            round(benchmark_v1$rmse, 2), " vs ", round(benchmark_v2$rmse, 2))
  } else {
    message("Wariant 2 (VSN) ma lepszą dokładność FC (niższe RMSE): ", 
            round(benchmark_v2$rmse, 2), " vs ", round(benchmark_v1$rmse, 2))
  }
}

# ==============================================================================
# 6. WYKRES EXPECTED VS OBSERVED
# ==============================================================================

message("\n=== Generowanie wykresu Expected vs Observed FC ===")

create_expected_vs_observed_plot <- function(benchmark_result, title) {
  
  if (is.null(benchmark_result$ups_comparison) || 
      nrow(benchmark_result$ups_comparison) == 0) {
    message("Brak danych UPS dla wykresu: ", title)
    return(NULL)
  }
  
  plot_data <- benchmark_result$ups_comparison
  
  # Kolory per grupa
  group_colors <- c(A = "#E41A1C", B = "#377EB8", C = "#4DAF4A", 
                    D = "#984EA3", E = "#FF7F00")
  
  p <- ggplot(plot_data, aes(x = Expected_log2FC, y = logFC, color = Group)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
    geom_point(size = 3, alpha = 0.7) +
    geom_text_repel(aes(label = UniProt), size = 2.5, max.overlaps = 15) +
    scale_color_manual(values = group_colors) +
    labs(
      title = title,
      subtitle = paste0("r = ", round(benchmark_result$correlation, 3), 
                        " | RMSE = ", round(benchmark_result$rmse, 2)),
      x = "Expected log2 Fold Change",
      y = "Observed log2 Fold Change",
      color = "UPS Group"
    ) +
    coord_fixed() +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )
  
  return(p)
}

# Wykres dla wariantu 1
if (!is.null(benchmark_v1$ups_comparison)) {
  p1 <- create_expected_vs_observed_plot(benchmark_v1, 
                                          "Expected vs Observed FC - Median")
  if (!is.null(p1)) {
    ggsave("results/figures/expected_vs_observed_median.png", p1, 
           width = 10, height = 8, dpi = 150)
    message("Zapisano: results/figures/expected_vs_observed_median.png")
  }
}

# Wykres dla wariantu 2
if (!is.null(benchmark_v2$ups_comparison)) {
  p2 <- create_expected_vs_observed_plot(benchmark_v2, 
                                          "Expected vs Observed FC - VSN")
  if (!is.null(p2)) {
    ggsave("results/figures/expected_vs_observed_vsn.png", p2, 
           width = 10, height = 8, dpi = 150)
    message("Zapisano: results/figures/expected_vs_observed_vsn.png")
  }
}

# ==============================================================================
# 7. ANALIZA ROC/AUC I PRECISION-RECALL
# ==============================================================================
# 
# W tym benchmarku:
# - TRUE POSITIVES: białka UPS które są wykryte jako istotne (powinny być różnicowe)
# - TRUE NEGATIVES: białka E. coli które NIE są istotne (powinny być stałe)
# - FALSE POSITIVES: białka E. coli wykryte jako istotne (błąd!)
# - FALSE NEGATIVES: białka UPS nie wykryte jako istotne (przegapione)
#
# Używamy -log10(p-value) lub -log10(adj.P.Val) jako score do ROC
# ==============================================================================

message("\n=== Analiza ROC/AUC i Precision-Recall ===")

run_roc_analysis <- function(results_table, variant_name) {
  
  # Przygotuj dane - potrzebujemy score i etykietę (UPS=1, ECOLI=0)
  roc_data <- results_table %>%
    filter(Species %in% c("UPS", "ECOLI")) %>%
    mutate(
      # Etykieta: UPS = pozytywne (powinny być różnicowe)
      is_ups = ifelse(Species == "UPS", 1, 0),
      # Score: -log10(p-value) - wyższy score = bardziej istotny
      score = -log10(P.Value),
      # Dodatkowy score z adj.P.Val
      adj_score = -log10(adj.P.Val)
    ) %>%
    filter(!is.na(score) & !is.infinite(score))
  
  n_ups <- sum(roc_data$is_ups == 1)
  n_ecoli <- sum(roc_data$is_ups == 0)
  
  message("\n--- ", variant_name, " ---")
  message("Białek UPS: ", n_ups)
  message("Białek E. coli: ", n_ecoli)
  
  if (n_ups < 5 || n_ecoli < 5) {
    message("Za mało danych do analizy ROC!")
    return(NULL)
  }
  
  # Oblicz podstawowe metryki przy ustalonym progu
  # Używamy kryterium: adj.P.Val < 0.05 AND |logFC| > 1
  significant <- results_table$adj.P.Val < 0.05 & abs(results_table$logFC) > 1
  
  roc_data <- roc_data %>%
    mutate(significant = adj.P.Val < 0.05 & abs(logFC) > 1)
  
  TP <- sum(roc_data$is_ups == 1 & roc_data$significant, na.rm = TRUE)
  FP <- sum(roc_data$is_ups == 0 & roc_data$significant, na.rm = TRUE)
  TN <- sum(roc_data$is_ups == 0 & !roc_data$significant, na.rm = TRUE)
  FN <- sum(roc_data$is_ups == 1 & !roc_data$significant, na.rm = TRUE)
  
  sensitivity <- TP / (TP + FN)  # True Positive Rate (Recall)
  specificity <- TN / (TN + FP)  # True Negative Rate
  precision <- ifelse(TP + FP > 0, TP / (TP + FP), 0)  # Positive Predictive Value
  fdr <- ifelse(TP + FP > 0, FP / (TP + FP), 0)  # False Discovery Rate
  
  message("\nPrzy progu adj.P.Val < 0.05 & |logFC| > 1:")
  message("  TP (UPS istotne): ", TP)
  message("  FP (E. coli istotne): ", FP)
  message("  TN (E. coli nieistotne): ", TN)
  message("  FN (UPS nieistotne): ", FN)
  message("  Sensitivity (Recall): ", round(sensitivity * 100, 1), "%")
  message("  Specificity: ", round(specificity * 100, 1), "%")
  message("  Precision: ", round(precision * 100, 1), "%")
  message("  FDR (Empirical): ", round(fdr * 100, 1), "%")
  
  # Oblicz ROC curve i AUC
  roc_result <- NULL
  auc_value <- NA
  
  if (has_pROC) {
    tryCatch({
      # ROC na podstawie raw p-value score
      roc_result <- roc(roc_data$is_ups, roc_data$score, quiet = TRUE)
      auc_value <- auc(roc_result)
      
      message("\n  AUC (Area Under ROC Curve): ", round(auc_value, 3))
      
    }, error = function(e) {
      message("Błąd w obliczeniu ROC: ", e$message)
    })
  }
  
  return(list(
    variant = variant_name,
    n_ups = n_ups,
    n_ecoli = n_ecoli,
    TP = TP, FP = FP, TN = TN, FN = FN,
    sensitivity = sensitivity,
    specificity = specificity,
    precision = precision,
    fdr = fdr,
    auc = as.numeric(auc_value),
    roc_curve = roc_result,
    roc_data = roc_data
  ))
}

# Uruchom analizę ROC dla obu wariantów
roc_v1 <- run_roc_analysis(results_v1, "Median (bez imputacji)")
roc_v2 <- run_roc_analysis(results_v2, "VSN")

# ==============================================================================
# 8. WYKRESY ROC
# ==============================================================================

if (has_pROC && !is.null(roc_v1) && !is.null(roc_v2) && 
    !is.null(roc_v1$roc_curve) && !is.null(roc_v2$roc_curve)) {
  
  message("\n=== Generowanie wykresów ROC ===")
  
  # Wykres ROC dla obu wariantów
  png("results/figures/roc_comparison.png", width = 800, height = 700, res = 100)
  
  # Rysuj pierwszą krzywą
  plot(roc_v1$roc_curve, 
       main = "ROC Curves - Comparison of Normalization Methods",
       col = "#E41A1C", lwd = 3, lty = 1,  # gruba linia ciągła
       print.auc = FALSE, legacy.axes = TRUE)
  
  # Druga krzywa z innym stylem - przerywana i cieńsza
  plot(roc_v2$roc_curve, 
       add = TRUE, col = "#377EB8", lwd = 2, lty = 2)  # linia przerywana
  
  # Linia diagonalna (random classifier)
  abline(a = 0, b = 1, lty = 3, col = "gray50")
  
  legend("bottomright", 
         legend = c(
           paste0("Median (AUC = ", round(roc_v1$auc, 3), ")"),
           paste0("VSN (AUC = ", round(roc_v2$auc, 3), ")")
         ),
         col = c("#E41A1C", "#377EB8"), 
         lwd = c(3, 2), 
         lty = c(1, 2),  # solid, dashed
         cex = 1.1)
  
  dev.off()
  message("Zapisano: results/figures/roc_comparison.png")
}

# ==============================================================================
# 9. PODSUMOWANIE ROC/AUC
# ==============================================================================

if (!is.null(roc_v1) && !is.null(roc_v2)) {
  roc_summary <- data.frame(
    Wariant = c("Median (bez imputacji)", "VSN"),
    N_UPS = c(roc_v1$n_ups, roc_v2$n_ups),
    N_Ecoli = c(roc_v1$n_ecoli, roc_v2$n_ecoli),
    TP = c(roc_v1$TP, roc_v2$TP),
    FP = c(roc_v1$FP, roc_v2$FP),
    Sensitivity_pct = round(c(roc_v1$sensitivity, roc_v2$sensitivity) * 100, 1),
    Specificity_pct = round(c(roc_v1$specificity, roc_v2$specificity) * 100, 1),
    Precision_pct = round(c(roc_v1$precision, roc_v2$precision) * 100, 1),
    AUC = round(c(roc_v1$auc, roc_v2$auc), 3)
  )
  
  message("\n=== PODSUMOWANIE ROC/AUC ===")
  print(roc_summary)
  
  write_csv(roc_summary, "results/tables/roc_summary.csv")
  message("\nZapisano: results/tables/roc_summary.csv")
}

# ==============================================================================
# 10. ZAPISZ WYNIKI BENCHMARK
# ==============================================================================

message("\n=== Zapisywanie wyników benchmark ===")

benchmark_results <- list(
  v1 = benchmark_v1,
  v2 = benchmark_v2,
  comparison = comparison_summary,
  ups_expected = ups_expected_fc,
  roc_v1 = if(exists("roc_v1")) roc_v1 else NULL,
  roc_v2 = if(exists("roc_v2")) roc_v2 else NULL,
  roc_summary = if(exists("roc_summary")) roc_summary else NULL
)

saveRDS(benchmark_results, "results/benchmark_results.rds")
write_csv(comparison_summary, "results/tables/benchmark_comparison.csv")

# EKSTRA TABELA UPS - wg sugestii
write_csv(ups_expected_fc, "results/tables/ups_expected_legend.csv")
message("Zapisano: results/tables/ups_expected_legend.csv (Legenda stężeń UPS)")

message("Zapisano: results/benchmark_results.rds")
message("Zapisano: results/tables/benchmark_comparison.csv")

# ==============================================================================
# 11. TABELA FALSE POSITIVES (E. coli)
# ==============================================================================

message("\n=== Generowanie tabeli False Positives ===")

# Podsumowanie FP
fp_summary <- data.frame(
  Wariant = c("Median", "VSN"),
  Ecoli_total = c(benchmark_v1$n_ecoli, benchmark_v2$n_ecoli),
  False_Positives = c(benchmark_v1$n_ecoli_fp, benchmark_v2$n_ecoli_fp),
  FP_Rate_pct = c(round(benchmark_v1$empirical_fdr, 2), round(benchmark_v2$empirical_fdr, 2))
)

message("=== FALSE POSITIVES SUMMARY ===")
message("Note: FP_Rate on E. coli serves as an empirical estimate of FDR (since E. coli should be null).")
print(fp_summary)

write_csv(fp_summary, "results/tables/false_positives_summary.csv")
message("Zapisano: results/tables/false_positives_summary.csv")

# Szczegóły FP dla każdego wariantu
if (benchmark_v1$n_ecoli_fp > 0) {
  fp_details_v1 <- benchmark_v1$ecoli_results %>%
    filter(Significant == TRUE) %>%
    select(Protein_ID, logFC, adj.P.Val) %>%
    mutate(logFC = round(logFC, 3), adj.P.Val = signif(adj.P.Val, 3))
  write_csv(fp_details_v1, "results/tables/false_positives_median_details.csv")
  message("FP Median details: ", nrow(fp_details_v1), " proteins")
}

if (benchmark_v2$n_ecoli_fp > 0) {
  fp_details_v2 <- benchmark_v2$ecoli_results %>%
    filter(Significant == TRUE) %>%
    select(Protein_ID, logFC, adj.P.Val) %>%
    mutate(logFC = round(logFC, 3), adj.P.Val = signif(adj.P.Val, 3))
  write_csv(fp_details_v2, "results/tables/false_positives_vsn_details.csv")
  message("FP VSN details: ", nrow(fp_details_v2), " proteins")
}

# ==============================================================================
# 12. WYKRESY PORÓWNAWCZE (SIDE-BY-SIDE)
# ==============================================================================

message("\n=== Generowanie wykresów porównawczych ===")

# Sprawdź czy pakiety są dostępne
if (requireNamespace("png", quietly = TRUE) && 
    requireNamespace("grid", quietly = TRUE) && 
    requireNamespace("gridExtra", quietly = TRUE)) {
  
  library(png)
  library(grid)
  library(gridExtra)
  
  load_png_as_grob <- function(path) {
    if (file.exists(path)) {
      img <- readPNG(path)
      return(rasterGrob(img, interpolate = TRUE))
    }
    return(NULL)
  }
  
  # PCA Comparison
  pca_median <- load_png_as_grob("results/figures/pca_median.png")
  pca_vsn <- load_png_as_grob("results/figures/pca_vsn.png")
  
  if (!is.null(pca_median) && !is.null(pca_vsn)) {
    png("results/figures/comparison_pca.png", width = 1600, height = 700, res = 100)
    grid.arrange(pca_median, pca_vsn, ncol = 2,
                 top = textGrob("Porownanie PCA: Median vs VSN", 
                               gp = gpar(fontsize = 16, fontface = "bold")))
    dev.off()
    message("Zapisano: results/figures/comparison_pca.png")
  }
  
  # Volcano Comparison
  volcano_median <- load_png_as_grob("results/figures/volcano_median.png")
  volcano_vsn <- load_png_as_grob("results/figures/volcano_vsn.png")
  
  if (!is.null(volcano_median) && !is.null(volcano_vsn)) {
    png("results/figures/comparison_volcano.png", width = 1600, height = 700, res = 100)
    grid.arrange(volcano_median, volcano_vsn, ncol = 2,
                 top = textGrob("Porownanie Volcano Plot: Median vs VSN", 
                               gp = gpar(fontsize = 16, fontface = "bold")))
    dev.off()
    message("Zapisano: results/figures/comparison_volcano.png")
  }
  
  # Expected vs Observed Comparison
  exp_median <- load_png_as_grob("results/figures/expected_vs_observed_median.png")
  exp_vsn <- load_png_as_grob("results/figures/expected_vs_observed_vsn.png")
  
  if (!is.null(exp_median) && !is.null(exp_vsn)) {
    png("results/figures/comparison_expected_vs_observed.png", width = 1600, height = 700, res = 100)
    grid.arrange(exp_median, exp_vsn, ncol = 2,
                 top = textGrob("Porownanie Expected vs Observed FC: Median vs VSN", 
                               gp = gpar(fontsize = 16, fontface = "bold")))
    dev.off()
    message("Zapisano: results/figures/comparison_expected_vs_observed.png")
  }
  
} else {
  message("Brak pakietów png/grid/gridExtra - pomijam wykresy porównawcze")
  message("Zainstaluj: install.packages(c('png', 'gridExtra'))")
}

message("\n=== Benchmark zakończony ===")
