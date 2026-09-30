# Differential Gene Expression
# Data Science Projektarbeit - Golub Dataset
# Jan


cat("===== DIFFERENTIAL GENE EXPRESSION =====\n")

# Ordner für Ergebnisse meines Teils anlegen
dir.create("differential_expression", showWarnings = FALSE)

# 1. Vorbereitete Daten laden
source("01_data_import_final.R")


# 2. Überblick über die vorbereiteten Daten
dim(analysis_matrix)
dim(metadata)


# 3. Verteilung der Diagnosen
table(metadata$Diagnosis)


# 4. Verteilung von Datensatz und Diagnose
table(metadata$Dataset, metadata$Diagnosis)


# 5. Initialen Datensatz auswählen
initial_samples <- metadata$Dataset == "initial"

analysis_matrix_initial <- analysis_matrix[initial_samples, ]
metadata_initial <- metadata[initial_samples, ]


# 6. Kontrolle des ausgewählten Datensatzes
cat("\n===== KONTROLLE INITIALER DATENSATZ =====\n")

print(dim(analysis_matrix_initial))
print(dim(metadata_initial))
print(table(metadata_initial$Diagnosis))

cat("===== ENDE KONTROLLE =====\n")

# 7. Kontrolle auf fehlende Werte
cat("\n===== KONTROLLE FEHLENDE WERTE =====\n")

print(sum(is.na(analysis_matrix_initial)))

cat("Anzahl Gene mit mindestens einem NA:\n")
print(sum(colSums(is.na(analysis_matrix_initial)) > 0))

cat("===== ENDE NA-KONTROLLE =====\n")

cat("\n===== GENAUERE NA-KONTROLLE =====\n")

# Anzahl NAs pro Gen
na_pro_gen <- colSums(is.na(analysis_matrix_initial))

# Nur Gene mit mindestens einem NA
na_gene_liste <- na_pro_gen[na_pro_gen > 0]

print(na_gene_liste)

cat("\nAnzahl betroffener Gene:\n")
print(length(na_gene_liste))

cat("\nMaximale Anzahl NAs in einem Gen:\n")
print(max(na_gene_liste))

cat("===== ENDE NA-KONTROLLE =====\n")

# 7.1. Test der Log2-Transformation

#cat("\n===== TEST LOG2-TRANSFORMATION =====\n")

#min_value <- min(analysis_matrix_initial, na.rm = TRUE)

#analysis_matrix_log <- log2(
 # analysis_matrix_initial - min_value + 1


#print(summary(as.vector(analysis_matrix_log)))
#print(range(analysis_matrix_log, na.rm = TRUE))

#cat("===== ENDE TEST =====\n")

# 7.2 Prüfen, ob die Daten bereits standardisiert sind

cat("\n===== STANDARDISIERUNGSCHECK =====\n")

gene_means <- rowMeans(analysis_matrix_initial, na.rm = TRUE)
gene_sds <- apply(analysis_matrix_initial, 1, sd, na.rm = TRUE)

print(summary(gene_means))
print(summary(gene_sds))

cat("===== ENDE STANDARDISIERUNGSCHECK =====\n")

# 8. Differential Expression mit limma

library(limma)

group <- factor(metadata_initial$Diagnosis,
                levels = c("ALL", "AML"))

design <- model.matrix(~ group)

fit <- lmFit(t(analysis_matrix_initial), design)
fit <- eBayes(fit)

results <- topTable(fit,
                    coef = "groupAML",
                    number = Inf,
                    adjust.method = "BH")

cat("\n===== TOP DIFFERENTIALLY EXPRESSED GENES =====\n")
print(head(results, 10))

cat("\n===== SIGNIFIKANTE GENE (FDR < 0.05) =====\n")
print(sum(results$adj.P.Val < 0.05))

# 9. Ergebnisse speichern

write.csv(results,
          "differential_expression/differential_expression_results.csv",
          row.names = TRUE)

cat("\n===== ERGEBNISSE GESPEICHERT =====\n")

# 9.1. Top differentiell exprimierte Gene

top_genes <- results[results$adj.P.Val < 0.05, ]

top_genes <- top_genes[
  order(top_genes$adj.P.Val),
]

top_genes <- head(top_genes, 20)

cat("\n===== TOP 20 DIFFERENTIELL EXPRIMIERTE GENE =====\n")
print(top_genes)

# 10. Differential-Expression-Plot

library(ggplot2)

# Signifikanz definieren
results$significant <- results$adj.P.Val < 0.05

# Differential-Expression-Plot
de_plot <- ggplot(results,
                  aes(x = logFC,
                      y = -log10(P.Value),
                      color = significant)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    title = "Differential Gene Expression: AML vs. ALL",
    x = "Expression difference (AML vs. ALL)",
    y = "-log10(p-value)",
    color = "FDR < 0.05"
  ) +
  theme_minimal()

print(de_plot)

# Plot speichern
ggsave(
  "differential_expression/differential_expression_plot.png",
  de_plot,
  width = 8,
  height = 6,
  dpi = 300
)

cat("\n===== DIFFERENTIAL EXPRESSION PLOT ERSTELLT =====\n")

# Check der extremen Werte

cat("\n===== EXTREME WERTE =====\n")

extreme <- which(abs(analysis_matrix_initial) > 10000,
                 arr.ind = TRUE)

print(nrow(extreme))

print(data.frame(
  Sample = rownames(analysis_matrix_initial)[extreme[,1]],
  Feature = colnames(analysis_matrix_initial)[extreme[,2]],
  Value = analysis_matrix_initial[extreme]
)[1:20, ])

cat("===== ENDE =====\n")

# Extremwerte genauer untersuchen

cat("\n===== EXTREMSTE WERTE =====\n")

max_value <- max(analysis_matrix_initial, na.rm = TRUE)

print(max_value)

max_pos <- which(
  analysis_matrix_initial == max_value,
  arr.ind = TRUE
)

print(data.frame(
  Sample = rownames(analysis_matrix_initial)[max_pos[, "row"]],
  Feature = colnames(analysis_matrix_initial)[max_pos[, "col"]],
  Value = analysis_matrix_initial[max_pos]
))

cat("\nTop 20 Werte:\n")
print(
  sort(
    as.vector(analysis_matrix_initial),
    decreasing = TRUE,
    na.last = NA
  )[1:20]
)

cat("===== ENDE =====\n")
