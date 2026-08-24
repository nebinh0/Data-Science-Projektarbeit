# PCA und Visualisierung
# Data Science Projektarbeit - Golub Dataset


# ============================================================
# 1. Daten laden
# ============================================================

# Vorbereitete Expressionsmatrix und Metadaten laden
source("01_data_import_final.R")


# ============================================================
# 2. Daten für PCA vorbereiten
# ============================================================

# Features mit fehlenden Werten identifizieren
features_mit_na <- colSums(is.na(analysis_matrix)) > 0

# Für die PCA nur Features ohne fehlende Werte verwenden
pca_matrix <- analysis_matrix[, !features_mit_na]

# Kontrolle
dim(pca_matrix)
sum(is.na(pca_matrix))

# Prüfen, ob Features ohne Variation vorhanden sind
feature_sd <- apply(pca_matrix, 2, sd)
sum(feature_sd == 0)


# ============================================================
# 3. PCA
# ============================================================

# PCA mit Skalierung der Features
pca <- prcomp(
  pca_matrix,
  scale. = TRUE
)

summary(pca)


# ============================================================
# 4. Erklärte Varianz
# ============================================================

# Anteil der erklärten Varianz in Prozent
varianz <- pca$sdev^2 / sum(pca$sdev^2) * 100

# Erklärte Varianz von PC1 und PC2
varianz[1:2]

# Gemeinsam erklärte Varianz von PC1 und PC2
sum(varianz[1:2])


# ============================================================
# 5. PCA-Koordinaten und Metadaten
# ============================================================

# Koordinaten der Proben in den Hauptkomponenten
pca_scores <- as.data.frame(pca$x)

# Prüfen, ob Reihenfolge der Proben mit Metadaten übereinstimmt
all(rownames(pca_scores) == metadata$Sample)

# Diagnose und Datensatz ergänzen
pca_scores$Diagnosis <- metadata$Diagnosis

pca_scores$Dataset <- factor(
  metadata$Dataset,
  levels = c("initial", "independent")
)


# ============================================================
# 6. PCA-Plot ALL vs. AML
# ============================================================

plot(
  pca_scores$PC1,
  pca_scores$PC2,
  col = pca_scores$Diagnosis,
  pch = 19,
  xlab = paste0("PC1 (", round(varianz[1], 2), " %)"),
  ylab = paste0("PC2 (", round(varianz[2], 2), " %)"),
  main = "PCA der Genexpressionsdaten"
)

legend(
  "topright",
  legend = levels(pca_scores$Diagnosis),
  col = 1:2,
  pch = 19
)


# ============================================================
# 7. Screeplot
# ============================================================

plot(
  varianz[1:10],
  type = "b",
  xlab = "Hauptkomponente",
  ylab = "Erklärte Varianz (%)",
  main = "Screeplot",
  xaxt = "n"
)

axis(
  1,
  at = 1:10,
  labels = paste0("PC", 1:10)
)


# ============================================================
# 8. Kontrolle: initialer vs. unabhängiger Datensatz
# ============================================================

plot(
  pca_scores$PC1,
  pca_scores$PC2,
  col = pca_scores$Dataset,
  pch = 19,
  xlab = paste0("PC1 (", round(varianz[1], 2), " %)"),
  ylab = paste0("PC2 (", round(varianz[2], 2), " %)"),
  main = "PCA nach Datensatz"
)

legend(
  "topright",
  legend = levels(pca_scores$Dataset),
  col = 1:2,
  pch = 19
)


# ============================================================
# 9. Loadings
# ============================================================

# Loadings der Hauptkomponenten
loadings <- pca$rotation

# Features mit den stärksten absoluten Loadings auf PC1 und PC2
pc1_order <- order(
  abs(loadings[, "PC1"]),
  decreasing = TRUE
)

pc2_order <- order(
  abs(loadings[, "PC2"]),
  decreasing = TRUE
)

head(loadings[pc1_order, "PC1"], 10)
head(loadings[pc2_order, "PC2"], 10)

# Fünf stärkste positive und negative PC1-Loadings
pc1_pos <- sort(
  loadings[, "PC1"],
  decreasing = TRUE
)[1:5]

pc1_neg <- sort(
  loadings[, "PC1"],
  decreasing = FALSE
)[1:5]

pc1_top <- c(pc1_neg, pc1_pos)

# Loading-Plot
par(mar = c(5, 6, 4, 2))

barplot(
  pc1_top,
  horiz = TRUE,
  las = 1,
  xlab = "Loading",
  main = "Stärkste Loadings auf PC1"
)


# ============================================================
# 10. Abbildungen speichern
# ============================================================

# PCA ALL vs. AML
png(
  "figures/PCA_ALL_AML.png",
  width = 1200,
  height = 900,
  res = 150
)

plot(
  pca_scores$PC1,
  pca_scores$PC2,
  col = pca_scores$Diagnosis,
  pch = 19,
  xlab = paste0("PC1 (", round(varianz[1], 2), " %)"),
  ylab = paste0("PC2 (", round(varianz[2], 2), " %)"),
  main = "PCA der Genexpressionsdaten"
)

legend(
  "topright",
  legend = levels(pca_scores$Diagnosis),
  col = 1:2,
  pch = 19
)

dev.off()


# Screeplot
png(
  "figures/Screeplot.png",
  width = 1200,
  height = 900,
  res = 150
)

plot(
  varianz[1:10],
  type = "b",
  xlab = "Hauptkomponente",
  ylab = "Erklärte Varianz (%)",
  main = "Screeplot",
  xaxt = "n"
)

axis(
  1,
  at = 1:10,
  labels = paste0("PC", 1:10)
)

dev.off()


# PC1-Loadings
png(
  "figures/PC1_Loadings.png",
  width = 1200,
  height = 900,
  res = 150
)

par(mar = c(5, 6, 4, 2))

barplot(
  pc1_top,
  horiz = TRUE,
  las = 1,
  xlab = "Loading",
  main = "Stärkste Loadings auf PC1"
)

dev.off()