# ============================================================
# Clustering
# Data Science Projekt - Golub Dataset
# Alena
# ============================================================

# Paket für Heatmap laden
library(pheatmap)

# Datenimport und Preprocessing
source("01_data_import_final.R")

# Ausgangsdaten kontrollieren
dim(analysis_matrix)
table(metadata$Diagnosis)
table(metadata$Dataset, metadata$Diagnosis)


# ============================================================
# 1. Vorbereitung der Daten für das Clustering
# ============================================================

# Features mit fehlenden Werten identifizieren
features_with_na <- colSums(is.na(analysis_matrix)) > 0

# Anzahl betroffener Features
sum(features_with_na)

# Features mit fehlenden Werten entfernen
clustering_matrix <- analysis_matrix[, !features_with_na]

# Dimensionen und fehlende Werte kontrollieren
dim(clustering_matrix)
sum(is.na(clustering_matrix))


# Varianz jedes Features über alle Proben berechnen
feature_variance <- apply(clustering_matrix, 2, var)

# Verteilung der Varianzen
summary(feature_variance)

# Die 25 % variabelsten Features auswählen
variance_cutoff <- quantile(feature_variance, 0.75)
variable_features <- feature_variance >= variance_cutoff

# Anzahl ausgewählter Features
sum(variable_features)

# Matrix mit den variabelsten Features erstellen
clustering_variable <- clustering_matrix[, variable_features]
dim(clustering_variable)


# Features standardisieren
clustering_scaled <- scale(clustering_variable)

# Standardisierung kontrollieren
dim(clustering_scaled)
colMeans(clustering_scaled)[1:5]
apply(clustering_scaled, 2, sd)[1:5]


# ============================================================
# 2. Hierarchisches Clustering
# ============================================================

# Euklidische Distanzen zwischen den Proben
sample_dist <- dist(clustering_scaled, method = "euclidean")

# Hierarchisches Clustering mit Ward-Methode
hc <- hclust(sample_dist, method = "ward.D2")

# Dendrogramm
plot(
  hc,
  main = "Hierarchisches Clustering der Proben",
  xlab = "Probe",
  ylab = "Distanz",
  sub = ""
)

# Dendrogramm in zwei Cluster aufteilen
hc_clusters <- cutree(hc, k = 2)

# Anzahl der Proben pro Cluster
table(hc_clusters)


# ============================================================
# 3. Vergleich der Cluster mit den Metadaten
# ============================================================

# 3.1 Diagnose

# ALL / AML
table(
  Cluster = hc_clusters,
  Diagnosis = metadata$Diagnosis
)

# Detaillierte Diagnose: B-ALL / T-ALL / AML
table(
  Cluster = hc_clusters,
  Cancer = metadata$Cancer
)


# 3.2 Datensatz

# Initial / Independent
table(
  Cluster = hc_clusters,
  Dataset = metadata$Dataset
)


# 3.3 Weitere Metadaten

# Probenmaterial: Knochenmark / peripheres Blut
table(
  Cluster = hc_clusters,
  BM_PB = metadata$BM_PB,
  useNA = "ifany"
)

# Herkunft der Proben
table(
  Cluster = hc_clusters,
  Source = metadata$Source,
  useNA = "ifany"
)

# Geschlecht
table(
  Cluster = hc_clusters,
  Gender = metadata$Gender,
  useNA = "ifany"
)

# Kombination aus Probenmaterial und Geschlecht
table(
  Cluster = hc_clusters,
  Tissue_MF = metadata$Tissue_MF,
  useNA = "ifany"
)


# ============================================================
# 4. Bewertung der Übereinstimmung von Cluster und Diagnose
# ============================================================

# Kontingenztabelle
cluster_diagnosis <- table(
  Cluster = hc_clusters,
  Diagnosis = metadata$Diagnosis
)

cluster_diagnosis

# Fisher-Exakt-Test
fisher.test(cluster_diagnosis)

# Ergebnis:
# Keine statistisch nachweisbare Assoziation zwischen
# Clusterzugehörigkeit und ALL/AML-Diagnose (p = 0.7752).


# ============================================================
# 5. Heatmap der Expressionsprofile
# ============================================================

# 50 Features mit der höchsten Varianz auswählen
top50_names <- names(sort(feature_variance, decreasing = TRUE)[1:50])

# Standardisierte Expressionswerte der Top-50-Features
heatmap_matrix <- clustering_scaled[, top50_names]

# Für Heatmap: Features in Zeilen, Proben in Spalten
heatmap_matrix_t <- t(heatmap_matrix)

# Diagnose als Annotation der Proben
annotation_col <- data.frame(
  Diagnosis = metadata$Diagnosis
)
rownames(annotation_col) <- metadata$Sample

# Übereinstimmung der Probenreihenfolge kontrollieren
identical(
  colnames(heatmap_matrix_t),
  as.character(metadata$Sample)
)

# Heatmap
pheatmap(
  heatmap_matrix_t,
  annotation_col = annotation_col,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D2",
  show_colnames = FALSE,
  main = "Expressionsprofile der 50 variabelsten Features"
)


# ============================================================
# 6. K-means-Clustering
# ============================================================

# Zufallsstart reproduzierbar machen
set.seed(123)

# K-means mit zwei Clustern
kmeans_result <- kmeans(
  clustering_scaled,
  centers = 2,
  nstart = 25
)

# Anzahl der Proben pro Cluster
table(kmeans_result$cluster)


# Vergleich mit ALL / AML
kmeans_diagnosis <- table(
  Cluster = kmeans_result$cluster,
  Diagnosis = metadata$Diagnosis
)

kmeans_diagnosis

# Fisher-Exakt-Test
fisher.test(kmeans_diagnosis)


# Vergleich mit dem hierarchischen Clustering
table(
  Hierarchical = hc_clusters,
  Kmeans = kmeans_result$cluster
)

# Beide Verfahren ergeben dieselbe Zuordnung:
# Cluster 1: 16 Proben
# Cluster 2: 56 Proben