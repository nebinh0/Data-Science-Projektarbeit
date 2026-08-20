# QA/QC & Exploration
# Data Science Projekt - Golub Dataset
# Tim 

# Datenimport / vorbereitete Matrix von Linda
source("01_data_import_final.R")


#=============================================================
#1. Qualitätskontrolle der Daten
#=============================================================

# Dimensionen der Analyse-Matrix
dim(analysis_matrix)
# Dimensionen der Metadaten
dim(metadata)
# Struktur der Analyse-Matrix
str(analysis_matrix)
# Struktur der Metadaten
str(metadata)



#=============================================================
#2. Fehlende Werte prüfen
#=============================================================
# Anzahl aller fehlender Werte
sum(is.na(analysis_matrix))

# Fehlende Werte pro Gen
colSums(is.na(analysis_matrix))
colSums(is.na(analysis_matrix))[colSums(is.na(analysis_matrix))>0]



#=============================================================


#Grundlegende Verteilung der Expressionswerte 
summary(as.vector(analysis_matrix))

hist(as.vector(analysis_matrix),
     breaks = 50,
     main = "Verteilung der Genexpressionswerte",
     xlab = "Expressionswert")

#hist(log2(as.vector(analysis_matrix)+1),
     #breaks = 50,
     #main = "Verteilung der Genexpressionswerte",
     #xlab = "log2(Expressionswert + 1)")

#log2-Transformation hier nicht verwendet, da negative Expressionswerte vorhanden sind


#Verteilung pro Probe: Sind manche Proben insgesamt deutlich anderts als die anderen?
sample_means <- rowMeans(analysis_matrix, na.rm=TRUE)
summary(sample_means)

boxplot(
  t(analysis_matrix),
  las =2,
  outline=FALSE,
  main = "Verteilung der Genexpressionswerte pro Probe",
  ylab = "Expressionswert")

dim(analysis_matrix)
head(rownames(analysis_matrix))
head(colnames(analysis_matrix))

# Mittelwerte der Expressionswerte pro Probe
sample_means <- rowMeans(analysis_matrix, na.rm=TRUE)
summary(sample_means)

plot (
  sample_means,
  type = "h",
  main = "Mittlere Expressionswerte pro Probe",
  xlab = "Probe",
  ylab = "mittlerer Expressionswert"
)

#Proben mit besonders hohen mittleren Expressionswerten
sort(sample_means, decreasing=TRUE)

# Mittlere Expressionswerte mit Probenamen:
names (sample_means) <- rownames (analysis_matrix)
sort (sample_means, decreasing = TRUE)

#Wie ähnlich sind die Proben insgesamt in ihrem Genexpressionsprofil=
#Ähnlichkeit der Proben untereinander
sample_cor <- cor(
  t(analysis_matrix),
  use= "pairwise.complete.obs"
)

#Verteilung der Korrelationswerte
summary(as.vector(sample_cor))

# Korrelationsmatrix als Heatmap
heatmap (
  sample_cor,
  main= "Korrelation zwischen den Proben",
  xlab = "Probe",
  ylab = "Probe"
)

# Proben mit besonders niedriger mittlerer Korrelation
mean_cor <- rowMeans(sample_cor, na.rm=TRUE)
names  (mean_cor) <- rownames(analysis_matrix)
sort (mean_cor)

# Korrelation der auffälligsten Probe mit allen anderen Proben
sort(sample_cor["39",])

# Mittlere Korrelation ohne Selbstkorrelation
sample_cor_no_diag <- sample_cor
diag(sample_cor_no_diag) <- NA

mean_cor <- rowMeans(sample_cor_no_diag, na.rm=TRUE)
names(mean_cor) <- rownames(analysis_matrix)

sort(mean_cor)

# Probe 6 weist mit mittleren Korrelationswert von ca.0,221 die geringste Ähnlichkeit zu den übrigen Proben auf. Danach folgen Proben 39 (0,272), Probe 3 (0,281), Probe 20 (0,294) und Probe 48 (0,295)

#Expressionswerte der einzelnen Proben hinsichtlich ihrer Verteilung wurden zur Qualitätskontrole untersucht. Anschließend wurden die mittleren Expressionswerte sowie die Korrelation zwischen Proben bestimmmt. Um Einfluss der Selbskorrelation zu vermeiden, wurde die Diagonale der Korrelationsmatrix auf NA gesetzt. Probe 6 zeigt mit einer mittleren Korrelation von 0,221 die geringste Ähnlichkeit zu den übrigen Proben und identifziert sichg somit als auffällige Probe. 

