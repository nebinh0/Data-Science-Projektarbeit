# 01_data_import.R
# Data import and preprocessing
# Project: Data Science Projektarbeit 
# Author: Linda

data_raw <- read.csv (
  "data/data .csv",
  header = FALSE,
  check.names = FALSE
)
  dim(data_raw)
  head(data_raw, 15)
  str(data_raw)
  
  names(data_raw)[1:30]
  data_raw[1:4, 1:10]
  str(data_raw)
  
  data_raw[1:4, 1:30]
  
  # Anzahl der Zeilen und Spalten
  dim(data_raw)
  #N Namen der ersten 10 Spalten
  names(data_raw) [1:10]
  # Namen der letzten 10 Spalten
  names (data_raw) [2965:2974]
  #Messdaten der ersten 6 Spalten ansehen
  data_raw[1:10, 1:6]
  # Untersuchung der ersten vier Zeilen
  data_raw[1:4, 1:15]
  data_raw[1:4, 2965:2974]  
  dim(data_raw)  
  names(data_raw)[1:15]  

  #Metadaten extrahieren
  metadata <- data_raw[4:nrow(data_raw), 1:6]
  head(metadata)
  
  data_raw[1:30, 1:6]  
  which(data_raw[, 2] %in% c("BM", "PB"))  
  data_raw[
    which(data_raw[, 2] %in% c("BM", "PB")),
    1:6
  ]  
  
  #Zeilen mit Sample-Metadaten identifizieren
  sample_rows <- which(
    data_raw[, 2] %in% c("BM", "PB")
  )
  sample_rows  
  #Anzahl der Samples überprüfen
  length(sample_rows)
  metadata <- data_raw[
    sample_rows,
    1:6
  ]  

  #Spalten benennen
  names(metadata) <- c(
    "Sample",
    "BM_PB",
    "Gender",
    "Source",
    "Tissue_MF",
    "Cancer"
  )
  head(metadata)  
  tail(metadata)  

  #Überprüfung der Diagnosen
  table(metadata$Cancer, useNA = "ifany")
  metadata$Diagnosis <- ifelse(
    tolower(metadata$Cancer) == "aml",
    "AML",
    ifelse(
      tolower(metadata$Cancer) %in% c("allb", "allt"),
      "ALL",
      NA
    )
  )  
  table(metadata$Diagnosis, useNA = "ifany")  

  #Fehlende Metadatenwerte
  colSums(is.na(metadata))

  #Darstellung Expressionsdaten
  data_raw[1:8, 7:20]
  data_raw[174:184, 1:20]  
  ncol(data_raw)  
  length(sample_rows)  
  sum(!is.na(unlist(data_raw[4:8, ])))  
  sum(!is.na(unlist(data_raw[9:13, ])))  

  #Gene entsprechenden Werten zuordnen
  data_raw[1:3, 2960:2974]
  data_raw[4:8, 2960:2974]  
  data_raw[4, 1:10]  
  data_raw[5, 1:10]  
  data_raw[6, 1:10]  
  data_raw[7, 1:10]  
  data_raw[8, 1:10]  

  #Rohdaten als einzelne CSV-Zeilen einlesen
  raw_lines <- readLines(
    "data/data .csv",
    encoding = "UTF-8",
    warn = FALSE
  )
  
  length(raw_lines)
  raw_split <- strsplit(
    raw_lines,
    ",",
    fixed = TRUE
  )  
  length(raw_split[[1]])  
  length(raw_split[[2]])  
  length(raw_split[[3]])  
  length(raw_split[[4]])
  length(raw_split[[5]])
  length(raw_split[[6]])
  length(raw_split[[7]])
  length(raw_split[[8]])  

  #Zusammenführen von Sample-Block 
  sample39 <- unlist(
    raw_split[4:8]
  )
  length(sample39)  
  sample39[1:10]  
  sample_start <- seq(
    from = 4,
    by = 5,
    length.out = 72
  )
  head(sample_start)  
  tail(sample_start)  
  samples_list <- lapply(
    sample_start,
    function(i) {
      unlist(raw_split[i:(i + 4)])
    }
  )  
  samples_list[[1]]  
  sapply(samples_list, length)  
  raw_lines <- readLines(
    "data/data .csv",
    encoding = "UTF-8",
    warn = FALSE
  )
  
  raw_split <- strsplit(
    raw_lines,
    ",",
    fixed = TRUE
  )
  
  length(raw_lines)
  
  length(raw_split[[1]])
  length(raw_split[[2]])
  length(raw_split[[3]])
  
  length(raw_split[[4]])
  length(raw_split[[5]])
  length(raw_split[[6]])
  length(raw_split[[7]])
  length(raw_split[[8]])
  
  sample39 <- unlist(raw_split[4:8])
  
  length(sample39)
  
  sample39[1:10]  

  #Datenrekonstruktion
  feature_header <- unlist(
    raw_split[1:3]
  )
  length(feature_header)  
  sample_info <- sample39[1:6]  
  expression_values <- sample39[7:length(sample39)]  
  sample_info  
  expression_values[1:10]  
  length(sample39) - 6  
  header1 <- raw_split[[1]]
  header2 <- raw_split[[2]]
  header3 <- raw_split[[3]]
  
  length(header1)
  length(header2)
  length(header3)  

  sample_start <- seq(
    from = 4,
    by = 5,
    length.out = 72
  )  
  head(sample_start)  
  tail(sample_start)  

  samples_list <- lapply(
    sample_start,
    function(i) {
      unlist(raw_split[i:(i + 4)])
    }
  )  
  lengths(samples_list)  

  #Vergleich der Sample-Blöcke
  table(lengths(samples_list))
  which(lengths(samples_list) != 7137)  
  metadata_temp <- data_raw[
    sample_rows,
    1:6
  ]
  
  metadata_temp[
    which(lengths(samples_list) != 7137),
    1:6
  ]  
  
  #Welche Blöcke sind unterschiedlich lang?
  block_lengths <- lapply(
    sample_start,
    function(i) {
      lengths(raw_split[i:(i + 4)])
    }
  )
  do.call(rbind, block_lengths)  
 
  #Kontrolle des vorherigen strsplitt Befehls
  field_counts <- count.fields(
    "data/data .csv",
    sep = ",",
    quote = "\""
  )
  
  field_counts[1:20]
  table(field_counts)  
  lengths(raw_split[4:8])  
  field_counts[4:8]  

  #Weitere Samples-Analyse
  table(lengths(samples_list))
  sample40 <- unlist(raw_split[9:13])
  
  length(sample40)  
  sample40[1:10]
  tail(sample40, 15)  
  tail(sample39, 15)  
 
  header <- c(
    raw_split[[1]],
    raw_split[[2]],
    raw_split[[3]]
  )
  
  length(header)  
  header[7120:7137]  
  sample39[7120:7137]  
  sample40[7120:7138]  

  sample40[7120:7138]  
  header[7134:7137]  
 
  #Identifikation der Position zusätzlicher Messwerte
  for (i in 9:13) {
    cat(
      "\nZeile", i,
      "Anzahl:", length(raw_split[[i]]),
      "\nErste Werte:",
      paste(head(raw_split[[i]], 6), collapse = " | "),
      "\nLetzte Werte:",
      paste(tail(raw_split[[i]], 6), collapse = " | "),
      "\n"
    )
  }
 
  which(sample40 == "") 
  sum(sample40 == "")  

  header[1:10]  
  length(sample39) - 6  
  length(sample40) - 6  
  non_numeric <- which(
    is.na(suppressWarnings(as.numeric(sample40[7:length(sample40)])))
  )
  
  length(non_numeric)  

  length(header) - 6  
  header[7:16]  
  tail(header, 10)  

  features <- header[7:length(header)]
  
  length(features)  
  length(unique(features))  
  features[duplicated(features)]  

  #Enthält Sample 40 einen Wert mehr oder sind die Werte gegneinander verschoben?
  expression39 <- sample39[7:length(sample39)]
  expression40 <- sample40[7:length(sample40)]
  length(expression39)
  length(expression40)  

  tail(expression39, 10)  
  tail(expression40, 10)  

  sum(expression40 == 0)  
  which(expression40 == 0)  

  block40 <- lengths(raw_split[9:13])
  block39 <- lengths(raw_split[4:8])
  
  rbind(
    Sample39 = block39,
    Sample40 = block40
  )  

  block40 - block39  

  #Enthält die Originaldatei absichtlich unterschiedlich lange Samples?
  data.frame(
    Feature = features[1:20],
    Sample39 = expression39[1:20],
    Sample40 = expression40[1:20]
  )

  field_counts <- count.fields(
    "data/data .csv",
    sep = ",",
    quote = "\""
  )
  
  field_counts[9:13]  
  lengths(raw_split[9:13])  

  #Wo beginnt die Abweichung von sample 40? 
  sample_lengths <- data.frame(
    Sample = metadata_temp$V1,
    TotalFields = lengths(samples_list)
  )
  
  head(sample_lengths, 20)
  table(sample_lengths$TotalFields)  

  #Bei der Rekonstruktion der Sample-Blöcke wurden unterschiedlich lange Datensätze festgestellt (7137–7139 Felder). Daher wurde die finale Analysematrix erst nach einer zusätzlichen Strukturprüfung erzeugt
  
  table(lengths(samples_list))
  sample_lengths <- data.frame(
    Sample = metadata_temp$V1,
    TotalFields = lengths(samples_list)
  )
  
  sample_lengths  
  
  table(sample_lengths$TotalFields)
  sample_lengths$Difference <-
    sample_lengths$TotalFields - 7137
  
  table(sample_lengths$Difference)  
  sample42 <- unlist(raw_split[14:18])
  
  length(sample42)
  
  sample42[1:10]
  
  sample42[7130:7139]  
  sample39[7120:7137]
  
  sample42[7120:7139]  
  
  # Samples mit exakt 7137 Feldern
  ok_samples <- which(
    sample_lengths$TotalFields == 7137
  )
  
  # Welche Samples sind das?
  sample_lengths[ok_samples, ]
  
  # Letzte drei Feature-Namen
  features[7129:7131]
  
  # Letzte Werte dieser Samples
  for (i in ok_samples) {
    cat(
      "\nSample:",
      sample_lengths$Sample[i],
      "\n"
    )
    
    print(
      tail(samples_list[[i]], 9)
    )
  }
  features[1:6]  
  features[7125:7131]  
  
  #Lösung der Problematik wird verschoben
  
  #Metadaten-Tabelle
  metadata <- data.frame(
    Sample = metadata_temp$V1,
    BM_PB = metadata_temp$V2,
    Gender = metadata_temp$V3,
    Source = metadata_temp$V4,
    Tissue_MF = metadata_temp$V5,
    Cancer = metadata_temp$V6,
    stringsAsFactors = FALSE
  )

  metadata$Diagnosis <- ifelse(
    tolower(metadata$Cancer) == "aml",
    "AML",
    ifelse(
      tolower(metadata$Cancer) %in% c("allb", "allt"),
      "ALL",
      NA
    )
  )  

  table(metadata$Diagnosis, useNA = "ifany") 

  #Fehlende Werte untersuchen
  colSums(is.na(metadata))
  table(metadata$Gender, useNA = "ifany")  

  #Vorläufige Expressionsmatrix
  expression_list <- lapply(
    samples_list,
    function(x) {
      x[7:length(x)]
    }
  )

  expression_lengths <- data.frame(
    Sample = metadata$Sample,
    Number_of_features = lengths(expression_list)
  )
  
  expression_lengths  

  valid_samples <- which(
    lengths(expression_list) == length(features)
  )  

  metadata[valid_samples, ]  

  analysis_matrix_temp <- do.call(
    rbind,
    expression_list[valid_samples]
  )  
  dim(analysis_matrix_temp)    

  colnames(analysis_matrix_temp) <- features  

  rownames(analysis_matrix_temp) <-
    metadata$Sample[valid_samples]  

  dim(analysis_matrix_temp)
  analysis_matrix_temp[1:5, 1:10]  

  # Datentyp prüfen
  str(analysis_matrix_temp)
  
  # In numerische Werte umwandeln
  analysis_matrix_temp <- apply(
    analysis_matrix_temp,
    2,
    as.numeric
  )
  
  # Namen wieder setzen
  colnames(analysis_matrix_temp) <- features
  rownames(analysis_matrix_temp) <- metadata$Sample[valid_samples]
  
  # Kontrolle
  is.numeric(analysis_matrix_temp)
  
  # Matrix ansehen
  analysis_matrix_temp[1:5, 1:10]
  
  # Fehlende Werte
  sum(is.na(analysis_matrix_temp))  
  
  colSums(is.na(metadata))
  round(
    colMeans(is.na(metadata)) * 100,
    2
  )  

  table(metadata$Cancer, useNA = "ifany")  

  table(metadata$Diagnosis, useNA = "ifany")  
  
  metadata$Diagnosis <- factor(
    ifelse(
      tolower(metadata$Cancer) == "aml",
      "AML",
      ifelse(
        tolower(metadata$Cancer) %in% c("allb", "allt"),
        "ALL",
        NA
      )
    ),
    levels = c("ALL", "AML")
  )

  table(metadata$Diagnosis, useNA = "ifany")  

  metadata$Gender <- factor(metadata$Gender)  

  table(metadata$Gender, useNA = "ifany")  

  #Missing-Value-Analyse:
 # Die Metadaten wurden auf fehlende Werte überprüft. Für die Variablen Sample, BM_PB, Source, Tissue_MF und Cancer wurden keine fehlenden Werte festgestellt. In der Variable Gender fehlten Angaben bei 23 von 72 Samples (31,94 %). Da es sich bei Gender um eine kategoriale Variable handelt und die Rekonstruktion so weiterhin möglich ist, wurden die Werte als NA beibehalten. Die Diagnosevariable Diagnosis wurde aus Cancer abgeleitet und enthält keine fehlenden Werte.
  
  #Prozentualer Anteil der Diagnosen
  prop.table(table(metadata$Diagnosis)) * 100

  #Zusammenführung Metadaten und Matrix
  metadata_valid <- metadata[valid_samples, ]
  metadata_valid$Sample  
  rownames(analysis_matrix_temp)  

  #Verknüpfung Metadaten und Expressionsmatrix über Samples
  metadata_valid
  analysis_matrix_temp  
  
  #Aufbereitung zur Trennung von Test- und Independent-Samples
  nrow(metadata)
  length(features)  
  anyDuplicated(features)  
  dim(analysis_matrix_temp)  

  #Definition der gemeinsamen Analysematrix
  metadata$Dataset <- c(
    rep("initial", 38),
    rep("independent", 34)
  )

  table(
    metadata$Dataset,
    metadata$Diagnosis,
    useNA = "ifany"
  )  

  metadata$Dataset <- NULL  

  metadata[, c("Sample", "Diagnosis")]  
 
  # Sample-IDs des initialen/training Datensatzes
  initial_ids <- c(
    1:38
  )
  
  # Alle übrigen Samples gehören zum unabhängigen Datensatz
  independent_ids <- c(
    39, 40, 42, 47, 48, 49, 41, 43, 44, 45, 46,
    70, 71, 72, 68, 69, 67, 55, 56, 59, 52, 53,
    51, 50, 54, 57, 58, 60, 61, 65, 66, 63, 64, 62
  ) 

  metadata$Dataset <- ifelse(
    metadata$Sample %in% initial_ids,
    "initial",
    "independent"
  )  

  table(
    metadata$Dataset,
    metadata$Diagnosis
  )  

  #Warum haben manche Samples 7132 bzw 7133 Werte?
  
  sample40 <- samples_list[[which(metadata$Sample == 40)]]
  
  
  length(sample40)

  tail(sample40, 15)  

  tail(features, 15)  
  
  sample40[1:10]
  
  sample40[7:20]
  sample40[7120:7138]  
  features[7120:7131]  

  data.frame(
    Position = 7126:7138,
    Value = sample40[7126:7138],
    Feature = c(
      features[7120:7131],
      "ZUSÄTZLICHER_WERT"
    )
  )  

  extra_values <- lapply(
    samples_list,
    function(x) {
      if (length(x) > 7137) {
        x[7138:length(x)]
      } else {
        character(0)
      }
    }
  )
  
  table(lengths(extra_values))  
 
  #Analysematrix
  analysis_matrix <- do.call(
    rbind,
    lapply(samples_list, function(x) {
      x[7:7137]
    })
  )
  dim(analysis_matrix)  
  colnames(analysis_matrix) <- features  
  rownames(analysis_matrix) <- metadata$Sample  

  analysis_matrix <- apply(
    analysis_matrix,
    2,
    as.numeric
  )  
  warnings()  

  #Betroffene features untersuchen
  sum(is.na(analysis_matrix))
  na_columns <- which(colSums(is.na(analysis_matrix)) > 0)
  length(na_columns)  
  features[na_columns]  

  expression_char <- do.call(
    rbind,
    lapply(samples_list, function(x) x[7:7137])
  )
  
  rownames(expression_char) <- metadata$Sample
  colnames(expression_char) <- features  

  problem_cells <- which(
    is.na(
      suppressWarnings(
        apply(expression_char, 2, as.numeric)
      )
    ),
    arr.ind = TRUE
  )  

  data.frame(
    Sample = rownames(expression_char)[problem_cells[, 1]],
    Feature = colnames(expression_char)[problem_cells[, 2]],
    OriginalValue = expression_char[problem_cells]
  )  

  # Expressionsmatrix als Character
  expression_char <- do.call(
    rbind,
    lapply(samples_list, function(x) {
      x[7:7137]
    })
  )
  
  rownames(expression_char) <- metadata$Sample
  colnames(expression_char) <- features  

  expression_char[expression_char == "-"] <- NA  

  analysis_matrix <- apply(
    expression_char,
    2,
    as.numeric
  )  

  colnames(analysis_matrix) <- features
  rownames(analysis_matrix) <- metadata$Sample  

  is.numeric(analysis_matrix)  

  dim(analysis_matrix)  

  sum(is.na(analysis_matrix))  

  
  #Fehlende Expressionswerte
  # Anzahl fehlender Expressionswerte pro Sample
  missing_per_sample <- rowSums(is.na(analysis_matrix))
  
  # Nur Samples mit fehlenden Werten anzeigen
  missing_per_sample[missing_per_sample > 0]

  sum(missing_per_sample > 0)  

  missing_per_feature <- colSums(is.na(analysis_matrix))
  
  missing_per_feature[missing_per_feature > 0]  

  18 / (72 * 7131) * 100  

  round(
    missing_per_feature[missing_per_feature > 0] / 72 * 100,
    2
  )  

  #Fehlende Werte:
  #Die Genexpressionsdaten wurden auf fehlende Werte überprüft. In der Rohdatei waren einzelne fehlende Expressionswerte durch das Zeichen - gekennzeichnet. Diese wurden bei der Aufbereitung als NA behandelt. Insgesamt wurden 18 fehlende Werte in 17 von 7.131 Features identifiziert. Die fehlenden Werte betreffen 16 der 72 Samples. Dies entspricht einem Anteil von lediglich 0,0035 % aller Expressionsmatrix-Einträge. Aufgrund des sehr geringen Anteils wurden die fehlenden Werte zunächst nicht imputiert und die betroffenen Samples bzw. Features nicht entfernt.
  
  total_missing <- sum(is.na(analysis_matrix))

  
  # Anteil fehlender Werte an der gesamten Matrix
  missing_percent <- total_missing /
    (nrow(analysis_matrix) * ncol(analysis_matrix)) * 100
  
  missing_percent  

  #Preprocessing
  
  # Die Rohdaten wurden zunächst in Metadaten und
  # Genexpressionsdaten aufgeteilt.
  #
  # Expressionsmatrix:
  # 72 Samples × 7131 Features
  #
  # Fehlende Werte:
  # "-" in der Rohdatei wurden als NA übernommen.
  #
  # Die 38 initialen Samples und 34 unabhängigen Samples
  # werden über metadata$Dataset unterschieden.
  #
  # Das Paper verwendet für den Class Predictor
  # logarithmierte Expressionswerte und anschließend
  # eine Standardisierung anhand des initialen Datensatzes.
  #
  # Da die vorliegende Expressionsmatrix bereits negative
  # Werte enthält, wird keine zusätzliche Log-Transformation
  # vorgenommen, bevor die Skala der bereitgestellten Daten
  # eindeutig geklärt ist.
  
  summary(as.vector(analysis_matrix))
  range(
    analysis_matrix,
    na.rm = TRUE
  )  
  mean(
    analysis_matrix,
    na.rm = TRUE
  )  
  sd(
    analysis_matrix,
    na.rm = TRUE
  )  

  #Ausreißer untersuchen
  max_value <- max(
    analysis_matrix,
    na.rm = TRUE
  )
  
  max_value
  
  which(
    analysis_matrix == max_value,
    arr.ind = TRUE
  )

  max_position <- which(
    analysis_matrix == max_value,
    arr.ind = TRUE
  )
  
  rownames(analysis_matrix)[max_position[1, "row"]]
  
  colnames(analysis_matrix)[max_position[1, "col"]]  

  analysis_matrix[
    max_position[1, "row"],
    1:20
  ]  

  sort(
    as.vector(analysis_matrix),
    decreasing = TRUE,
    na.last = NA
  )[1:20]  
 
  #Prüfung Feature Z70222_at 
  summary(analysis_matrix[, "Z70222_at"])

  colnames(analysis_matrix)[5197]  
  grep("Z7022", colnames(analysis_matrix), value = TRUE)  
  grep("Z70222", features, value = TRUE)  

  summary(analysis_matrix[, 5197])  
  sample35 <- samples_list[[which(metadata$Sample == 35)]]
  
  sample35[5203]  
  sort(
    analysis_matrix[, 5197],
    decreasing = TRUE,
    na.last = NA
  )[1:10]  
  
  large_pos <- which(
    abs(analysis_matrix) > 10000,
    arr.ind = TRUE
  )
  
  large_values <- data.frame(
    Sample = rownames(analysis_matrix)[large_pos[, "row"]],
    Feature = colnames(analysis_matrix)[large_pos[, "col"]],
    Value = analysis_matrix[large_pos]
  )
  
  large_values

  nrow(large_values)  

  #Suche nach extremen Werten
  extreme_pos <- which(
    abs(analysis_matrix) > 1e6,
    arr.ind = TRUE
  )
  
  extreme_values <- data.frame(
    Sample = rownames(analysis_matrix)[extreme_pos[, "row"]],
    Feature = colnames(analysis_matrix)[extreme_pos[, "col"]],
    Value = analysis_matrix[extreme_pos]
  )
  
  nrow(extreme_values)
  extreme_values  
  
  table(extreme_values$Feature)
  length(unique(extreme_values$Feature))  

  #Sind die Extremwerte ausschließlich positiv?
  sum(analysis_matrix > 1e6, na.rm = TRUE)
  sum(analysis_matrix < -1e6, na.rm = TRUE)  
  #Die extremen Werte treten ausschließlich auf der positiven Seite auf.
  
  #64 betroffene Features prüfen
  extreme_features <- unique(extreme_values$Feature)
  
  extreme_features
  
  extreme_summary <- data.frame(
    Feature = extreme_features,
    Median = sapply(
      extreme_features,
      function(f) median(analysis_matrix[, f], na.rm = TRUE)
    ),
    Maximum = sapply(
      extreme_features,
      function(f) max(analysis_matrix[, f], na.rm = TRUE)
    )
  )
  
  extreme_summary  
  extreme_summary[
    order(extreme_summary$Maximum, decreasing = TRUE),
  ]  
  
  #Gibt es Fehler in der Rohdatenrepräsentation?
  # Alle Werte von Z7022_at
  z7022 <- analysis_matrix[, 5197]
  
  z7022
  # Die größten Werte
  sort(z7022, decreasing = TRUE, na.last = NA)  
  sample35[5203]  
  sample49 <- samples_list[[which(metadata$Sample == 49)]]
  
  sample49[5203]  
  sample8 <- samples_list[[which(metadata$Sample == 8)]]
  
  sample8[5203]  

  #Fehlen Dezimalstellen systematisch?
  test_value <- as.numeric("95985426204822")
  
  test_value / 10^11
  
  grepl("\\.", sample35[5203])
  grepl("\\.", sample49[5203])  

  integer_extreme <- do.call(
    rbind,
    lapply(seq_along(samples_list), function(i) {
      
      x <- samples_list[[i]][7:length(samples_list[[i]])]
      
      is_integer_extreme <- grepl(
        "^-?[0-9]+$",
        x
      ) &
        abs(as.numeric(x)) > 1e6
      
      if (any(is_integer_extreme)) {
        
        data.frame(
          Sample = metadata$Sample[i],
          Feature = features[which(is_integer_extreme)],
          OriginalValue = x[is_integer_extreme],
          stringsAsFactors = FALSE
        )
        
      } else {
        NULL
      }
    })
  )
  
  nrow(integer_extreme)
  head(integer_extreme, 20)  
 
  #Positionen der Dezimalstellen bestimmen
  integer_extreme$CorrectedValue <- sapply(
    integer_extreme$OriginalValue,
    function(x) {
      x_num <- as.numeric(x)
      digits <- nchar(gsub("-", "", x))
      
      x_num / 10^(digits - 3)
    }
  )
  
  head(integer_extreme, 20)

  integer_extreme$FeatureMedian <- sapply(
    integer_extreme$Feature,
    function(f) {
      median(
        analysis_matrix[, f],
        na.rm = TRUE
      )
    }
  )
  
  integer_extreme[
    ,
    c("Sample", "Feature", "OriginalValue",
      "CorrectedValue", "FeatureMedian")
  ]  
 
  integer_extreme$FeatureMin <- sapply(
    seq_len(nrow(integer_extreme)),
    function(i) {
      
      f <- integer_extreme$Feature[i]
      s <- integer_extreme$Sample[i]
      
      values <- analysis_matrix[, f]
      
      # Den aktuell problematischen Extremwert ausschließen
      values <- values[
        seq_along(values) != which(metadata$Sample == s)
      ]
      
      min(values, na.rm = TRUE)
    }
  )

  integer_extreme$FeatureMax <- sapply(
    seq_len(nrow(integer_extreme)),
    function(i) {
      
      f <- integer_extreme$Feature[i]
      s <- integer_extreme$Sample[i]
      
      values <- analysis_matrix[, f]
      
      values <- values[
        seq_along(values) != which(metadata$Sample == s)
      ]
      
      max(values, na.rm = TRUE)
    }
  )  

  integer_extreme[
    ,
    c(
      "Sample",
      "Feature",
      "OriginalValue",
      "CorrectedValue",
      "FeatureMin",
      "FeatureMax"
    )
  ]  

  integer_extreme$Q01 <- sapply(
    seq_len(nrow(integer_extreme)),
    function(i) {
      
      f <- integer_extreme$Feature[i]
      s <- integer_extreme$Sample[i]
      
      values <- analysis_matrix[, f]
      values <- values[
        seq_along(values) != which(metadata$Sample == s)
      ]
      
      quantile(values, 0.01, na.rm = TRUE)
    }
  )

  integer_extreme$Q99 <- sapply(
    seq_len(nrow(integer_extreme)),
    function(i) {
      
      f <- integer_extreme$Feature[i]
      s <- integer_extreme$Sample[i]
      
      values <- analysis_matrix[, f]
      values <- values[
        seq_along(values) != which(metadata$Sample == s)
      ]
      
      quantile(values, 0.99, na.rm = TRUE)
    }
  )  

  integer_extreme[
    ,
    c(
      "Sample",
      "Feature",
      "CorrectedValue",
      "Q01",
      "Q99"
    )
  ]  

  integer_extreme$WithinQRange <-
    integer_extreme$CorrectedValue >= integer_extreme$Q01 &
    integer_extreme$CorrectedValue <= integer_extreme$Q99  

  table(integer_extreme$WithinQRange)  

  #Bei 100 extrem großen Werten wurde geprüft, ob eine Rekonstruktion des Dezimalpunkts nach den ersten drei Ziffern zu plausiblen Werten führt. Bei 91 % der Werte lag der rekonstruierte Wert innerhalb des 1.–99. Perzentilbereichs des jeweiligen Features. Bei 9 % lag er außerhalb dieses Bereichs. Die Rekonstruktion wird daher als plausibel eingestuft, wobei die 9 verbleibenden Werte gesondert geprüft werden.
  
  #Überprüfung der 9 Werte
  integer_extreme[
    integer_extreme$WithinQRange == FALSE,
    c(
      "Sample",
      "Feature",
      "OriginalValue",
      "CorrectedValue",
      "Q01",
      "Q99",
      "FeatureMedian"
    )
  ]

  #Erstellen einer Rangposition
  extreme_check <- integer_extreme[
    integer_extreme$WithinQRange == FALSE,
    c("Sample", "Feature", "OriginalValue", "CorrectedValue")
  ]

  extreme_check$Min <- NA
  extreme_check$Max <- NA
  extreme_check$Median <- NA
  extreme_check$RankPercent <- NA  

  for (i in 1:nrow(extreme_check)) {
    
    feature <- extreme_check$Feature[i]
    sample <- extreme_check$Sample[i]
    
    # Werte des Features
    values <- analysis_matrix[, feature]
    
    # problematischen Originalwert entfernen
    sample_position <- which(metadata$Sample == sample)
    values <- values[-sample_position]
    
    values <- values[!is.na(values)]
    
    corrected <- extreme_check$CorrectedValue[i]
    
    extreme_check$Min[i] <- min(values)
    extreme_check$Max[i] <- max(values)
    extreme_check$Median[i] <- median(values)
    
    # Rang des korrigierten Wertes
    extreme_check$RankPercent[i] <-
      mean(values <= corrected) * 100
  }

  extreme_check  

  #Gibt es eine Dezimalposition, bei der der Wert in den beobachteten Bereich des Features fällt?
  test_decimal_positions <- function(x, min_value, max_value) {
    
    x_clean <- gsub("-", "", x)
    
    results <- data.frame(
      Position = integer(),
      Value = numeric(),
      InsideRange = logical()
    )
    
    for (p in 1:(nchar(x_clean) - 1)) {
      
      value <- as.numeric(x) / 10^(
        nchar(x_clean) - p
      )
      
      results <- rbind(
        results,
        data.frame(
          Position = p,
          Value = value,
          InsideRange = value >= min_value &
            value <= max_value
        )
      )
    }
    
    results
  }

  test_decimal_positions(
    "95985426204822",
    -2132.5879,
    764.8682
  )  

  test_decimal_positions(
    "012976373338",
    -4322.6894,
    -170.2612
  )  
  
  #Vorherige Hypothese kann nicht vollständig belegt werden
  
  u60319 <- analysis_matrix[, "U60319_at"]
  summary(u60319)
  sort(u60319)  
  z70222 <- analysis_matrix[, "Z70222_at"]
  summary(z70222)  
  sort(z70222)  
  
  #Die Dezimalstellen werden nicht mehr automatisch rekonstruiert, da keine zuverlässige Aussage über die ursprüngliche Position des Dezimalpunktes getroffen werden kann.
  
  #Bereinigung der Analysematrix
  analysis_matrix_clean <- analysis_matrix

  for (i in 1:nrow(integer_extreme)) {
    
    sample_id <- integer_extreme$Sample[i]
    feature_id <- integer_extreme$Feature[i]
    
    row <- which(metadata$Sample == sample_id)
    col <- which(colnames(analysis_matrix_clean) == feature_id)
    
    analysis_matrix_clean[row, col] <- NA
  }  
  sum(is.na(analysis_matrix_clean))  

  sum(analysis_matrix_clean > 1e6, na.rm = TRUE)  
  sum(analysis_matrix_clean < -1e6, na.rm = TRUE)  

  dim(analysis_matrix_clean)  

  analysis_matrix_clean[
    which(metadata$Sample == 35),
    "Z70222_at"
  ]  

  analysis_matrix_clean[
    which(metadata$Sample == 5),
    "U60319_at"
  ]  

  #Test
  summary(as.vector(analysis_matrix_clean))
  sum(is.na(analysis_matrix_clean))  

  #Die bereinigte Matrix enthält: 72 Samples, 7.131 Features, 118 fehlende Werte (NA), keine Werte > 1e6, keine Werte < -1e6 
  
  max_clean <- max(
    analysis_matrix_clean,
    na.rm = TRUE
  )
  
  max_clean

  which(
    analysis_matrix_clean == max_clean,
    arr.ind = TRUE
  )  

  #Auffälligen Wert kontrollieren
  colnames(analysis_matrix_clean)[1723]
  analysis_matrix_clean[
    which(metadata$Sample == 48),
    1723
  ]  
  summary(analysis_matrix_clean[, 1723])  

  sort(
    analysis_matrix_clean[, "M14539_at"],
    decreasing = TRUE,
    na.last = NA
  )[1:10]  

  samples_list[[which(metadata$Sample == 48)]][1723]  

  #Korrekten Wert übernehmen
  analysis_matrix_clean[
    which(metadata$Sample == 48),
    "M14539_at"
  ] <- as.numeric(
    samples_list[[which(metadata$Sample == 48)]][1723]
  )

  #Kontrolle
  analysis_matrix_clean[
    which(metadata$Sample == 48),
    "M14539_at"
  ]

  summary(analysis_matrix_clean[, "M14539_at"])  

  sort(
    analysis_matrix_clean[, "M14539_at"],
    decreasing = TRUE,
    na.last = NA
  )[1:10]  

  which(
    analysis_matrix_clean[, "M14539_at"] == 439437,
    arr.ind = TRUE
  )  
  metadata$Sample[15]  
  colnames(analysis_matrix_clean)[49]  

  sample_id <- metadata$Sample[15]
  samples_list[[which(metadata$Sample == sample_id)]][49]  

  analysis_matrix_clean[
    which(metadata$Sample == 68),
    "AFFX-HUMTFRR/M11507_M_at"
  ] <- as.numeric(
    samples_list[[which(metadata$Sample == 68)]][49]
  )  

  analysis_matrix_clean[
    which(metadata$Sample == 68),
    "AFFX-HUMTFRR/M11507_M_at"
  ]  

  max(
    analysis_matrix_clean,
    na.rm = TRUE
  )  

  sum(
    analysis_matrix_clean > 1e6,
    na.rm = TRUE
  )  

  max_position <- which(
    analysis_matrix_clean == 939111,
    arr.ind = TRUE
  )
  
  max_position  

  metadata$Sample[max_position[1, "row"]]  
  colnames(analysis_matrix_clean)[max_position[1, "col"]]  

  sample_id <- metadata$Sample[max_position[1, "row"]]
  
  samples_list[[which(metadata$Sample == sample_id)]][
    max_position[1, "col"]
  ]  

  analysis_matrix_clean[
    which(metadata$Sample == 4),
    "U59325_at"
  ]  

  analysis_matrix_clean[
    which(metadata$Sample == 4),
    "U59325_at"
  ] <- as.numeric(
    samples_list[[which(metadata$Sample == 4)]][3464]
  )  

  analysis_matrix_clean[
    which(metadata$Sample == 4),
    "U59325_at"
  ]  

  colnames(analysis_matrix_clean)[3464]  
  features[3464]  
  samples_list[[which(metadata$Sample == 4)]][
    which(header == features[3464])
  ]  

  sample_id <- 48
  feature_id <- "M14539_at"
  
  samples_list[[which(metadata$Sample == sample_id)]][
    which(header == feature_id)
  ]  

  sample_id <- 68
  feature_id <- "AFFX-HUMTFRR/M11507_M_at"
  
  samples_list[[which(metadata$Sample == sample_id)]][
    which(header == feature_id)
  ]  
  
  analysis_matrix_clean[
    which(metadata$Sample == 68),
    "AFFX-HUMTFRR/M11507_M_at"
  ] <- as.numeric(
    samples_list[[which(metadata$Sample == 68)]][
      which(header == "AFFX-HUMTFRR/M11507_M_at")
    ]
  )

  analysis_matrix_clean[
    which(metadata$Sample == 48),
    "M14539_at"
  ] <- as.numeric(
    samples_list[[which(metadata$Sample == 48)]][
      which(header == "M14539_at")
    ]
  )  

  analysis_matrix_clean[
    which(metadata$Sample == 4),
    "U59325_at"
  ]  

  analysis_matrix_clean[
    which(metadata$Sample == 4),
    "U59325_at"
  ] <- as.numeric(
    samples_list[[which(metadata$Sample == 4)]][
      which(header == "U59325_at")
    ]
  )  

  analysis_matrix_clean[
    which(metadata$Sample == 4),
    "U59325_at"
  ]  

  analysis_matrix_clean[
    which(metadata$Sample == 68),
    "AFFX-HUMTFRR/M11507_M_at"
  ] <- as.numeric(
    samples_list[[which(metadata$Sample == 68)]][
      which(header == "AFFX-HUMTFRR/M11507_M_at")
    ]
  )  

  analysis_matrix_clean[
    which(metadata$Sample == 48),
    "M14539_at"
  ] <- as.numeric(
    samples_list[[which(metadata$Sample == 48)]][
      which(header == "M14539_at")
    ]
  )  

  analysis_matrix_clean[
    which(metadata$Sample == 68),
    "AFFX-HUMTFRR/M11507_M_at"
  ]
  
  analysis_matrix_clean[
    which(metadata$Sample == 48),
    "M14539_at"
  ]  

  #Vergleich analysis_matrix_clean mit Originaldaten
  analysis_matrix_original <- do.call(
    rbind,
    lapply(samples_list, function(x) {
      as.numeric(x[7:7137])
    })
  )

  colnames(analysis_matrix_original) <- features
  rownames(analysis_matrix_original) <- metadata$Sample  
  dim(analysis_matrix_original)  
  sum(is.na(analysis_matrix_original))  

  warnings()  

  #nicht-numerische Expressionswerte suchen
  non_numeric_values <- list()
  
  for (i in seq_along(samples_list)) {
    
    values <- samples_list[[i]][7:7137]
    numeric_values <- suppressWarnings(as.numeric(values))
    
    bad <- which(is.na(numeric_values) & !is.na(values) & values != "")
    
    if (length(bad) > 0) {
      non_numeric_values[[i]] <- data.frame(
        Sample = metadata$Sample[i],
        Position = bad,
        Feature = features[bad],
        OriginalValue = values[bad],
        stringsAsFactors = FALSE
      )
    }
  }

  non_numeric_values  
  
  #Kontrolle der 18 NAs
  table(
    is.na(analysis_matrix_original)
  )
  sum(is.na(analysis_matrix_original))  
  sum(is.na(analysis_matrix_clean))  
  
  #Originalwerte identifizieren, die als NA dargestellt sind
  additional_na <- is.na(analysis_matrix_clean) &
    !is.na(analysis_matrix_original)
  
  sum(additional_na)

  additional_na_positions <- which(
    additional_na,
    arr.ind = TRUE
  )
  
  additional_na_table <- data.frame(
    Sample = rownames(analysis_matrix_clean)[additional_na_positions[, "row"]],
    Feature = colnames(analysis_matrix_clean)[additional_na_positions[, "col"]],
    OriginalValue = analysis_matrix_original[additional_na],
    stringsAsFactors = FALSE
  )
  
  head(additional_na_table, 20)  

  summary(additional_na_table$OriginalValue)  

  range(
    additional_na_table$OriginalValue,
    na.rm = TRUE
  )  

  table(additional_na_table$Feature)  
  length(unique(additional_na_table$Feature))  
  table(additional_na_table$Sample)  

  additional_na_table$Digits <- nchar(
    format(
      additional_na_table$OriginalValue,
      scientific = FALSE,
      trim = TRUE
    )
  )
  
  table(additional_na_table$Digits)  

  head(
    additional_na_table[
      order(additional_na_table$OriginalValue, decreasing = TRUE),
    ],
    20
  ) 
  
  #Vergleich innerhalb desselben Features
  z70222 <- analysis_matrix_original[, "Z70222_at"]
  sort(z70222, decreasing = TRUE, na.last = NA)
 
  z70222_clean <- z70222[z70222 < 10000]
  summary(z70222_clean) 
  quantile(
    z70222_clean,
    probs = c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99),
    na.rm = TRUE
  )  

  mean(z70222_clean < 95.985426204822)  
  mean(z70222_clean <= 95.985426204822)  
  sum(z70222_clean >= 95.985426204822)  

  feature_medians <- apply(
    analysis_matrix_original,
    2,
    median,
    na.rm = TRUE
  )  
  additional_na_table$FeatureMedian <- feature_medians[
    match(additional_na_table$Feature,
          names(feature_medians))
  ]  
  additional_na_table$RatioToMedian <- 
    abs(additional_na_table$OriginalValue) /
    abs(additional_na_table$FeatureMedian)  
  summary(additional_na_table$RatioToMedian)  

  additional_na_table$Corrected_12 <- 
    additional_na_table$OriginalValue / 10^12
  
  head(
    additional_na_table[
      order(additional_na_table$OriginalValue, decreasing = TRUE),
    ],
    20
  )  
  
  #Wertebereich bestimmen, die nicht als Extremwerte gelten
  feature_q01 <- apply(
    analysis_matrix_original,
    2,
    quantile,
    probs = 0.01,
    na.rm = TRUE
  )
  
  feature_q99 <- apply(
    analysis_matrix_original,
    2,
    quantile,
    probs = 0.99,
    na.rm = TRUE
  )

  #Mögliche Korrekturen
  additional_na_table$Q01 <- feature_q01[
    match(additional_na_table$Feature,
          names(feature_q01))
  ]
  
  additional_na_table$Q99 <- feature_q99[
    match(additional_na_table$Feature,
          names(feature_q99))
  ]

  #Testen der Dezimalverschiebungen
  additional_na_table$Corrected_9  <- additional_na_table$OriginalValue / 10^9
  additional_na_table$Corrected_10 <- additional_na_table$OriginalValue / 10^10
  additional_na_table$Corrected_11 <- additional_na_table$OriginalValue / 10^11
  additional_na_table$Corrected_12 <- additional_na_table$OriginalValue / 10^12
  additional_na_table$Corrected_13 <- additional_na_table$OriginalValue / 10^13

  additional_na_table$InRange_9 <-
    additional_na_table$Corrected_9 >= additional_na_table$Q01 &
    additional_na_table$Corrected_9 <= additional_na_table$Q99
  
  additional_na_table$InRange_10 <-
    additional_na_table$Corrected_10 >= additional_na_table$Q01 &
    additional_na_table$Corrected_10 <= additional_na_table$Q99
  
  additional_na_table$InRange_11 <-
    additional_na_table$Corrected_11 >= additional_na_table$Q01 &
    additional_na_table$Corrected_11 <= additional_na_table$Q99
  
  additional_na_table$InRange_12 <-
    additional_na_table$Corrected_12 >= additional_na_table$Q01 &
    additional_na_table$Corrected_12 <= additional_na_table$Q99
  
  additional_na_table$InRange_13 <-
    additional_na_table$Corrected_13 >= additional_na_table$Q01 &
    additional_na_table$Corrected_13 <= additional_na_table$Q99  

  table(additional_na_table$InRange_9)
  table(additional_na_table$InRange_10)
  table(additional_na_table$InRange_11)
  table(additional_na_table$InRange_12)
  table(additional_na_table$InRange_13)  

  feature_q01_clean <- apply(
    analysis_matrix_clean,
    2,
    quantile,
    probs = 0.01,
    na.rm = TRUE
  )
  
  feature_q99_clean <- apply(
    analysis_matrix_clean,
    2,
    quantile,
    probs = 0.99,
    na.rm = TRUE
  )  

  additional_na_table$Q01_clean <- feature_q01_clean[
    match(
      additional_na_table$Feature,
      names(feature_q01_clean)
    )
  ]
  
  additional_na_table$Q99_clean <- feature_q99_clean[
    match(
      additional_na_table$Feature,
      names(feature_q99_clean)
    )
  ]  

  additional_na_table$InRange_9_clean <-
    additional_na_table$Corrected_9 >= additional_na_table$Q01_clean &
    additional_na_table$Corrected_9 <= additional_na_table$Q99_clean
  
  additional_na_table$InRange_10_clean <-
    additional_na_table$Corrected_10 >= additional_na_table$Q01_clean &
    additional_na_table$Corrected_10 <= additional_na_table$Q99_clean
  
  additional_na_table$InRange_11_clean <-
    additional_na_table$Corrected_11 >= additional_na_table$Q01_clean &
    additional_na_table$Corrected_11 <= additional_na_table$Q99_clean
  
  additional_na_table$InRange_12_clean <-
    additional_na_table$Corrected_12 >= additional_na_table$Q01_clean &
    additional_na_table$Corrected_12 <= additional_na_table$Q99_clean
  
  additional_na_table$InRange_13_clean <-
    additional_na_table$Corrected_13 >= additional_na_table$Q01_clean &
    additional_na_table$Corrected_13 <= additional_na_table$Q99_clean  

  table(additional_na_table$InRange_9_clean)
  
  table(additional_na_table$InRange_10_clean)
  
  table(additional_na_table$InRange_11_clean)
  
  table(additional_na_table$InRange_12_clean)
  
  table(additional_na_table$InRange_13_clean)  

  additional_na_table[
    additional_na_table$InRange_12_clean == FALSE,
    c(
      "Sample",
      "Feature",
      "OriginalValue",
      "Corrected_9",
      "Corrected_10",
      "Corrected_11",
      "Corrected_12",
      "Corrected_13",
      "Q01_clean",
      "Q99_clean"
    )
  ]  

  table(
    additional_na_table$InRange_12_clean &
      !additional_na_table$InRange_13_clean
  )  
  table(
    additional_na_table$InRange_13_clean &
      !additional_na_table$InRange_12_clean
  )  

  additional_na_table[
    additional_na_table$InRange_12_clean == TRUE &
      additional_na_table$InRange_13_clean == TRUE,
    c(
      "Sample",
      "Feature",
      "OriginalValue",
      "Corrected_11",
      "Corrected_12",
      "Corrected_13",
      "Q01_clean",
      "Q99_clean"
    )
  ]  

  feature_medians_clean <- apply(
    analysis_matrix_clean,
    2,
    median,
    na.rm = TRUE
  )  

  additional_na_table$MedianClean <- feature_medians_clean[
    match(
      additional_na_table$Feature,
      names(feature_medians_clean)
    )
  ]  

  additional_na_table$Dist_9 <- abs(
    additional_na_table$Corrected_9 -
      additional_na_table$MedianClean
  )
  
  additional_na_table$Dist_10 <- abs(
    additional_na_table$Corrected_10 -
      additional_na_table$MedianClean
  )
  
  additional_na_table$Dist_11 <- abs(
    additional_na_table$Corrected_11 -
      additional_na_table$MedianClean
  )
  
  additional_na_table$Dist_12 <- abs(
    additional_na_table$Corrected_12 -
      additional_na_table$MedianClean
  )
  
  additional_na_table$Dist_13 <- abs(
    additional_na_table$Corrected_13 -
      additional_na_table$MedianClean
  )  
  
  additional_na_table$Dist_13 <- abs(
    additional_na_table$Corrected_13 -
      additional_na_table$MedianClean
  )
  
  additional_na_table[
    additional_na_table$Feature == "Z70222_at",
    c(
      "OriginalValue",
      "MedianClean",
      "Corrected_9",
      "Corrected_10",
      "Corrected_11",
      "Corrected_12",
      "Corrected_13",
      "Dist_9",
      "Dist_10",
      "Dist_11",
      "Dist_12",
      "Dist_13"
    )
  ]

  sum(is.na(analysis_matrix_clean))  
  sum(analysis_matrix_clean > 1e6, na.rm = TRUE)  
  sum(analysis_matrix_clean < -1e6, na.rm = TRUE)  
  summary(as.vector(analysis_matrix_clean))  

  table(
    metadata$Diagnosis,
    useNA = "ifany"
  )  
  table(
    metadata$Dataset,
    metadata$Diagnosis,
    useNA = "ifany"
  )  
  all(rownames(analysis_matrix_clean) == metadata$Sample)  

  #Finale Matrix
  analysis_matrix <- analysis_matrix_clean
  dim(analysis_matrix)  
  sum(is.na(analysis_matrix))  
  stopifnot(
    nrow(analysis_matrix) == nrow(metadata),
    ncol(analysis_matrix) == length(features),
    all(rownames(analysis_matrix) == metadata$Sample),
    all(colnames(analysis_matrix) == features)
  )  
  