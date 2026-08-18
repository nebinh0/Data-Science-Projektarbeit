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
  