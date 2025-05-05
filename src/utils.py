import os
import pandas as pd
import numpy as np
from scipy.stats import zscore

# === PARAMETER ===
input_dir = "data/Ceyssens_directional_full_raw_counts.tsv"           # Ordner mit Count-Tabellen (*.tsv)
output_dir = "output/Ceyssens_cleaned.csv"        # Ergebnisordner für bereinigte Dateien und Bericht
normalize = True             # True = Normalisierung aktivieren (log1p)

# === HILFSFUNKTIONEN ===
def detect_outliers_iqr(df):
    Q1 = df.quantile(0.25)
    Q3 = df.quantile(0.75)
    IQR = Q3 - Q1
    lower = Q1 - 1.5 * IQR
    upper = Q3 + 1.5 * IQR
    return (df < lower) | (df > upper)

def normalize_log1p(df):
    return np.log1p(df)


# === DATEN EINLESEN ===
df = pd.read_csv(input_file, sep='\t', index_col=0)
numeric_df = df.select_dtypes(include=[np.number])
gene_metadata = df.drop(columns=numeric_df.columns, errors='ignore')

# === AUSREIßER ERKENNEN UND ENTFERNEN ===
outliers = detect_outliers_iqr(numeric_df)
cleaned_df = numeric_df.mask(outliers)
n_outliers = outliers.sum().sum()

# === NORMALISIERUNG ===
normalized_df = normalize_log1p(cleaned_df)

# === ERGEBNIS SPEICHERN ===
final_df = pd.concat([gene_metadata, normalized_df], axis=1)
final_df.to_csv(output_file)

print(" Datei verarbeitet:", input_file)
print(" Ausreißer entfernt:", n_outliers)
print(" Gespeichert als:", output_file)