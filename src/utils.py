from pathlib import Path
import os
import pandas as pd
import numpy as np


# ORDNER
project_root = Path(__file__).resolve().parents[1]
input_dir = project_root / "data"
output_dir = project_root / "output"
output_dir.mkdir(exist_ok=True)

# UNSERE DATEIEN
my_team_files = [
    "ceyssens_2014/Ceyssens_directional_full_raw_counts.tsv",
    "finstrlova_2022/Finstrlova_SH1000_full_raw_counts.tsv",
    "guegler_2021/Guegler_T4_plusToxIN_full_raw_counts.tsv",
    "kuptsov_2022/Kuptsov_full_raw_counts.tsv",
    "meaden_2021/Meaden_BIM_full_raw_counts.tsv",
    "sprenger_2024/Sprenger_VC_WT_VP882_WT_full_raw_counts.tsv",
    "zhong_2020/Zhong_full_raw_counts.tsv"
]

# HILFSFUNKTIONEN FÜR ERKENNUNG VON AUSREIßER UND NORMALISIERUNG
def detect_outliers_iqr(df):
    Q1 = df.quantile(0.25)
    Q3 = df.quantile(0.75)
    IQR = Q3 - Q1
    lower = Q1 - 1.5 * IQR
    upper = Q3 + 1.5 * IQR
    return (df < lower) | (df > upper)

def normalize_log1p(df):
    return np.log1p(df)

# NUR UNSERE DATEIEN VERARBEITEN
# === DATEIEN VERARBEITEN ===
for file in my_team_files:
    input_file = input_dir / file
    relative_path = input_file.relative_to(input_dir)
    output_file = output_dir / relative_path.with_name(relative_path.stem + "_cleaned.csv")
    output_file.parent.mkdir(parents=True, exist_ok=True)

    if not input_file.exists():
        print(f"Datei nicht gefunden: {file}")
        continue  

    # DATEN EINLESEN 
    df = pd.read_csv(input_file, sep='\t', index_col=0)
    numeric_df = df.select_dtypes(include=[np.number])
    gene_metadata = df.drop(columns=numeric_df.columns, errors='ignore')

    # AUSREIßER ERKENNEN UND ENTFERNEN 
    outliers = detect_outliers_iqr(numeric_df)
    cleaned_df = numeric_df.mask(outliers)
    n_outliers = outliers.sum().sum()

    # NORMALISIERUNG (To be edited eventually for User Story 3)
    normalized_df = normalize_log1p(cleaned_df)

    # ERGEBNIS SPEICHERN UND AUSGEBEN ===
    final_df = pd.concat([gene_metadata, normalized_df], axis=1)
    final_df.to_csv(output_file)

    print(" Datei verarbeitet:", input_file)
    print(" Ausreißer entfernt:", n_outliers)
    print(" Gespeichert als:", output_file)