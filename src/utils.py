from pathlib import Path
import os
import pandas as pd
import numpy as np
import seaborn as sb
import matplotlib.pyplot as plt


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

# ordnet jedem Gen ein Label zu, arbeitet mit den normalisierten, bereinigten Daten
def classify_temporal_expression(row, time_cols):
    max_timepoint_idx = row[time_cols].values.argmax() # Index des höchsten Werts der Zeile
    n_timepoints = len(time_cols) # Anzahl der Werte
    if max_timepoint_idx < n_timepoints / 3: # wenn Index im ersten Drittel liegt -> early
        return 'early'
    elif max_timepoint_idx < 2 * n_timepoints / 3: # wenn Index im zweiten Drittel liegt -> middle
        return 'middle'
    else: # wenn Index im dritten Drittel liegt -> late
        return 'late'

# NUR UNSERE DATEIEN VERARBEITEN
# === DATEIEN VERARBEITEN ===
for file in my_team_files:
    input_file = input_dir / file
    relative_path = input_file.relative_to(input_dir)
    output_cleaned_file = output_dir / relative_path.with_name(relative_path.stem + "_cleaned.csv")
    output_normalized_file = output_dir / relative_path.with_name(relative_path.stem + "_normalized.csv")
    output_cleaned_file.parent.mkdir(parents=True, exist_ok=True)
    output_normalized_file.parent.mkdir(parents=True, exist_ok=True)
    output_dir = Path(__file__).resolve().parents[1] / "output"
    boxplot_output_dir = output_dir / "boxplots"
    boxplot_output_dir.mkdir(exist_ok=True) 

    if not input_file.exists():
        print(f"Datei nicht gefunden: {file}")
        continue  

    # DATEN EINLESEN 
    df = pd.read_csv(input_file, sep='\t', index_col=0)
    numeric_df = df.select_dtypes(include=[np.number])
    gene_metadata = df.drop(columns=numeric_df.columns, errors='ignore')

    # NORMALISIERUNG (To be edited eventually for User Story 3)
    normalized_df = normalize_log1p(numeric_df)
    
     # AUSREIßER ERKENNEN UND ENTFERNEN 
    outliers = detect_outliers_iqr(normalized_df)
    cleaned_df = normalized_df.mask(outliers)
    n_outliers = outliers.sum().sum()

    # ERGEBNIS SPEICHERN UND AUSGEBEN ===
    
    # 1. Bereinigte Daten speichern
    cleaned_final_df = pd.concat([gene_metadata, cleaned_df], axis=1)
    cleaned_final_df.to_csv(output_cleaned_file)
    
    # 2. Normalisierte Daten speichern
    normalized_final_df = pd.concat([gene_metadata, normalized_df], axis=1)
    normalized_final_df.to_csv(output_normalized_file)
    
    # Einteilung der Gene
    time_cols = list(normalized_final_df.columns[3:]) #betrachtet erste Zeile ab Spalte 4
    normalized_final_df["Temporal_Class"] = normalized_final_df.apply(lambda row: classify_temporal_expression(row, time_cols), axis=1) # ordnet jedem Gen ein Label zu
    
    # BOXPLOTS ERSTELLEN
    
    # Extrahiere nur die numerischen Spalten (Genexpressionswerte)
    cleaned_data = cleaned_final_df.select_dtypes(include=[float, int])
    normalized_data = normalized_final_df.select_dtypes(include=[float, int])

    # Schmelze die Daten für die Boxplot-Darstellung
    cleaned_data_melted = cleaned_data.melt(var_name='Gene', value_name='Expression')
    cleaned_data_melted['Condition'] = 'Cleaned'

    normalized_data_melted = normalized_data.melt(var_name='Gene', value_name='Expression')
    normalized_data_melted['Condition'] = 'Normalized'

    # Kombiniere beide DataFrames
    combined_data = pd.concat([cleaned_data_melted, normalized_data_melted])

    # Boxplot erstellen
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))  # Zwei Subplots nebeneinander

    # Boxplot für bereinigte Daten
    sb.boxplot(x='Condition', y='Expression', data=cleaned_data_melted, ax=axes[0], color='#B0E0E6')
    axes[0].set_title(f'{file} - Vor der Normalisierung (bereinigt)')
    axes[0].set_ylabel('Genexpression')

    # Boxplot für normalisierte Daten
    sb.boxplot(x='Condition', y='Expression', data=normalized_data_melted, ax=axes[1], color='#B0E0E6')
    axes[1].set_title(f'{file} - Nach der Normalisierung')
    axes[1].set_ylabel('Genexpression')

    # Layout anpassen
    plt.tight_layout()

    # Ausgabe-Dateiname generieren
    output_file = boxplot_output_dir / f"{Path(file).stem}_boxplot.png"
    
    # Speichern der Boxplots als PNG-Datei
    plt.savefig(output_file)

    # Speichere Dateien mit Labels 'early', 'middle', 'late' in neue Datei
    output_labeled_file = output_dir / relative_path.with_name(relative_path.stem + "_labeled.csv")
    normalized_final_df.to_csv(output_labeled_file)
    
    print(" Datei verarbeitet:", input_file)
    print(" Ausreißer entfernt:", n_outliers)
    print(" Vor Normalisierung:", output_cleaned_file)
    print(" Nach Normalisierung:", output_normalized_file)
    print(f" Boxplot für {file} gespeichert als {output_file}")