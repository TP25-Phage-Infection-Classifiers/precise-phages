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

# Hilfsfunktionen für Erkennung der Ausreißer, Normalisierung und Klassifizierung
def detect_outliers_iqr(df):
    Q1 = df.quantile(0.25)
    Q3 = df.quantile(0.75)
    IQR = Q3 - Q1
    lower = Q1 - 1.5 * IQR
    upper = Q3 + 1.5 * IQR
    return (df < lower) | (df > upper)

def normalize_log1p(df):
    return np.log1p(df)

# ordnet jedem Gen ein Label zu, arbeitet mit den normalisierten, bereinigten Phagengendaten
def classify_temporal_expression(row, time_cols):
    max_timepoint_idx = row[time_cols].values.argmax() # Index des höchsten Werts der Zeile
    n_timepoints = len(time_cols) # Anzahl der Werte
    if max_timepoint_idx < n_timepoints / 3: # wenn Index im ersten Drittel liegt -> early
        return 'early'
    elif max_timepoint_idx < 2 * n_timepoints / 3: # wenn Index im zweiten Drittel liegt -> middle
        return 'middle'
    else: # wenn Index im dritten Drittel liegt -> late
        return 'late'

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

    # Normalisierung
    normalized_df = normalize_log1p(numeric_df)
    
    # Normalisierte Daten speichern
    normalized_final_df = pd.concat([gene_metadata, normalized_df], axis=1)
    normalized_final_df.to_csv(output_normalized_file)
    
    # Nur Phagengene rausfiltern
    phage_df = normalized_final_df[normalized_final_df["Entity"] == "phage"]
    numeric_phage_df = phage_df.select_dtypes(include=[np.number])
    phage_metadata = phage_df.drop(columns=numeric_phage_df.columns, errors='ignore')
    
    # Nur normalisierte Phagendaten speichern
    phage_output_file = output_normalized_file.with_name(output_normalized_file.stem + "_phage_only.csv")
    phage_df.to_csv(phage_output_file, index=True)

    # Ausreißer erkennen und entfernen
    outliers = detect_outliers_iqr(numeric_phage_df)
    cleaned_df = numeric_phage_df.mask(outliers)
    n_outliers = outliers.sum().sum()
    
    # Bereinigte Daten speichern
    cleaned_final_df = pd.concat([phage_metadata, cleaned_df], axis=1)
    cleaned_final_df.to_csv(output_cleaned_file)
    
    # Einteilung der Gene
    time_cols = list(cleaned_final_df.columns[3:]) #betrachtet erste Zeile ab Spalte 4
    cleaned_final_df["Temporal_Class"] = cleaned_final_df.apply(lambda row: classify_temporal_expression(row, time_cols), axis=1) # ordnet jedem Gen ein Label zu
    
    # Speichere Dateien mit Labels 'early', 'middle', 'late' in neue Datei
    output_labeled_file = output_dir / relative_path.with_name(relative_path.stem + "_labeled.csv")
    cleaned_final_df.to_csv(output_labeled_file)
    
    # BOXPLOTS ERSTELLEN - es werden normalisierte Daten mit rohen Daten verglichen
    # Extrahiere nur die numerischen Spalten (Genexpressionswerte)
    raw_data = df.select_dtypes(include=[float, int])
    #cleaned_data = cleaned_final_df.select_dtypes(include=[float, int])
    normalized_data = normalized_final_df.select_dtypes(include=[float, int])

    # Schmelze die Daten für die Boxplot-Darstellung
    raw_data_melted = raw_data.melt(var_name='Gene', value_name='Expression')
    raw_data_melted['Condition'] = 'Raw'
    #cleaned_data_melted = cleaned_data.melt(var_name='Gene', value_name='Expression')
    #cleaned_data_melted['Condition'] = 'Cleaned'

    normalized_data_melted = normalized_data.melt(var_name='Gene', value_name='Expression')
    normalized_data_melted['Condition'] = 'Normalized'

    # Kombiniere beide DataFrames
    combined_data = pd.concat([raw_data_melted, normalized_data_melted])

    # Boxplot erstellen
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))  # Zwei Subplots nebeneinander

    # Boxplot für rohe Daten
    sb.boxplot(x='Condition', y='Expression', data=raw_data_melted, ax=axes[0], color='#B0E0E6')
    axes[0].set_title(f'{file} - Vor der Normalisierung')
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


    print(" Datei verarbeitet:", input_file)
    print(" Nach Normalisierung:", output_normalized_file)
    print(f" Boxplot für {file} gespeichert als {output_file}")
    print(" Normalisierte Phagengene: ", phage_output_file)
    print(" Bereinigte Phagengene: ", output_cleaned_file)
    print(" Gelabelte Phagengene: ", output_labeled_file)
    print(" Ausreißer entfernt:", n_outliers)




    # Testen, ob wo der Branch landet
    #########



   
