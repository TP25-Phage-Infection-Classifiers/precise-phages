import pandas as pd
import numpy as np
from pathlib import Path
import seaborn as sb


import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


import preprocess as pre
import label
import utils






def draw_piechart(label_counts, file):
    colors = ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3']
    fig, ax = plt.subplots(figsize=(5, 5)) # Pie Plot mit matplotlib direkt, bessere Kontrolle
    # Prozentwerte berechnen
    total = label_counts.sum() # Anzahl aller gelabelten Gene
    percentages = (label_counts / total * 100).round(1) # Verteilung in Prozent auf eine Nachkommastelle gerundet
    legend_labels = [f"{label}  {pct}%" for label, pct in zip(label_counts.index, percentages)]
    
    wedges, _ = ax.pie( 
        label_counts,
        startangle=90, # startet oben
        counterclock=False, # im Uhrzeigersinn
        colors=colors,
        textprops={'fontsize': 8}
        )   
    # Legende mit Label + Prozent
    ax.legend(wedges, legend_labels, title="Labels", loc="center left", bbox_to_anchor=(1, 0, 0.5, 1))
    
    ax.set_title(f'Label-Verteilung Diagramm\n{file}', fontsize=10) # Titel
    ax.axis('equal') # Kreis statt Oval
    
    plt.tight_layout() # Layout anpassen
    output_pie_file = pie_output_dir / f"{Path(file).stem}_pie.png" # Ausgabe-Datei generieren
    plt.savefig(output_pie_file) # Kuchendiagramm speichern als PNG-Datei








# === DATEIEN VERARBEITEN ===
boxplot_output_dir = utils.output_dir / "boxplots"
pie_output_dir = utils.output_dir / "pie chart"
boxplot_output_dir.mkdir(exist_ok=True) 
pie_output_dir.mkdir(exist_ok=True)



def generate_output():
    for file in utils.my_team_files:
        input_file = utils.input_dir / file
        relative_path = input_file.relative_to(utils.input_dir)
   
        output_cleaned_file = utils.output_dir / relative_path.with_name(relative_path.stem + "_cleaned.csv")
        output_normalized_file = utils.output_dir / relative_path.with_name(relative_path.stem + "_normalized.csv")
        output_cleaned_file.parent.mkdir(parents=True, exist_ok=True)
        output_normalized_file.parent.mkdir(parents=True, exist_ok=True)

    
    
    
        if not input_file.exists():
            print(f"Datei nicht gefunden: {file}")
            continue  

        # DATEN EINLESEN 
        df = pd.read_csv(input_file, sep='\t', index_col=0)
        numeric_df = df.select_dtypes(include=[np.number])
        gene_metadata = df.drop(columns=numeric_df.columns, errors='ignore')

        # Normalisierung
        normalized_df = pre.normalize_log1p(numeric_df)
    
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
        outliers = pre.detect_outliers_iqr(numeric_phage_df)
        cleaned_df = numeric_phage_df.mask(outliers)
        n_outliers = outliers.sum().sum()
    
        # Bereinigte Daten speichern
        cleaned_final_df = pd.concat([phage_metadata, cleaned_df], axis=1)
        cleaned_final_df.to_csv(output_cleaned_file)
    
        # Einteilung der Gene
        time_cols = list(cleaned_final_df.columns[3:]) #betrachtet erste Zeile ab Spalte 4
        cleaned_final_df["Temporal_Class"] = cleaned_final_df.apply(lambda row: label.classify_temporal_expression(row, time_cols), axis=1) # ordnet jedem Gen ein Label zu
    
        # Speichere Dateien mit Labels 'early', 'middle', 'late' in neue Datei
        output_labeled_file = utils.output_dir / relative_path.with_name(relative_path.stem + "_labeled.csv")
        cleaned_final_df.to_csv(output_labeled_file)
    
        # Speichern der Daten im Exportformat mit nur GenID und Label
        # Erstelle DataFrame mit nur GenID und zugehörigem Label
        export_df = cleaned_final_df[["Temporal_Class"]].copy()
        export_df.index.name = "GeneID" 
        export_df.reset_index(inplace=True)
        output_export_file = utils.output_dir / relative_path.with_name(relative_path.stem + "_export.csv")
        export_df.to_csv(output_export_file, index=False)
    
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

        # Analyse Label-Distribution
        label_counts = label.count_labels(export_df) # Anzahl Labels pro Datei
        label.total_label_counts += label_counts # Anzahl Labels aller Dateien

        draw_piechart(label_counts, file)

        # Kuchendiagramm Ergebnis Dokumentation und Interpretation:
        # Verteilung der Labels zwischen den Datensätzen variiert stark.
        # Innerhalb eines Datensatzes ist die Verteilung oft unausgeglichen:
        # in der Regel dominieren 1 bis 2 Label deutlich, während das dritte Label nur selten vorkommt.

        print(" Datei verarbeitet:", input_file)
        print(" Nach Normalisierung:", output_normalized_file)
        print(f" Boxplot für {file} gespeichert als {output_file}")
        print(" Normalisierte Phagengene: ", phage_output_file)
        print(" Bereinigte Phagengene: ", output_cleaned_file)
        print(" Gelabelte Phagengene: ", output_labeled_file)
        print(" Ausreißer entfernt:", n_outliers)
        print(" Exportfile: ", output_export_file)
        print(" Label Verteilung: ", label_counts)
