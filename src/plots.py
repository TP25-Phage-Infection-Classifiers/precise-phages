import pandas as pd
from pathlib import Path
import seaborn as sb
import matplotlib.pyplot as plt
import files

boxplot_output_dir = files.output_dir / "boxplots"
pie_output_dir = files.output_dir / "pie chart"
boxplot_output_dir.mkdir(exist_ok=True) 
pie_output_dir.mkdir(exist_ok=True)


def draw_piechart(label_counts, file: str):
    if label_counts.empty or label_counts.sum() == 0:
        print(f"⚠️  Keine Daten für Kuchendiagramm vorhanden für Datei: {file}")
        return
    
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


def create_boxplots(file: str, df_raw: pd.DataFrame, df_normalized: pd.DataFrame):
    """Erstellt und speichert Boxplots für rohe und normalisierte Daten."""
    # Schmelze die Daten
    raw_data_melted = df_raw.melt(var_name='Gene', value_name='Expression')
    raw_data_melted['Condition'] = 'Raw'

    normalized_data_melted = df_normalized.melt(var_name='Gene', value_name='Expression')
    normalized_data_melted['Condition'] = 'Normalized'

    # Subplots erstellen
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    sb.boxplot(x='Condition', y='Expression', data=raw_data_melted, ax=axes[0], color='#B0E0E6')
    axes[0].set_title(f'{file} - Vor der Normalisierung')
    axes[0].set_ylabel('Genexpression')

    sb.boxplot(x='Condition', y='Expression', data=normalized_data_melted, ax=axes[1], color='#B0E0E6')
    axes[1].set_title(f'{file} - Nach der Normalisierung')
    axes[1].set_ylabel('Genexpression')

    plt.tight_layout()

    output_file = boxplot_output_dir / f"{Path(file).stem}_boxplot.png"
    plt.savefig(output_file)
    plt.close()

    return output_file
