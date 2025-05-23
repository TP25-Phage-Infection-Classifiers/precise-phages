from pathlib import Path
import os
import pandas as pd
import numpy as np
import seaborn as sb
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from PIL import Image
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord 
# pip3 install bcbio-gff
from BCBio import GFF



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
my_gff3_files=[
    "ceyssens_2014/Pseudomonas_phage_phiKZ.gff3",
    "finstrlova_2022/Staphylococcus_phage_K.gff3",
    "guegler_2021/Enterobacteria_phage_T4.gff3",
    "kuptsov_2022/Staphylococcus_phage_vB_SauM-515A1.gff3",
    "meaden_2021/Pseudomonas_phage_DMS3.gff3",
    "sprenger_2024/Vibrio_phage_VP882.gff3",
    "zhong_2020/Pseudomonas_phage_phiYY_complete.gff3"
]
my_genome_flies=[
    "ceyssens_2014/Pseudomonas_phage_phiKZ.fasta",
    "finstrlova_2022/Staphylococcus_phage_K.fasta",
    "guegler_2021/Enterobacteria_phage_T4.fasta",
    "kuptsov_2022/Staphylococcus_phage_vB_SauM-515A1.fasta",
    "meaden_2021/Pseudomonas_phage_DMS3.fasta",
    "sprenger_2024/Vibrio_phage_VP882.fasta",
    "zhong_2020/Pseudomonas_phage_phiYY_complete.fasta",
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
# hierfür wird geschaut, in welchem drittel der timeslots der maximale Wert an counts liegt. 
# liegt dieser im ersten Drittel wird das Gen als 'early' gelabled, liegt er im mittleren Drittel als 'middle' und im hinteren als 'late'
# Sind die normalisierten Counts alle unter einem Wert von 2.5 wird das Gen als 'undefined gelabled und nicht für das ML verwendet
def classify_temporal_expression(row, time_cols):
    max_value = row[time_cols].max()
    max_timepoint_idx = row[time_cols].values.argmax() # Index des höchsten Werts der Zeile
    n_timepoints = len(time_cols) # Anzahl der Werte
    if max_value < 2.5: # durchgehend niedrig exprimierte Gene werden als undefined definiert 
        return 'undefined'
    elif max_timepoint_idx < n_timepoints / 3: # wenn Index im ersten Drittel liegt -> early
        return 'early'
    elif max_timepoint_idx < 2 * n_timepoints / 3: # wenn Index im zweiten Drittel liegt -> middle
        return 'middle'
    else:
        return 'late'
    
# Label-Distribution Analyse
label_order = ['undefined', 'early', 'middle','late'] # Reihenfolge, in der counts ausgegeben
total_label_counts = pd.Series(0, index=label_order) # Counts aller Dateien, wird später berechnet

# zählt, wie oft jedes Label in einer Datei vorkommt
def count_labels(df):
    counts = df['Temporal_Class'].value_counts()
    return counts.reindex(label_order, fill_value=0) 

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
boxplot_output_dir = output_dir / "boxplots"
pie_output_dir = output_dir / "pie chart"
boxplot_output_dir.mkdir(exist_ok=True) 
pie_output_dir.mkdir(exist_ok=True)


for file in my_team_files:
    input_file = input_dir / file
    relative_path = input_file.relative_to(input_dir)
   
    output_cleaned_file = output_dir / relative_path.with_name(relative_path.stem + "_cleaned.csv")
    output_normalized_file = output_dir / relative_path.with_name(relative_path.stem + "_normalized.csv")
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
    
    # Speichern der Daten im Exportformat mit nur GenID und Label
    # Erstelle DataFrame mit nur GenID und zugehörigem Label
    export_df = cleaned_final_df[["Temporal_Class"]].copy()
    export_df.index.name = "GeneID" 
    export_df.reset_index(inplace=True)
    output_export_file = output_dir / relative_path.with_name(relative_path.stem + "_export.csv")
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
    label_counts = count_labels(export_df) # Anzahl Labels pro Datei
    total_label_counts += label_counts # Anzahl Labels aller Dateien

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

print("Gesamte Label Verteilung: ")
print(total_label_counts) # in zwei Zeilen geprintet damit Formatierung in Terminal schöner
draw_piechart(total_label_counts, "all_files") # pie chart für alle Daten

# Die Verteilung der Labels 'early', 'middle' und 'late' über alle Datensets ist relativ ausgeglichen. 
# Nur ca. 3% der Gene waren durchgehen so niedrig exprimiert, dass sie als undefined gelabled wurden. 

png_files = ["output/pie chart/Ceyssens_directional_full_raw_counts_pie.png", 
             "output/pie chart/Finstrlova_SH1000_full_raw_counts_pie.png",
             "output/pie chart/Guegler_T4_plusToxIN_full_raw_counts_pie.png",
             "output/pie chart/Kuptsov_full_raw_counts_pie.png",
             "output/pie chart/Meaden_BIM_full_raw_counts_pie.png",
             "output/pie chart/Sprenger_VC_WT_VP882_WT_full_raw_counts_pie.png",
             "output/pie chart/Zhong_full_raw_counts_pie.png",
             "output/pie chart/all_files_pie.png"]

images = [Image.open(f).convert("RGB") for f in png_files]
images[0].save("output.pdf", save_all=True, append_images=images[1:])

## Genomkarte erstellen

def read_gff3(file_path):
    
    ##Liest eine GFF3 Datei ein und extrahiert die Gene mit deren Start, End und GeneID.
    
    gff3 = pd.read_csv(file_path, sep='\t', comment='#', header=None,
                       names=['seqid','source','type','start','end','score','strand','phase','attributes'])
    genes = gff3[gff3['type'] == 'gene'].copy()
    genes['GeneID'] = genes['attributes'].str.extract(r'ID=([^;]+)')
    return genes[['GeneID','start','end']]

def read_labels(file_path):
    
    #Liest die Labels ein und stellt sicher, dass 'GeneID' und 'Temporal_Class' vorhanden sind.

    labels = pd.read_csv(file_path)
    # Prüfen ob die Spalten vorhanden sind
    required_cols = ['GeneID', 'Temporal_Class']
    for col in required_cols:
        if col not in labels.columns:
            raise ValueError(f"Spalte '{col}' nicht in Label-Datei gefunden. Vorhandene Spalten: {labels.columns.tolist()}")
    return labels

def read_fasta_length(fasta_path):
    
    ##Liest eine FASTA-Datei ein und gibt die Länge des Genoms zurück.
    record = next(SeqIO.parse(fasta_path, "fasta"))
    return len(record.seq)

def plot_genome_map(genes_df, genome_length, title, output_path):
    
    ##Zeichnet eine Genomkarte, farblich nach Temporal_Class, auf Grundlage der Positionen.
    colors = {'early':'#fc8d62','middle':'#8da0cb','late':'#e78ac3','undefined':'#cccccc'}
    fig, ax = plt.subplots(figsize=(12,2))

    # Falls Temporal_Class nicht definiert, setze 'undefined'
    if 'Temporal_Class' not in genes_df.columns:
        genes_df['Temporal_Class'] = 'undefined'
    else:
        genes_df['Temporal_Class'] = genes_df['Temporal_Class'].fillna('undefined')

    for _, row in genes_df.iterrows():
        color = colors.get(row['Temporal_Class'], '#999999')
        start = row['start']
        width = row['end'] - row['start']
        ax.add_patch(mpatches.Rectangle((start, 0), width, 1, color=color))

    ax.set_xlim(0, genome_length + 1000)
    ax.set_ylim(0,1)
    ax.set_yticks([])
    ax.set_xlabel('Genomposition (bp)')
    ax.set_title(title)
    legend = [mpatches.Patch(color=clr, label=lbl) for lbl, clr in colors.items()]
    ax.legend(handles=legend, loc='upper right')
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


# Genome Maps für alle Datensätze

genome_map_dir = Path("output/genome_maps")
genome_map_dir.mkdir(parents=True, exist_ok=True)

for gff_rel, fasta_rel, count_rel in zip(my_gff3_files,
                                        my_genome_flies,
                                         my_team_files ):
    #Absoluten Pfad für jede Eingabedatei zusammenbauen
    gff_path   = input_dir / gff_rel     #Gen Annotation (GFF3)
    fasta_path = input_dir / fasta_rel   #Gesamtgenom (Fasta)

    # Label Datei finden in 'output'
    label_path = (output_dir / Path(count_rel).parent /
                f"{Path(count_rel).stem}_export.csv")
    # Falls Label Datei nicht existiert, nächsten Datensatz nehmen
    if not label_path.exists():
        print(f"Überspringe {count_rel}: Label-Datei fehlt")
        continue
   

    # Daten laden
    gff_genes = read_gff3(gff_path)   #Gene + Positio aus GFF3
    labels = read_labels(label_path)

    #Gene und Labels zusammenführen
    # Merge auf GeneID (Vorsicht: Duplizierte IDs oder fehlende könnten Probleme machen)
    merged = pd.merge(gff_genes, labels, on='GeneID', how='left')
   
    # Alle Gene ohne Temporal_Class als 'undefined' kennzeichnen (wird in Plot gefixt)
    merged["Temporal_Class"] = merged["Temporal_Class"].fillna("undefined")
    
    # Genomlänge aus Fasta bestimmen
    genome_length = read_fasta_length(fasta_path)
    
    # Pfad für die Ausgabegrafik
    out_img = genome_map_dir / f"{gff_path.stem}_genome_map.png"

    # Plot erstellen
    plot_genome_map(merged,
                    genome_length,
                    f"{gff_path.stem} Genomkarte",  # Titel = Dateiname + Genomkarte
                    out_img)

    print(f"Genomkarte gespeichert: {out_img}")



# DNA + Protein extrahieren
    genome = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))

    dna_records = []
    protein_records = []

    with open(gff_path) as gff_handle:
        for rec in GFF.parse(gff_handle, base_dict=genome):
            for feature in rec.features:
                if feature.type == "gene":
                    gene_id = feature.id if feature.id else "unknown_gene"
                    seq = feature.extract(rec.seq)

                    # DNA speichern
                    dna_record = SeqRecord(seq, id=gene_id, description="")
                    dna_records.append(dna_record)

                    # Protein speichern
                    protein_seq = seq.translate(to_stop=True)
                    protein_record = SeqRecord(protein_seq, id=gene_id, description="")
                    protein_records.append(protein_record)


    # Ausgabepfade für DNA und Protein (pro Datensatz)
    output_dna_file = output_dir / f"database_DNA/{gff_path.stem}_genes.fasta"
    output_protein_file = output_dir / f"database_Protein/{gff_path.stem}_genes_protein.fasta"

    output_dna_file.parent.mkdir(parents=True, exist_ok=True)
    output_protein_file.parent.mkdir(parents=True, exist_ok=True)

    with open(output_dna_file, "w") as out_handle:
        SeqIO.write(dna_records, out_handle, "fasta")

    with open(output_protein_file, "w") as out_handle:
        SeqIO.write(protein_records, out_handle, "fasta")

    print(f"DNA-Sequenzen gespeichert: {output_dna_file}")
    print(f"Protein-Sequenzen gespeichert: {output_protein_file}")