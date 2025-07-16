from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from BCBio import GFF
import files
import numpy as np

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

    # Start und Ende in Prozent umrechnen
    genes_df['start_pct'] = (genes_df['start'] / genome_length) * 100
    genes_df['end_pct'] = (genes_df['end'] / genome_length) * 100    

    for _, row in genes_df.iterrows():
        color = colors.get(row['Temporal_Class'], '#999999')
        start = row['start_pct']
        width = row['end_pct'] - row['start_pct']
        ax.add_patch(mpatches.Rectangle((start, 0), width, 1, color=color))

    ax.set_xlim(0, 100)
    ax.set_ylim(0,1)
    ax.set_yticks([])
    ax.set_xlabel('Genomposition (%)')
    ax.set_title(title)
    legend = [mpatches.Patch(color=clr, label=lbl) for lbl, clr in colors.items()]
    ax.legend(handles=legend, loc='upper right')
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

# Genome Maps für alle Datensätze
def generate_genome_map():

    genome_map_dir = Path("output/genome_maps")
    genome_map_dir.mkdir(parents=True, exist_ok=True)

    for gff_rel, fasta_rel, count_rel in zip(files.my_gff3_files,
                                            files.my_genome_files,
                                            files.my_team_files ):
        #Absoluten Pfad für jede Eingabedatei zusammenbauen
        gff_path   = files.input_dir / gff_rel     #Gen Annotation (GFF3)
        fasta_path = files.input_dir / fasta_rel   #Gesamtgenom (Fasta)

        # Label Datei finden in 'output'
        label_path = (files.output_dir / Path(count_rel).parent /
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

        # print(f"Genomkarte gespeichert: {out_img}")


def plot_temporal_class_position_distribution(genes_df, genome_length, output_path):

    # Umrechnung auf Prozentposition (z.B. Startposition)
    genes_df['start_pct'] = (genes_df['start'] / genome_length) * 100

    # Temporal Classes definieren und Farben
    classes = ['early', 'middle', 'late', 'undefined']
    colors = {'early':'#fc8d62', 'middle':'#8da0cb', 'late':'#e78ac3', 'undefined':'#cccccc'}

    plt.figure(figsize=(12,4))

    bins = np.linspace(0, 100, 50)  # 50 bins von 0 bis 100%

    for tc in classes:
        subset = genes_df[genes_df['Temporal_Class'] == tc]
        if subset.empty:
            continue
        plt.hist(subset['start_pct'], bins=bins, alpha=0.6, label=tc, color=colors[tc])

    plt.xlabel('Genomposition (%)')
    plt.ylabel('Anzahl Gene')
    plt.title('Verteilung der Gene entlang des Genoms nach Temporal Class')
    plt.legend(title='Temporal Class')
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
def generate_temporal_class_position_distributions():

    output_dir = Path("output/temporal_class_position_distributions")
    output_dir.mkdir(parents=True, exist_ok=True)

    for gff_rel, fasta_rel, count_rel in zip(files.my_gff3_files,
                                            files.my_genome_files,
                                            files.my_team_files):
        gff_path = files.input_dir / gff_rel
        fasta_path = files.input_dir / fasta_rel

        label_path = (files.output_dir / Path(count_rel).parent /
                      f"{Path(count_rel).stem}_export.csv")
        if not label_path.exists():
            print(f"Überspringe {count_rel}: Label-Datei fehlt")
            continue

        # Daten laden
        gff_genes = read_gff3(gff_path)
        labels = read_labels(label_path)

        merged = pd.merge(gff_genes, labels, on='GeneID', how='left')
        merged["Temporal_Class"] = merged["Temporal_Class"].fillna("undefined")

        genome_length = read_fasta_length(fasta_path)

        out_img = output_dir / f"{gff_path.stem}_temporal_class_position_distribution.png"

        # Plot-Funktion intern verwenden
        plot_temporal_class_position_distribution(merged, genome_length, out_img)
        print(f"Positionsverteilung gespeichert: {out_img}")