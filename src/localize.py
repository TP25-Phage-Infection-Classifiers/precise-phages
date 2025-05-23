# Datei für die user stories zur Genomischen Lokalisierung der Klassen

from pathlib import Path
import pandas as pd
import numpy as np
import seaborn as sb
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from BCBio import GFF


import utils




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


def generate_genome_map():

    genome_map_dir = Path("output/genome_maps")
    genome_map_dir.mkdir(parents=True, exist_ok=True)

    for gff_rel, fasta_rel, count_rel in zip(utils.my_gff3_files,
                                            utils.my_genome_flies,
                                            utils.my_team_files ):
        #Absoluten Pfad für jede Eingabedatei zusammenbauen
        gff_path   = utils.input_dir / gff_rel     #Gen Annotation (GFF3)
        fasta_path = utils.input_dir / fasta_rel   #Gesamtgenom (Fasta)

        # Label Datei finden in 'output'
        label_path = (utils.output_dir / Path(count_rel).parent /
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
        output_dna_file = utils.output_dir / f"database_DNA/{gff_path.stem}_genes.fasta"
        output_protein_file = utils.output_dir / f"database_Protein/{gff_path.stem}_genes_protein.fasta"

        output_dna_file.parent.mkdir(parents=True, exist_ok=True)
        output_protein_file.parent.mkdir(parents=True, exist_ok=True)

        with open(output_dna_file, "w") as out_handle:
            SeqIO.write(dna_records, out_handle, "fasta")

        with open(output_protein_file, "w") as out_handle:
            SeqIO.write(protein_records, out_handle, "fasta")

        print(f"DNA-Sequenzen gespeichert: {output_dna_file}")
        print(f"Protein-Sequenzen gespeichert: {output_protein_file}")