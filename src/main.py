import pandas as pd
import numpy as np
import glob
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path
from PIL import Image
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord 
# pip3 install bcbio-gff
from BCBio import GFF
import prepare_data as pre 
import label  
import plots
import files
import genome_maps
import dna_aa_seqs
import feature_engineering as fe


def generate_output():
    for file in files.my_team_files:
        input_file = files.input_dir / file
        relative_path = input_file.relative_to(files.input_dir)
   
        output_cleaned_file = files.output_dir / relative_path.with_name(relative_path.stem + "_cleaned.csv")
        output_normalized_file = files.output_dir / relative_path.with_name(relative_path.stem + "_normalized.csv")
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
        output_labeled_file = files.output_dir / relative_path.with_name(relative_path.stem + "_labeled.csv")
        cleaned_final_df.to_csv(output_labeled_file)
    
        # Speichern der Daten im Exportformat mit nur GenID und Label
        # Erstelle DataFrame mit nur GenID und zugehörigem Label
        export_df = cleaned_final_df[["Temporal_Class"]].copy()
        export_df.index.name = "GeneID" 
        export_df.reset_index(inplace=True)
        output_export_file = files.output_dir / relative_path.with_name(relative_path.stem + "_export.csv")
        export_df.to_csv(output_export_file, index=False)
    
        # Extrahiere nur die numerischen Spalten (Genexpressionswerte)
        raw_data = df.select_dtypes(include=[float, int])
        normalized_data = normalized_final_df.select_dtypes(include=[float, int])

        # Boxplot erstellen
        plots.create_boxplots(file, raw_data, normalized_data)

        # Analyse Label-Distribution
        label_counts = label.count_labels(export_df) # Anzahl Labels pro Datei
        label.total_label_counts += label_counts # Anzahl Labels aller Dateien

        # Piechart erstellen
        plots.draw_piechart(label_counts, file)

        # Kuchendiagramm Ergebnis Dokumentation und Interpretation:
        # Verteilung der Labels zwischen den Datensätzen variiert stark.
        # Innerhalb eines Datensatzes ist die Verteilung oft unausgeglichen:
        # in der Regel dominieren 1 bis 2 Label deutlich, während das dritte Label nur selten vorkommt.

        print("Datei verarbeitet:", input_file)
        #print(" Nach Normalisierung:", output_normalized_file)
        #print(" Normalisierte Phagengene: ", phage_output_file)
        #print(" Bereinigte Phagengene: ", output_cleaned_file)
        #print(" Gelabelte Phagengene: ", output_labeled_file)
        print("Ausreißer entfernt:", n_outliers)
        print("Exportfile: ", output_export_file)
        print("Label Verteilung:")
        print(label_counts)
        print()

generate_output()
print("Gesamte Label Verteilung: ")
print(label.total_label_counts) # in zwei Zeilen geprintet damit Formatierung in Terminal schöner
plots.draw_piechart(label.total_label_counts, "all_files") # pie chart für alle Daten

# Die Verteilung der Labels 'early', 'middle' und 'late' über alle Datensets ist relativ ausgeglichen. 
# Nur ca. 3% der Gene waren durchgehen so niedrig exprimiert, dass sie als undefined gelabled wurden. 

images = [Image.open(f).convert("RGB") for f in files.png_files]
images[0].save("output.pdf", save_all=True, append_images=images[1:])

# Alle Dateien, die auf '_export.csv' enden
dateien = glob.glob('output/**/*_export.csv', recursive=True)

# Zusamenfügen aller _export.csv Dateien
df_liste = [pd.read_csv(datei) for datei in dateien]
zusammen = pd.concat(df_liste, ignore_index=True)
zusammen = zusammen[zusammen['Temporal_Class'] != 'undefined']
zusammen.to_csv('output/zusammengefuegt.csv', index=False)

# Datei einlesen
df = pd.read_csv('output/zusammengefuegt.csv')

# Transponieren, sodass GeneIDs Spalten werden und Temporal_Class die Werte sind
df_transponiert = df.set_index('GeneID').T
df_transponiert.to_csv('output/zusammengefuegt_transponiert.csv')

genome_maps.generate_genome_map()
genome_maps.generate_temporal_class_position_distributions()
dna_aa_seqs.write_sequences()


# Eingabe-Ordner mit FASTA-Dateien
fasta_dir = Path("output/database_DNA")

# Ziel-Ordner für GFF-Dateien
gff_dir = Path("output/database_GFF")
gff_dir.mkdir(parents=True, exist_ok=True)  # Ordner anlegen, falls nicht vorhanden

# Schleife durch alle FASTA-Dateien im Projekt
for fasta_path in fasta_dir.glob("*.fasta"):
    # Zielpfad der zugehörigen GFF-Datei
    gff_path = gff_dir / (fasta_path.stem + ".gff3")

    # Falls GFF-Datei schon existiert → nichts tun
    if gff_path.exists():
        continue

    # FASTA-Sequenzen einlesen
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if not records:
        continue  # Datei enthält keine Gene

    # GFF3-Headerzeile
    lines = ["##gff-version 3"]

    # Startposition für das erste Gen im GFF
    start_pos = 1

    # Für jeden Gen-Record eine GFF-Zeile erzeugen
    for record in records:
        gene_id = record.id  # ID aus FASTA-Header
        end_pos = start_pos + len(record.seq) - 1  # Ende = Start + Länge - 1

        # GFF-Zeile zusammenbauen: fester Strang (+), keine Phase, keine Scores
        line = f"seq1\tRefSeq\tgene\t{start_pos}\t{end_pos}\t.\t+\t.\tID={gene_id}"
        lines.append(line)

        # Startposition für das nächste Gen
        start_pos = end_pos + 100

    # GFF-Datei speichern
    with open(gff_path, "w") as f:
        f.write("\n".join(lines))

    print(f"GFF-Datei automatisch erzeugt: {gff_path.name}")


fe.extract_all_features()
fe.merge_csvs("output/feature_engineering")