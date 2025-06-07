"""

Extrahiert Features  aus FASTA Einzelndateien oder Batch.

Features
---------
* Länge des Gens
* GC-Gehalt (%)
* Motiv-Vorkommen (TATA-Box, Sigma-70-Motife)
* K-Mer-Frequenzen (variable k-Liste)
* Codon Usage Bias

"""

import importlib.util
import re
from collections import Counter
from itertools import product
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from Bio.Data import CodonTable
from collections import defaultdict, Counter

# #Versuche normalen Import
# try:
#     from src.files import input_dir, my_genome_files

# #Wenn das Paket nicht erkannt wird:
# except ImportError:

# #src/files.py manuell laden
#     spec = importlib.util.spec_from_file_location("files", Path(__file__).parent / "files.py")
#     files = importlib.util.module_from_spec(spec)
#     spec.loader.exec_module(files)  
#     #Zugriff auf die Variablen
#     input_dir = files.input_dir
#     my_genome_files = files.my_genome_files

# -----------------------------------------------------------------------------

#Reguläre Ausdrücke zur Erkennung bakterielle Promotor Motive
#Sigma 70 abhängige Transkriptionsstartstellen

MOTIFS = {
    "TATA": re.compile(r"TATA[AT]A[AT]"), #TATA-Box (mögliche Matches:TATAAA, TATAAT, TATATA, TATATT)
    "sigma35": re.compile(r"TTGACA"),     # -35 Region bakterielle Promotoren
    "sigma10": re.compile(r"TATAAT"),     # Pribnow Box oder -10 Region
}


#GC Verhätlnis
def compute_gc(seq: str) -> float:
    
    seq = seq.upper()

    #Anzahl der G und C Basen
    gc = seq.count("G") + seq.count("C")

    #GC Gehalt in Prozent berechnen, auf 2 Nachkommastellen gerundet
    return round(100 * gc / len(seq), 2) if seq else 0.0


#Erzeugt alle möglichen Kombinationen von A, C, G, T der Länge k
# Beispiel: k=2 →'AA', 'AC', 'AG' 
def all_kmers(k: int):
    return ["".join(p) for p in product("ACGT", repeat=k)]


def kmer_freqs(seq: str, k: int):
    
    seq = seq.upper()

    #Gesamtzahl möglicher K-Mere
    total = max(len(seq) - k + 1, 1)

    #Zähle alle K-Mere in der Sequenz
    counts = Counter(seq[i : i + k] for i in range(len(seq) - k + 1))

    #Normalisiere K-Mer Häufigkeiten auf (0,1), auf 4 Nachkommastellen runden
    return {
        f"k{k}_{mer}": round(counts.get(mer, 0) / total, 4) 
        for mer in all_kmers(k)
        }


def motif_counts(seq: str):
    
    seq = seq.upper()

    #Durchlaufe alle Motive und zähle, wie oft jedes Motiv vorkommt
    #Als Ergebnis hat man Dictionary mit Namen und Vorkommenshäufigkeit
    return {
        f"motif_{name}": len(p.findall(seq))
          for name, p in MOTIFS.items()
          }

def codon_usage_bias(seq: str):
    
    seq = seq.upper()
    
    # Codons in gegebener Reihenfolge aus Sequenz extrahieren
    codons = [seq[i:i+3] for i in range(0, len(seq) - 2, 3) if len(seq[i:i+3]) == 3]
    table = CodonTable.unambiguous_dna_by_id[1]
    
    # Codon -> Aminosäure Mapping
    codon_to_aa = {codon: aa for codon, aa in table.forward_table.items()}
    
    # Falls vorhanden Stop-Codons entfernen
    for stop in table.stop_codons:
        codon_to_aa.pop(stop, None)
    
    # Aminosäure -> Liste aller zugehörigen Codons    
    aa_to_codons = defaultdict(list)
    for codon, aa in codon_to_aa.items():
        aa_to_codons[aa].append(codon)
    
    # Codonhäufigkeit wird gezählt    
    codon_counts = Counter(codons)  
    
    rscu = {}
    for aa, codons_for_aa in aa_to_codons.items():
        total = sum(codon_counts[c] for c in codons_for_aa)
        n = len(codons_for_aa)
        for c in codons_for_aa:
            val = codon_counts[c] / (total / n) if total > 0 else 0
            rscu[f"rscu_{c}"] = round(val, 4)  
    
    return rscu

# Gen ID (Header) extrahieren, um ID in .gff3 zu finden 
def get_gene_ids_from_fasta(fasta_path):
    return [record.id for record in SeqIO.parse(fasta_path, "fasta")]

#Liest Start, Ende und Strangrichtung für Gene aus einer GFF3-Datei
# Diese Funktion wird für jede FASTA-Datei aufgerufen, wenn die passende GFF3 vorhanden ist
# Die Gen-IDs aus dem FASTA werden mit denen aus der GFF-Datei abgeg
def extract_positions_for_genes(gff3_path, gene_ids):
    positions = {}   # speichert: Gen-ID -> (start, end, strand, seqid)
    gene_id_set = set(gene_ids)

    with open(gff3_path, "r") as f:
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue
            cols = line.strip().split("\t")
            if len(cols) != 9:
                continue  # ungültige Zeile überspringen

            seqid, source, feature, start, end, score, strand, phase, attributes = cols

            # Nur Gene berücksichtigen
            if feature.lower() != "gene":
                continue

            attr_dict = {}
            for attr in attributes.split(";"):
                if "=" in attr:
                    key, value = attr.split("=", 1)
                    attr_dict[key] = value

            gene_id = attr_dict.get("ID")
            if gene_id and gene_id in gene_id_set:
                positions[gene_id] = (int(start), int(end), strand, seqid)

    return positions       

# Labels aus transponierten export-file laden    
def load_labels(label_file):
    df_labels = pd.read_csv(label_file, header=None)
    # Erste Zeile (ab der 2. Spalte) sind Gen-IDs
    gene_ids = df_labels.iloc[0, 1:].tolist()
    # Zweite Zeile (ab der 2. Spalte) sind die Labels
    labels = df_labels.iloc[1, 1:].tolist()
    # Erstelle ein Dictionary: {Gen-ID: Label}
    return dict(zip(gene_ids, labels))
    
# Labels als weitere Zeile in output .csv-Datei laden
def add_label_row(df, label_dict):
    labels = []
    for gene in df.columns:
        base_gene = gene.split("_")[0]  # "gene-T4p161_0" -> "gene-T4p161"
        label = label_dict.get(base_gene, "unknown")
        labels.append(label)
    df.loc["Temporal_Class"] = labels
    return df
    
def extract_features(records, ks, positions=None): # positions optional übergeben
    #Initialisiere leeres Dictionary für alle Sequenzen
    data = {}

    #Durchlaufe alle FASTA Records (Liste von Bio.SeqRecord-Objekten aus FASTA)
    for i, r in enumerate(records):

        # Länge und GC Gehalt
        row = {
            "length": len(r.seq),
            "GC_content": compute_gc(str(r.seq)),
        }

        gene_id = r.id
        # Positionen (start, end, strand) als Feature hinzufügen
        if positions and gene_id in positions:
            start, end, strand, _ = positions[gene_id]
            row.update({
                "start" : start,
                "end" : end,
                "strand" : 1 if strand == "+" else -1 # Plus-Strang (5' -> 3'), Minus-Strang (3' <- 5')
            })

        #Zähle Promotor Motive
        row.update(motif_counts(str(r.seq)))

        #Berechne K-Mer Frequenzen für alle k-Werte
        for k in ks:
            row.update(kmer_freqs(str(r.seq), k))

        #Berechne Codon Usage Bias
        row.update(codon_usage_bias(str(r.seq)))
        
        #Weise die Feature Zeile der  Sequenz-ID zu
        uid = f"{r.id}_{i}"  # verhindert Überschreiben
        data[uid] = row

  #Erstelle ein DataFrame: Zeilen = Gene, Spalten = Features
    return pd.DataFrame.from_dict(data, orient="index")


def extract_all_features():
    input_dir = Path("output/database_DNA")
    gff_dir = Path("output/database_GFF")  
    output_dir = Path("output/feature_engineering")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Label_file laden
    label_file = Path("output/zusammengefuegt_transponiert.csv")
    label_dict = load_labels(label_file)

    ks = [3, 4]
    transpose = True

    for fasta_path in input_dir.glob("*.fasta"):
        out_path = output_dir / f"{fasta_path.stem}_features.csv"

        records = list(SeqIO.parse(fasta_path, "fasta"))
        if not records:
            continue

        # GFF3-Datei mit gleichem Namen wie FASTA laden
        gff3_path = gff_dir / f"{fasta_path.stem}.gff3"
        if not gff3_path.exists():
            print(f"Keine GFF3-Datei gefunden für {fasta_path.name}")
            continue

        gene_ids = [r.id for r in records]
        #Extrahiere Positionen für bekannte GEN IDs
        positions = extract_positions_for_genes(gff3_path, gene_ids)

        #Übergabe Positionen am Feature Funktion
        df = extract_features(records, ks, positions)
        if transpose:
            df = df.T

        # Labelzeile hinzufügen
        df = add_label_row(df, label_dict)
        df.to_csv(out_path, index=True, index_label="GeneID")