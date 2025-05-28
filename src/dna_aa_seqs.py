from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from BCBio import GFF
import files
import genome_maps

# DNA + Protein extrahieren
def write_sequences(): 
    
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
        labels = genome_maps.read_labels(label_path)
        valid_genes = set(gene_id.lower() for gene_id in labels[labels["Temporal_Class"].isin(["early", "middle", "late"])]["GeneID"])
        genome = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))

        dna_records = []
        protein_records = []

        with open(gff_path) as gff_handle:
            for rec in GFF.parse(gff_handle, base_dict=genome):
                for feature in rec.features:
                    if feature.type == "gene":
                        gene_id = feature.id if feature.id else "unknown_gene"
                        # Nur Gene extrahieren, die in valid_genes sind
                        if gene_id.lower() in valid_genes:
                            seq = feature.extract(rec.seq)

                            # DNA speichern
                            dna_record = SeqRecord(seq, id=gene_id, description="")
                            dna_records.append(dna_record)

                            # Protein speichern
                            protein_seq = seq.translate(to_stop=True)
                            protein_record = SeqRecord(protein_seq, id=gene_id, description="")
                            protein_records.append(protein_record)

        # Ausgabepfade für DNA und Protein (pro Datensatz)
        output_dna_file = files.output_dir / f"database_DNA/{gff_path.stem}_genes.fasta"
        output_protein_file = files.output_dir / f"database_Protein/{gff_path.stem}_genes_protein.fasta"

        output_dna_file.parent.mkdir(parents=True, exist_ok=True)
        output_protein_file.parent.mkdir(parents=True, exist_ok=True)

        with open(output_dna_file, "w") as out_handle:
            SeqIO.write(dna_records, out_handle, "fasta")

        with open(output_protein_file, "w") as out_handle:
            SeqIO.write(protein_records, out_handle, "fasta")

        #print(f"DNA-Sequenzen gespeichert: {output_dna_file}")
        #print(f"Protein-Sequenzen gespeichert: {output_protein_file}")