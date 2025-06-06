from Bio import SeqIO
from pathlib import Path

# Ordner definieren
fasta_dir = Path("output/database_DNA")
gff_dir = Path("output/database_GFF")
gff_dir.mkdir(parents=True, exist_ok=True)

#Durch alle FASTA Dateien gehen
for fasta_path in fasta_dir.glob("*.fasta"):
    records = list(SeqIO.parse(fasta_path, "fasta"))
    lines = ["##gff-version 3"]
    start_pos = 100

    for record in records:
        gene_id = record.id
        end_pos = start_pos + len(record.seq) - 1
        line = f"seq1\tRefSeq\tgene\t{start_pos}\t{end_pos}\t.\t+\t.\tID={gene_id}"
        lines.append(line)
        start_pos = end_pos + 100 # Abstand zwischen Genen

    gff_name = gff_dir / (fasta_path.stem + ".gff3")
    with open(gff_name, "w") as f:
        f.write("\n".join(lines))

    print(f"GFF Datei erzeugt: {gff_name}") #Ausgabe für jede Datei