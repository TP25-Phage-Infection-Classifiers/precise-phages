"""
split_fasta_per_protein.py
---------------------------------
• liest alle *.fasta im Ordner output/cleaned_Protein
• legt zu jeder Datei einen Unterordner *_split an
• schreibt dort pro Sequenz eine neue FASTA-Datei
"""

from pathlib import Path
from Bio import SeqIO   

BASE = Path(r"C:\Users\sanja\OneDrive\Documents\GitHub\precise-phages\output")
IN_DIR = BASE / "cleaned_Protein"

for fasta in IN_DIR.glob("*.fasta"):
    split_dir = IN_DIR / (fasta.stem + "_split")
    split_dir.mkdir(exist_ok=True)
    count = 0
    for rec in SeqIO.parse(fasta, "fasta"):
        out_fasta = split_dir / f"{rec.id}.fasta"
        SeqIO.write(rec, out_fasta, "fasta")
        count += 1
    print(f" {count:>4} Proteine → {split_dir}")
