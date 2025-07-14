"""
Split multi-FASTA files into single-sequence FASTA files
and write a ColabFold-batch command list.

Input-Dir :  output/cleaned_Protein
Output-Dir:  output/cleaned_Protein_split
Cmd-File  :  colabfold_batch_cmd.txt   (in Output-Dir)

Run with:  python split_fasta_and_prepare_batch.py
"""

from pathlib import Path
from Bio import SeqIO   

# --- Pfade 
BASE = Path(r"C:\Users\sanja\OneDrive\Documents\GitHub\precise-phages\output")
INPUT_DIR  = BASE / "cleaned_Protein"
SPLIT_DIR  = BASE / "cleaned_Protein_split"
CMD_FILE   = SPLIT_DIR / "colabfold_batch_cmd.txt"
# -----------------------------------------------------------------------------

SPLIT_DIR.mkdir(parents=True, exist_ok=True)
cmd_parts = []

for fasta in INPUT_DIR.glob("*.fasta"):
    for rec in SeqIO.parse(fasta, "fasta"):
        # Dateiname = FASTA-Header
        out_path = SPLIT_DIR / f"{rec.id}.fasta"
        SeqIO.write(rec, out_path, "fasta")
        cmd_parts.append(str(out_path))
print(f"{len(cmd_parts)} Sequenzen geschrieben nach {SPLIT_DIR}")

