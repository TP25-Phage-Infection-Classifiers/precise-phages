from pathlib import Path

input_dir = Path("output/database_Protein")


def clean_fasta_file(input_path: Path, output_path: Path):
    with input_path.open('r') as infile, output_path.open('w') as outfile:
        sequence = ''
        for line in infile:
            line = line.strip()
            if line.startswith('>'):
                if sequence:
                    for i in range(0, len(sequence), 60):
                        outfile.write(sequence[i:i+60] + '\n')
                    sequence = ''
                outfile.write(line + '\n')
            else:
                sequence += line.replace(" ", "").replace("\t", "")
        if sequence:
            for i in range(0, len(sequence), 60):
                outfile.write(sequence[i:i+60] + '\n')


def clean_all_fastas(input_dir: Path, output_dir: Path):
    output_dir.mkdir(parents=True, exist_ok=True)
    for fasta_file in input_dir.glob("*.fasta"):
        clean_fasta_file(fasta_file, output_dir / fasta_file.name)
    print(f"Bereinigt: {len(list(output_dir.glob('*.fasta')))} Dateien gespeichert in {output_dir}")


if __name__ == "__main__":
    output_dir = Path("output/cleaned_Protein")  
    clean_all_fastas(input_dir, output_dir)