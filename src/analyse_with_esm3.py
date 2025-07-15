import os
from pathlib import Path
from huggingface_hub import login
from esm.models.esm3 import ESM3
from esm.sdk.api import ESM3InferenceClient, ESMProtein, GenerationConfig

# KONFIGURATION 
TOKEN = "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"  # Wird ggf. als Parameter ersetzt
INPUT_DIR = Path("outout/cleaned_Protein")
OUTPUT_DIR = Path("output/output_structures")
DEVICE = "cpu"
STEPS = 8
TEMP = 0.7

# LOGIN UND MODELL LADEN
login(token=TOKEN)
model: ESM3InferenceClient = ESM3.from_pretrained("esm3-open").to(DEVICE)

# PROTEINORDNER DURCHGEHEN 
for phage_folder in INPUT_DIR.iterdir():
    if not phage_folder.is_dir():
        continue

    out_folder = OUTPUT_DIR / phage_folder.name
    out_folder.mkdir(parents=True, exist_ok=True)

    for file in phage_folder.glob("*.txt"):
        with open(file, "r") as f:
            seq = f.read().strip()
        if not seq:
            print(f"Leere Datei: {file}")
            continue

        protein = ESMProtein(sequence=seq)
        try:
            protein = model.generate(protein, GenerationConfig(track="structure", num_steps=STEPS, temperature=TEMP))
            out_path = out_folder / (file.stem + ".pdb")
            protein.to_pdb(out_path)
            print(f" {file.name} -> {out_path.name}")
        except Exception as e:
            print(f"Fehler bei {file.name}: {e}")
