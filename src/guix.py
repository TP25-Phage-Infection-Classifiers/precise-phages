
import streamlit as st
import pandas as pd
import os
import shutil
from pathlib import Path
import subprocess

st.set_page_config(page_title="Precise Phages GUI", layout="centered")

st.title("🧬 Precise Phages: Analyse und Vorhersage")

# == 1. Datei-Upload ==
st.header(" Datei hochladen")
uploaded_file = st.file_uploader("Lade eine Genom- oder CSV-Datei hoch", type=["csv", "fasta", "fa", "txt"])

input_dir = Path("input")
input_dir.mkdir(exist_ok=True)

if uploaded_file:
    file_path = input_dir / uploaded_file.name
    with open(file_path, "wb") as f:
        f.write(uploaded_file.getbuffer())
    st.success(f"Datei erfolgreich hochgeladen: {uploaded_file.name}")

# == 2. Daten aufbereiten ==
st.header(" Datenaufbereitung")

if st.button("Daten vorbereiten"):
    try:
        from src.main import generate_output
        generate_output()
        st.success("👩‍🔬 Datenaufbereitung abgeschlossen.")
    except Exception as e:
        st.error(f"🧐 Fehler bei der Aufbereitung: {e}")

# == 3. Vorhersage starten ==
st.header(" Vorhersage mit XGBoost")

if st.button("Vorhersage ausführen"):
    try:
        subprocess.run(["python3", "src/predict.py"], check=True)
        st.success(" Vorhersage abgeschlossen.")
    except subprocess.CalledProcessError as e:
        st.error(f" Fehler bei der Vorhersage: {e}")

# == 4. Ergebnisse anzeigen ==
st.header(" Vorhersagen anzeigen")

pred_path = input_dir / "predictions.csv"
if pred_path.exists():
    df = pd.read_csv(pred_path)
    st.dataframe(df)
    st.download_button(" Vorhersagen herunterladen", df.to_csv(index=False), file_name="predictions.csv")
else:
    st.info("Noch keine Vorhersagen gefunden.")
