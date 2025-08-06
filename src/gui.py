import streamlit as st
from Bio import SeqIO
from BCBio import GFF
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd
import joblib
import tempfile
import os
import feature_engineering as fe

# Modell laden
model = joblib.load("output/models/xgboost_model.pkl")
label_encoder = joblib.load("output/models/label_encoder.joblib")

# Erwartete Feature-Reihenfolge laden
expected_columns = pd.read_csv("output/feature_matrix_with_structure_reduced.csv")
expected_columns = expected_columns.drop(columns=["GeneID", "Temporal_Class"]).columns

st.set_page_config(page_title="Precise-Phages", layout="centered")
st.title("🧬 Precise-Phages: Klassifikation mit FASTA + GFF3")

st.markdown("Lade eine FASTA- und die zugehörige GFF3-Datei hoch, um die temporale Expression deiner Gene vorherzusagen.")

fasta_file = st.file_uploader("📄 FASTA-Datei hochladen", type=["fasta", "fa"])
gff_file = st.file_uploader("📄 GFF3-Datei hochladen", type=["gff3"])

if fasta_file and gff_file:
    if st.button("🔍 Vorhersagen starten"):
        with st.spinner("Verarbeite Sequenzen & GFF3…"):
            # Temporäre Dateien erzeugen (Streamlit braucht Pfade für GFF Parser)
            with tempfile.TemporaryDirectory() as tmpdir:
                fasta_path = os.path.join(tmpdir, "input.fasta")
                gff_path = os.path.join(tmpdir, "input.gff3")

                with open(fasta_path, "wb") as f:
                    f.write(fasta_file.read())
                with open(gff_path, "wb") as f:
                    f.write(gff_file.read())

                # FASTA einlesen
                fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))

                # GFF3 einlesen -> extrahierte Gene als SeqRecord-Liste
                gene_records = []
                with open(gff_path) as gff_handle:
                    for rec in GFF.parse(gff_handle, base_dict=fasta_dict):
                        for feature in rec.features:
                            if feature.type != "gene":
                                continue
                            gene_seq = feature.extract(rec.seq)
                            gene_id = feature.qualifiers.get("ID", ["unknown_gene"])[0]
                            gene_record = SeqRecord(gene_seq, id=gene_id, description="")
                            gene_records.append(gene_record)

                if not gene_records:
                    st.error("Keine Gene im GFF3 gefunden.")
                else:
                    # Feature Engineering
                    features = fe.extract_features(gene_records, ks=[3, 4], positions=None)

                    # Feature-Alignment (wie beim Training)
                    features = features.reindex(columns=expected_columns, fill_value=0)

                    # Prediction
                    prediction = model.predict(features)
                    predicted_labels = label_encoder.inverse_transform(prediction)

                    # Ergebnis anzeigen
                    results = pd.DataFrame({
                        "GeneID": [r.id for r in gene_records],
                        "Temporal_Class": predicted_labels
                    })

                    st.success("Vorhersage abgeschlossen!")
                    st.dataframe(results)

                    # Download
                    csv = results.to_csv(index=False).encode("utf-8")
                    st.download_button("📥 Ergebnisse als CSV herunterladen", data=csv, file_name="vorhersage.csv", mime="text/csv")
