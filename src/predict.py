"""
predict.py

This script loads the trained model and makes predictions on new data.
"""
import pandas as pd
import joblib
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import feature_engineering as fe

# Modell laden
model = joblib.load("output/models/xgboost_model.pkl")
label_encoder = joblib.load("output/models/label_encoder.joblib")

# Die Spalten (Features), die das Modell erwartet:
expected_columns = pd.read_csv("output/feature_matrix_with_structure_reduced.csv").drop(columns=["GeneID", "Temporal_Class"]).columns

def predict_temporal_class(gene: str):
    record = SeqRecord(Seq(gene), id="geneX")
    features = fe.extract_features([record], ks=[3,4], positions=None)

    # Setze dieselbe Feature-Reihenfolge wie beim Training
    features = features.reindex(columns=expected_columns, fill_value=0)

    prediction = model.predict(features)
    decoded = label_encoder.inverse_transform(prediction)
    return decoded[0]

# Beispiel-Gene
gene1 = "ATGAAATCATATAAAGTAAATTTAGAACTTTTTGATAAAGCAGTTCATCGAGAATATAGAATCATTCAACGCTTTTTCGATATGGGAGAAGCCGAAGAATTTAAAACCCGCTTTAAAGATATTAGAGATAAAATTCAATCCGACACCGCAACTAAAGATGAACTACTAGAAGTTGCTGAAGTTATTAAGCGTAATATGAATTAA"
print("Vorhersage:", predict_temporal_class(gene1))