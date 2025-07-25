"""
predict.py

This script loads the trained model and makes predictions on new data.
"""

import pandas as pd
import joblib
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import feature_engineering as fe

# Modelle vorbereiten
model_paths = {
    "Random Forest": "output/models/random_forest.joblib",
    "SVM": "output/models/svm.joblib",
    "XGBoost": "output/models/xgboost.joblib"
}

# Daten vorbereiten
feature_data = pd.read_csv("output/feature_engineering_reduced.csv", index_col=0)
labels = feature_data.iloc[-1]
features = feature_data.iloc[:-1].transpose()
features = features.apply(pd.to_numeric)

# LabelEncoder laden oder neu erzeugen
label_encoder = None
try:
    label_encoder = joblib.load("output/models/label_encoder.joblib")
except:
    from sklearn.preprocessing import LabelEncoder
    label_encoder = LabelEncoder()
    label_encoder.fit(labels)
    joblib.dump(label_encoder, "output/models/label_encoder.joblib")

# Prediction Methode
def predict_temporal_class(gene: str, model_name: str):
    try:
        model = joblib.load(model_paths[model_name])
        record = SeqRecord(Seq(gene), id="geneX")
        records = [record]
        ks = [3, 4]
        features = fe.extract_features(records, ks, positions = None)
        prediction = model.predict(features)
        # Dekodiere falls XGBoost (numerisches Label)
        if "xgboost" in model_paths[model_name]:
            prediction = label_encoder.inverse_transform(prediction)
        return prediction[0]
    except KeyError:
        return "Gen nicht gefunden"
    except Exception as e:
        return f"Fehler: {e}"

# # Save predictions
# output = pd.DataFrame({'Id': new_data.index, 'Prediction': predictions})
# output.to_csv('../input/predictions.csv', index=False)
# print("Predictions saved successfully.")

gene1 = "ATGAAATCATATAAAGTAAATTTAGAACTTTTTGATAAAGCAGTTCATCGAGAATATAGAATCATTCAACGCTTTTTCGATATGGGAGAAGCCGAAGAATTTAAAACCCGCTTTAAAGATATTAGAGATAAAATTCAATCCGACACCGCAACTAAAGATGAACTACTAGAAGTTGCTGAAGTTATTAAGCGTAATATGAATTAA"

print(predict_temporal_class(gene1, "Random Forest"))