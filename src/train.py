"""
train.py

This script trains the machine learning model using the dataset and saves the trained model to the models directory.
"""

import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
import joblib
import os
from sklearn.metrics import classification_report, accuracy_score
from xgboost import XGBClassifier
from sklearn.preprocessing import LabelEncoder


# Load data
data = pd.read_csv('output/feature_matrix_with_structure_reduced.csv')

# Features und Zielspalte
X = data.drop(columns=['GeneID', 'Temporal_Class'])
y = data['Temporal_Class']


# Split data into training and validation sets
X_train, X_val, y_train, y_val = train_test_split(X, y, test_size=0.2, random_state=42)

# Train model
#model = RandomForestClassifier(n_estimators=100, random_state=42)

model = XGBClassifier(use_label_encoder=False, eval_metric='mlogloss', random_state=42)

le = LabelEncoder()
y_train_encoded = le.fit_transform(y_train)
y_val_encoded = le.transform(y_val)

model.fit(X_train, y_train_encoded)

# Save trained model
os.makedirs('output/models', exist_ok=True)
#joblib.dump(model, 'output/models/random_forest_model.pkl')
#print("Model trained and saved successfully: output/models/random_forest_model.pkl.")
joblib.dump(model, 'output/models/xgboost_model.pkl')
print("Model trained and saved successfully: output/models/xgboost_model.pkl")

y_pred = model.predict(X_val)
y_pred_labels = le.inverse_transform(y_pred)
print("\n Modellbewertung auf dem Validierungsset:")
print(f"Accuracy: {accuracy_score(y_val, y_pred):.4f}")
print(classification_report(y_val_encoded, y_pred, target_names=le.classes_))
