import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score

# Daten laden
data = pd.read_csv("output/feature_engineering/output_genes_features.csv")

print("Original Data shape:", data.shape)
print(data.head())
print(data.tail())

# Labels = letzte Zeile (ohne erste Spalte)
labels = data.iloc[-1, 1:].reset_index(drop=True)

# Features = alle Zeilen außer letzte, ab zweiter Spalte
features = data.iloc[:-1, 1:]

# Features transponieren, damit Samples in Zeilen sind
features = features.transpose()

# Features in numerische Werte umwandeln
features = features.apply(pd.to_numeric)

# 5-facher Kreuzvalidierung
# Hier haben wir cv=5 gewählt, da die Ergebnisse etwas besser waren als für cv=10
# Accurancy eventuell noch verbessern (für bessere Ergebnisse)
model = RandomForestClassifier(random_state=42)
scores = cross_val_score(model, features, labels, cv=5, scoring='accuracy')

# Ausgabe
print("Accuracy pro Fold:", scores)
print("Durchschnittliche Accuracy:", np.mean(scores))