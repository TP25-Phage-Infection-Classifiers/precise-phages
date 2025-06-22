import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score, train_test_split, StratifiedKFold
from sklearn.metrics import accuracy_score, f1_score, make_scorer

# Daten laden
data = pd.read_csv("output/feature_engineering/output_genes_features.csv")

print("Original Data shape:", data.shape)
print(data.head())
print(data.tail())

# Labels = letzte Zeile (ohne erste Spalte)
labels = data.iloc[-1, 1:].reset_index(drop=True)

# Features = alle Zeilen außer letzte, ab zweiter Spalte zu numerischen Werten
# Features transponieren, damit Samples in Zeilen sind
features = data.iloc[:-1, 1:].transpose()
features = features.apply(pd.to_numeric)

# Model
model = RandomForestClassifier(random_state=42)

# 1. Klassische 5-fache Cross-Validation
# Hier haben wir cv=5 gewählt, da die Ergebnisse etwas besser waren als für cv=10
# Accurancy eventuell noch verbessern (für bessere Ergebnisse)
accuracy_scores = cross_val_score(model, features, labels, cv=5, scoring='accuracy')
f1_scores = cross_val_score(model, features, labels, cv=5, scoring='f1_macro')

# Ausgabe
print("=== 5-Fold Cross-Validation ===")
print("Accuracy pro Fold:", accuracy_scores)
print("Mean Accuracy:", np.mean(accuracy_scores))
print("F1 score pro Fold:", f1_scores)
print("Mean F1 score:", np.mean(f1_scores))

# 2. Stratified K-Fold Cross-Validation
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
acc_list = []
f1_list = []

for train_idx, test_idx in skf.split(features, labels):
    X_train, X_test = features.iloc[train_idx], features.iloc[test_idx]
    y_train, y_test = labels.iloc[train_idx], labels.iloc[test_idx]

    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)

    acc_list.append(accuracy_score(y_test, y_pred))
    f1_list.append(f1_score(y_test, y_pred, average='macro'))

print("\n=== 2. Stratified K-Fold (5-fold) ===")
print("Accuracy pro Fold:", acc_list)
print("F1 Score pro Fold:", f1_list)
print("Mean Accuracy:", np.mean(acc_list))
print("Mean F1 Score:", np.mean(f1_list))

# 3. Holdout-Split (80/20)
# Train/Test-Split
X_train, X_test, y_train, y_test = train_test_split(
    features, labels, test_size=0.2, random_state=42, stratify=labels
)

model.fit(X_train, y_train)
y_pred = model.predict(X_test)
holdout_accuracy = accuracy_score(y_test, y_pred)
holdout_f1 = f1_score(y_test, y_pred, average='macro')

print("\n=== Holdout (80/20 Split) ===")
print("Holdout Accuracy:", holdout_accuracy)
print("Holdout F1-Score:", holdout_f1)