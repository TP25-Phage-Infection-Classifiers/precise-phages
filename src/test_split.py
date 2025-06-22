import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score, train_test_split, StratifiedKFold
from sklearn.metrics import accuracy_score, f1_score

import pprint


# Vergleich von verschiedenen Split Strategien um Daten in Trainings- und Testdaten einzuteilen
# Strategien: 
# - Klassische 5-fache Cross-Validation
# - Stratified K-Fold Cross-Validation
# - Holdout-Split (80/20)
# Metriken zum Vergleich der Strategien:
# - Mean Accuracy
# - F1 Score

# Daten laden
data = pd.read_csv("output/feature_engineering_merged.csv")

# Liste aller Phagennamen
phagenames = ["PHIKZ", "CPT_phageK", "T4", "phage515_", "DMS3", "VPVV882", "phiYY"]

# Phagen Gene-IDS einlesen
# gene_IDs = Liste von allen GeneIDs
gene_IDs = data.columns.tolist()[1:]

# Dictionary mit allen Phagennamen und leeren Listen (hier kommen die GeneIDS dann rein)
phage_to_genes = {phage: [] for phage in phagenames}

# Dictionary befüllen
for gene in gene_IDs:
    for phage in phagenames:
        if gene.startswith(f"gene-{phage}"):
            phage_to_genes[phage].append(gene)
            break  

# pprint.pprint(phage_to_genes)

# phage_to_genes = Pro Phage eine Liste, die alle Gene der Phage speichert
# --> Für saubere Trennung von Test- und Trainingsdaten

# Labels = letzte Zeile (ohne erste Spalte), Index setzen auf Gen-Namen
labels = data.iloc[-1, 1:]
labels.index = data.columns[1:]

# Features = alle Zeilen außer letzte, ab zweiter Spalte zu numerischen Werten
# Features transponieren, damit Samples in Zeilen sind
features = data.iloc[:-1, 1:].transpose()
features = features.apply(pd.to_numeric)

# Model
model = RandomForestClassifier(random_state=42)

# 1. Klassische 5-fache Cross-Validation
accuracy_scores = cross_val_score(model, features, labels, cv=10, scoring='accuracy')
f1_scores = cross_val_score(model, features, labels, cv=10, scoring='f1_macro')

print("=== 5-Fold Cross-Validation ===")
print("Accuracy pro Fold:", accuracy_scores)
print("F1 score pro Fold:", f1_scores)
print("Mean Accuracy:", np.mean(accuracy_scores))
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

print("\n=== Stratified K-Fold (5-fold) ===")
print("Accuracy pro Fold:", acc_list)
print("F1 Score pro Fold:", f1_list)
print("Mean Accuracy:", np.mean(acc_list))
print("Mean F1 Score:", np.mean(f1_list))

# 3. Holdout-Split (80/20)
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

# 4. Leave-One-Phage-Out Cross-Validation (jede Phage einmal als Testset)
print("\n=== Leave-One-Phage-Out Cross-Validation ===")

results = {}

for test_phage in phagenames:
    test_genes = phage_to_genes[test_phage]
    train_genes = [gene for gene in features.index if gene not in test_genes]
    
    X_train = features.loc[train_genes]
    y_train = labels.loc[train_genes]
    
    X_test = features.loc[test_genes]
    y_test = labels.loc[test_genes]
    
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)
    
    acc = accuracy_score(y_test, y_pred)
    f1 = f1_score(y_test, y_pred, average='macro')
    
    results[test_phage] = {"accuracy": acc, "f1_score": f1}

for phage, metrics in results.items():
    print(f"Test-Phage: {phage}")
    print(f"  Accuracy: {metrics['accuracy']:.4f}")
    print(f"  F1 Score: {metrics['f1_score']:.4f}")
    print()