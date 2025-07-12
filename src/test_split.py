import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from xgboost import XGBClassifier
from sklearn.model_selection import cross_val_score, train_test_split, StratifiedKFold
from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score
import matplotlib.pyplot as plt
from sklearn.preprocessing import LabelEncoder
import pprint


# Vergleich von verschiedenen Split Strategien um Daten in Trainings- und Testdaten einzuteilen
# Strategien: 
# - Klassische 5-fache Cross-Validation
# - Stratified K-Fold Cross-Validation
# - Holdout-Split (80/20)
# Metriken zum Vergleich der Strategien:
# - Mean Accuracy
# - F1 Score (macro)
# - Precision
# - Recall


# Daten laden
data = pd.read_csv("output/feature_matrix_with_structure.csv")
phagenames = ["PHIKZ", "CPT_phageK", "T4", "phage515_", "DMS3", "VPVV882", "phiYY"]
gene_IDs = data.columns.tolist()[1:] #Liste von allen GeneIDs
# Dictionary mit allen Phagennamen und leeren Listen (hier kommen die GeneIDS dann rein)
phage_to_genes = {phage: [] for phage in phagenames}

# Dictionary befüllen
for gene in gene_IDs:
    for phage in phagenames:
        if gene.startswith(f"gene-{phage}"):
            phage_to_genes[phage].append(gene)
            break  

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
models = {
    "Random Forest": RandomForestClassifier(random_state=42),
    "SVM": SVC(kernel='rbf', probability=True, random_state=42),
    "XGBoost": XGBClassifier(eval_metric='mlogloss', random_state=42)
}

def evaluate_model(name, model, features, labels):
    results = {}
    # Labels encoden für XGBoost (braucht numerische Labels)
    le = LabelEncoder()
    y_encoded = pd.Series(le.fit_transform(labels), index=labels.index)
    if "XGBoost" in name:
        y_to_use = y_encoded
    else:
        y_to_use = labels

    # 1. Klassische 5-fache Cross-Validation
    acc_scores = cross_val_score(model, features, y_to_use, cv=10, scoring='accuracy')
    f1_scores = cross_val_score(model, features, y_to_use, cv=10, scoring='f1_macro')
    precision_scores = cross_val_score(model, features, y_to_use, cv=10, scoring='precision_macro')
    recall_scores = cross_val_score(model, features, y_to_use, cv=10, scoring='recall_macro')
    results['5-Fold Mean Accuracy'] = np.mean(acc_scores)
    results['5-Fold Mean F1'] = np.mean(f1_scores)
    results['5-Fold Mean Precision'] = np.mean(precision_scores)
    results['5-Fold Mean Recall'] = np.mean(recall_scores)

    # 2. Stratified K-Fold Cross-Validation
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    acc_list = []
    f1_list = []
    precision_list = []
    recall_list = []

    for train_idx, test_idx in skf.split(features, y_to_use):
        X_train, X_test = features.iloc[train_idx], features.iloc[test_idx]
        y_train, y_test = y_to_use.iloc[train_idx], y_to_use.iloc[test_idx]
        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)
        acc_list.append(accuracy_score(y_test, y_pred))
        f1_list.append(f1_score(y_test, y_pred, average='macro'))
        precision_list.append(precision_score(y_test, y_pred, average='macro'))
        recall_list.append(recall_score(y_test, y_pred, average='macro'))
        
    results['Stratified K-Fold Mean Accuracy'] = np.mean(acc_list)
    results['Stratified K-Fold Mean F1'] = np.mean(f1_list)
    results['Stratified K-Fold Mean Precision'] = np.mean(precision_list)
    results['Stratified K-Fold Mean Recall'] = np.mean(recall_list)

    # 3. Holdout-Split (80/20)
    X_train, X_test, y_train, y_test = train_test_split(
        features, y_to_use, test_size=0.2, random_state=42, stratify=y_to_use
    )
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)
    results['Holdout Accuracy'] = accuracy_score(y_test, y_pred)
    results['Holdout F1'] = f1_score(y_test, y_pred, average='macro')
    results['Holdout Precision'] = precision_score(y_test, y_pred, average='macro')
    results['Holdout Recall'] = recall_score(y_test, y_pred, average='macro')

    # 4. Leave-One-Phage-Out Cross-Validation (jede Phage einmal als Testset)
    lopo_acc = []
    lopo_f1 = []
    lopo_precision = []
    lopo_recall = []
    for test_phage in phagenames:
        test_genes = phage_to_genes[test_phage]
        train_genes = [gene for gene in features.index if gene not in test_genes]
        X_train = features.loc[train_genes]
        y_train = y_to_use.loc[train_genes]
        X_test = features.loc[test_genes]
        y_test = y_to_use.loc[test_genes]
        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)
        lopo_acc.append(accuracy_score(y_test, y_pred))
        lopo_f1.append(f1_score(y_test, y_pred, average='macro'))
        lopo_precision.append(precision_score(y_test, y_pred, average='macro'))
        lopo_recall.append(recall_score(y_test, y_pred, average='macro'))
        # Optional: Pro Phage anzeigen
        print(f"[{name}] Test-Phage: {test_phage}, Accuracy={lopo_acc[-1]:.4f}, F1={lopo_f1[-1]:.4f}, Precision={lopo_precision[-1]:.4f}, Recall={lopo_recall[-1]:.4f}")

    results['Leave-One-Phage-Out Mean Accuracy'] = np.mean(lopo_acc)
    results['Leave-One-Phage-Out Mean F1'] = np.mean(lopo_f1)
    results['Leave-One-Phage-Out Mean Precision'] = np.mean(lopo_precision)
    results['Leave-One-Phage-Out Mean Recall'] = np.mean(lopo_recall)
        
    return results

comparison_results = {}

for name, model in models.items():
    print(f"\nEvaluating {name}...")
    comparison_results[name] = evaluate_model(name, model, features, labels)

# Ergebnisse als DataFrame
results_df = pd.DataFrame(comparison_results).T
results_df.to_csv("output/models/Vergleichstabelle.csv", index=True)
print("\n=== Vergleichstabelle ===")
print(results_df)

# Ergebnisse plotten
results_df[['5-Fold Mean Accuracy', 'Stratified K-Fold Mean Accuracy', 'Holdout Accuracy', 'Leave-One-Phage-Out Mean Accuracy']].plot(
    kind='bar', figsize=(10, 6), color=['#FF9999', '#99FF99', '#9999FF', '#FFD699'])
plt.title("Accuracy Vergleich der Modelle")
plt.ylabel("Accuracy")
plt.xticks(rotation=45)
plt.ylim(0, 1)
plt.tight_layout()
plt.savefig("output/models/accuracy_comparison.png")

results_df[['5-Fold Mean F1', 'Stratified K-Fold Mean F1', 'Holdout F1', 'Leave-One-Phage-Out Mean F1']].plot(
    kind='bar', figsize=(10, 6), color=['#FF9999', '#99FF99', '#9999FF', '#FFD699'])
plt.title("F1-Score Vergleich der Modelle")
plt.ylabel("F1 Score")
plt.xticks(rotation=45)
plt.ylim(0, 1)
plt.tight_layout()
plt.savefig("output/models/f1_comparison.png")

results_df[['5-Fold Mean Precision', 'Stratified K-Fold Mean Precision', 'Holdout Precision', 'Leave-One-Phage-Out Mean Precision']].plot(
    kind='bar', figsize=(10, 6), color=['#FF9999', '#99FF99', '#9999FF', '#FFD699'])
plt.title("Precision-Score Vergleich der Modelle")
plt.ylabel("Precision Score")
plt.xticks(rotation=45)
plt.ylim(0, 1)
plt.tight_layout()
plt.savefig("output/models/precision_comparison.png")

results_df[['5-Fold Mean Recall', 'Stratified K-Fold Mean Recall', 'Holdout Recall', 'Leave-One-Phage-Out Mean Recall']].plot(
    kind='bar', figsize=(10, 6), color=['#FF9999', '#99FF99', '#9999FF', '#FFD699'])
plt.title("Recall-Score Vergleich der Modelle")
plt.ylabel("Recall Score")
plt.xticks(rotation=45)
plt.ylim(0, 1)
plt.tight_layout()
plt.savefig("output/models/recall_comparison.png")

