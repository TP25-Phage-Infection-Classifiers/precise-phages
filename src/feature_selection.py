import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.inspection import permutation_importance
import shap 
import matplotlib.pyplot as plt

# feature Datei wird eingelesen
df = pd.read_csv("output/feature_matrix_with_structure.csv", index_col=0)

y = df.loc["Temporal_Class"]
X = df.drop("Temporal_Class", axis=0).T

# Train und Testsplit
X_train, X_test, y_train, y_test = train_test_split(
    X, y, stratify=y, random_state=42
)

model = RandomForestClassifier(n_estimators=100, random_state=42)
model.fit(X_train, y_train)

#Permutation Importance berechnen
perm = permutation_importance(
    model, X_test, y_test,
    n_repeats=5, random_state=42, scoring="accuracy" # repeats von 10 auf 5 für laufzeit
)
perm_imp = pd.Series(perm.importances_mean, index=X.columns) \
             .sort_values(ascending=False)

#SHAP-Werte berechnen
explainer = shap.TreeExplainer(model)
raw_shap_vals = explainer.shap_values(X_test)
n_feat = X_test.shape[1]
shap_vals = [
    arr.T if arr.shape[0] == n_feat else arr
    for arr in raw_shap_vals
]
shap_df  = pd.DataFrame(shap_vals[1], columns=X_test.columns)
shap_imp = shap_df.abs().mean().sort_values(ascending=False)

#Irrelevante Features erkennen dafür (5%-Perzentil)
quantile = 0.05
perm_thresh = perm_imp.quantile(quantile)
shap_thresh = shap_imp.quantile(quantile)

irrelevant_perm = set(perm_imp[perm_imp < perm_thresh].index)
irrelevant_shap = set(shap_imp[shap_imp < shap_thresh].index)

to_drop = sorted(irrelevant_perm | irrelevant_shap)
##print(f"Features zum Entfernen ({len(to_drop)}):", to_drop)

X_reduced = X.drop(columns=to_drop)

df_reduced = pd.concat([
    X_reduced.T,                           # transponierte Features
    pd.DataFrame([y], index=["Temporal_Class"])  # Klasse als neue Zeile
])

df_reduced.to_csv("output/feature_engineering_reduced.csv")
print("Neue Datei geschrieben: output/feature_engineering_reduced.csv")

