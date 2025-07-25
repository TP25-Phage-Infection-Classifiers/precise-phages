from pathlib import Path
import pandas as pd


# Hilfsfunktion: Phage-Gruppe aus Struktur-GeneID ableiten

def extract_phage_group(name: str) -> str | None:
    s = name.lower()
    if "dms3" in s:
        return "DMS3"
    if "vpvv882" in s or "vp882" in s:
        return "VPVV882"
    if "phikz" in s:
        return "PHIKZ"
    if "phiyy" in s:
        return "phiYY"
    if "phage_k" in s and "sau" not in s:
        return "CPT_phageK"
    if "sau" in s or "515a1" in s:
        return "phage515_"
    if "t4" in s:
        return "T4"
    return None



# Merge-Funktion

def merge_feature_and_structure(feature_csv: Path, structure_csv: Path, out_csv: Path):
    # Lade Feature-Matrix (Gene als Spalten, letzte Zeile = Labels)
    df = pd.read_csv(feature_csv, index_col=0)
    labels = df.loc["Temporal_Class"].copy()
    df = df.drop(index="Temporal_Class")

    # Transponieren: Gene in Zeilen, Features in Spalten
    feature_df = df.T.copy()
    feature_df["Temporal_Class"] = labels
    feature_df.reset_index(inplace=True)
    feature_df = feature_df.rename(columns={"index": "GeneID"})

    # Lade Strukturdaten und berechne Mittelwerte pro Phage-Gruppe
    struct = pd.read_csv(structure_csv)
    struct["PhageGroup"] = struct["GeneID"].apply(extract_phage_group)
    struct_grouped = struct.groupby("PhageGroup").mean(numeric_only=True).reset_index()

    # Mappe Phage-Gruppe zu jedem Gen in Feature-Matrix
    feature_df["PhageGroup"] = feature_df["GeneID"].apply(
        lambda x: next((g for g in struct_grouped["PhageGroup"] if g and g.lower() in x.lower()), None)
    )

    # Mergen per PhageGroup
    merged = feature_df.merge(struct_grouped, how="left", left_on="PhageGroup", right_on="PhageGroup")

    # Aufräumen und speichern
    merged.drop(columns=["PhageGroup"], inplace=True)
    merged.set_index("GeneID", inplace=True)
    merged = merged.astype(float, errors="ignore")
    merged.to_csv(out_csv)
    print(f" Feature-Matrix mit Struktur gespeichert: {out_csv.resolve()}")



# Ausführung

if __name__ == "__main__":
    merge_feature_and_structure(
        feature_csv=Path("output/feature_engineering_merged.csv"),
        structure_csv=Path("output/structure_features.csv"),
        out_csv=Path("output/feature_matrix_with_structure.csv")
    )
