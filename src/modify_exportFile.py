import pandas as pd
import glob

# Alle Dateien, die auf '_export.csv' enden
dateien = glob.glob('output/**/*_export.csv', recursive=True)

# Zusamenfügen aller _export.csv Dateien
df_liste = [pd.read_csv(datei) for datei in dateien]
zusammen = pd.concat(df_liste, ignore_index=True)
zusammen.to_csv('output/zusammengefuegt.csv', index=False)

# Transponieren der zusammengefuegt.csv Datei

# Datei einlesen
df = pd.read_csv('zusammengefuegt.csv')

# Transponieren, sodass GeneIDs Spalten werden und Temporal_Class die Werte sind
df_transponiert = df.set_index('GeneID').T

# Optional: In eine neue CSV speichern
df_transponiert.to_csv('output/zusammengefuegt_transponiert.csv')