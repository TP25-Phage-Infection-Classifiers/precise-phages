"""
utils.py

This module contains utility functions for data preprocessing and evaluation.
"""

import pandas as pd
from sklearn.metrics import accuracy_score
from scipy.stats import zscore
import numpy as np


def tsv_to_csv(tsv_filepath, csv_filepath):
    # load tsv file in dataframe
    df = pd.read_csv(tsv_filepath, sep='\t')  # 'sep' gibt an, dass die Trennung durch Tabulatoren erfolgt
    
    # save dataframe as csv file
    df.to_csv(csv_filepath, index=False)  # 'index=False' entfernt die Index-Spalte beim Speichern
    
def load_data(filepath):
    """Load dataset from a CSV file."""
    return pd.read_csv(filepath, sep='\t')

def identify_timepoint_columns(data):
    """Identifiziere Zeitreihen-Spalten (alles zwischen erster und den letzten zwei Spalten)."""
    columns = data.columns.tolist()
    return columns[1:-2]  # alle Spalten außer erster (Geneid) und den letzten beiden (Entity, Symbol)

def evaluate_model(model, X, y):
    """Evaluate the model and print accuracy."""
    predictions = model.predict(X)
    accuracy = accuracy_score(y, predictions)
    print(f"Model Accuracy: {accuracy:.2f}")

filepath = "data/brandao_2021/Brandao_LB_full_raw_counts.tsv"
output_filepath = "data/brandao_2021/Brandao_with_outliers_marked.csv"  # Name der neuen Datei


def detectOutliers(filepath, output_filepath):

    # 1. Daten einlesen
    data = load_data(filepath)

    # 2. Identifizierung der Zeitpunkts-Spalten
    time_columns = identify_timepoint_columns(data)

    # 3. Berechnung der IQR für jede Zeitreihe
    outliers_iqr = pd.DataFrame(False, index=data.index, columns=time_columns)  # DataFrame für Ausreißer

    for col in time_columns:
        # Berechne Q1 (25. Perzentil) und Q3 (75. Perzentil)
        Q1 = data[col].quantile(0.25)
        Q3 = data[col].quantile(0.75)
        
        # Berechne den IQR
        IQR = Q3 - Q1
        
        # Bestimme die Grenzen für Ausreißer
        lower_bound = Q1 - 1.5 * IQR
        upper_bound = Q3 + 1.5 * IQR
        
        # Identifiziere die Ausreißer: Werte unter der unteren Grenze oder über der oberen Grenze
        outliers_iqr[col] = (data[col] < lower_bound) | (data[col] > upper_bound)

    # 4. Entfernen der Ausreißer aus den Zeitspalten (ersetzen durch NaN)
    for idx, col in enumerate(time_columns):
        # Setze alle Ausreißer (True-Werte) in den entsprechenden Spalten auf NaN
        data.loc[outliers_iqr[col], col] = np.nan

    # 6. Speichern der bereinigten Datei in einer neuen CSV-Datei
    data.to_csv(output_filepath, index=False)  # Speichert das aktualisierte DataFrame ohne Ausreißer
    
    print(f"Gefundene Ausreißer insgesamt: {outliers_iqr.sum().sum()}")
    print(data.head())
    # Rückgabe des bereinigten DataFrames
    return data

if __name__ == "__main__":
    print(detectOutliers(filepath,output_filepath))



