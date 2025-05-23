import pandas as pd
import numpy as np



# ordnet jedem Gen ein Label zu, arbeitet mit den normalisierten, bereinigten Phagengendaten
# hierfür wird geschaut, in welchem drittel der timeslots der maximale Wert an counts liegt. 
# liegt dieser im ersten Drittel wird das Gen als 'early' gelabled, liegt er im mittleren Drittel als 'middle' und im hinteren als 'late'
# Sind die normalisierten Counts alle unter einem Wert von 2.5 wird das Gen als 'undefined gelabled und nicht für das ML verwendet
def classify_temporal_expression(row, time_cols):
    max_value = row[time_cols].max()
    max_timepoint_idx = row[time_cols].values.argmax() # Index des höchsten Werts der Zeile
    n_timepoints = len(time_cols) # Anzahl der Werte
    if max_value < 2.5: # durchgehend niedrig exprimierte Gene werden als undefined definiert 
        return 'undefined'
    elif max_timepoint_idx < n_timepoints / 3: # wenn Index im ersten Drittel liegt -> early
        return 'early'
    elif max_timepoint_idx < 2 * n_timepoints / 3: # wenn Index im zweiten Drittel liegt -> middle
        return 'middle'
    else:
        return 'late'
    
# Label-Distribution Analyse
label_order = ['undefined', 'early', 'middle','late'] # Reihenfolge, in der counts ausgegeben
total_label_counts = pd.Series(0, index=label_order) # Counts aller Dateien, wird später berechnet

# zählt, wie oft jedes Label in einer Datei vorkommt
def count_labels(df):
    counts = df['Temporal_Class'].value_counts()
    return counts.reindex(label_order, fill_value=0)