#### DNA-basierte Features
Vorteile:

- robust, numerisch -> GC-Gehalt, Genlänge
- Direkter Bezug zur Genexpression
- Regulatorische Informationen enthalten
- Einfach aus FASTA-Dateien extrahierbar
- Gut für Klassifikation nach Expressionszeitpunkt

Nachteile:

- keine direkten Informationen über Funktionalität
- Begrenzte Information über Proteinfaltung oder Interaktion

#### Protein-basierte Features
Vorteile:

- Rückschlüsse auf Funktion des Proteins möglich
- Biologisch interpretierbar

Nachteile:

- Berechnung oft komplexer
- Kein direkter Bezug zur Transkriptionsregulation

-> Wir haben uns für die Verwendung von DNA-basierten Features entschieden.

Gründe:

- größere Auswahl an möglichen Features
- Großteil unserer bevorzugetn Features funktioniert besser auf DNA-Ebene
- Simplere Features wie Genlänge oder GC-Gehalt gut geeignet für ML und leichter umsetzbar
- Promotormuster sind sehr aussagekräftig über Expressionszeitpunkt -> auf DNA-Ebene
    
    
Features:

- Genlänge
- GC-Gehalt
- K-Mer-Frequenz
- Position im Genom (Gencluster)
- Promotor-Regionen -> Motiverkennung
- Codon Usage Bias
