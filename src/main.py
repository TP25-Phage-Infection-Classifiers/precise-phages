from pathlib import Path
import os
import pandas as pd
import numpy as np
import seaborn as sb
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from PIL import Image
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord 
# pip3 install bcbio-gff
from BCBio import GFF


import preprocess as pre 
import label  
import plots
import utils
import localize

 

# === DATEIEN VERARBEITEN ===
boxplot_output_dir = utils.output_dir / "boxplots"
pie_output_dir = utils.output_dir / "pie chart"
boxplot_output_dir.mkdir(exist_ok=True) 
pie_output_dir.mkdir(exist_ok=True)


plots.generate_output()

print("Gesamte Label Verteilung: ")
print(label.total_label_counts) # in zwei Zeilen geprintet damit Formatierung in Terminal schöner
plots.draw_piechart(label.total_label_counts, "all_files") # pie chart für alle Daten

# Die Verteilung der Labels 'early', 'middle' und 'late' über alle Datensets ist relativ ausgeglichen. 
# Nur ca. 3% der Gene waren durchgehen so niedrig exprimiert, dass sie als undefined gelabled wurden. 

png_files = ["output/pie chart/Ceyssens_directional_full_raw_counts_pie.png", 
             "output/pie chart/Finstrlova_SH1000_full_raw_counts_pie.png",
             "output/pie chart/Guegler_T4_plusToxIN_full_raw_counts_pie.png",
             "output/pie chart/Kuptsov_full_raw_counts_pie.png",
             "output/pie chart/Meaden_BIM_full_raw_counts_pie.png",
             "output/pie chart/Sprenger_VC_WT_VP882_WT_full_raw_counts_pie.png",
             "output/pie chart/Zhong_full_raw_counts_pie.png",
             "output/pie chart/all_files_pie.png"]

images = [Image.open(f).convert("RGB") for f in png_files]
images[0].save("output.pdf", save_all=True, append_images=images[1:])

localize.generate_genome_map()
