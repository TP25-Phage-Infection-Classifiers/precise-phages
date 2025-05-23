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


# ORDNER
project_root = Path(__file__).resolve().parents[1]
input_dir = project_root / "data"
output_dir = project_root / "output"
output_dir.mkdir(exist_ok=True)

# UNSERE DATEIEN
my_team_files = [
    "ceyssens_2014/Ceyssens_directional_full_raw_counts.tsv",
    "finstrlova_2022/Finstrlova_SH1000_full_raw_counts.tsv",
    "guegler_2021/Guegler_T4_plusToxIN_full_raw_counts.tsv",
    "kuptsov_2022/Kuptsov_full_raw_counts.tsv",
    "meaden_2021/Meaden_BIM_full_raw_counts.tsv",
    "sprenger_2024/Sprenger_VC_WT_VP882_WT_full_raw_counts.tsv",
    "zhong_2020/Zhong_full_raw_counts.tsv"
]
my_gff3_files=[
    "ceyssens_2014/Pseudomonas_phage_phiKZ.gff3",
    "finstrlova_2022/Staphylococcus_phage_K.gff3",
    "guegler_2021/Enterobacteria_phage_T4.gff3",
    "kuptsov_2022/Staphylococcus_phage_vB_SauM-515A1.gff3",
    "meaden_2021/Pseudomonas_phage_DMS3.gff3",
    "sprenger_2024/Vibrio_phage_VP882.gff3",
    "zhong_2020/Pseudomonas_phage_phiYY_complete.gff3"
]
my_genome_flies=[
    "ceyssens_2014/Pseudomonas_phage_phiKZ.fasta",
    "finstrlova_2022/Staphylococcus_phage_K.fasta",
    "guegler_2021/Enterobacteria_phage_T4.fasta",
    "kuptsov_2022/Staphylococcus_phage_vB_SauM-515A1.fasta",
    "meaden_2021/Pseudomonas_phage_DMS3.fasta",
    "sprenger_2024/Vibrio_phage_VP882.fasta",
    "zhong_2020/Pseudomonas_phage_phiYY_complete.fasta",
]