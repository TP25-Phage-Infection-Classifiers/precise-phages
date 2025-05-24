from Bio import SeqIO

def fasta_parser(filepath):
    return list(SeqIO.parse(filepath, "fasta"))

# Genlänge
def gene_length(filepath):
    lengths = []
    for gene in fasta_parser(filepath):
        lengths.append(len(gene.seq))
    return lengths
        
        
# GC-Gehalt 
def gc_content(filepath):
    GCs = []
    for gene in fasta_parser(filepath):
        count = 0
        for base in gene.seq.upper():
            if base == 'G' or base == 'C':
                count += 1
        GCs.append(count / len(gene.seq))
    return GCs