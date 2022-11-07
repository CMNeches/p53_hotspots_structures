#!/usr/bin/env python3

#import packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import kstest
import argparse

#input arguments
parser = argparse.ArgumentParser()
parser.add_argument("short_topo", help = "the column in IARC with tissue types in all capital letters (case sensitive).  Write the full list for classes including multiple tissue types")
parser.add_argument("mutation_type", help = 'the "Effect" column in IARC')
parser.add_argument("bio_sex", help = "M, F, or all")
parser.add_argument("output_dir", help = "Output directory")
args = parser.parse_args()
short_topo = args.short_topo
mutation_type = args.mutation_type
bio_sex = args.bio_sex
output_dir = args.output_dir

#load in the data
iarc_data = pd.read_csv('../data/Germline.tsv')

#filtering by biological sex
if bio_sex != "all":
    iarc_data = iarc_data[iarc_data["Sex"] == bio_sex]
else:
    pass

df = iarc_data[iarc_data["Short_topo"] == short_topo]

#we are defining the codon array
codon_array = np.arange(1,394)
all_codons = df["Codon_number"].to_numpy()

mutation_counts = np.array([(all_codons == c).sum() for c in codon_array])

temporary_df = df[df["Effect"] == t]
codon_array = np.arange(1,394)
all_codons = temporary_df["Codon_number"].to_numpy()
mutation_counts = np.array([(all_codons == c).sum() for c in codon_array])
    
#discovering novel hotspot locations
#mutation counts is the number of mutations at each codon, and then codon array is an array from 1 to 393, argsort sorts an array in ascending order.  That is why we sort on the negative of the mutation count array (we take the top 10 mutations)
for k,v in domain_dict.items():
    start, end, _ = domain_dict[k]
    new_codon_array = codon_array[start : end]
    new_mutation_counts = mutation_counts[start : end]
    hotspots = new_codon_array[(-new_mutation_counts).argsort()][:10]
    freqs = (new_mutation_counts[(-new_mutation_counts).argsort()]/new_mutation_counts.sum())[:10]

with open(f"{output_dir}/hotspots.txt", "w") as f:
    title = [f"{short_topo}", f"{mutation_type}", f"{bio_sex}"]
    to_write = "\n".join([str(hotspots), str(freqs)]
    f.write(title + "\n" + to_write + "\n")
