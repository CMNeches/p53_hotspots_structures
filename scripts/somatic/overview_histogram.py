#!usr/bin/env/python3

# import packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import kstest

# load in the data
df = pd.read_csv('../../data/TVDr20.tsv', sep='\t')

# prepare to make the histograms
codon_array = np.arange(1,394)
all_codons = df["Codon_number"].to_numpy()
mutation_counts = np.array([(all_codons==c).sum() for c in codon_array])

plt.figure()
plt.title("Codon Position vs. Mutation Counts in Somatic Mutations of TP53")
plt.bar(codon_array, mutation_counts)
plt.xlabel("Codon Postion")
plt.ylabel("Mutatation Counts")
plt.show()
