#!/usr/bin/env python3

#import packages 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import kstest
import argparse

#input arguments
parser = argparse.ArgumentParser()
parser.add_argument("short_topo", help = 'Use the IARC "Short_topo" column, and remember to use all capital letters')
parser.add_argument("mutation_type", help = 'the "Effect" column in IARC')
parser.add_argument("bio_sex", help = "M, F, or all")
args = parser.parse_args()
short_topo = args.short_topo
mutation_type = args.mutation_type
bio_sex = args.bio_sex
output_dir = args.output_dir

#import the data
iarc_data = pd.read_csv('../../data/Germline.tsv', sep="\t")

#filtering by biological sex
if bio_sex != "all":
    iarc_data = iarc_data[iarc_data["Sex"] == bio_sex]
else:
    pass

df = iarc_data[iarc_data["Short_topo"] == short_topo]

#make all of the graphs for the different types
codon_array = np.arange(1,394)
all_codons = df["Codon_number"].to_numpy()
mutation_counts = np.array([(all_codons == c).sum() for c in codon_array])
temporary_df = df[df["Effect"] == t]
all_codons = temporary_df["Codon_number"].to_numpy()
mutation_counts = np.array([(all_codons == c).sum() for c in codon_array])
plt.figure()

#whole gene
plt.title("Codon Position vs. Mutation Counts in " + mutation_type + " Germline Mutations in " + short_topo + " " + bio_sex + " of TP53")
plt.bar(codon_array[1:41], mutation_counts[1:41], color = 'b')
plt.bar(codon_array[41:62], mutation_counts[41:62], color = 'r')
plt.bar(codon_array[65:93], mutation_counts[65:93], color = 'g')
plt.bar(codon_array[101:301], mutation_counts[101:301], color = 'c')
plt.bar(codon_array[324:356], mutation_counts[324:356], color = 'm')
plt.bar(codon_array[365:394], mutation_counts[365:394], color = 'y')
plt.xlabel("Codon Position")
plt.ylabel("Mutation Counts")
plt.yscale("log")
counts = codon_array.sum()
kstest = kstest(codon_array, 'uniform')

#normalized frequency
plt.bar(codon_array[1:41], mutation_counts[1:41]/mutation_counts[1:41].sum(), color = 'b')
plt.bar(codon_array[41:62], mutation_counts[41:62]/mutation_counts[41:62].sum(), color = 'r')
plt.bar(codon_array[65:93], mutation_counts[65:93]/mutation_counts[65:93].sum(), color = 'g')
plt.bar(codon_array[101:301], mutation_counts[101:301]/mutation_counts[101:301].sum(), color = 'c')
plt.bar(codon_array[324:356], mutation_counts[324:356]/mutation_counts[324:356].sum(), color = 'm')
plt.bar(codon_array[365:394], mutation_counts[365:394]/mutation_counts[365:394].sum(), color = 'y')
plt.title("Codon Position vs. Frequency of Mutation in " + mutation_type + " Germline Mutations in " + short_topo + " " + bio_sex + " of TP53")
plt.xlabel("Codon Position")
plt.ylabel("Frequency")

#per domain graphs
fig, axs = plt.subplots(2,3)
axs[0,0].bar(codon_array[1:41], mutation_counts[1:41]/mutation_counts[1:41].sum(), color = 'b')
axs[0,0].set_title("TAD1")
axs[0,0].set_xlabel("Codon Position")
axs[0,0].set_ylabel("Mutation Frequency (%)")
counts_TAD1 = codon_array[1:41].sum()
kstest_TAD1 = kstest(codon_array[1:41], 'uniform')

axs[0,1].bar(codon_array[41:62], mutation_counts[41:62]/mutation_counts[41:62].sum(), color = 'r')
axs[0,1].set_title("TAD2")
axs[0,1].set_xlabel("Codon Position")
axs[0,1].set_ylabel("Mutation Frequency (%)")
counts_TAD2 = codon_array[41:62].sum()
kstest_TAD2 = kstest(codon_array[41:62], 'uniform')

axs[0,2].bar(codon_array[65:93], mutation_counts[65:93]/mutation_counts[65:93].sum(), color = 'g')
axs[0,2].set_title("PR")
axs[0,2].set_xlabel("Codon Position")
axs[0,2].set_ylabel("Mutation Frequency (%)")
counts_PR = codon_array[65:93].sum()
kstest_PR = kstest(codon_array[65:93], 'uniform')

axs[1,0].bar(codon_array[101:301], mutation_counts[101:301]/mutation_counts[101:301].sum(), color = 'c')
axs[1,0].set_title("DBD")
axs[1,0].set_xlabel("Codon Position")
axs[1,0].set_ylabel("Mutation Frequency (%)")
counts_DBD = codon_array[101:301].sum()
kstest_DBD = kstest(codon_array[101:301], 'uniform')

axs[1,1].bar(codon_array[323:355], mutation_counts[323:355]/mutation_counts[323:355].sum(), color = 'm')
axs[1,1].set_title("OD")
axs[1,1].set_xlabel("Codon Position")
axs[1,1].set_ylabel("Mutation Frequency (%)")
counts_OD = codon_array[323:355].sum()
kstest_OD = kstest(codon_array[323:355], 'uniform')

axs[1,2].bar(codon_array[364:394], mutation_counts[364:394]/mutation_counts[364:394].sum(),color = 'y')
axs[1,2].set_title("C-terminus")
axs[1,2].set_xlabel("Codon Position")
axs[1,2].set_ylabel("Mutation Frequency (%)")
counts_C_terminus = codon_array[364:394].sum()
kstest_C_terminus = kstest(codon_array[364:394], 'uniform')

fig.tight_layout(rect = [0,0.03, 1, 0.95])
fig.suptitle("Frequency of " + t + " Germline Mutations in " + short_topo + " " + bio_sex + " of TP53")
plt.show()

#statistics
counts_list = [counts, counts_TAD1, counts_TAD2, counts_PR, counts_DBD, counts_OD, counts_C_terminus]
kstest_list = [kstest, kstest_TAD1, kstest_TAD2, kstest_PR, kstest_DBD, kstest_OD, kstest_C_terminus]
names_list = ["whole protein", "TAD1", "TAD2", "PR", "DBD", "OD", "C-terminus"]

with open(f"{output_dir}/statistics.txt", 'w') as f:
    to_write = []
    title = "\n".join(["Somatic", f"Mutation_type : {mutation_type}", f"Tissue_type : {short_topo}", f"Biological_sex: {bio_sex}"])
    for n, c, k in zip(names_list, counts_list, kstest_list):
        to_write.append(f"{n} Mutation Counts: {c}")
        to_write.append(f"{n} KS test: {k}")
    f.write(title + "\n" + "\n".join(to_write))
