#!usr/bin/env python3

# import packages
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from statannot import add_stat_annotation
import itertools

# import data
df = pd.read_csv("../../data/foldX/SNV_free_energies.tsv", sep="\t")
df = df.drop_duplicates()
df = df[~df["Mutation"].str.contains("STOP")]

iarc_df = pd.read_csv("../../data/TVDr20.tsv", sep="\t")
#hotspot_df = pd.read_csv("../../results/hotspots/somatic/all/all/all/amino_acid_mutations.txt")
iarc_df = iarc_df[iarc_df["Effect"] == "missense"]
iarc_muts = np.array(iarc_df["ProtDescription"].unique())
iarc_muts_counts = np.array([iarc_df[iarc_df["ProtDescription"] == m].shape[0] for m in iarc_muts])

df = df[df["Mutation"].str.contains("|".join(iarc_muts))]

# domain dictionary
domain_dict = {"TAD1":[1,40],
        "TAD2": [41, 61],
        "PR":[64,92],
        "DBD":[100,300],
        "OD":[323, 355],
        "C-terminus":[364, 393]}

domains = []
for mut in df["Mutation"].to_list():
    codon = int(mut[3:-1])
    for k, v in domain_dict.items():
        start, end = v
        if start <= codon <= end:
            domain = k
            break
        else:
            domain = "Other" 
            pass
    domains.append(domain)

df["Domain"] = domains
df["Free Energy Difference (MT-WT)"] = df["Free_Energy_(kcal/mol)"]-338.445

# make the figure
order = ["TAD1", "TAD2", "PR", "DBD", "OD", "C-terminus", "Other"]
pairs = list(itertools.combinations(order, 2))

fig, ax = plt.subplots()
ax = sns.violinplot(x="Domain", y="Free Energy Difference (MT-WT)", data = df, order=order)
#add_stat_annotation(ax=ax, x="Domain", y="Free Energy Difference (MT-WT)", data = df, order=order,
                    #box_pairs = pairs, test="Mann-Whitney", text_format="star", loc="inside",
                    #verbose=2, comparisons_correction=None)
plt.show()
