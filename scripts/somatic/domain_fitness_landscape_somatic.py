#!/usr/bin/env python3

#import packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import kstest
import argparse
import os

#input arguments
parser = argparse.ArgumentParser()
parser.add_argument("bio_sex", help = "M, F, or all")
args = parser.parse_args()
bio_sex = args.bio_sex

#read in fitness database
fit_df = pd.read_csv("../../data/fitness/transactivation/systematicFunctionalAssessmentIARC_TP53_Database_R20.txt", sep="\t")
cols = ["ProtDescription","WAF1nWT", "MDM2nWT", "BAXnWT", "h1433snWT", "AIP1nWT", "GADD45nWT", "NOXAnWT", "P53R2nWT"]
fit_df = fit_df[cols]
fit_df["Fitness"] = fit_df[cols[1:]].median(axis=1)

domain_dict = {"TAD1": [1,40],
                "TAD2": [41,61],
                "PR" : [64,92],
                "DBD" : [100, 300],
                "OD" : [323, 355],
                "C-terminus": [364, 393]}

domains = ["TAD1", "TAD2", "PR", "DBD", "OD", "C-terminus"]

# go across all tissues and read in domain statistics
def get_domain_stats(stats_file):
    stats_df = pd.read_csv(stats_file, sep="\t")
    stats_df = stats_df[stats_df["Amino Acid Mutation"].str.contains("|".join(fit_df["ProtDescription"].to_list()))]
    mutations = stats_df["Amino Acid Mutation"].to_numpy()
    codons = np.array([int(m[3:-1]) for m in mutations])
    counts = stats_df["Mutation Counts"].to_numpy()
    weights = counts/counts.sum()
    domain_fits = np.zeros(6)
    for i,domain in enumerate(domains):
        domain_range = domain_dict[domain]
        indices = np.where(np.logical_and((domain_range[0] <= codons), (codons <= domain_range[1])))
        domain_muts = mutations[indices]
        domain_weights = weights[indices]
        domain_fit = []
        for m in domain_muts:
            f = fit_df[fit_df["ProtDescription"] == m]["Fitness"].to_list()[0]
            domain_fit.append(f)
        domain_fit = (domain_weights * domain_fit).sum()
        domain_fits[i] = domain_fit
        
        if domain_fits[i] == 0:
            domain_fits[i] = 100
    return domain_fits

prefix = "../../results/hotspots/somatic/"
df = []
tissues = []
for tissue in os.listdir(prefix):
    stats_file = f"{prefix}/{tissue}/missense/{bio_sex}/amino_acid_mutations.txt"
    if os.path.isfile(stats_file):
        domain_array = get_domain_stats(stats_file)
        df.append(list(domain_array))
        tissues.append(tissue)

df = pd.DataFrame(df)
df.columns = ["TAD1", "TAD2", "PR", "DBD", "OD", "C-term"]
df.index = tissues
df = df.fillna(100)

print(df)

ax = sns.clustermap(df, yticklabels=True)
ax.ax_heatmap.set_yticklabels(ax.ax_heatmap.get_ymajorticklabels(), fontsize=8)
plt.savefig(f"../../results/fitness/missense_{bio_sex}_domain_fit.png", bbox_inches="tight")
plt.show()
