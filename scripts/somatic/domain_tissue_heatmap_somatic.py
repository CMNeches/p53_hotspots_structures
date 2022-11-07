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
parser.add_argument("mutation_type", help = 'the "Effect" column in IARC')
parser.add_argument("bio_sex", help = "M, F, or all")
args = parser.parse_args()
mutation_type = args.mutation_type
bio_sex = args.bio_sex

# go across all tissues and read in domain statistics
def get_domain_stats(stats_file):
    stats_df = pd.read_csv(stats_file, sep="\t")
    domain_array = stats_df["Domain Cross Section"].to_numpy()
    normalization_constant = stats_df["Number of Mutations in Domain"].sum()
    domain_array = domain_array/normalization_constant
    return domain_array

prefix = "../../results/hotspots/somatic/"
df = []
tissues = []
for tissue in os.listdir(prefix):
    stats_file = f"{prefix}/{tissue}/{mutation_type}/{bio_sex}/domain_statistics.txt"
    if os.path.isfile(stats_file):
        domain_array = get_domain_stats(stats_file)
        df.append(list(domain_array))
        tissues.append(tissue)

df = pd.DataFrame(df)
df.columns = ["TAD1", "TAD2", "PR", "DBD", "OD", "C-term"]
df.index = tissues
df = df.fillna(0)

print(df)

ax = sns.clustermap(df, yticklabels=True)
ax.ax_heatmap.set_yticklabels(ax.ax_heatmap.get_ymajorticklabels(), fontsize=8)
plt.show()
#plt.savefig(f"../../results/heatmaps/{mutation_type}_{bio_sex}_domain_dist_heatmap.png", bbox_inches="tight")
