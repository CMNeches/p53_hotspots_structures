#!/usr/bin/env python3

#import packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import chisquare
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
    chi_square_test = chisquare(domain_array)
    return chi_square_test

prefix = "../../results/hotspots/germline/"
df = []
tissues = []
for tissue in os.listdir(prefix):
    stats_file = f"{prefix}/{tissue}/{mutation_type}/{bio_sex}/domain_statistics.txt"
    if os.path.isfile(stats_file):
        chi_square_test = get_domain_stats(stats_file)
        df.append(list(chi_square_test))
        tissues.append(tissue)

df = pd.DataFrame(df)
df.columns = ["Chi-square Statistic", "Chi-square p-value"]
df["Tissue"] = tissues
df.to_csv("../../results/heatmaps/chi_square_domain_tissue_statistics_germline.tsv", sep="\t", index = False)

# make the figure
plt.figure()
plt.bar(range(len(tissues)), -np.log10(df["Chi-square p-value"].to_numpy()))
plt.axhline(-np.log10(0.05)/len(tissues), ls = '--', c='r', label = "Bonferroni-corrected significance cutoff")
plt.xticks(range(len(tissues)), tissues, rotation=90, fontsize=6)
plt.yscale("log")
plt.ylabel("-log10(Chi-square p-value)")
plt.legend()
plt.savefig("../../results/heatmaps/chi_square_domain_tissue_statitics_germline.png", bbox_inches="tight")
plt.show()
