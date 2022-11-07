#!/usr/bin/env python3

#import packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import kstest
import argparse

#input arguments
parser = argparse.ArgumentParser()
parser.add_argument("short_topo", help = "the column in IARC with tissue types in all capital letters (case sensitive).  Input the full name for classes containing multiple tissue types")
parser.add_argument("mutation_type", help = 'the "Effect" column in IARC')
parser.add_argument("bio_sex", help = "M, F, or all")
parser.add_argument("output_dir", help = "Output directory")
args = parser.parse_args()
short_topo = args.short_topo
mutation_type = args.mutation_type
bio_sex = args.bio_sex
output_dir = args.output_dir

# read in codon mutation count information
df = pd.read_csv(f"{output_dir}/codon_mutations.txt", sep="\t")
codon_array = df["Codon"].to_numpy()
mutation_counts = df["Codon Counts"].to_numpy()

# dictionary defining domain regions
domain_dict = {"TAD1" : [1, 40],
               "TAD2" : [41, 61],
               "PR" : [64, 92],
               "DBD" : [100, 300],
               "OD": [323, 355],
               "C-terminus" : [364, 393]}

# domains
domains = ["TAD1", "TAD2", "PR", "DBD", "OD", "C-terminus"]
colors = ['b', 'r', 'g', 'c', 'm', 'y']

# plot frequency normalized by domain

# function to plot frequency distribution for a domain
def make_domain_dist(domain, color, axs_subset):
    # make bar graph
    start, end = domain_dict[domain]
    x_data = codon_array[start - 1 : end]
    y_data = mutation_counts[start - 1 : end]

    if y_data.sum() == 0:
        ks = ("Not Applicable", "Not Applicable")
    else:
        y_data = y_data/y_data.sum()
        axs_subset.bar(x_data, y_data, color=color)

        # label hotspots
        # threshold needs to exclude zeros
        thresh = np.percentile(y_data[y_data != 0], 90)
        indices = (y_data >= thresh)
        # if there are a few points, annotate all of them
        if indices.sum() <= 5:
            hotspots = x_data[y_data != 0]
            hotspot_freqs = y_data[y_data != 0]
        else:
            hotspots = x_data[indices]
            hotspot_freqs = y_data[indices]

        ks = kstest(y_data, 'uniform')

        for h,f in zip(hotspots, hotspot_freqs):
            axs_subset.text(h, f + y_data.max()/10, h, ha='center', 
                            bbox=dict(facecolor='wheat', alpha=0.8))
        
            if not np.isnan(y_data.max()):
                axs_subset.set_ylim(None, 1.2*y_data.max())

    # annotation
    axs_subset.set_title(domain)
    axs_subset.set_xlabel("Codon Position")
    axs_subset.set_ylabel("Mutation Frequency (%)")
    return ks

fig, axs = plt.subplots(2, 3)
ks_tad1 = make_domain_dist("TAD1", 'b', axs[0,0])
ks_tad2 = make_domain_dist("TAD2", 'r', axs[0,1])
ks_PR = make_domain_dist("PR", 'g', axs[0,2])
ks_DBD = make_domain_dist("DBD", 'c', axs[1,0])
ks_OD = make_domain_dist("OD", 'm', axs[1,1])
ks_C = make_domain_dist("C-terminus", 'c', axs[1,2])

fig.tight_layout(rect = [0, 0.03, 1, 0.95])
fig.suptitle("Frequency of " + mutation_type + " in Somatic TP53 mutations " + short_topo + " " + bio_sex)

plt.savefig(f"{output_dir}/domain_mut_dist.png", bbox_inches="tight")

df_stats = pd.DataFrame()
df_stats["Domain"] = domains
df_stats["KS Test Statistic"], df_stats["KS Test p-value"] = zip(*[ks_tad1, ks_tad2, ks_PR, ks_DBD, ks_OD, ks_C])
df_stats.to_csv(f"{output_dir}/domain_mut_dist_stats.tsv", sep="\t", index=False)
