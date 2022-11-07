#!/usr/bin/env python3

#import packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import kstest
import argparse

# input arguments
parser = argparse.ArgumentParser()
parser.add_argument("short_topo", help = "the column in IARC with tissue types in all capital letters (case sensitive).  Write the full list for classes including multiple tissue types")
parser.add_argument("bio_sex", help = "M, F, or all")
parser.add_argument("output_dir", help = "Output directory")
args = parser.parse_args()
short_topo = args.short_topo
bio_sex = args.bio_sex
output_dir = args.output_dir

#import the data
iarc_data = pd.read_csv('../../data/TVDr20.tsv', sep="\t")

#filtering by biological sex
if bio_sex != "all":
    iarc_data = iarc_data[iarc_data["Sex"] == bio_sex]
else:
    pass

#filtering by tissue type
if short_topo != "all tissues":
    iarc_data = iarc_data[iarc_data["Short_topo"] == short_topo]
else:
    pass

#deletions
deletions = iarc_data[iarc_data["Type"] == "del"]
descriptions = deletions["Description"]
description_nums = descriptions.str.split("del", expand = True)
data = description_nums[1].dropna()
data = data[~data.str.contains("|".join(["exon", "\?"]))].to_numpy().astype("int")

data = np.array(data)
hist, bins = np.histogram(data)
plt.figure()
logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]))
plt.hist(data, bins=logbins, ec = 'k')
plt.title("Base Pair Deletions in Somatic Mutations of TP53" + " " + short_topo + " "  + bio_sex)
plt.xlabel("Length of Deletion")
plt.ylabel("Frequency")
plt.xscale("log")
plt.yscale("log")
plt.savefig(f"{output_dir}/{short_topo}_{bio_sex}_basepair_deletion_dist.png", bbox_inches="tight")

#insertions
insertions = iarc_data[iarc_data["Type"] == "ins"]
descriptions = insertions["Description"].to_list()

data = []
for mut in descriptions:
    if ("dup Promoter" not in mut) and ("exon" not in mut) and ("del" not in mut):
        mut = mut.split("ins")[1]
        if mut != "":
            mut = int(mut)
            data.append(mut)

hist, bins = np.histogram(data)
plt.figure()
logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]))
plt.hist(data, bins=logbins, ec = 'k')
plt.title("Base Pair Insertions in Somatic Mutations of TP53" + " " + short_topo + " " + bio_sex)
plt.xlabel("Length of Insertion")
plt.ylabel("Frequency")
plt.xscale("log")
plt.yscale("log")
plt.savefig(f"{output_dir}/{short_topo}_{bio_sex}_basepair_insertion_dist.png", bbox_inches="tight")

#frameshift
frameshifts = iarc_data[iarc_data["Effect"] == "FS"]
c_description = frameshifts["c_description"].to_list()

if len(c_description) == 0:
    pass
else:
    data = []
    for mut in c_description:
        if ("?" not in mut) and ("+" not in mut) and ("-" not in mut):
            mut = mut.split(".")[1].split("del")[0].split("ins")[0].split("_")
            if len(mut) == 2:
                start,end = mut
                length_of_fs = int(end) - int(start) + 1
            else:
                length_of_fs = 1
            data.append(length_of_fs)

    hist,bins = np.histogram(data)
    plt.figure()
    logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]))
    plt.hist(data, bins = logbins, ec = 'k')
    plt.title("Frameshift Mutations in Somatic Mutations of TP53" + " " + short_topo + " " + bio_sex)
    plt.xlabel("Number of Base Pairs")
    plt.xscale("log")
    plt.yscale("log")
    plt.savefig(f"{output_dir}/{short_topo}_{bio_sex}_frameshift_dist.png", bbox_inches="tight")
