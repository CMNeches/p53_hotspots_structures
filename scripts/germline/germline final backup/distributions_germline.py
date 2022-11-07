#!/usr/bin/env python3

#import packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import kstest
import argparse

#input arguments
parser = argparse.ArgumentParser()
parser.add_argument("short_topo", help = 'Use the IARC "Short_topo" column, and remember to use all capital letters, or "all tissues" to not filter')
parser.add_argument("bio_sex", help = 'M, F, or all')
parser.add_argument("distribution", help = 'insertion, deletion, frameshift (available in the "Type" column in IARC')
args = parser.parse_args()
short_topo = args.short_topo
bio_sex = args.bio_sex
distribution = args.distribution

#import the data
iarc_data = pd.read_csv('../../data/Germline.tsv', sep="\t")

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
deletions = df[df["Type"] == "del"]
descriptions = deletions["Description"]
description_nums = descriptions.str.split("del", expand = True)
data = descritpion_num[1].dropna()
data = data[~data.str.contains("|".join(["exon", "\?"]))].to_numpy().astype("int")

data = np.array(data)
hist, bins, _ = plt.hist(data)
plt.figure()
logbins = np.logspace(np.log10(bins[0])), np.log10(bins[-1])))
plt.hist(data, bins=logbins, ec='k')
plt.title("Base Pair Deletions in Germline Mutations of TP53" + (mutation_type + " " + short_topo + " " + bio_sex))
plt.xlabel("Number of Base Pairs")
plt.xscale("log")
plt.yscale("log")
plt.show()

#insertions
insertions = df[df["Type"] == "ins"]
descriptions = insertons["Description"].to_list()

data = []
for mut in descriptions:
    if ("dup Promoter" not in mut) and ("exon" not in mut) and ("del" not in mut):
        mut = int(mut.split("ins")[1])
        print(mut)
        data.append(mut)

plt.figure()
plt.hist(data, bins=10, ec = 'k')
plt.xlabel("Length of Insertion")
plt.ylabel("Frequency")
plt.yscale("log")
plt.show()

#frameshift
framshifts = df[df["Effect"] == "FS"]
c_description = frameshifts["c_description"].to_list()
data = []

plt.hist(data, bins = 10, ec = 'k')
plt.title("Frameshift Mutations in Germline Mutations of TP53" + (mutation_type + " " + short_topo + " " + bio+sex))
plt.xlabel("Number of Base Pairs")
plt.xscale("log")
plt.yscale("log")
plt.show()
