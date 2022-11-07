#!/usr/bin/env python3

#import packages 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#import the data
df = pd.read_csv('../data/Germline.csv')
#print(df)
df = df[df["Sex"] == "M"]

#pie chart
piedf = df["Effect"].value_counts()
print(piedf)
labels = ["Missense", "FS", "Nonsense", "Splice", "Other", "Large Del", "Intronic", "Silent"]
sizes = [983, 124, 117, 11, 24, 22, 7, 2]
plt.figure()
plt.pie(sizes, labels=labels, autopct = '%1.1f%%', startangle=90)
plt.show()

#make all of the graphs for the different types
codon_array = np.arange(1,394)
all_codons = df["Codon_number"].to_numpy()
#print(all_codons)

mutation_counts = np.array([(all_codons == c).sum() for c in codon_array])

types = ["missense", "FS", "nonsense", "silent", "splice", "other", "intronic", "large del"]
for t in types:
    temporary_df = df[df["Effect"] == t]
    codon_array = np.arange(1, 394)
    all_codons = temporary_df["Codon_number"].to_numpy()
    mutation_counts = np.array([(all_codons == c).sum() for c in codon_array])
    plt.figure()

    #task 2 a from hw
    plt.title("Codon Position vs. Mutation Counts in " + t + " Germline Mutations in TP53")
    plt.bar(codon_array[1:41], mutation_counts[1:41], color = 'b')
    plt.bar(codon_array[41:62], mutation_counts[41:62], color = 'r')
    plt.bar(codon_array[65:93], mutation_counts[65:93], color = 'g')
    plt.bar(codon_array[101:301], mutation_counts[101:301], color = 'c')
    plt.bar(codon_array[324:356], mutation_counts[324:356], color = 'm')
    plt.bar(codon_array[365:394], mutation_counts[365:394], color = 'y')
    plt.xlabel("Codon Position")
    plt.ylabel("Mutation Counts")
    plt.yscale("log")
    plt.show()
    plt.close()

    #task 2 b from hw
    plt.bar(codon_array[1:41], mutation_counts[1:41]/mutation_counts[1:41].sum(), color = 'b')
    plt.bar(codon_array[41:62], mutation_counts[41:62]/mutation_counts[41:62].sum(), color = 'r')
    plt.bar(codon_array[65:93], mutation_counts[65:93]/mutation_counts[65:93].sum(), color = 'g')
    plt.bar(codon_array[101:301], mutation_counts[101:301]/mutation_counts[101:301].sum(), color = 'c')
    plt.bar(codon_array[324:356], mutation_counts[324:356]/mutation_counts[324:356].sum(), color = 'm')
    plt.bar(codon_array[365:394], mutation_counts[365:394]/mutation_counts[365:394].sum(), color = 'y')
    plt.title("Codon Position vs. Frequency of Mutation in " + t + " Germline Mutations of TP53")
    plt.xlabel("Codon Position")
    plt.ylabel("Frequency")
    plt.show()
    plt.close()

    #per domain graphs
    fig, axs = plt.subplots(2,3)
    axs[0,0].bar(codon_array[1:41], mutation_counts[1:41]/mutation_counts[1:41].sum(), color = 'b')
    axs[0,0].set_title("TAD1")
    axs[0,0].set_xlabel("Codon Position")
    axs[0,0].set_ylabel("Mutation Frequency (%)")

    axs[0,1].bar(codon_array[41:62], mutation_counts[41:62]/mutation_counts[41:62].sum(), color = 'r')
    axs[0,1].set_title("TAD2")
    axs[0,1].set_xlabel("Codon Position")
    axs[0,1].set_ylabel("Mutation Frequency (%)")

    axs[0,2].bar(codon_array[65:93], mutation_counts[65:93]/mutation_counts[65:93].sum(), color = 'g')
    axs[0,2].set_title("PR")
    axs[0,2].set_xlabel("Codon Position")
    axs[0,2].set_ylabel("Mutation Frequency (%)")

    axs[1,0].bar(codon_array[101:301], mutation_counts[101:301]/mutation_counts[101:301].sum(), color = 'c')
    axs[1,0].set_title("DBD")
    axs[1,0].set_xlabel("Codon Position")
    axs[1,0].set_ylabel("Mutation Frequency (%)")

    axs[1,1].bar(codon_array[323:355], mutation_counts[323:355]/mutation_counts[323:355].sum(), color = 'm')
    axs[1,1].set_title("OD")
    axs[1,1].set_xlabel("Codon Position")
    axs[1,1].set_ylabel("Mutation Frequency (%)")

    axs[1,2].bar(codon_array[364:393], mutation_counts[364:393]/mutation_counts[364:393].sum(),
color = 'y')
    axs[1,2].set_title("C-terminus")
    axs[1,2].set_xlabel("Codon Position")
    axs[1,2].set_ylabel("Mutation Frequency (%)")

    fig.tight_layout(rect = [0,0.03, 1, 0.95])
    fig.suptitle("Frequency of " + t + " Germline Mutations in TP53")
    plt.show()
    plt.close()
    d = codon_array[1:41]
    f = mutation_counts[1:41]
    print(d)
    print(f)
    print("Max of TAD1", d[f == f.max()])

#deletions
deletions = df[df["Type"] == "del"]
descriptions = deletions["Description"]
description_nums = descriptions.str.split("del", expand = True)
data = description_nums[1].dropna()
data = data[~data.str.contains("|".join(["exon","\?"]))].to_numpy().astype("int")
#print(data)

data = np.array(data)
mean = np.mean(data)
median = np.median(data)
hist, bins, _ = plt.hist(data)
plt.figure()
logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]))
plt.hist(data, bins=logbins, ec = 'k')
plt.title("Base Pair Deletions in Germline Mutations of TP53")
plt.xlabel("Number of Base Pairs")
plt.annotate(mean, (mean, 150))
plt.annotate(median, (median,150))
plt.xscale("log")
plt.yscale("log")
plt.show()
plt.close()

#insertions
insertions  = df[df["Type"] == "ins"]
print(insertions)
print()
descriptions = insertions["Description"].to_list()
print(descriptions)
print()

data = []
for mut in descriptions:
    if ("dup Promoter" not in mut) and ("exon" not in mut) and ("del" not in mut):
        mut = int(mut.split("ins")[1])
        print(mut)
        data.append(mut)

plt.figure()
plt.hist(data, bins = 10, ec = 'k')
plt.title("Base Pair Insertions in Germline Mutations of TP53")
plt.xlabel("Length of Insertion")
plt.ylabel("Frequency")
plt.yscale("log")
plt.show()

#frameshift
frameshifts = df[df["Effect"] == "FS"]
c_description = frameshifts["c_description"].to_list()
data = []
for mut in c_description:
    if ("?" not in mut) and ("exon" not in mut) and ("NNNN" not in mut) and ("AGACC" not in mut)and ("+" not in mut):
        mut = mut.split(".")[1].split("del")[0].split("ins")[0].split("_")
        if len(mut) == 2:
            start, end = mut
            length_of_fs = int(end) - int(start) +1
        else:
            length_of_fs = 1
        data.append(length_of_fs)

plt.figure()
plt.hist(data, bins = 20, ec = 'k')
plt.title("Length of Frameshifts in Somatic Mutations of the TP53 Gene")
plt.xlabel("Length of Frameshift")
plt.ylabel("Frequency")
plt.yscale("log")
plt.show()

