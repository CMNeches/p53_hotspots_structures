#!/usr/bin/env python3

#import packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#load in the data
df = pd.read_csv('../data/TVDr20.csv')
df = df[df["Topo_code"] == "C18"]

#pie chart
piedf = df["Effect"].value_counts()
labels = ["Missense", "Nonsense", "FS", "Silent", "Other", "Splice", "Intronic", "Large Del"]
sizes = [878, 107, 84, 30, 26, 18, 11, 1]
plt.figure()
plt.pie(sizes, labels = labels, autopct= '%1.1f%%', startangle=90)
plt.show()


#we are defining the codon array (Python starts counting at 0)
codon_array = np.arange(1,394)

#df is the dataset that we have (ProtDescription is a column, we have already filtered for missense), and we are converting that column (or series) into a numpy array.  This shows the slow way
#all_mutations = df["ProtDescription"].to_numpy()
 
#here is the fast way
all_codons = df["Codon_number"].to_numpy()
print(all_codons)

#this is called a list comprehension (everything within the square brackets), which is a one line for loop.  For example, when c=1, the code is all_codons ==1, which means that we are looking for all of the places within all of the codons where the value is equal to one.  .sum() is the sum of all of the comparisons
mutation_counts = np.array([(all_codons == c).sum() for c in codon_array])

#here it is step by step
#print()
#print([c for c in codon_array])
#print()
#print(all_codons == 1)
#print()
#print(all_codons == 143)
#print()
#print(False == 0)
#print()
#print((all_codons == 42).sum())

#here is how you loop through the mutation types
types = ["missense", "FS", "nonsense", "silent", "splice", "other", "intronic", "large del"]
for t in types:
    temporary_df = df[df["Effect"] == t]
    codon_array = np.arange(1,394)
    all_codons = temporary_df["Codon_number"].to_numpy()
    mutation_counts = np.array([(all_codons == c).sum() for c in codon_array])
    #fig, ax = plt.subplots()
    
    #plt.figtext(20.5, -0.1, "TAD1", transform = ax.transAxes)
    #ax.text(20.5, -0.1, "TAD1", transform = ax.transAxes)

    #plt.show()
    domain_dict = {"TAD1" : [1, 41, 'b'], 
                   "TAD2" : [41, 62, 'r'],
                   "PR" : [65, 93, 'g'],
                   "DBD" : [100, 300, 'c'],
                   "OD": [323, 355, 'm'],
                   "C-terminus" : [364, 393, 'y']}
    #task 2a from hw 
    plt.figure()
    plt.title("Codon Position vs. Mutation Counts " + t + " Somatic Mutations of TP53")
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

    #for k,v in domain_dict.items():
        #start_codon, end_codon, color = v
        #plt.bar(codon_array[start_codon:end_codon], mutation_counts[start_codon:end_codon], color = color)
        #plt.title("Codon Position vs. Mutation Counts " + t + " Somatic Mutations of TP53")
        #plt.xlabel("Codon Position")
        #plt.ylabel("Mutation Counts")
        #plt.yscale("log")
        #plt.show()

    #task 2b from hw
    plt.figure()
    plt.title("Codon Position vs. Frequency in " + t + " Somatic Mutations in TP53")
    plt.bar(codon_array[1:41], mutation_counts[1:41]/mutation_counts[1:41].sum(), color = 'b')
    plt.bar(codon_array[41:62], mutation_counts[41:62]/mutation_counts[41:62].sum(), color = 'r')
    plt.bar(codon_array[65:93], mutation_counts[65:93]/mutation_counts[65:93].sum(), color = 'g')
    plt.bar(codon_array[101:301], mutation_counts[101:301]/mutation_counts[101:301].sum(), color = 'c')
    plt.bar(codon_array[324:356], mutation_counts[324:356]/mutation_counts[324:356].sum(), color = 'm')
    plt.bar(codon_array[365:394], mutation_counts[365:394]/mutation_counts[365:394].sum(), color = 'y')
    plt.xlabel("Codon Position")
    plt.ylabel("Frequency")
    plt.show()
    plt.close()

    #for k, v in domain_dict.items():
        #start_codon, end_codon, color = v
        #plt.bar(codon_array[start_codon:end_codon], mutation_counts[start_codon:end_codon]/mutation_counts.sum(), color = color)
        #plt.title("Codon Position vs. Frequency of Mutations" + t + " Somatic Mutations of TP53")
        #plt.xlabel("Codon Position")
        #plt.ylabel("Frequency")
        #plt.yscale("log")
        #plt.show()

    #per domain graphs
    #how do I put titles on these guys?  Thank you
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

    axs[1,0].bar(codon_array[101:301], mutation_counts[100:300]/mutation_counts[101:301].sum(), color = 'c')
    axs[1,0].set_title("DBD")
    axs[1,0].set_xlabel("Codon Position")
    axs[1,0].set_ylabel("Mutation Frequency (%)")

    axs[1,1].bar(codon_array[323:355], mutation_counts[323:355]/mutation_counts[323:355].sum(), color = 'm')
    axs[1,1].set_title("OD")
    axs[1,1].set_xlabel("Codon Position")
    axs[1,1].set_ylabel("Mutation Frequency (%)")

    axs[1,2].bar(codon_array[364:393], mutation_counts[364:393]/mutation_counts[364:393].sum(), color = 'y')
    axs[1,2].set_title("C-terminus")
    axs[1,2].set_xlabel("Codon Position")
    axs[1,2].set_ylabel("Mutation Frequency (%)")

    fig.tight_layout(rect = [0, 0.03, 1, 0.95])
    plt.suptitle("Frequency of " + t + " Somatic Mutations in TP53")
    plt.show()
    plt.close()


#we are pausing on the labels below thing because it is a pain
#task 4 from hw
deletions = df[df["Type"] == "del"]
description = deletions["Description"]
description_nums= description.str.split("del", expand = True)
data = description_nums[1].dropna()#.to_list()
data = data[~data.str.contains("|".join(["exon", "\?", "\)","\(","-"]))].to_numpy().astype("int")

data = np.array(data)
mean = np.mean(data)
median = np.median(data)
print(mean)
print(median)
print()
hist, bins, _ = plt.hist(data)
plt.figure()
logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]))
plt.hist(data, bins=logbins, ec = 'k')
plt.title("Base Pair Deletions in TP53")
plt.xlabel("Number of Base Pairs")
plt.annotate(mean, (mean, 500))
plt.annotate(median, (median, 500))
plt.xscale("log")
plt.yscale("log")
plt.show()
plt.close()

#task 5 from hw
insertions = df[df["Type"] == "ins"]
description = insertions["Description"]
description_nums = description.str.split("ins", expand = True)
data = description_nums[1].replace('', np.nan, inplace = True)
data = description_nums[1].dropna()#.to_list()
data = data[~data.str.contains("|".join(["InFrame", "\?"]))].to_numpy().astype("int")
data = np.array(data)
mean = np.mean(data)
median = np.median(data)
print(mean)
print(median)
print()
hist, bins, _ = plt.hist(data)
plt.figure()
logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]))
plt.hist(data, bins = logbins, ec = 'k')
plt.title("Length of Insertions in TP53 Mutations")
plt.xlabel("Length of Insertion")
plt.annotate(mean, (mean, 500))
plt.annotate(median, (median, 500))
plt.xscale("log")
plt.yscale("log")
plt.show()
plt.close()

#task 6 from hw
frameshifts = df[df["Effect"] == "FS"]
print(frameshifts)
print()
c_description = frameshifts["c_description"].to_list()
print(c_description)
print()

data = []
for mut in c_description:
    if ("?" not in mut) and ("+" not in mut) and ("-" not in mut):
        mut = mut.split(".")[1].split("del")[0].split("ins")[0].split("_")
        if len(mut) == 2:
            start, end = mut
            length_of_fs = int(end) - int(start) + 1
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
