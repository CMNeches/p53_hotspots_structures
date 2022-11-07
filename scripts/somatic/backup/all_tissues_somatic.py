#!/usr/bin/env python3

#import packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import kstest
import argparse

#input arguments
parser = argparse.ArgumentParser()
parser.add_argument("topo_code", help = "the column in IARC with C and a number")
parser.add_argument("mutation_type", help = 'the "effects" column in IARC')
args = parser.parse_args()
topo_code = args.topo_code
mutation_type = args.mutation_type

#load in the data
iarc_data = pd.read_csv('../../data/TVDr20.csv')

tissues = ['C00', 'C00-C14', 'C01', 'C02', 'C03', 'C04', 'C05', 'C06', 'C07', 'C08', 'C09', 'C10', 'C11', 'C12', 'C13', 'C14', 'C15', 'C15-C26', 'C16', 'C17', 'C18', 'C18-C20', 'C19', 'C20', 'C21', 'C22', 'C23', 'C24', 'C25', 'C26', 'C30', 'C30-C39', 'C31', 'C32', 'C33', 'C34', 'C37', 'C38', 'C39', 'C40', 'C41', 'C42', 'C44', 'C47', 'C48', 'C49', 'C50', 'C51', 'C51-58', 'C52', 'C53', 'C54', 'C55', 'C56', 'C57', 'C58', 'C60', 'C60-C63', 'C61', 'C62', 'C63', 'C64', 'C64-C68', 'C65', 'C66', 'C67', 'C68', 'C69', 'C70', 'C71', 'C72', 'C73', 'C73-C75', 'C74', 'C75', 'C76', 'C77', 'C80']

for topo_code in tissues:
    df = iarc_data[iarc_data["Topo_code"] == topo_code]

    #pie chart
    piedf = df["Effect"].value_counts()
    plt.figure()
    plt.title("Types of Somatic Mutations " + topo_code)
    plt.pie(piedf)
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

    #here is how you loop through the mutation types
    types = ["missense", "FS", "nonsense", "silent", "splice", "other", "intronic", "large del"]
    for t in types:
        temporary_df = df[df["Effect"] == t]
        codon_array = np.arange(1,394)
        all_codons = temporary_df["Codon_number"].to_numpy()
        mutation_counts = np.array([(all_codons == c).sum() for c in codon_array])

        #mutations in each region
        plt.figure()
        plt.title("Codon Position vs. Mutation Counts " + t + " Somatic Mutations of TP53 in " + topo_code)
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

        #normalized frequency
        plt.figure()
        plt.title("Codon Position vs. Frequency in " + t + " Somatic Mutations in TP53 in " + topo_code)
        plt.bar(codon_array[1:41], mutation_counts[1:41]/mutation_counts[1:41].sum(), color = 'b')
        plt.bar(codon_array[41:62], mutation_counts[41:62]/mutation_counts[41:62].sum(), color = 'r')
        plt.bar(codon_array[65:93], mutation_counts[65:93]/mutation_counts[65:93].sum(), color = 'g')
        plt.bar(codon_array[101:301], mutation_counts[101:301]/mutation_counts[101:301].sum(), color = 'c')
        plt.bar(codon_array[324:356], mutation_counts[324:356]/mutation_counts[324:356].sum(), color = 'm')
        plt.bar(codon_array[365:394], mutation_counts[365:394]/mutation_counts[365:394].sum(), color = 'y')
        plt.xlabel("Codon Position")
        plt.ylabel("Frequency")
        plt.show()

        #per domain graphs normalized by region
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
        plt.suptitle("Frequency of " + t + " Somatic Mutations in TP53 in " + topo_code)
        plt.show()

        #statistics by type
        print(kstest(codon_array, 'uniform'))

    #deletions
    deletions = df[df["Type"] == "del"]
    description = deletions["Description"]
    if description.size == 0:
        pass
    else:
        description_nums= description.str.split("del", expand = True)
        data = description_nums[1].dropna()#.to_list()
        data = data[~data.str.contains("|".join(["exon", "\?", "\)","\(","-"]))].to_numpy().astype("int")

        data = np.array(data)
        mean = np.mean(data)
        median = np.median(data)
    
        hist, bins, _ = plt.hist(data)
        plt.figure()
        logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]))
        plt.hist(data, bins=logbins, ec = 'k')
        plt.title("Somatic Deletions in TP53 in" + topo_code)
        plt.xlabel("Number of Base Pairs")
        plt.annotate(mean, (mean, 500))
        plt.annotate(median, (median, 500))
        plt.xscale("log")
        plt.yscale("log")
        plt.show()

    #insertions
    insertions = df[df["Type"] == "ins"]
    description = insertions["Description"]
    if description.size == 0:
        pass
    else:
        description_nums = description.str.split("ins", expand = True)
        data = description_nums[1].replace('', np.nan, inplace = True)
        data = description_nums[1].dropna()#.to_list()
        data = data[~data.str.contains("|".join(["InFrame", "\?"]))].to_numpy().astype("int")
        data = np.array(data)
        mean = np.mean(data)
        median = np.median(data)

        hist, bins, _ = plt.hist(data)
        plt.figure()
        logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]))
        plt.hist(data, bins = logbins, ec = 'k')
        plt.title("Length of Insertions in TP53 Mutations in " + topo_code)
        plt.xlabel("Length of Insertion")
        plt.annotate(mean, (mean, 500))
        plt.annotate(median, (median, 500))
        plt.xscale("log")
        plt.yscale("log")
        plt.show()

    #frameshifts
    frameshifts = df[df["Effect"] == "FS"]
    c_description = frameshifts["c_description"]#.to_list()
    if c_description.size == 0:
        pass
    else:
        c_description = c_description.to_list()
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
                plt.title("Length of Frameshifts in Somatic Mutations of the TP53 Gene in " + topo_code)
                plt.xlabel("Length of Frameshift")
                plt.ylabel("Frequency")
                plt.yscale("log")
                plt.show()
