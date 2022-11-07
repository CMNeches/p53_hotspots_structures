#!/usr/bin/env python3

# import packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import kstest
import argparse

# input arguments
parser = argparse.ArgumentParser()
parser.add_argument("short_topo", help = "the column in IARC with tissue types in all capital letters (case sensitive).  Write the full list for classes including multiple tissue types")
parser.add_argument("mutation_type", help = 'the "Effect" column in IARC')
parser.add_argument("bio_sex", help = "M, F, or all")
parser.add_argument("output_dir", help = "Output directory")
args = parser.parse_args()
short_topo = args.short_topo
mutation_type = args.mutation_type
bio_sex = args.bio_sex
output_dir = args.output_dir

# load in the IARC data
df = pd.read_csv('../../data/TVDr20.tsv', sep = '\t')

# filtering by biological sex
if bio_sex != "all":
    df = df[df["Sex"] == bio_sex]
else:
    pass

# filtering by tissue type
if short_topo != "all":
    df = df[df["Short_topo"] == short_topo]
else: 
    pass

#filtering by mutation type
if mutation_type != "all":
    df = df[df["Effect"] == mutation_type]
else:
    pass

# define the codon array (positions on the gene, the numbers 1 through 393)
codon_array = np.arange(1,394)

#mutated codons in the dataset
all_codons = df["Codon_number"].to_numpy()

# mutation counts per codon in dataset
mutation_counts = np.array([(all_codons == c).sum() for c in codon_array])

#keep track of amino acid mutations
aa_muts = df["ProtDescription"].to_list()

temporary_df = df[df["Effect"] == mutation_type]
all_codons = temporary_df["Codon_number"].to_numpy()

domain_dict = {"TAD1": [1, 40],
                "TAD2": [41, 61],
                "PR": [64, 92],
                "DBD": [100, 300],
                "OD": [323, 355],
                "C-terminus": [364, 393]}

                
# discovering novel hotspot locations in different domains
domain_list, num_domain_mutations, cross_domain_mutations = [], [], [] 
for domain, interval in domain_dict.items():
    start, end = interval

    # start - 1 due to Python counting, whcih starts at 0, end is not inclusive in Python
    new_codon_array = codon_array[start -1 : end]
    new_mutation_counts = mutation_counts[start -1 : end]

    #defining hotspots as mutations within 90th frequency percentile
    threshold = np.percentile(new_mutation_counts, 90)
    hotspots = new_codon_array[new_mutation_counts >= threshold]
    freqs = np.array([new_mutation_counts[new_codon_array == h].sum() for h in hotspots])/new_mutation_counts.sum()

    # 95% confidence intervals
    ci_offset = 1.96*np.std(freqs)/np.sqrt(new_mutation_counts.sum()) 

    # sort the indices
    sort_indices = (-freqs).argsort()
    hotspots = hotspots[sort_indices]
    freqs = freqs[sort_indices]

    # save the number of mutations within the domain
    domain_list.append(domain)
    num_domain_mutations.append(new_mutation_counts.sum())
    cross_domain_mutations.append(new_mutation_counts.sum()/(end - start + 1))

    # write to file
    with open(f"{output_dir}/{domain}_hotspot_statistics.txt", "w") as f:
        title = "\t".join(["Mutation", "Domain Frequency +/- 95% Confidence Interval"])
        to_write = "\n".join(["{}\t{}+/-{}".format(h,f,ci_offset) for h,f in zip(hotspots,freqs)])
        f.write(title + "\n" + to_write + "\n")

# write the domain statistics to file
with open(f"{output_dir}/domain_statistics.txt", "w") as f:
    title = "\t".join(["Domain", "Number of Mutations in Domain", "Domain Cross Section"])
    to_write = "\n".join(["{}\t{}\t{}".format(d, m, c) for d,m,c in zip(domain_list, num_domain_mutations,cross_domain_mutations)])
    f.write(title + "\n" + to_write + "\n")

#write codon mutation to file
with open(f"{output_dir}/codon_mutations.txt", "w") as f:
    title = "Codon\tCodon Counts"
    to_write = "\n".join(["{}\t{}".format(m,c) for m,c in zip(codon_array, mutation_counts)])
    f.write(title + "\n" + to_write + "\n")

#write amino acid mutations to file
unique_aa_muts = np.array(list(set(aa_muts)))
unique_aa_muts_counts = np.array([aa_muts.count(m) for m in unique_aa_muts])
sort_indices = (-unique_aa_muts_counts).argsort()
unique_aa_muts = unique_aa_muts[sort_indices]
unique_aa_muts_counts = unique_aa_muts_counts[sort_indices]

with open(f"{output_dir}/amino_acid_mutations.txt", "w") as f:
    title = "Amino Acid Mutation\tMutation Counts"
    to_write = "\n".join(["{}\t{}".format(m,c) for m,c in zip(unique_aa_muts, unique_aa_muts_counts)])
    f.write(title + "\n" + to_write + "\n")
