#!/usr/bin/env python3

# import packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# gathering data
prefix = "../../results/hotspots/somatic/"
consistency_count = []
tissue_df = []
all_tissues = []
for tissue in os.listdir(prefix):
    data_file = f"{prefix}/{tissue}/missense/all/codon_mutations.txt"
    if os.path.isfile(data_file):
        df = pd.read_csv(data_file, sep="\t")
        df["Tissue"] = tissue
        df["Notes"] = "Passes 90th percentile threshold"
        codons = df["Codon"].to_numpy()
        counts = df["Codon Counts"].to_numpy()
        threshold = np.percentile(counts[counts>0], 90)
        if threshold > 0:
            all_tissues.append(tissue)
            tissue_df.append(df[df["Codon Counts"] >= threshold])
            hotspots = codons[counts >= threshold]
            consistency_count += list(hotspots)
    else:
        pass

tissue_df = pd.concat(tissue_df)
tissue_df.columns = ["Codon", "Codon Counts", "Tissues", "Notes"]

# plot hotspot consistency 
residues = np.arange(1, 394)
y_data = []
for res in residues:
    y_data.append(consistency_count.count(res))

plt.figure()
plt.bar(residues, y_data)
plt.xlabel("Residue")
plt.ylabel("Hotspot Consistency Count")
plt.savefig("../../results/hotspots/somatic/tissue_hotspot_consistency.png", bbox_inches="tight")

print(tissue_df)
tissue_df.to_csv("../../results/hotspots/somatic/tissue_hotspots.tsv", sep="\t", index=False)
