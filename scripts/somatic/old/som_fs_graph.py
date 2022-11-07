#!/usr/bin/env python3

#import packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#load in the data
df = pd.read_csv("../data/TVDr20.csv")

#get the frameshifts
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
