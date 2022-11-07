#!/usr/bin/env python3

#import packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#import the data
df = pd.read_csv('../data/TVDr20.csv')

#make a pie chart
piedf = df["Effect"].value_counts()
print(piedf)
labels = ["Missense", "FS", "Nonsense", "Silent", "Splice", "Other", "Intronic", "Large Del"]
#plt.pie(piedata, labels = labels)
#plt.show()
sizes = [21781, 2682, 2417, 1238, 730, 658, 217, 49]
plt.figure()
plt.pie(sizes, labels=labels, autopct = '%1.1f%%', startangle=90)
plt.show()
