#!/usr/bin/env python3

#import packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#import data
df = pd.read_csv("../../data/TVDr20.csv")
df_germline = pd.read_csv("../../data/Germline.csv")

#make the insertions
somatic_insertions = df[df["Type"] == "ins"][""]
germline_insertions = df[df["Type"] == "ins"][""]

#set automatically uniquifies data
insertions = list(set.intersection(set(somatic_insertions), set(germline_insertions)))
print(insertions)

