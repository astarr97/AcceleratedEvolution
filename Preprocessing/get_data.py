import sys
import pandas as pd
import numpy as np
import os
from math import floor

#Name of the .txt file (always the same)
file = sys.argv[1]

#Whether to filter out the sites that might be polymorphic (1 is yes, 0 is no)
filt_wgs = int(sys.argv[2])

#Variant effect prediction effect size cutoff
p_cut = float(sys.argv[3])

#Species support cutoff
spec_sup = float(sys.argv[4])

#Whether to control for the type of mutation (1 is yes, 0 is no)
cont_mut = int(sys.argv[5])

#Whether to bin by variant effect prediction when doing the assignment
cont_var = int(sys.argv[6])

#Percentile cutoff to consider an element "accessible"
percent_cut = float(sys.argv[7])

#Must be NC, 3UTR, 5UTR, or Mis
variant_cat = sys.argv[8]

#Folder to write the matrices out to
out_folder = sys.argv[9]

#Either the name of a deep learning prediction set or "PhyloP447"
var_pred = sys.argv[10]

#Variant effect prediction metric to use, can be "abs logfc" or "logfc" if it is deep learning, can only be "PhyloP447" otherwise
metric = sys.argv[11].replace("_", " ")

assert(metric in ["abs logfc", "logfc", "PhyloP447"])
assert(variant_cat in ["NC", "Mis", "3UTR", "5UTR"])

if out_folder.split("/")[-1] not in os.listdir("/".join(out_folder.split("/")[0:-1])):
    os.mkdir(out_folder)

#Read in the sites
v = pd.read_csv(file, sep = "\t")

#Remove some extraneous information
v = v.drop(["Hum|Chp|Gor", "PhastCons447", "NearestDist"], axis = 1).copy()

print(v.shape)

#Replace mutation information that is WW or SS with WW_SS
v["MutCat"] = v["MutCat"].replace("WW", "WW_SS")
v["MutCat"] = v["MutCat"].replace("SS", "WW_SS")


#Do some preliminary filtering
v = v[v["PhyloP447"] != "."]
v["PhyloP447"] = v["PhyloP447"].astype(float)

print(v.shape)

v = v[v["SpecSup447"] != "."]
v["SpecSup447"] = v["SpecSup447"].astype(float)

v = v[v["NearestGene"] != "."]

print(v.shape)

v = v[(v["Derived"] == "H") | (v["Derived"] == "C")]

print(v.shape)

v = v[v["Category"] != "."]
v = v[v["KeptAfterFilt"] != "."].copy()
v.index = v["Position"]

print(v.shape)

v.to_csv("HumChp_All_PhyloP.txt", sep = "\t", index = False)
