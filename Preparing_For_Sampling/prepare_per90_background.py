import pandas as pd
import numpy as np
from collections import Counter
from numpy.random import choice
from scipy.stats import mannwhitneyu as mwu
from scipy.stats import ttest_ind
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import fdrcorrection
import sys
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

#Whether we are controlling for the predictions in the background
back_cont_var = int(sys.argv[7])

#Percentile cutoff to consider an element "accessible"
percent_cut = float(sys.argv[8])

#Must be NC, 3UTR, 5UTR, or Mis
variant_cat = sys.argv[9]

#Folder to write the matrices out to
out_folder = sys.argv[10]

var_pred = sys.argv[11]

assert(variant_cat in ["NC", "Mis", "3UTR", "5UTR"])

if out_folder.split("/")[-1] not in os.listdir("/".join(out_folder.split("/")[0:-1])):
    os.mkdir(out_folder)

#Read in the sites
v = pd.read_csv(file, sep = "\t").drop_duplicates("Position")

#Remove some extraneous information
v = v.drop(["Hum|Chp|Gor", "PhastCons447", "NearestDist"], axis = 1).copy()

print(v.shape)

#Replace mutation information that is WW or SS with WW_SS
v["MutCat"] = v["MutCat"].replace("WW", "WW_SS")
v["MutCat"] = v["MutCat"].replace("SS", "WW_SS")


#Do some preliminary filtering
v = v[v["PhyloP447"] != "."]
v["PhyloP447"] = v["PhyloP447"].astype(float)
v = v[v["SpecSup447"] != "."]
v["SpecSup447"] = v["SpecSup447"].astype(float)
v = v[(v["Derived"] == "H") | (v["Derived"] == "C")]
v = v[v["Category"] != "."]
v = v[v["KeptAfterFilt"] != "."].copy()
v.index = v["Position"]
v = v[v["NearestGene"] != "."]

print(v.shape)

#If it is a deep learning model, then we need to read in the deep learning predictions
#We will use the chimp-referenced for chimp-derived and human-referenced for human-derived
#We only did predictions for our 3 WGS filtered, so that parameter is irrelevant
dl_path = "/oak/stanford/groups/hbfraser/astarr/ForMikeChromBPNet/Variants_Grouped_AccelEvol/"

cpred = pd.read_csv(dl_path + var_pred + "_ChimpDerived_Chpreffed.txt", sep = "\t")
cpred = cpred[cpred["variant_id"] != "variant_id"]
cpred["allele1_pred_counts"] = cpred["allele1_pred_counts"].astype(float)
cpred["allele2_pred_counts"] = cpred["allele2_pred_counts"].astype(float)
cpred["logfc"] = cpred["logfc"].astype(float)
cpred["jsd"] = cpred["jsd"].astype(float)

hpred = pd.read_csv(dl_path + var_pred + "_HumanDerived_Humreffed.txt", sep = "\t")
hpred = hpred[hpred["variant_id"] != "variant_id"]
hpred["allele1_pred_counts"] = hpred["allele1_pred_counts"].astype(float)
hpred["allele2_pred_counts"] = hpred["allele2_pred_counts"].astype(float)
hpred["logfc"] = hpred["logfc"].astype(float)
hpred["jsd"] = hpred["jsd"].astype(float)

dl_pred = pd.concat([cpred, hpred]).set_index("variant_id")

if percent_cut != 0:
    cpred_to_cut = pd.DataFrame(np.max(cpred[["allele1_pred_counts", "allele2_pred_counts"]].astype(float), axis = 1)).sort_values(0)
    hpred_to_cut = pd.DataFrame(np.max(hpred[["allele1_pred_counts", "allele2_pred_counts"]].astype(float), axis = 1)).sort_values(0)
    cutoff_cpred = list(cpred_to_cut[0])[floor(percent_cut*hpred_to_cut.shape[0]//100)]
    cutoff_hpred = list(hpred_to_cut[0])[floor(percent_cut*cpred_to_cut.shape[0]//100)]
    
    cut = (cutoff_cpred + cutoff_hpred)/2

cpred = 0
hpred = 0
v = v.join(dl_pred).dropna()
dl_pred = 0

#Now filter based on our arguments

#Filter on what category of variant
v = v[v["Category"].isin([variant_cat])]

#Filter out those filtered out by 3 WGS filtering if desired
if filt_wgs:
    v = v[v["KeptAfterFilt"] == "Y"].copy()
#Filter on spec sup
v = v[v["SpecSup447"] > spec_sup]

v["abs logfc"] = np.abs(v["logfc"])

print(v.shape)
#Filter based on percentile accessibility
if percent_cut != 0:
    v = v[(v["allele1_pred_counts"] > cut) | (v["allele2_pred_counts"] > cut)]

#Done filtering

#Remove any with N
alphabet = ["A", "T", "C", "G"]
toss = ["NNN"]
for i in alphabet:
    toss.append(i + "NN")
    toss.append("N" + i + "N")
    toss.append("NN" + i)
    for j in alphabet:
        toss.append(i + j + "N")
        toss.append(i + "N" + j)
        toss.append("N" + i + j)
        toss.append(j + i + "N")
        toss.append(j + "N" + i)
        toss.append("N" + j + i)
toss = list(set(toss))

#Write out the 1000 base pair region flanking each qualifying site so that we can construct background matrices for only the sites in the 90th percentile of accessibility
v = v[~v["AncTrinuc"].isin(toss)].copy()
v = v[~v["DerTrinuc"].isin(toss)].copy()
print(v)
v["Chrom"] = [x.split(":")[0] for x in v["Position"]]
v["Pos1"] = [int(x.split(":")[1]) - 501 for x in v["Position"]]
v["Pos2"] = [int(x.split(":")[1]) + 500 for x in v["Position"]]
print(v)
v[["Chrom", "Pos1", "Pos2"]].to_csv(out_folder + "/Per90_Regions.bed", sep = "\t", header = False, index = False)
