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


#If it is a deep learning model, then we need to read in the deep learning predictions
#We will use the chimp-referenced for chimp-derived and human-referenced for human-derived
#We only did predictions for our 3 WGS filtered, so that parameter is irrelevant
dl_path = "/oak/stanford/groups/hbfraser/astarr/ForMikeChromBPNet/Variants_Grouped_AccelEvol/"

if var_pred != "PhyloP447":
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

print(v)

#Now filter based on our arguments

#Filter on what category of variant
v = v[v["Category"].isin([variant_cat])]

#Filter out those filtered out by 3 WGS filtering if desired
if filt_wgs:
    v = v[v["KeptAfterFilt"] == "Y"].copy()
    print(v.shape)
#Filter on spec sup
v = v[v["SpecSup447"] > spec_sup]
print(v.shape)
#Filter on metric cutoff
if metric == "PhyloP447" or p_cut == 1:
    print(v.shape[0])
    v = v[v["PhyloP447"] >= p_cut].copy()
    print(v.shape[0])
elif metric == "abs logfc" or metric == "logfc":
    v = v[np.abs(v["logfc"]) > p_cut].copy()

if metric == "abs logfc":
    v["abs logfc"] = np.abs(v["logfc"])

print(v.shape)
#Filter based on percentile accessibility
if percent_cut != 0:
    v = v[(v["allele1_pred_counts"] > cut) | (v["allele2_pred_counts"] > cut)]

print(v.shape)
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

v = v[~v["AncTrinuc"].isin(toss)].copy()
v = v[~v["DerTrinuc"].isin(toss)].copy()

d_comp = {"A":"T", "C":"G", "G":"C", "T":"A"}
def revcomp(s):
    new_s = ""
    for i in s[::-1]:
        new_s = new_s + d_comp[i]
    return new_s

#Make trinuc equivalence classes
tris = list(set(v["AncTrinuc"]))
print(tris)
tris.sort()

d_equiv = {}
for i in tris:
    if i in d_equiv.values():
        d_equiv[i] = i
    else:
        d_equiv[i] = revcomp(i)

print(d_equiv)

def equiv(s):
    return d_equiv[s]
    
for i in d_equiv.keys():
    assert(d_equiv[i] == i or d_equiv[i] == revcomp(i))

v["AncTrinuc"] = v["AncTrinuc"].apply(equiv)
v["DerTrinuc"] = v["DerTrinuc"].apply(equiv)

#If we want to control for the type of mutation, append that to the trinucleotide context
#Increases baseline number of bins from 64 to 192
#This was not used
if cont_mut:
    v["AncTrinuc"] = v["AncTrinuc"] + "-" + v["MutCat"]

tris = list(set(v["AncTrinuc"]))
tris.sort()

v = v.drop_duplicates("Position")

print(v)
#If we aren't controlling for variant effect prediction effect size, then we can just output the proportion in each trinuc that is human vs chimp
if not cont_var:
    out = []
    for i in tris:
        human_tri = v[(v["AncTrinuc"].isin([i])) & (v["Derived"] == "H")].shape[0]
        chimp_tri = v[(v["AncTrinuc"].isin([i])) & (v["Derived"] == "C")].shape[0]
        out.append(human_tri/(human_tri + chimp_tri))
    df = pd.DataFrame(out)
    df.index = tris
    df.columns = ["All_Scores"]
    df.to_csv(out_folder + "/Species_TriProbs.txt", sep = "\t")

#If we are controlling for the metric, then we need to make our bins
else:
    hist_bins = np.array([np.float64(x) for x in np.histogram(v[metric], bins = 100)[1]])
    
    print(list(hist_bins))
    
    #Bin the sites by the metric
    v["Digitized"] = np.digitize(v[metric], hist_bins) - 1
    
    out = []
    for i in tris:
        human_tri = v[(v["AncTrinuc"].isin([i])) & (v["Derived"] == "H")]
        chimp_tri = v[(v["AncTrinuc"].isin([i])) & (v["Derived"] == "C")]
        #print(human_tri)
        human_hist = np.histogram(human_tri[metric], bins = hist_bins)
        #print(human_hist)
        #Make a histogram with the same bin edges for chimp
        chimp_hist = np.histogram(chimp_tri[metric], bins = hist_bins)
        #Compute probability in each bin of the histogram
        #When there are zero total mutations, it will output zero
        #But that is fine as that will never be encountered later on
        species_prob_hist = human_hist[0]/(np.maximum(human_hist[0] + chimp_hist[0], 1))
        out.append(species_prob_hist)
    df = pd.DataFrame(out)
    df.index = tris
    df.columns = hist_bins[0:len(hist_bins)-1]
    df.to_csv(out_folder + "/Species_TriProbs.txt", sep = "\t")
