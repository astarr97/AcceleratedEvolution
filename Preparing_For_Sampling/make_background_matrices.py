import sys
import pandas as pd
import numpy as np
import os
from math import floor
from collections import Counter

#E.g. python make_background_matrices.py Background_AccelEvolInput.5UTR.Final.sort.txt 5UTR 1 -100 0 0 PhyloP/BackgroundMatrix_PhyloP_5UTR_PhyloP-100_SpecSup0_NCM_Per0 HumChp_AccelEvolInput.Final.txt
#Name of the .txt file (not always the same)
#MUST BE SORTED BY GENE NAME
file = sys.argv[1]

#The type of variant
variant_cat = sys.argv[2]

#Whether to filter human-chimp variants by WGS
filt_wgs = int(sys.argv[3])

#ONLY APPLICABLE FOR PHYLOP
#Variant effect prediction effect size cutoff
p_cut = float(sys.argv[4])

#Species support cutoff
spec_sup = float(sys.argv[5])

#ONLY APPLICABLE FOR PHYLOP
#Whether to bin by variant effect prediction when doing the assignment
cont_var = int(sys.argv[6])

#Folder to write the matrices out to
out_folder = sys.argv[7]

#Either the name of a deep learning prediction set or "PhyloP447"
var_pred = "PhyloP447"

#File containing the human-chimp variants
var_file = sys.argv[8]

print(file)

o = open(file)

print(out_folder.split("/"))
if out_folder.split("/")[-1] not in os.listdir("/".join(out_folder.split("/")[0:-1])):
    os.mkdir(out_folder)
    
#Make trinuc equivalence classes
tris = []
for i in ["A", "C", "G", "T"]:
    for j in ["A", "C", "G", "T"]:
        for k in ["A", "C", "G", "T"]:
            tris.append(i + j + k)
print(tris)
tris.sort()

d_comp = {"A":"T", "C":"G", "G":"C", "T":"A"}
def revcomp(s):
    new_s = ""
    for i in s[::-1]:
        new_s = new_s + d_comp[i]
    return new_s

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

v = pd.read_csv(var_file, sep = "\t").drop_duplicates("Position")

#First, need to get the bins for each trinucleotide using the human-chimp variants
#This just reproduces the code in make make_species_matrices.py
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
v = v[v["PhyloP447"] >= p_cut].copy()

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

v["AncTrinuc"] = v["AncTrinuc"].apply(equiv)
v["DerTrinuc"] = v["DerTrinuc"].apply(equiv)

keep_genes = np.unique(v["NearestGene"])

v = v.drop_duplicates("Position").copy()

#Now we actually get the background matrices
#If we are not controlling for the variant effect prediction, then we can go through and make our genes x trinucleotides table

if not cont_var:
    order = list(set(d_equiv.values()))
    order.sort()
    prev_gene = 0
    c = Counter(list(set(d_equiv.values())))
    out = []

    for line in o:
        
        l = line.replace("\n", "").split("\t")
        #If it is a qualifying site according to our cutoffs
        if float(l[1]) > p_cut and float(l[2]) > spec_sup and "N" not in l[-1]:
            new_tri = d_equiv[l[-1]]
            if not prev_gene:
                prev_gene = l[3]
                c[new_tri] += 1
            elif prev_gene != l[3]:
                #add the accumulated information to our list
                to_app = []
                for key in order:
                    to_app.append(c[key] - 1)
                out.append([prev_gene] + to_app)
                c = Counter(list(set(d_equiv.values())))
                c[new_tri] += 1
                prev_gene = l[3]
            else:
                #Increment counter if it is the same gene
                c[new_tri] += 1
    out.append([prev_gene] + to_app)
    df_back = pd.DataFrame(out)
    print(df_back)
    df_back.columns = ["Gene"] + order
    df_back = df_back[df_back["Gene"].isin(keep_genes)]
    df_back = df_back.set_index("Gene")
    print(df_back)
    #Add back for any missing genes so it plays nicely with the other scripts
    for gene in list(np.setdiff1d(keep_genes, df_back.index)):
        df_back.loc[gene] = np.repeat(0, df_back.shape[1])
    print(df_back)
    for i in order:
        df_back[i] = df_back[i]/np.sum(df_back[i])
    #Write out our gene trinucleotide probabilities
    df_back.sort_index().to_csv(out_folder + "/Gene_TriProbs.txt", sep = "\t", header = True, index = True)
#THIS WAS NOT USED, PLEASE IGNORE
else:

    hist_bins = np.array([np.float64(x) for x in np.histogram(v[var_pred], bins = 100)[1]])
    
    #Read in the whole thing as a dataframe, but only use the PhyloP column
    #Then can get the histogram and go from there
    for key in d_equiv.keys():
        if d_equiv[key] != key:
            to_keep = [key, d_equiv[key]]
            df = pd.read_csv(file, sep = "\t", header = None, usecols = [0, 1, 2, 3, 5], chunksize = 10000000)
            new_df = pd.DataFrame()
            for chunk in df:
                chunk2 = chunk[chunk[5].isin(to_keep)]
                new_df = pd.concat([new_df, chunk2])
            df = 0
            df = new_df.copy().drop_duplicates(0).drop([0], axis = 1)
            df_back = df[(df[1] > p_cut) & (df[2] > spec_sup)].copy()
            df = 0

            df_back = df_back.sort_values(3)
            df_back[1] = df_back[1].astype(np.float64)
            
            df2 = pd.DataFrame(df_back.groupby(3).agg({1: lambda x: np.histogram(x, hist_bins)}))
            genes = list(df2.index)
            to_df2 = []
            [to_df2.append(x[0]) for x in list(df2[1])]
            df2 = pd.DataFrame(to_df2)
            df2["Gene"] = genes
            df2.index = df2["Gene"]
            df2 = df2.drop(["Gene"], axis = 1).copy()
            df2.columns = hist_bins[0:len(hist_bins)-1]
            
            #Adding this to restrict to only genes of interest and add back any missing genes
            df2 = df2.loc[np.intersect1d(df2.index, keep_genes)]
            for gene in list(np.setdiff1d(keep_genes, df2.index)):
                df2.loc[gene] = np.repeat(0, df2.shape[1])
            
            df2 = df2.sort_index()
            df2 = df2.divide(df2.sum(axis = 0), axis = 1)
            df2 = df2.replace(np.nan, np.float64(1/df2.shape[0]))
            df2.sort_index().to_csv(out_folder + "/" + d_equiv[key] + "_Gene_TriProbs.txt", sep = "\t", header = True, index = True)
            
o.close()
