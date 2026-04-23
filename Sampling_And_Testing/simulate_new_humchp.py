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

#Either the name of a deep learning prediction set or "PhyloP447"
var_pred = sys.argv[11]

#Variant effect prediction metric to use, can be "abs logfc" or "logfc" if it is deep learning, can only be "PhyloP447" otherwise
metric = sys.argv[12].replace("_", " ")

#Folder where the species probabilities are stored
spec_prob_folder = sys.argv[13]

#Folder where the background probabilities are stored
back_prob_folder = sys.argv[14]

#Permutation start and end, comma-delimited
pse = sys.argv[15]
start_iter = int(pse.split(",")[0])
end_iter = int(pse.split(",")[1])

#Whether we are doing this at the gene set level
#Input should be of the form file,maximum number of genes per category,minimum number of genes per category
try:
    gene_set_string = sys.argv[16]
    gene_set = gene_set_string.split(",")[0]
    min_genes = int(gene_set_string.split(",")[1])
    max_genes = int(gene_set_string.split(",")[2])
except:
    gene_set = 0
print(gene_set)
#Only write out the permutations if we want to
if variant_cat == "NC":
    write_perm = 1
else:
    write_perm = 0

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

#Now filter based on our arguments

#Filter on what category of variant
v = v[v["Category"].isin([variant_cat])]

#Filter out those filtered out by 3 WGS filtering if desired
if filt_wgs:
    v = v[v["KeptAfterFilt"] == "Y"].copy()
#Filter on spec sup
v = v[v["SpecSup447"] > spec_sup]

#Filter on metric cutoff
if metric == "PhyloP447" or p_cut == 1:
    v = v[v["PhyloP447"] >= p_cut].copy()
elif metric == "abs logfc" or metric == "logfc":
    v = v[np.abs(v["logfc"]) > p_cut].copy()

if metric == "abs logfc" or metric == "logfc":
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
tris.sort()

d_equiv = {}
for i in tris:
    if i in d_equiv.values():
        d_equiv[i] = i
    else:
        d_equiv[i] = revcomp(i)

def equiv(s):
    return d_equiv[s]
    
for i in d_equiv.keys():
    assert(d_equiv[i] == i or d_equiv[i] == revcomp(i))

v["AncTrinuc"] = v["AncTrinuc"].apply(equiv)
v["DerTrinuc"] = v["DerTrinuc"].apply(equiv)

#If we are controlling for the type of mutation, append it to the trinucleotide
#THIS WAS NOT USED
if cont_mut:
    v["AncTrinuc"] = v["AncTrinuc"] + "-" + v["MutCat"]

v = v.drop(["DerTrinuc", "SpecSup447", "Category", "KeptAfterFilt"], axis = 1).copy()
v = v.drop_duplicates("Position").copy()

#Create non-negative PhyloP if we are using PhyloP
if var_pred == "PhyloP447":
    v["PhyloP NonNeg"] = np.maximum(v["PhyloP447"], 0)

all_genes = np.unique(v["NearestGene"])

species1 = v[v["Derived"] == "H"].copy()
species2 = v[v["Derived"] == "C"].copy()

species1["Total_Vars"] = np.repeat(1, species1.shape[0])
species2["Total_Vars"] = np.repeat(1, species2.shape[0])

species1_sum = species1.groupby(["NearestGene"]).sum(numeric_only=1)
species2_sum = species2.groupby(["NearestGene"]).sum(numeric_only=1)
if metric == "PhyloP447":
    species1_sum.columns = ["Species1 Sum PhyloP447", "Species1 Sum PhyloP NonNeg", "Species1 Sum Total_Vars"]
    species2_sum.columns = ["Species2 Sum PhyloP447", "Species2 Sum PhyloP NonNeg", "Species2 Sum Total_Vars"]
elif metric == "abs logfc" or metric == "logfc" or metric == "abs_logfc":
    species1_sum.columns = ["Species1 Sum PhyloP447", "Species1 Sum allele1_pred_counts", "Species1 Sum allele2_pred_counts", "Species1 Sum logfc", "Species1 Sum jsd", "Species1 Sum abs logfc", "Species1 Sum Total_Vars"]
    species2_sum.columns = ["Species2 Sum PhyloP447", "Species2 Sum allele1_pred_counts", "Species2 Sum allele2_pred_counts", "Species2 Sum logfc", "Species2 Sum jsd", "Species2 Sum abs logfc", "Species2 Sum Total_Vars"]

new_sum = species1_sum.join(species2_sum, how = "outer").fillna(0)

toss_gene = []
for i in np.unique(new_sum.index):
    if i.startswith("PCDHGA") or i.startswith("PCDHGB") or i.startswith("PCDHGC") or i.startswith("MAPKBP"):
        toss_gene.append(i)

def sum_gs(genes, df_metric):
    new_sum_genes = df_metric.loc[np.intersect1d(df_metric.index, np.setdiff1d(genes.split(";"), toss_gene))]
    new_sum_genes = new_sum_genes.sum(axis = 0)
    return new_sum_genes

#Write out Mutation category count for checking that everything worked
result = v.groupby('NearestGene')['MutCat'].value_counts().reset_index(name='Count')

mut_cat_count = []
for i in np.unique(result["NearestGene"]):
    result2 = result[result["NearestGene"] == i]
    try:
        sw = str(int(result2[result2["MutCat"] == "SW"]["Count"].iloc[0]))
    except:
        sw = "NA"
    try:
        ws = str(int(result2[result2["MutCat"] == "WS"]["Count"].iloc[0]))
    except:
        ws = "NA"
    try:
        ww_ss = str(int(result2[result2["MutCat"] == "WW_SS"]["Count"].iloc[0]))
    except:
        ww_ss = "NA"
    mut_cat_count.append([i, ";".join([sw, ws, ww_ss])])
mut_cat_df = pd.DataFrame(mut_cat_count)
mut_cat_df.columns = ["Gene", "SW;WS;WW_SS"]
mut_cat_df.to_csv(out_folder + "/MutCatCount.txt", sep = "\t", index = False)

print(v)

#If it is a gene set, do it at the gene set level
if gene_set:
    #Create a dataframe for the gene set
    gs_df = pd.read_csv(gene_set, sep = "\t")
    keep_terms = []
    for index, row in gs_df.iterrows():
        genes_oi = np.intersect1d(row["Genes"].split(";"), new_sum.index)
        num_genes = len(genes_oi)
        if num_genes >= min_genes and num_genes <= max_genes:
            keep_terms.append([row["Term"], ";".join(list(genes_oi))])
    gs_df = pd.DataFrame(keep_terms)
    gs_df.columns = ["Term", "Genes"]
    
    #Sum things per gene set
    new_sum_gs = gs_df["Genes"].apply(sum_gs, args=(new_sum, ))
    new_sum_gs.index = gs_df["Term"]
    
    #Compute our statistics at the gene set level
    if metric == "PhyloP447":
        new_sum_gs["PhyloP Difference"] = new_sum_gs["Species1 Sum PhyloP NonNeg"] - new_sum_gs["Species2 Sum PhyloP NonNeg"]
        new_sum_gs["PhyloP L2FC"] = np.log2((new_sum_gs["Species1 Sum PhyloP NonNeg"] + 10)/(new_sum_gs["Species2 Sum PhyloP NonNeg"] + 10))
        new_sum_gs.to_csv(out_folder + "/AccelEvol_Init.txt", sep = "\t")
        new_sum_gs["Total_Sims"] = np.repeat(0, new_sum_gs.shape[0])
    elif metric == "abs logfc":
        new_sum_gs["Abs Logfc Difference"] = new_sum_gs["Species1 Sum abs logfc"] - new_sum_gs["Species2 Sum abs logfc"]
        new_sum_gs.to_csv(out_folder + "/AccelEvol_Init.txt", sep = "\t")
        new_sum_gs["Total_Sims"] = np.repeat(0, new_sum_gs.shape[0])
    elif metric == "logfc":
        new_sum_gs["Logfc Difference"] = new_sum_gs["Species1 Sum logfc"] - new_sum_gs["Species2 Sum logfc"]
        new_sum_gs.to_csv(out_folder + "/AccelEvol_Init.txt", sep = "\t")
        new_sum_gs["Total_Sims"] = np.repeat(0, new_sum_gs.shape[0])
    new_sum_gs["Better_Sims_L2FC"] = np.repeat(0, new_sum_gs.shape[0])
    new_sum_gs["Better_Sims_Sum"] = np.repeat(0, new_sum_gs.shape[0])
#Otherwise, just compute our statistics
else:
    if metric == "PhyloP447":
        new_sum["PhyloP Difference"] = new_sum["Species1 Sum PhyloP NonNeg"] - new_sum["Species2 Sum PhyloP NonNeg"]
        new_sum["PhyloP L2FC"] = np.log2((new_sum["Species1 Sum PhyloP NonNeg"] + 10)/(new_sum["Species2 Sum PhyloP NonNeg"] + 10))
        new_sum.to_csv(out_folder + "/AccelEvol_Init.txt", sep = "\t")
        new_sum["Total_Sims"] = np.repeat(0, new_sum.shape[0])
    elif metric == "abs logfc":
        new_sum["Abs Logfc Difference"] = new_sum["Species1 Sum abs logfc"] - new_sum["Species2 Sum abs logfc"]
        new_sum.to_csv(out_folder + "/AccelEvol_Init.txt", sep = "\t")
        new_sum["Total_Sims"] = np.repeat(0, new_sum.shape[0])
    elif metric == "logfc":
        new_sum["Logfc Difference"] = new_sum["Species1 Sum logfc"] - new_sum["Species2 Sum logfc"]
        new_sum.to_csv(out_folder + "/AccelEvol_Init.txt", sep = "\t")
        new_sum["Total_Sims"] = np.repeat(0, new_sum.shape[0])
    new_sum["Better_Sims_L2FC"] = np.repeat(0, new_sum.shape[0])
    new_sum["Better_Sims_Sum"] = np.repeat(0, new_sum.shape[0])

#Read in the background probabilities
if not back_cont_var:
    gene_probs = pd.read_csv(back_prob_folder + "/Gene_TriProbs.txt", sep = "\t").set_index("Gene")
#THIS ELSE WAS NOT USED
else:
    d_gene_probs = {}
    for file in os.listdir(back_prob_folder):
        v = pd.read_csv(back_prob_folder + "/" + file, sep = "\t").set_index("Gene")
        v.columns = list(range(v.shape[1]))
        d_gene_probs[file.split("_")[0]] = v.copy()
#Bin sites if we are controlling for PhyloP or abs logfc, this was used
if cont_var:    
    species_probs = pd.read_csv(spec_prob_folder + "/Species_TriProbs.txt", sep = "\t").set_index("Unnamed: 0")
    hist_bins = np.array([np.float64(jj) for jj in species_probs.columns])
    species_probs.columns = list(range(species_probs.shape[1]))
    #Add bin indices
    print(metric)
    if metric == "PhyloP447":
        species1["Bin index"] = np.digitize(species1["PhyloP447"], bins = hist_bins) - 1
        species2["Bin index"] = np.digitize(species2["PhyloP447"], bins = hist_bins) - 1
    elif metric == "abs logfc" or metric == "abs_logfc":
        species1["Bin index"] = np.digitize(species1["abs logfc"], bins = hist_bins) - 1
        species2["Bin index"] = np.digitize(species2["abs logfc"], bins = hist_bins) - 1
    elif metric == "logfc":
        species1["Bin index"] = np.digitize(species1["logfc"], bins = hist_bins) - 1
        species2["Bin index"] = np.digitize(species2["logfc"], bins = hist_bins) - 1
else:
    species_probs = pd.read_csv(spec_prob_folder + "/Species_TriProbs.txt", sep = "\t").set_index("Unnamed: 0")

#Filtering to make sure gene lists match exactly
if not back_cont_var:
    gene_probs = gene_probs.loc[np.intersect1d(new_sum.index, [str(jj) for jj in gene_probs.index])].copy()
    new_sum = new_sum.loc[np.intersect1d(new_sum.index, [str(jj) for jj in gene_probs.index])].copy()
    new_sum = new_sum.sort_index()
    gene_probs = gene_probs.sort_index()
    assert(list(new_sum.index) == list(gene_probs.index))
#This was not sued
else:
    for key in d_gene_probs.keys():
        gene_probs = d_gene_probs[key]
        gene_probs = gene_probs.loc[np.intersect1d(new_sum.index, [str(jj) for jj in gene_probs.index])].copy()
        new_sum = new_sum.loc[np.intersect1d(new_sum.index, [str(jj) for jj in gene_probs.index])].copy()
        new_sum = new_sum.sort_index()
    for key in d_gene_probs.keys():
        gene_probs = d_gene_probs[key]
        gene_probs = gene_probs.loc[np.intersect1d(new_sum.index, [str(jj) for jj in gene_probs.index])].copy()
        new_sum = new_sum.loc[np.intersect1d(new_sum.index, [str(jj) for jj in gene_probs.index])].copy()
        new_sum = new_sum.sort_index()
        gene_probs = gene_probs.sort_index()
        d_gene_probs[key] = gene_probs
        assert(list(new_sum.index) == list(d_gene_probs[key].index))
#Prepare for permutations
input_sim = pd.concat([species1, species2]).sort_values("AncTrinuc")

print("Input to simulation:")
print(input_sim)

print(gene_probs)
tris = np.unique(input_sim["AncTrinuc"])

#For each permutation
tot_sim = 0
c = 0
for j in range(start_iter, end_iter):
    if j % 1000 == 0:
        print(j)
    
    c += 1
    tot_sim += 1
    np.random.seed(j)
    #print(j)
    #Go through the sorted list of trinucleotides and assign each mutation to a new gene and new species
    new_gene = []
    new_species = []
    
    #Iterate through each trinucleotide, assigning substiutions to genes and a lineage
    for tri in tris:
        input_sim_tri = input_sim[input_sim["AncTrinuc"].isin([tri])].copy()
        
        #Newer version, much faster!
        if not cont_var and not back_cont_var:
            #Little bit of normalization just to be sure
            gene_prob_dist = np.array(gene_probs[tri])/np.sum(np.array(gene_probs[tri]))
            input_sim_tri["New_Gene"] = list(choice(np.array(gene_probs.index), p=gene_prob_dist, size = len(list(input_sim_tri.index))))

            s1prob = np.float64(species_probs.loc[tri].iloc[0])
            input_sim_tri["New_Species"] = list(choice(["S1Der", "S2Der"], p = [s1prob, 1-s1prob], size = len(list(input_sim_tri.index))))
        elif cont_var and not back_cont_var:
            #Sort things by bin_index
            input_sim_tri = input_sim_tri.sort_values("Bin index")
            bin_indices = np.unique(input_sim_tri["Bin index"])
            #Identical to just_tri for New_Gene
            gene_prob_dist = np.array(gene_probs[tri])/np.sum(np.array(gene_probs[tri]))
            input_sim_tri["New_Gene"] = list(choice(np.array(gene_probs.index), p=gene_prob_dist, size = len(list(input_sim_tri.index))))
            
            #Iterate through the bin_indices
            new_species = []
            prob_list = np.array(species_probs.loc[tri])
        
            for bin_index in bin_indices:
                #Get the size of the sample (number of species differences with the bin_index and trinucleotide)
                size_to_sample = len(list(input_sim_tri[input_sim_tri["Bin index"].isin([bin_index])].index))
                s1prob = np.float64(prob_list[bin_index])
                #print(s1prob, bin_index, size_to_sample)
                new_species = new_species + list(choice(["S1Der", "S2Der"], p = [s1prob, 1-s1prob], size = size_to_sample))
            from collections import Counter
            #print(Counter(new_species), tri)
            input_sim_tri["New_Species"] = new_species
        elif cont_var and back_cont_var:
            #Sort things by bin_index
            input_sim_tri = input_sim_tri.sort_values("Bin index")
            bin_indices = np.unique(input_sim_tri["Bin index"])

            #Iterate through the bin_indices and assign genes
            new_species = []
            new_gene = []
            prob_list = np.array(species_probs.loc[tri])
            gene_probs = d_gene_probs[tri]

            for bin_index in bin_indices:
                #Get the size of the sample (number of species differences with the bin_index and trinucleotide)
                size_to_sample = len(list(input_sim_tri[input_sim_tri["Bin index"].isin([bin_index])].index))
                s1prob = np.float64(prob_list[bin_index])
                new_species = new_species + list(choice(["S1Der", "S2Der"], p = [s1prob, 1-s1prob], size = size_to_sample))
                new_gene = new_gene + list(choice(np.array(gene_probs.index), p = np.array(gene_probs[bin_index]), size = size_to_sample))
            input_sim_tri["New_Gene"] = new_gene
            input_sim_tri["New_Species"] = new_species
        if tri == tris[0]:
            input_sim_new = input_sim_tri
        else:
            input_sim_new = pd.concat([input_sim_new, input_sim_tri])
    input_sim_species1 = input_sim_new[input_sim_new["New_Species"].isin(["S1Der"])]
    input_sim_species2 = input_sim_new[input_sim_new["New_Species"].isin(["S2Der"])]

    if cont_var or back_cont_var:
        input_sim_species1 = input_sim_species1.drop(["Bin index"], axis = 1).copy()
        input_sim_species2 = input_sim_species2.drop(["Bin index"], axis = 1).copy()
    
    #Do the same computation as done for the original
    species1_sum = input_sim_species1.groupby(["New_Gene"]).sum(numeric_only=1)
    species2_sum = input_sim_species2.groupby(["New_Gene"]).sum(numeric_only=1)
    #print(species1_sum)
    if metric == "PhyloP447":
        species1_sum.columns = ["Species1 Sum PhyloP447 Sim", "Species1 Sum PhyloP NonNeg Sim" , "Species1 Sum Total_Vars Sim"]
        species2_sum.columns = ["Species2 Sum PhyloP447 Sim", "Species2 Sum PhyloP NonNeg Sim" , "Species2 Sum Total_Vars Sim"]
    elif metric == "abs logfc" or metric == "logfc":
        species1_sum.columns = ["Species1 Sum PhyloP447 Sim", "Species1 Sum allele1_pred_counts Sim", "Species1 Sum allele2_pred_counts Sim", "Species1 Sum logfc Sim", "Species1 Sum jsd Sim", "Species1 Sum abs logfc Sim", "Species1 Sum Total_Vars Sim"]
        species2_sum.columns = ["Species2 Sum PhyloP447 Sim", "Species2 Sum allele1_pred_counts Sim", "Species2 Sum allele2_pred_counts Sim", "Species2 Sum logfc Sim", "Species2 Sum jsd Sim", "Species2 Sum abs logfc Sim", "Species2 Sum Total_Vars Sim"]

    new_sim = species1_sum.join(species2_sum, how = "outer").fillna(0)
    #print(new_sim)
    #if metric == "abs logfc" or metric == "logfc":
    #    new_sum = new_sum.drop(["Species1 Sum allele1_pred_counts", "Species1 Sum allele2_pred_counts", "Species2 Sum allele1_pred_counts", "Species2 Sum allele2_pred_counts", "Species1 Sum jsd", "Species2 Sum jsd"], axis = 1).copy()
    
    #If it is a gene set, compute our statistics for this sampling
    if gene_set:
        if not write_perm:
            new_sim_gs = gs_df["Genes"].apply(sum_gs, args=(new_sim, ))
            new_sim_gs.index = gs_df["Term"]
            new_sum_gs = new_sum_gs.join(new_sim_gs).fillna(0)
        else:
            new_sim_gs = gs_df["Genes"].apply(sum_gs, args=(new_sim, ))
            new_sim_gs.index = gs_df["Term"]
            
        #print(new_sum_gs)
        if metric == "PhyloP447":
            if not write_perm:
                new_sum_gs["PhyloP Difference Sim"] = new_sum_gs["Species1 Sum PhyloP NonNeg Sim"] - new_sum_gs["Species2 Sum PhyloP NonNeg Sim"]
                new_sum_gs["PhyloP L2FC Sim"] = np.log2((new_sum_gs["Species1 Sum PhyloP NonNeg Sim"] + 10)/(new_sum_gs["Species2 Sum PhyloP NonNeg Sim"] + 10))
            
                new_sum_gs["Better_Sims_L2FC"] = new_sum_gs["Better_Sims_L2FC"] + pd.DataFrame(((new_sum_gs["PhyloP L2FC"] <= new_sum_gs["PhyloP L2FC Sim"]) & (new_sum_gs["PhyloP L2FC"] >= 0)) | ((new_sum_gs["PhyloP L2FC"] >= new_sum_gs["PhyloP L2FC Sim"]) & (new_sum_gs["PhyloP L2FC"] < 0))).replace({True:1, False:0})[0]
                new_sum_gs["Better_Sims_Sum"] = new_sum_gs["Better_Sims_Sum"] + pd.DataFrame(((new_sum_gs["PhyloP Difference"] <= new_sum_gs["PhyloP Difference Sim"]) & (new_sum_gs["PhyloP Difference"] >= 0)) | ((new_sum_gs["PhyloP Difference"] >= new_sum_gs["PhyloP Difference Sim"]) & (new_sum_gs["PhyloP Difference"] < 0))).replace({True:1, False:0})[0]
            else:
                new_sim_gs["PhyloP Difference Sim"] = new_sim_gs["Species1 Sum PhyloP NonNeg Sim"] - new_sim_gs["Species2 Sum PhyloP NonNeg Sim"]
                new_sim_gs["PhyloP L2FC Sim"] = np.log2((new_sim_gs["Species1 Sum PhyloP NonNeg Sim"] + 10)/(new_sim_gs["Species2 Sum PhyloP NonNeg Sim"] + 10))

        elif metric == "abs logfc":
            print(new_sum_gs.columns)
            new_sim_gs["Abs Logfc Difference Sim"] = new_sim_gs["Species1 Sum abs logfc Sim"] - new_sim_gs["Species2 Sum abs logfc Sim"]

        elif metric == "logfc":
            new_sim_gs["Logfc Difference Sim"] = new_sim_gs["Species1 Sum logfc Sim"] - new_sim_gs["Species2 Sum logfc Sim"]

        new_sum_gs["Total_Sims"] = np.repeat(tot_sim, new_sum_gs.shape[0])
        #If desired, write out the new_sim because we are going to be using the normal approximation for this one
        #Otherwise, we just count whether this sample beats the test statistic for the actual one
        if write_perm:
            new_sim_gs.to_csv(out_folder + "/Simulation" + str(j) + ".txt", sep = "\t")

        if metric == "PhyloP447":
            if not write_perm:
                new_sum_gs = new_sum_gs.drop(["Species1 Sum PhyloP447 Sim", "Species1 Sum PhyloP NonNeg Sim", "Species1 Sum Total_Vars Sim", "Species2 Sum PhyloP447 Sim", "Species2 Sum PhyloP NonNeg Sim", "Species2 Sum Total_Vars Sim", "PhyloP Difference Sim", "PhyloP L2FC Sim"], axis = 1)
            else:
                new_sim_gs = new_sim_gs.drop(["Species1 Sum PhyloP447 Sim", "Species1 Sum PhyloP NonNeg Sim", "Species1 Sum Total_Vars Sim", "Species2 Sum PhyloP447 Sim", "Species2 Sum PhyloP NonNeg Sim", "Species2 Sum Total_Vars Sim", "PhyloP Difference Sim", "PhyloP L2FC Sim"], axis = 1)
        elif metric == "abs logfc":
            new_sim_gs = new_sim_gs.drop(['Species1 Sum PhyloP447 Sim', 'Species1 Sum allele1_pred_counts Sim', 'Species1 Sum allele2_pred_counts Sim', 'Species1 Sum logfc Sim', 'Species1 Sum jsd Sim', 'Species1 Sum abs logfc Sim', 'Species1 Sum Total_Vars Sim', 'Species2 Sum PhyloP447 Sim', 'Species2 Sum allele1_pred_counts Sim', 'Species2 Sum allele2_pred_counts Sim', 'Species2 Sum logfc Sim', 'Species2 Sum jsd Sim', 'Species2 Sum abs logfc Sim', 'Species2 Sum Total_Vars Sim', 'Abs Logfc Difference Sim',], axis = 1)
        elif metric == "logfc":
            new_sim_gs = new_sim_gs.drop(['Species1 Sum PhyloP447 Sim', 'Species1 Sum allele1_pred_counts Sim', 'Species1 Sum allele2_pred_counts Sim', 'Species1 Sum logfc Sim', 'Species1 Sum jsd Sim', 'Species1 Sum abs logfc Sim', 'Species1 Sum Total_Vars Sim', 'Species2 Sum PhyloP447 Sim', 'Species2 Sum allele1_pred_counts Sim', 'Species2 Sum allele2_pred_counts Sim', 'Species2 Sum logfc Sim', 'Species2 Sum jsd Sim', 'Species2 Sum abs logfc Sim', 'Species2 Sum Total_Vars Sim', 'Logfc Difference Sim',], axis = 1)
        #print(new_sum_gs)
    else:
        if not write_perm:
            new_sum = new_sum.join(new_sim).fillna(0)
            #print(new_sum)
        else:
            new_sum = new_sim.copy()
        if metric == "PhyloP447":
            new_sum["PhyloP Difference Sim"] = new_sum["Species1 Sum PhyloP NonNeg Sim"] - new_sum["Species2 Sum PhyloP NonNeg Sim"]
            new_sum["PhyloP L2FC Sim"] = np.log2((new_sum["Species1 Sum PhyloP NonNeg Sim"] + 10)/(new_sum["Species2 Sum PhyloP NonNeg Sim"] + 10))
            if not write_perm:
                new_sum["Better_Sims_L2FC"] = new_sum["Better_Sims_L2FC"] + pd.DataFrame(((new_sum["PhyloP L2FC"] <= new_sum["PhyloP L2FC Sim"]) & (new_sum["PhyloP L2FC"] >= 0)) | ((new_sum["PhyloP L2FC"] >= new_sum["PhyloP L2FC Sim"]) & (new_sum["PhyloP L2FC"] < 0))).replace({True:1, False:0})[0]
                new_sum["Better_Sims_Sum"] = new_sum["Better_Sims_Sum"] + pd.DataFrame(((new_sum["PhyloP Difference"] <= new_sum["PhyloP Difference Sim"]) & (new_sum["PhyloP Difference"] >= 0)) | ((new_sum["PhyloP Difference"] >= new_sum["PhyloP Difference Sim"]) & (new_sum["PhyloP Difference"] < 0))).replace({True:1, False:0})[0]

        elif metric == "abs logfc":
            new_sum["Abs Logfc Difference Sim"] = new_sum["Species1 Sum abs logfc Sim"] - new_sum["Species2 Sum abs logfc Sim"]

        elif metric == "logfc":
            new_sum["Logfc Difference Sim"] = new_sum["Species1 Sum logfc Sim"] - new_sum["Species2 Sum logfc Sim"]
        
        if write_perm:
            new_sum.to_csv(out_folder + "/Simulation" + str(j) + ".txt", sep = "\t")
        new_sum["Total_Sims"] = np.repeat(tot_sim, new_sum.shape[0])
        if metric == "PhyloP447":
            new_sum = new_sum.drop(["Species1 Sum PhyloP447 Sim", "Species1 Sum PhyloP NonNeg Sim", "Species1 Sum Total_Vars Sim","Species2 Sum PhyloP447 Sim",  "Species2 Sum PhyloP NonNeg Sim", "Species2 Sum Total_Vars Sim", "PhyloP Difference Sim", "PhyloP L2FC Sim"], axis = 1)
        elif metric == "abs logfc":
            new_sum = new_sum.drop(['Species1 Sum PhyloP447 Sim', 'Species1 Sum allele1_pred_counts Sim', 'Species1 Sum allele2_pred_counts Sim', 'Species1 Sum logfc Sim', 'Species1 Sum jsd Sim', 'Species1 Sum abs logfc Sim', 'Species1 Sum Total_Vars Sim', 'Species2 Sum PhyloP447 Sim', 'Species2 Sum allele1_pred_counts Sim', 'Species2 Sum allele2_pred_counts Sim', 'Species2 Sum logfc Sim', 'Species2 Sum jsd Sim', 'Species2 Sum abs logfc Sim', 'Species2 Sum Total_Vars Sim', 'Abs Logfc Difference Sim',], axis = 1)
        elif metric == "logfc":
            new_sum = new_sum.drop(['Species1 Sum PhyloP447 Sim', 'Species1 Sum allele1_pred_counts Sim', 'Species1 Sum allele2_pred_counts Sim', 'Species1 Sum logfc Sim', 'Species1 Sum jsd Sim', 'Species1 Sum abs logfc Sim', 'Species1 Sum Total_Vars Sim', 'Species2 Sum PhyloP447 Sim', 'Species2 Sum allele1_pred_counts Sim', 'Species2 Sum allele2_pred_counts Sim', 'Species2 Sum logfc Sim', 'Species2 Sum jsd Sim', 'Species2 Sum abs logfc Sim', 'Species2 Sum Total_Vars Sim', 'Logfc Difference Sim',], axis = 1)
    #if c % 110000 == 0:
    #    new_sum.to_csv(out_folder + "/AccelEvol_Sim_" + str(start_iter) + "-" + str(end_iter) + "Set" + str(c) + ".txt", sep = "\t")
    #    print(c)

if not gene_set:
    new_sum = new_sum.loc[np.setdiff1d(new_sum.index, toss_gene)]
    new_sum.to_csv(out_folder + "/AccelEvol_" + str(start_iter) + "-" + str(end_iter) + ".txt", sep = "\t")
else:
    new_sum_gs.to_csv(out_folder + "/AccelEvol_" + str(start_iter) + "-" + str(end_iter) + ".txt", sep = "\t")
    
