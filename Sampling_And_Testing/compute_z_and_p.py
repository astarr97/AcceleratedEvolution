import sys
import pandas as pd
import numpy as np
import os
from scipy.stats import norm
from statsmodels.stats.multitest import fdrcorrection

folder_name = sys.argv[1]
prefix = sys.argv[2]
metric = " ".join(sys.argv[3].split("_"))
out_folder = sys.argv[4]
min_var = int(sys.argv[5])
try:
    os.mkdir(out_folder)
except:
    pass

folder = os.listdir(folder_name)

folder_new = []
for fold in folder:
    if "Simulation" in fold:
        folder_new.append(fold)
        
folder_new.sort(key = lambda x: int(x.split("_")[0].replace("Simulation", "").replace(".txt", "")))
folder = ["AccelEvol_Init.txt"] + folder_new

#start_cols = ["Species1 Sum PhyloP447", "Num Species1 Var", "Species2 Sum PhyloP447", "Num Species2 Var", "PhyloP447 Difference", "PhyloP447 L2FC", "Median_PhyloP447_Species1", "Median_PhyloP447_Species2", "MWU p-value", "Centered PhyloP447 L2FC", "Centered PhyloP447 Difference"]
#start_cols = ["Species1 Sum PhyloP447", "Num Species1 Var", "Species2 Sum PhyloP447", "Num Species2 Var", "PhyloP447 Difference", "PhyloP447 L2FC", "Centered PhyloP447 L2FC", "Centered PhyloP447 Difference"]
if metric == "PhyloP":
    start_cols = ["Species1 Sum PhyloP NonNeg", "Species1 Sum Total_Vars", "Species2 Sum PhyloP NonNeg", "Species2 Sum Total_Vars", "PhyloP Difference", "PhyloP L2FC"]
    add_cols = ['PhyloP Difference Sim', 'PhyloP L2FC Sim', 'Species1 Sum Total_Vars Sim', 'Species2 Sum Total_Vars Sim', "Species1 Sum PhyloP NonNeg Sim", "Species2 Sum PhyloP NonNeg Sim"]

elif metric == "logfc":
    start_cols = ["Species1 Sum allele1_pred_counts", "Species1 Sum allele2_pred_counts", "Species1 Sum logfc", "Species1 Sum jsd", "Species1 Sum abs logfc", "Species1 Sum Total_Vars", \
    "Species2 Sum allele1_pred_counts", "Species2 Sum allele2_pred_counts", "Species2 Sum logfc", "Species2 Sum jsd", "Species2 Sum abs logfc", "Species2 Sum Total_Vars", "Logfc Difference"]
    add_cols = ["Species1 Sum logfc Sim", "Species1 Sum Total_Vars Sim", "Species2 Sum logfc Sim", "Species2 Sum Total_Vars Sim", "Logfc Difference Sim"]
elif metric == "abs logfc":
    start_cols = ["Species1 Sum allele1_pred_counts", "Species1 Sum allele2_pred_counts", "Species1 Sum logfc", "Species1 Sum jsd", "Species1 Sum abs logfc", "Species1 Sum Total_Vars", \
    "Species2 Sum allele1_pred_counts", "Species2 Sum allele2_pred_counts", "Species2 Sum logfc", "Species2 Sum jsd", "Species2 Sum abs logfc", "Species2 Sum Total_Vars", "Abs Logfc Difference"]
    add_cols = ["Species1 Sum abs logfc Sim", "Species1 Sum Total_Vars Sim", "Species2 Sum abs logfc Sim", "Species2 Sum Total_Vars Sim", "Abs Logfc Difference Sim"]
ind = 1
v_start = 0
v_temp = 0
seen = []
c = 0
for file in folder:
    c += 1
    if c % 100 == 0:
        print(c)
    v = pd.read_csv(folder_name + "/" + file, sep = "\t")
    v = v.set_index(v.columns[0])

    if ind:
        v_start = v[start_cols].copy()
        v_temp = v[[]]
        ind = 0
    else:
        v = v[add_cols].copy()
        for col in add_cols:
            if "Var" not in col:
                v[col] = np.round(v[col], 3)
        seen.append(file.split("_")[0].replace("Simulation", "").replace(".txt", ""))
        v.columns = [x + " " + file.split("_")[0].replace("Simulation", "").replace(".txt", "") for x in add_cols]
        v_temp = v_temp.join(v)

v_temp = v_temp.astype(str)
for col in add_cols:
    out_col = []
    v_col = v_temp[[col + " " + x for x in seen]]
    for index, row in v_col.iterrows():
        out_col.append([index, ";".join(list(row))])
    to_join = pd.DataFrame(out_col).set_index(0)
    to_join.columns = [col]
    v_start = v_start.join(to_join)
v_start.to_csv(out_folder + "/" + prefix + "_AllSims.txt", sep = "\t")

v_start = v_start[(v_start["Species1 Sum Total_Vars"] > min_var) & (v_start["Species2 Sum Total_Vars"] > min_var)]

y_linked = list(pd.read_csv("../../HumanYGenes.txt.gz", sep = "\t", header = None)[0])

v_start = v_start.loc[np.setdiff1d(v_start.index, y_linked)].copy()
print(v_start)
zscore_dif = []
zscore_l2fc = []
pvalue_dif = []
pvalue_l2fc = []

zscore_dif_corr_tot = []
zscore_l2fc_corr_tot = []
pvalue_dif_corr_tot = []
pvalue_l2fc_corr_tot = []

zscore_dif_corr_s = []
zscore_l2fc_corr_s = []
pvalue_dif_corr_s = []
pvalue_l2fc_corr_s = []
for index, row in v_start.iterrows():
    if metric == "PhyloP":
        
        #Get the actual number of variants
        sum_s1_s2 = row["Species1 Sum Total_Vars"] + row["Species2 Sum Total_Vars"]
        sum_s1 = row["Species1 Sum Total_Vars"]
        sum_s2 = row["Species2 Sum Total_Vars"]
        
        #Get the permutation number of variants
        s1_nv = np.array([np.float64(x) for x in row['Species1 Sum Total_Vars Sim'].split(";")])
        s2_nv = np.array([np.float64(x) for x in row['Species2 Sum Total_Vars Sim'].split(";")])
        
        #Get the statistic
        s1_stat = np.array([np.float64(x) for x in row['Species1 Sum PhyloP NonNeg Sim'].split(";")])
        s2_stat = np.array([np.float64(x) for x in row['Species2 Sum PhyloP NonNeg Sim'].split(";")])
        
        #Correct for the total number of variants assigned to the gene
        s1_stat_corr_tot = np.divide(s1_stat, s1_nv/(sum_s1_s2/2))
        s2_stat_corr_tot = np.divide(s2_stat, s2_nv/(sum_s1_s2/2))
        
        #Correct instead for the actual number of variants assigned to each allele in the process, effectively testing for shifts in the distribution
        s1_stat_corr_s1 = np.divide(s1_stat, s1_nv/sum_s1)
        s2_stat_corr_s2 = np.divide(s2_stat, s2_nv/sum_s2)
        
        #Compute the corrected differences
        difs_corr_tot = list(s1_stat_corr_tot - s2_stat_corr_tot)
        difs_corr_s = list(s1_stat_corr_s1 - s2_stat_corr_s2)
        
        l2fcs_corr_tot = list(np.log2((s1_stat_corr_tot + 10)/(s2_stat_corr_tot + 10)))
        l2fcs_corr_s = list(np.log2((s1_stat_corr_s1 + 10)/(s2_stat_corr_s2 + 10)))
        
        difs = [np.float64(x) for x in row['PhyloP Difference Sim'].split(";")]
        l2fcs = [np.float64(x) for x in row['PhyloP L2FC Sim'].split(";")]
        stdev_dif = np.std(difs)
        mean_dif = np.mean(difs)
        stdev_l2fc = np.std(l2fcs)
        mean_l2fc = np.mean(l2fcs)
    
        z_dif = (row["PhyloP Difference"] - mean_dif)/stdev_dif
        z_l2fc = (row["PhyloP L2FC"] - mean_l2fc)/stdev_l2fc
    
        zscore_dif.append(z_dif)
        zscore_l2fc.append(z_l2fc)
        p_dif = norm.sf(abs(z_dif))*2
        p_l2fc = norm.sf(abs(z_l2fc))*2
        pvalue_dif.append(p_dif)
        pvalue_l2fc.append(p_l2fc)
        
        #Repeat correcting for total vars assigned
        stdev_dif = np.std(difs_corr_tot)
        mean_dif = np.mean(difs_corr_tot)
        z_dif = (row["PhyloP Difference"] - mean_dif)/stdev_dif
        zscore_dif_corr_tot.append(z_dif)
        p_dif = norm.sf(abs(z_dif))*2
        pvalue_dif_corr_tot.append(p_dif)
        
        #Now correcting for allele assignment instead, effectively testing for shifts in the distribution
        stdev_dif = np.std(difs_corr_s)
        mean_dif = np.mean(difs_corr_s)
        z_dif = (row["PhyloP Difference"] - mean_dif)/stdev_dif
        zscore_dif_corr_s.append(z_dif)
        p_dif = norm.sf(abs(z_dif))*2
        pvalue_dif_corr_s.append(p_dif)
        
        #Repeat correcting for total vars assigned
        stdev_l2fc = np.std(l2fcs_corr_tot)
        mean_l2fc = np.mean(l2fcs_corr_tot)
        z_l2fc = (row["PhyloP L2FC"] - mean_l2fc)/stdev_l2fc
        zscore_l2fc_corr_tot.append(z_l2fc)
        p_l2fc = norm.sf(abs(z_l2fc))*2
        pvalue_l2fc_corr_tot.append(p_l2fc)
        
        #Now correcting for allele assignment instead, effectively testing for shifts in the distribution
        stdev_l2fc = np.std(l2fcs_corr_s)
        mean_l2fc = np.mean(l2fcs_corr_s)
        z_l2fc = (row["PhyloP L2FC"] - mean_l2fc)/stdev_l2fc
        zscore_l2fc_corr_s.append(z_l2fc)
        p_l2fc = norm.sf(abs(z_l2fc))*2
        pvalue_l2fc_corr_s.append(p_l2fc)
        
    elif metric == "abs logfc":
        #Get the actual number of variants
        sum_s1_s2 = row["Species1 Sum Total_Vars"] + row["Species2 Sum Total_Vars"]
        sum_s1 = row["Species1 Sum Total_Vars"]
        sum_s2 = row["Species2 Sum Total_Vars"]
        
        #Get the permutation number of variants
        s1_nv = np.array([np.float64(x) for x in row['Species1 Sum Total_Vars Sim'].split(";")])
        s2_nv = np.array([np.float64(x) for x in row['Species2 Sum Total_Vars Sim'].split(";")])
        
        #Get the statistic
        s1_stat = np.array([np.float64(x) for x in row['Species1 Sum abs logfc Sim'].split(";")])
        s2_stat = np.array([np.float64(x) for x in row['Species2 Sum abs logfc Sim'].split(";")])
        
        #Correct for the total number of variants assigned to the gene
        s1_stat_corr_tot = np.divide(s1_stat, s1_nv/(sum_s1_s2/2))
        s2_stat_corr_tot = np.divide(s2_stat, s2_nv/(sum_s1_s2/2))
        
        #Correct instead for the actual number of variants assigned to each allele in the process, effectively testing for shifts in the distribution
        s1_stat_corr_s1 = np.divide(s1_stat, s1_nv/sum_s1)
        s2_stat_corr_s2 = np.divide(s2_stat, s2_nv/sum_s2)
        
        #Compute the corrected differences
        difs_corr_tot = list(s1_stat_corr_tot - s2_stat_corr_tot)
        difs_corr_s = list(s1_stat_corr_s1 - s2_stat_corr_s2)
        
        difs = [np.float64(x) for x in row['Abs Logfc Difference Sim'].split(";")]
        stdev_dif = np.std(difs)
        mean_dif = np.mean(difs)
        z_dif = (row["Abs Logfc Difference"] - mean_dif)/stdev_dif
        zscore_dif.append(z_dif)
        p_dif = norm.sf(abs(z_dif))*2
        pvalue_dif.append(p_dif)
        
        #Repeat correcting for total vars assigned
        stdev_dif = np.std(difs_corr_tot)
        mean_dif = np.mean(difs_corr_tot)
        z_dif = (row["Abs Logfc Difference"] - mean_dif)/stdev_dif
        zscore_dif_corr_tot.append(z_dif)
        p_dif = norm.sf(abs(z_dif))*2
        pvalue_dif_corr_tot.append(p_dif)
        
        #Now correcting for allele assignment instead, effectively testing for shifts in the distribution
        stdev_dif = np.std(difs_corr_s)
        mean_dif = np.mean(difs_corr_s)
        z_dif = (row["Abs Logfc Difference"] - mean_dif)/stdev_dif
        zscore_dif_corr_s.append(z_dif)
        p_dif = norm.sf(abs(z_dif))*2
        pvalue_dif_corr_s.append(p_dif)
    elif metric == "logfc":
        #Get the actual number of variants
        sum_s1_s2 = row["Species1 Sum Total_Vars"] + row["Species2 Sum Total_Vars"]
        sum_s1 = row["Species1 Sum Total_Vars"]
        sum_s2 = row["Species2 Sum Total_Vars"]
        
        #Get the permutation number of variants
        s1_nv = np.array([np.float64(x) for x in row['Species1 Sum Total_Vars Sim'].split(";")])
        s2_nv = np.array([np.float64(x) for x in row['Species2 Sum Total_Vars Sim'].split(";")])
        
        #Get the statistic
        s1_stat = np.array([np.float64(x) for x in row['Species1 Sum logfc Sim'].split(";")])
        s2_stat = np.array([np.float64(x) for x in row['Species2 Sum logfc Sim'].split(";")])
        
        #Correct for the total number of variants assigned to the gene
        s1_stat_corr_tot = np.divide(s1_stat, s1_nv/(sum_s1_s2/2))
        s2_stat_corr_tot = np.divide(s2_stat, s2_nv/(sum_s1_s2/2))
        
        #Correct instead for the actual number of variants assigned to each allele in the process, effectively testing for shifts in the distribution
        s1_stat_corr_s1 = np.divide(s1_stat, s1_nv/sum_s1)
        s2_stat_corr_s2 = np.divide(s2_stat, s2_nv/sum_s2)
        
        #Compute the corrected differences
        difs_corr_tot = list(s1_stat_corr_tot - s2_stat_corr_tot)
        difs_corr_s = list(s1_stat_corr_s1 - s2_stat_corr_s2)
        
        #Don't correct for potential issues with the number of sites assigned to genes in the permutation
        difs = [np.float64(x) for x in row['Logfc Difference Sim'].split(";")]
        stdev_dif = np.std(difs)
        mean_dif = np.mean(difs)
        z_dif = (row["Logfc Difference"] - mean_dif)/stdev_dif
        zscore_dif.append(z_dif)
        p_dif = norm.sf(abs(z_dif))*2
        pvalue_dif.append(p_dif)
        
        #Repeat correcting for total vars assigned
        stdev_dif = np.std(difs_corr_tot)
        mean_dif = np.mean(difs_corr_tot)
        z_dif = (row["Logfc Difference"] - mean_dif)/stdev_dif
        zscore_dif_corr_tot.append(z_dif)
        p_dif = norm.sf(abs(z_dif))*2
        pvalue_dif_corr_tot.append(p_dif)
        
        #Now correcting for allele assignment instead, effectively testing for shifts in the distribution
        stdev_dif = np.std(difs_corr_s)
        mean_dif = np.mean(difs_corr_s)
        z_dif = (row["Logfc Difference"] - mean_dif)/stdev_dif
        zscore_dif_corr_s.append(z_dif)
        p_dif = norm.sf(abs(z_dif))*2
        pvalue_dif_corr_s.append(p_dif)
        

v_write = v_start[start_cols].copy()

v_write["Z-score Difference"] = zscore_dif
v_write["p-value Difference"] = pvalue_dif

v_write["Z-score Difference Corr Tot"] = zscore_dif_corr_tot
v_write["p-value Difference Corr Tot"] = pvalue_dif_corr_tot

v_write["Z-score Difference Corr S"] = zscore_dif_corr_s
v_write["p-value Difference Corr S"] = pvalue_dif_corr_s

if metric == "PhyloP":
    v_write["Z-score L2FC"] = zscore_l2fc
    v_write["p-value L2FC"] = pvalue_l2fc
    
v_write = v_write.dropna()
v_write["FDR Difference"] = fdrcorrection(v_write["p-value Difference"])[1]
v_write["FDR Difference Corr Tot"] = fdrcorrection(v_write["p-value Difference Corr Tot"])[1]
v_write["FDR Difference Corr S"] = fdrcorrection(v_write["p-value Difference Corr S"])[1]

if metric == "PhyloP":
    v_write["FDR L2FC"] = fdrcorrection(v_write["p-value L2FC"])[1]
v_write = v_write.sort_values("FDR Difference")
v_write.to_csv(out_folder + "/" + prefix + "_Results.txt", sep = "\t")

