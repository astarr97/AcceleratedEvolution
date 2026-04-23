import pandas as pd
import numpy as np
from collections import Counter
from numpy.random import choice

#Read in
hum = pd.read_csv("HumChp_Mis_Final_Tri447.bed", sep = "\t", header = None).drop_duplicates()
chp = pd.read_csv("HumChp_Mis_Final_ChpDer_Tri447.bed", sep = "\t", header = None).drop_duplicates()
print(hum.sort_values(5))
#Remove those that do not have PhyloP and convert to a float
hum = hum[~hum[5].isin(["."])]
chp = chp[~chp[5].isin(["."])]
hum[5] = hum[5].astype(float)
chp[5] = chp[5].astype(float)

#Zero out the negative values, then remove unneeded columns
hum[5] = np.maximum(hum[5], 0)
chp[5] = np.maximum(chp[5], 0)
hum = hum[[3, 5, 8]]
chp = chp[[3, 5, 8]]
hum["TotalMissense"] = np.repeat(1, hum.shape[0])
chp["TotalMissense"] = np.repeat(1, chp.shape[0])

#Compute the actual values for acceleration, both sum and l2fc
hsum = hum.groupby([3]).sum(numeric_only=1)
csum = chp.groupby([3]).sum(numeric_only=1)
hsum.columns = ["Hum Sum PhyloP", "Num Hum Var"]
csum.columns = ["Chp Sum PhyloP", "Num Chp Var"]
new_sum = hsum.join(csum, how = "outer").fillna(0)
new_sum["PhyloP Difference"] = new_sum["Hum Sum PhyloP"] - new_sum["Chp Sum PhyloP"]
new_sum["PhyloP L2FC"] = np.log2((new_sum["Hum Sum PhyloP"] + 10)/(new_sum["Chp Sum PhyloP"] + 10))
new_sum.sort_values("PhyloP Difference")

#Filter to only those that have greater than or equal to 3 missense mutations per gene
new_sum = new_sum[(new_sum["Num Hum Var"] >= 3) | (new_sum["Num Chp Var"] >= 3)]

#Iterate through the canonical isoforms from uniprot
o = open("uniprotkb_taxonomy_id_9606_AND_reviewed_2023_09_04.fasta")
out = []
cur_gene = ""
cur_seq = ""
all_aa_props = {}
ind = 1
for line in o:
    if ">" in line:
        #Make sure that it is in the geneset
        if cur_gene in list(hum[3]) + list(chp[3]):
            x = Counter(cur_seq)
            hum_g = hum[hum[3].isin([cur_gene])]
            #Correct for human derived mutations (effectively reconstructs the ancestral sequence for the background)
            for index, row in hum_g.iterrows():
                x[row[8].split("/")[0]] = max(x[row[8].split("/")[0]] - 1, 0)
                x[row[8].split("/")[1]] = x[row[8].split("/")[1]] + 1
            out.append([cur_gene, len(cur_seq), x])
            for key in x.keys():
                if key not in all_aa_props.keys():
                    all_aa_props[key] = x[key]
                else:
                    all_aa_props[key] = all_aa_props[key] + x[key]
        try:
            cur_gene = line.split("GN=")[1].split(" ")[0]
            cur_seq = ""
        except:
            cur_gene = line.split("|")[2].split(" ")[0].split("_")[1]
            cur_seq = ""
        
    else:
        cur_seq = cur_seq + line.replace("\n", "")
        
#Create a dataframe
df = pd.DataFrame(out)
df = df.drop_duplicates(0)
df = df.set_index(0)
df.columns = ["Gene_Length", "AminoAcid_Counts"]

#Convert to proportions of amino acids rather than counts
x = []
for index, row in df.iterrows():
    d = {}
    prop = row["AminoAcid_Counts"]
    for key in prop.keys():
        d[key] = prop[key]/all_aa_props[key]
    x.append(d)
df["AminoAcid_Probs"] = x

#Convert to amino acid proportions per gene
dpre = {}
for index, row in df.iterrows():
    r = row["AminoAcid_Probs"]
    for key in r.keys():
        if key not in dpre.keys():
            dpre[key] = [[index], [r[key]]]
        else:
            dpre[key][0].append(index)
            dpre[key][1].append(r[key])
d = {}
for key in dpre.keys():
    x = dpre[key]
    x[1] = np.array(x[1])/np.sum(x[1])
    d[key] = x
np.median(df["Gene_Length"])

#Write out the preliminary file with no added permutation information
new_sum = new_sum.loc[np.intersect1d(df.index, new_sum.index)]
new_sum.to_csv("AccelEvol_Missense.txt", sep = "\t")

#Prepare for simulations and set a random seed
new_sum["Total_Sims"] = np.repeat(0, new_sum.shape[0])
new_sum["Better_Sims_L2FC"] = np.repeat(0, new_sum.shape[0])
new_sum["Better_Sims_Sum"] = np.repeat(0, new_sum.shape[0])

np.random.seed(6)

hum = hum[hum[3].isin(np.intersect1d(df.index, new_sum.index))]
chp = chp[chp[3].isin(np.intersect1d(df.index, new_sum.index))]

#Get the ancestral amino acid
hum["Anc"] = [x.split("/")[1] for x in list(hum[8])]
chp["Anc"] = [x.split("/")[0] for x in list(chp[8])]
input_sim = pd.concat([hum, chp]).sort_values("Anc")
hprob = hum.shape[0]/(hum.shape[0] + chp.shape[0])
aa = list(set(list(input_sim["Anc"])))
aa.sort()
c = 0

#Perform 1,000,000 simulations
for j in range(1000000):
    new_gene = []
    if c % 100000 == 0:
        new_sum.to_csv("AccelEvol_Missense_Sim" + str(c) + ".txt", sep = "\t")
        print(c)
    #Reassign each mutation to a new gene
    for i in aa:
        new_gene = new_gene + list(choice(d[i][0], p=d[i][1], size = len(list(input_sim[input_sim['Anc'].isin([i])].index))))
    
    #Reassign to human or chimp derived with probability equal to the true probability from the data
    input_sim["New_Gene"] = new_gene
    input_sim["New_Species"] = choice(["HumDer", "ChpDer"], p = [hprob, 1-hprob], size = input_sim.shape[0])
    input_sim_hum = input_sim[input_sim["New_Species"].isin(["HumDer"])]
    input_sim_chp = input_sim[input_sim["New_Species"].isin(["ChpDer"])]
    
    #Compute our statistics for the permuted data
    hsum = input_sim_hum.groupby(["New_Gene"]).sum(numeric_only=1)
    csum = input_sim_chp.groupby(["New_Gene"]).sum(numeric_only=1)
    hsum.columns = ["Hum Sum PhyloP Sim", "Num Hum Var Sim"]
    csum.columns = ["Chp Sum PhyloP Sim", "Chp Hum Var Sim"]
    new_sim = hsum.join(csum, how = "outer").fillna(0)

    new_sum = new_sum.join(new_sim)
    new_sum = new_sum.fillna(0)

    new_sum["PhyloP Difference Sim"] = new_sum["Hum Sum PhyloP Sim"] - new_sum["Chp Sum PhyloP Sim"]
    new_sum["PhyloP L2FC Sim"] = np.log2((new_sum["Hum Sum PhyloP Sim"] + 10)/(new_sum["Chp Sum PhyloP Sim"] + 10))

    #See if it is better or worse than the previous one, two-sided
    new_sum["Total_Sims"] = new_sum["Total_Sims"] + 1
    new_sum["Better_Sims_L2FC"] = new_sum["Better_Sims_L2FC"] + pd.DataFrame(((new_sum["PhyloP L2FC"] <= new_sum["PhyloP L2FC Sim"]) & (new_sum["PhyloP L2FC"] >= 0)) | ((new_sum["PhyloP L2FC"] >= new_sum["PhyloP L2FC Sim"]) & (new_sum["PhyloP L2FC"] < 0))).replace({True:1, False:0})[0]
    new_sum["Better_Sims_Sum"] = new_sum["Better_Sims_Sum"] + pd.DataFrame(((new_sum["PhyloP Difference"] <= new_sum["PhyloP Difference Sim"]) & (new_sum["PhyloP Difference"] >= 0)) | ((new_sum["PhyloP Difference"] >= new_sum["PhyloP Difference Sim"]) & (new_sum["PhyloP Difference"] < 0))).replace({True:1, False:0})[0]
    
    #Remove unnecessary columns
    new_sum = new_sum.drop(["Hum Sum PhyloP Sim", "Num Hum Var Sim", "Chp Sum PhyloP Sim", "Chp Hum Var Sim", "PhyloP Difference Sim", "PhyloP L2FC Sim"], axis = 1)
    c += 1

new_sum["L2FC p-value"] = new_sum["Better_Sims_L2FC"] / new_sum["Total_Sims"]
new_sum["Sum p-value"] = new_sum["Better_Sims_Sum"] / new_sum["Total_Sims"]
new_sum.sort_values("L2FC p-value")

new_sum.to_csv("AccelEvol_Missense_Sim" + str(c) + ".txt", sep = "\t")