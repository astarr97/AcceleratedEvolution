import pandas as pd
import numpy as np
from collections import Counter
from numpy.random import choice
import sys

geneset = sys.argv[1]

#Read in HPO.
if geneset == "HPO":
    HPO = pd.read_csv("/oak/stanford/groups/hbfraser/astarr/AccelEvol/AccelEvol_New/GeneSets/HPO_AccelEvol_Input.txt", sep = "\t")
    d_GS = {}
    for index, row in HPO.iterrows():
        d_GS[row["Term"]] = row["Genes"].split(";")
elif geneset == "GOBP":
    GOBP = pd.read_csv("/oak/stanford/groups/hbfraser/astarr/AccelEvol/AccelEvol_New/GeneSets/GOBP_AccelEvol_Input.txt", sep = "\t")
    d_GS = {}
    for index, row in GOBP.iterrows():
        d_GS[row["Term"]] = row["Genes"].split(";")
else:
    print("Unknown gene set")
    assert(False)

#Read in
hum = pd.read_csv("HumChp_Mis_Final_Tri447.bed", sep = "\t", header = None).drop_duplicates()
chp = pd.read_csv("HumChp_Mis_Final_ChpDer_Tri447.bed", sep = "\t", header = None).drop_duplicates()
hum = hum[~hum[5].isin(["."])]
chp = chp[~chp[5].isin(["."])]
hum[5] = hum[5].astype(float)
chp[5] = chp[5].astype(float)
hum[5] = np.maximum(hum[5], 0)
chp[5] = np.maximum(chp[5], 0)
hum = hum[[3, 5, 8]]
chp = chp[[3, 5, 8]]
hum["TotalMissense"] = np.repeat(1, hum.shape[0])
chp["TotalMissense"] = np.repeat(1, chp.shape[0])

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
            #Correct for human derived mutations
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
df = pd.DataFrame(out)
df = df.drop_duplicates(0)
df = df.set_index(0)
df.columns = ["Gene_Length", "AminoAcid_Counts"]

out_h = pd.DataFrame()
out_c = pd.DataFrame()
tot_h = 0
for key in d_GS.keys():
    genes = d_GS[key]
    hum_s = pd.DataFrame(hum[hum[3].isin(np.intersect1d(genes, df.index))])
    chp_s = pd.DataFrame(chp[chp[3].isin(np.intersect1d(genes, df.index))])
    tot = len(list(set(list(hum_s[3]) + list(chp_s[3]))))
    if tot >= 10 and tot <= 200:
        tot_h += len(list(hum_s.index))
        hum_s["Geneset"] = np.repeat(key, hum_s.shape[0])
        out_h = pd.concat([out_h, hum_s])
        chp_s["Geneset"] = np.repeat(key, chp_s.shape[0])
        out_c = pd.concat([out_c, chp_s])
print(tot_h)

#Compute the actual values
hsum = out_h.groupby(["Geneset"]).sum(numeric_only=1)
csum = out_c.groupby(["Geneset"]).sum(numeric_only=1)
hsum.columns = ["Hum Sum PhyloP", "Num Hum Var"]
csum.columns = ["Chp Sum PhyloP", "Num Chp Var"]
new_sum = hsum.join(csum, how = "outer").fillna(0)
new_sum["PhyloP Difference"] = new_sum["Hum Sum PhyloP"] - new_sum["Chp Sum PhyloP"]
new_sum["PhyloP L2FC"] = np.log2((new_sum["Hum Sum PhyloP"] + 10)/(new_sum["Chp Sum PhyloP"] + 10))
new_sum.sort_values("PhyloP L2FC")

out = []
hpo_tot_counts = Counter()
tot_muts_hum = 0
tot_muts_chp = 0
for key in list(set(list(out_h["Geneset"]) + list(out_c["Geneset"]))):
    out_hk = out_h[out_h["Geneset"].isin([key])]
    out_ck = out_c[out_c["Geneset"].isin([key])]
    tot_muts_hum += out_hk.shape[0]
    tot_muts_chp += out_ck.shape[0]
    df2 = df.loc[list(set(list(out_hk[3]) + list(out_ck[3])))]
    nc = Counter()
    length = 0
    for index, row in df2.iterrows():
        nc = nc + row["AminoAcid_Counts"]
        length += row["Gene_Length"]
    hpo_tot_counts = hpo_tot_counts + nc
    out.append([key, length, nc])
df_hpo = pd.DataFrame(out)
df_hpo.columns = ["Geneset_Name", "Geneset_Length", "AminoAcid_Counts"]

#Convert to probabilities
x = []
for index, row in df_hpo.iterrows():
    d = {}
    prop = row["AminoAcid_Counts"]
    for key in prop.keys():
        d[key] = prop[key]/hpo_tot_counts[key]
    x.append(d)
df_hpo["AminoAcid_Probs"] = x

df_hpo.index = df_hpo["Geneset_Name"]
dpre = {}
for index, row in df_hpo.iterrows():
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
np.median(df_hpo["Geneset_Length"])

new_sum.to_csv("AccelEvol_Missense_" + geneset + ".txt", sep = "\t")
new_sum["Total_Sims"] = np.repeat(0, new_sum.shape[0])
new_sum["Better_Sims_L2FC"] = np.repeat(0, new_sum.shape[0])
new_sum["Better_Sims_Sum"] = np.repeat(0, new_sum.shape[0])

c = 900000
out_h.index = out_h["Geneset"]
out_c.index = out_c["Geneset"]
out_h["Anc"] = [x.split("/")[1] for x in list(out_h[8])]
out_c["Anc"] = [x.split("/")[0] for x in list(out_c[8])]

input_sim = pd.concat([out_h, out_c]).sort_values("Anc")
hprob = out_h.shape[0]/(out_h.shape[0] + out_c.shape[0])
aa = list(set(list(input_sim["Anc"])))
aa.sort()

for j in range(900000, 1100000):
    np.random.seed(j)
    new_gene_set = []
    if c % 100000 == 0 and c != 900000:
        new_sum.to_csv("AccelEvol_Missense_" + geneset + "_Sim" + str(c) + ".txt", sep = "\t")
        print(c)
    for i in aa:
        new_gene_set = new_gene_set + list(choice(d[i][0], p=d[i][1], size = len(list(input_sim[input_sim['Anc'].isin([i])].index))))

    input_sim["New_Geneset"] = new_gene_set
    input_sim["New_Species"] = choice(["HumDer", "ChpDer"], p = [hprob, 1-hprob], size = input_sim.shape[0])
    input_sim_hum = input_sim[input_sim["New_Species"].isin(["HumDer"])]
    input_sim_chp = input_sim[input_sim["New_Species"].isin(["ChpDer"])]

    hsum = input_sim_hum.groupby(["New_Geneset"]).sum(numeric_only=1)
    csum = input_sim_chp.groupby(["New_Geneset"]).sum(numeric_only=1)
    hsum.columns = ["Hum Sum PhyloP Sim", "Num Hum Var Sim"]
    csum.columns = ["Chp Sum PhyloP Sim", "Chp Hum Var Sim"]
    new_sim = hsum.join(csum, how = "outer").fillna(0)

    new_sum = new_sum.join(new_sim)
    new_sum = new_sum.fillna(0)

    new_sum["PhyloP Difference Sim"] = new_sum["Hum Sum PhyloP Sim"] - new_sum["Chp Sum PhyloP Sim"]
    new_sum["PhyloP L2FC Sim"] = np.log2((new_sum["Hum Sum PhyloP Sim"] + 10)/(new_sum["Chp Sum PhyloP Sim"] + 10))


    new_sum["Total_Sims"] = new_sum["Total_Sims"] + 1
    new_sum["Better_Sims_L2FC"] = new_sum["Better_Sims_L2FC"] + pd.DataFrame(((new_sum["PhyloP L2FC"] <= new_sum["PhyloP L2FC Sim"]) & (new_sum["PhyloP L2FC"] >= 0)) | ((new_sum["PhyloP L2FC"] >= new_sum["PhyloP L2FC Sim"]) & (new_sum["PhyloP L2FC"] < 0))).replace({True:1, False:0})[0]
    new_sum["Better_Sims_Sum"] = new_sum["Better_Sims_Sum"] + pd.DataFrame(((new_sum["PhyloP Difference"] <= new_sum["PhyloP Difference Sim"]) & (new_sum["PhyloP Difference"] >= 0)) | ((new_sum["PhyloP Difference"] >= new_sum["PhyloP Difference Sim"]) & (new_sum["PhyloP Difference"] < 0))).replace({True:1, False:0})[0]
    new_sum = new_sum.drop(["Hum Sum PhyloP Sim", "Num Hum Var Sim", "Chp Sum PhyloP Sim", "Chp Hum Var Sim", "PhyloP Difference Sim", "PhyloP L2FC Sim"], axis = 1)
    c += 1
    
new_sum["L2FC p-val"] = np.maximum(new_sum["Better_Sims_L2FC"]/new_sum["Total_Sims"], 1/1000000)
new_sum["Sum p-val"] = np.maximum(new_sum["Better_Sims_Sum"]/new_sum["Total_Sims"], 1/1000000)

new_sum.sort_values("L2FC p-value")

new_sum.to_csv("AccelEvol_Missense_" + geneset + "_Sim" + str(c) + ".txt", sep = "\t")

