import pandas as pd
import numpy as np
from scipy.stats import fisher_exact,binomtest,ttest_ind

o = open("Config_ToRun.txt")
out = []
for line in o:
    ct = line.replace("\n", "")
    if ct not in ["FetalHeartCardiacFibroblasts", "FetalHeartPericytes", "FetalHeartEndocardium", 'AdultHeartSmoothMuscle']:
        print(ct)
        vv = pd.read_csv("/oak/stanford/groups/hbfraser/astarr/PosSelect_New/ToDownload_Vars/Poly_MAF0.25_" + ct + ".txt.gz", sep = "\t")
        te_blacklist = pd.read_csv("BlacklistTE_Variants.txt", sep = "\t")
        vv = vv[~vv["Position"].isin(te_blacklist["Position"])]
        v = pd.read_csv(ct + "_AllSitesToDownload.txt.gz", sep = "\t")
        
        #Filter
        vx = v[v["SpecSup447"] > 248]
        vvx = vv[vv["SpecSup447"] > 250]
        vx["Chrom"] = [x.split(":")[0] for x in vx["Position"]]
        vx["Pos"] = [int(x.split(":")[1]) for x in vx["Position"]]
        vvx["Chrom"] = [x.split(":")[0] for x in vvx["Position"]]
        vvx["Pos"] = [int(x.split(":")[1]) for x in vvx["Position"]]
        inter = np.intersect1d(v["Position"], vv["Position"])
        vx = vx[~vx["Position"].isin(inter)]
        vxh = vx[vx["Derived"] == "H"]
        vxc = vx[vx["Derived"] == "C"]
        vxh = vxh[vxh["KeptAfterFilt"] == "Y"]
        vxc = vxc[vxc["KeptAfterFilt"] == "Y"]
        
        vxhu = vxh[(vxh["PhyloP447"] > 3) & (vxh["logfc"] < -0.5)]
        khu = []
        for index, row in vxhu.iterrows():
            khu = khu + [row["Chrom"] + ":" + str(x) for x in range(row["Pos"] - 100, row["Pos"] + 100)]
        
        vxhd = vxh[(vxh["PhyloP447"] > 3) & (vxh["logfc"] > 0.5)]
        khd = []
        for index, row in vxhd.iterrows():
            khd = khd + [row["Chrom"] + ":" + str(x) for x in range(row["Pos"] - 100, row["Pos"] + 100)]

        vxcu = vxc[(vxc["PhyloP447"] > 3) & (vxc["logfc"] < -0.5)]
        kcu = []
        for index, row in vxcu.iterrows():
            kcu = kcu + [row["Chrom"] + ":" + str(x) for x in range(row["Pos"] - 100, row["Pos"] + 100)]
        
        vxcd = vxc[(vxc["PhyloP447"] > 3) & (vxc["logfc"] > 0.5)]
        kcd = []
        for index, row in vxcd.iterrows():
            kcd = kcd + [row["Chrom"] + ":" + str(x) for x in range(row["Pos"] - 100, row["Pos"] + 100)]

        vvxhu = vvx[vvx["Position"].isin(khu)]
        vvxhd = vvx[vvx["Position"].isin(khd)]
        vvxcu = vvx[vvx["Position"].isin(kcu)]
        vvxcd = vvx[vvx["Position"].isin(kcd)]
        
        bh = vvxhd[vvxhd["PhyloP447"] > 3].shape[0]
        bc = vvxcd[vvxcd["PhyloP447"] > 3].shape[0]
        sh = vvxhd[vvxhd["PhyloP447"] < 3].shape[0]
        sc = vvxcd[vvxcd["PhyloP447"] < 3].shape[0]
        fe = fisher_exact([[bh, bc], [sh, sc]])
        out.append([ct, "PhyloP", 0.5, "Near down", bh, bc, sh, sc, fe[0], fe[1]])
        
        bh = vvxhu[vvxhu["PhyloP447"] > 3].shape[0]
        bc = vvxcu[vvxcu["PhyloP447"] > 3].shape[0]
        sh = vvxhu[vvxhu["PhyloP447"] < 3].shape[0]
        sc = vvxcu[vvxcu["PhyloP447"] < 3].shape[0]
        fe = fisher_exact([[bh, bc], [sh, sc]])
        out.append([ct, "PhyloP", 0.5, "Near up", bh, bc, sh, sc, fe[0], fe[1]])
        
        fh = vxhd.shape[0]
        fc = vxcd.shape[0]
        ph = vvxhd.shape[0]
        pc = vvxcd.shape[0]
        fe = fisher_exact([[ph, pc], [fh, fc]])
        out.append([ct, "Total SNPs", 0.5, "Near down", ph, pc, fh, fc, fe[0], fe[1]])
        
        fh = vxhu.shape[0]
        fc = vxcu.shape[0]
        ph = vvxhu.shape[0]
        pc = vvxcu.shape[0]
        fe = fisher_exact([[ph, pc], [fh, fc]])
        out.append([ct, "Total SNPs", 0.5, "Near up", ph, pc, fh, fc, fe[0], fe[1]])
        
        ### Repeating for 0.25-5
        vxhu = vxh[(vxh["PhyloP447"] > 3) & (vxh["logfc"] < -0.25) & (vxh["logfc"] > -0.5)]
        khu = []
        for index, row in vxhu.iterrows():
            khu = khu + [row["Chrom"] + ":" + str(x) for x in range(row["Pos"] - 100, row["Pos"] + 100)]
        
        vxhd = vxh[(vxh["PhyloP447"] > 3) & (vxh["logfc"] > 0.25) & (vxh["logfc"] < 0.5)]
        khd = []
        for index, row in vxhd.iterrows():
            khd = khd + [row["Chrom"] + ":" + str(x) for x in range(row["Pos"] - 100, row["Pos"] + 100)]

        vxcu = vxc[(vxc["PhyloP447"] > 3) & (vxc["logfc"] < -0.25) & (vxc["logfc"] > -0.5)]
        kcu = []
        for index, row in vxcu.iterrows():
            kcu = kcu + [row["Chrom"] + ":" + str(x) for x in range(row["Pos"] - 100, row["Pos"] + 100)]
        
        vxcd = vxc[(vxc["PhyloP447"] > 3) & (vxc["logfc"] > 0.25) & (vxc["logfc"] < 0.5)]
        kcd = []
        for index, row in vxcd.iterrows():
            kcd = kcd + [row["Chrom"] + ":" + str(x) for x in range(row["Pos"] - 100, row["Pos"] + 100)]

        vvxhu = vvx[vvx["Position"].isin(khu)]
        vvxhd = vvx[vvx["Position"].isin(khd)]
        vvxcu = vvx[vvx["Position"].isin(kcu)]
        vvxcd = vvx[vvx["Position"].isin(kcd)]
        
        bh = vvxhd[vvxhd["PhyloP447"] > 3].shape[0]
        bc = vvxcd[vvxcd["PhyloP447"] > 3].shape[0]
        sh = vvxhd[vvxhd["PhyloP447"] < 3].shape[0]
        sc = vvxcd[vvxcd["PhyloP447"] < 3].shape[0]
        fe = fisher_exact([[bh, bc], [sh, sc]])
        out.append([ct, "PhyloP", 0.25, "Near down", bh, bc, sh, sc, fe[0], fe[1]])
        
        bh = vvxhu[vvxhu["PhyloP447"] > 3].shape[0]
        bc = vvxcu[vvxcu["PhyloP447"] > 3].shape[0]
        sh = vvxhu[vvxhu["PhyloP447"] < 3].shape[0]
        sc = vvxcu[vvxcu["PhyloP447"] < 3].shape[0]
        fe = fisher_exact([[bh, bc], [sh, sc]])
        out.append([ct, "PhyloP", 0.25, "Near up", bh, bc, sh, sc, fe[0], fe[1]])
        
        fh = vxhd.shape[0]
        fc = vxcd.shape[0]
        ph = vvxhd.shape[0]
        pc = vvxcd.shape[0]
        fe = fisher_exact([[ph, pc], [fh, fc]])
        out.append([ct, "Total SNPs", 0.25, "Near down", ph, pc, fh, fc, fe[0], fe[1]])
        
        fh = vxhu.shape[0]
        fc = vxcu.shape[0]
        ph = vvxhu.shape[0]
        pc = vvxcu.shape[0]
        fe = fisher_exact([[ph, pc], [fh, fc]])
        out.append([ct, "Total SNPs", 0.25, "Near up", ph, pc, fh, fc, fe[0], fe[1]])
        
df = pd.DataFrame(out)
df.columns = ["Cell type", "Metric", "Cutoff for fixed abs logfc", "Var type", "Large effect human", "Large effect chimp", "Small effect human", "Small effect chimp", "Odds ratio", "p-value"]
df.to_csv("Checking_Constraint2_ML_Poly.txt", sep = "\t", index = False)
