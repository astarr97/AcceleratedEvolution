import numpy as np
import pandas as pd
from math import floor
from scipy.stats import fisher_exact
import os

out = []
for file in os.listdir():
    if file.endswith("_AllSitesToDownload.txt.gz"):
        ct = file.replace("_AllSitesToDownload.txt.gz", "")
        print(ct)
        z = pd.read_csv(file, sep = "\t")
        z = z.drop(["Position.1"], axis = 1)
        z = z[z["SpecSup447"] > 250]
        z = z[z["KeptAfterFilt"] == "Y"]
        z = z[z["Derived"] == "H"]
        z["abs logfc"] = np.abs(z["logfc"])
        z["Max"] = np.max(z[["allele1_pred_counts", "allele2_pred_counts"]], axis = 1)
        z = z.sort_values("Max")

        zt = z.tail(floor(0.1*z.shape[0]))
        znt = z.head(floor(0.9*z.shape[0]))
        
        ztt = zt.sort_values("abs logfc").tail(floor(0.1*zt.shape[0]))
        ztnt = zt.sort_values("abs logfc").head(floor(0.9*zt.shape[0]))
        
        for b in [(0, 3), (3, 6), (6, 9), (9, 12)]:
            cons_h = znt[(znt["PhyloP447"] > b[0]) & (znt["PhyloP447"] < b[1])].shape[0]
            ncons_h = znt[znt["PhyloP447"] < 0].shape[0]
            cons_t = zt[(zt["PhyloP447"] > b[0]) & (zt["PhyloP447"] < b[1])].shape[0]
            ncons_t = zt[zt["PhyloP447"] < 0].shape[0]
            fe = fisher_exact([[cons_t, ncons_t], [cons_h, ncons_h]])
            out.append([ct, cons_t, ncons_t, cons_h, ncons_h, fe[0], fe[1], "Total CA", b])
        
            cons_h = ztnt[(ztnt["PhyloP447"] > b[0]) & (ztnt["PhyloP447"] < b[1])].shape[0]
            ncons_h = ztnt[ztnt["PhyloP447"] < 0].shape[0]
            cons_t = ztt[(ztt["PhyloP447"] > b[0]) & (ztt["PhyloP447"] < b[1])].shape[0]
            ncons_t = ztt[ztt["PhyloP447"] < 0].shape[0]
            fe = fisher_exact([[cons_t, ncons_t], [cons_h, ncons_h]])
            out.append([ct, cons_t, ncons_t, cons_h, ncons_h, fe[0], fe[1], "Effect on CA", b])
        
df = pd.DataFrame(out)
df.columns = ["Cell type", "Conserved_Tail", "Unconserved_Tail", "Conserved_Head", "Unconserved_Head", "p-value", "Odds ratio", "Test", "Bin"]
df.to_csv("TestingConsEnrichment.csv", index = False)