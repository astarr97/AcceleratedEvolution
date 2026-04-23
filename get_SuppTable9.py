import os
import numpy as np
import pandas as pd

d_abrev = {"LiangSteinNeuron":"FC exc. neur.", "FetalChondrocytes":"F chond.", "SertoliMale":"FG sertoli", "preGC_IIaFemale":"FG preGC IIa",\
          "NeuralFemale":"FG neur.", "FetalGonadImmuneFemale":"FG immune", "VIP":"AC VIP inh. neur.", "LiangSteinProgenitor":"FC prog.",\
          "AdultHeartVentricularCardiomyocyte":"AH cardiomyo.", "AdultLoopOfHenle":"AK loop of henle", "FetalBrainNeurGlioblast_CB_VZ":"FCB glioblast",\
         "AdultProximalTubule":"AK prox. tub.", "FetalLeydigMale":"FG leydig", "SST":"AC SST inh neur.", "KosoyRoussosControlMicroglia":"AC microglia",\
         "FetalBrainFloorPlate":"FB fl. plate", "FetalArterialECs":"FH endoth.", "ASCT":"AC astro.", "FetalBrainCOP":"FB COP",\
         "AMY":"AA neur.", "PVALB":"AC PVALB inh neur.", "ITL23":"AC L2-3 IT neur.", "FetalBrainNeurCB_GNP_IPC_1":"FB inter. prog.", "FetalBrainNeurDAergic":"FB DA neur.",\
          "OGC":"AC Oligo.", "D1Pu":"AP D1 inh neur.", "FetalBrainNeurSerotonergic":"FB 5-HT neur.", "FetalBrainNeurDRG_2":"FS DRG neur.",\
          "FetalHeartPericytes":"FH peri.", "FetalHeartEndocardium":"FH endocard.", "FetalHeartCardiacFibroblasts":"FH fibro.", "FetalBrainNeurPurkinje_6":"FCB Purk. inh neur.",\
          "AdultHeartSmoothMuscle":"AH smooth musc.", "FetalBrainRoofPlate":"FB ro. plate"}

ind = 1
for file in os.listdir():
    if file.endswith("_AllSitesToDownload.txt.gz"):
        ct = file.replace("_AllSitesToDownload.txt.gz", "")
        print(ct)

        z = pd.read_csv(file, sep = "\t")
        z = z[z["SpecSup447"] > 250]
        z["Max acc."] = np.max(z[["allele1_pred_counts", "allele2_pred_counts"]], axis = 1)
        z = z.sort_values("Max acc.", ascending = False)
        z["Accessibility percentile"] = [x/z.shape[0] for x in list(range(1, 1 + z.shape[0]))[::-1]]
        z = z[(z["KeptAfterFilt"] == "Y") & (z["PhyloP447"] > 3) & (z["Derived"] == "H")]
        
        if ind:
            z = z[["Position", "PhyloP447", "SpecSup447", "NearestGene", "NearestDist", "logfc", "Accessibility percentile"]].copy()
            z.columns = ["Position, hg38", "PhyloP", "Gene symbol", "Distance to TSS of nearest gene", "Number of species in alignment", d_abrev[ct] + " l2fc", d_abrev[ct] + " accessibility percentile"]
            z = z.set_index("Position, hg38")
            df = z.copy()
            ind = 0
        else:
            z = z[["Position", "logfc", "Accessibility percentile"]].copy()
            z.columns = ["Position, hg38", d_abrev[ct] + " l2fc", d_abrev[ct] + " accessibility percentile"]
            df = df.join(z.set_index("Position, hg38"))
        print(df)
df.to_csv("Supplemental_Table9.csv")
