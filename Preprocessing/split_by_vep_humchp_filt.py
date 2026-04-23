o = open("All_Summarized_HumChp_GW_EnsConsol.bed")
out = open("HumChp_KeptBy3WGS_Final_VEP.bed", 'w')

for line in o:
    l = line.replace("\n", "").split("\t")
    anno = ""
    if "missense" not in line:
        if "UTR" in line:
            if "3_prime_UTR_variant" in line:
                anno = "3UTR"
            else:
                anno = "5UTR"
        else:
            anno = "NC"
    if "missense_variant" in l[5].split(";"):
        anno = "Mis"
    if "missense_variant" not in l[5].split(";") and "synonymous_variant" in l[5].split(";") and "TF_binding_site_variant" not in l[5].split(";") and "CTCF_binding_site" not in l[5].split(";") and "InPromoter" not in line:
        anno = anno + "_Syn"
    out.write("\t".join(l[0:5] + [anno]) + "\n")
out.close()
