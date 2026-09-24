o = open("All_Summarized_GorDer_Test_GW_EnsConsol.bed")
out = open("GorDer_Mis_UTR_ForFix.bed", "w")

for line in o:
    l = line.rstrip("\n").split("\t")

    # Consequences and genes are semicolon-separated
    consequences = l[5].split(";")
    genes = l[6].split(";")

    anno = ""
    gene = ""

    # Find missense or UTR annotation and corresponding gene
    for i, consequence in enumerate(consequences):

        if "missense_variant" in consequence:
            anno = "Mis"
            gene = genes[i] if i < len(genes) else ""
            break

        elif "3_prime_UTR_variant" in consequence:
            anno = "3UTR"
            gene = genes[i] if i < len(genes) else ""
            break

        elif "5_prime_UTR_variant" in consequence:
            anno = "5UTR"
            gene = genes[i] if i < len(genes) else ""
            break

    # Only write missense and UTR variants
    if anno:
        out.write("\t".join(l[0:5] + [anno, gene]) + "\n")

out.close()
o.close()
