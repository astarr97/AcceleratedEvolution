import os

c = 0
out = open("All_Summarized_GorDer_Test_GW_Ens.txt", 'w')
out.write("\t".join(["Chrom", "Pos", "ID", "Human", "Chimp", "Variant_type", "Gene", "Gene_type", "protein_change", "APPRIS", "SIFT", "PolyPhen"]) + "\n")
for file in os.listdir():
    if file.endswith(".vcf") and "Filt" not in file:
        print(file)
        o = open(file)
        for line in o:
            if "#" not in line:
                c += 1
                l = line.replace("\n", "").split("\t")
                if "," in l[3] or "," in l[4]:
                    print(l[0:6])
                vep_list = []
                vep = l[7].split(",")
                for j in vep:
                    jl = j.split("|")
                    #This try prevents errors due to premature stop codons which annoyingly include comma separated things
                    type_var = jl[1]
                    gene = jl[3]
                    gene_type = jl[7]
                    prot_pos = jl[14]
                    aas = jl[15]
                    isoform = jl[26]
                    sift = jl[27]
                    polyphen = jl[28]
                    vep_list = [type_var, gene, gene_type, prot_pos, aas, isoform, sift, polyphen]
                    out.write("\t".join([l[0].replace("chrchr", "chr"), l[1], l[2], l[3], l[4]] + vep_list) + "\n")
                if c % 100000 == 0:
                    print(c)
        o.close()
