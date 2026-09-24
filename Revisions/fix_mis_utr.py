# File containing the correct gene assignments
fix_file = "GorDer_Mis_UTR_ForFix.bed"

# File you want to modify
input_file = "GorDer_AllSubs_Input_ToFixUTRMissense.bed"

# Output
output_file = "GorDer_AllSubs_Input_Final.bed"


# Build lookup: (chromosome, position) -> gene
gene_lookup = {}

with open(fix_file) as f:
    for line in f:
        if not line.strip():
            continue

        l = line.rstrip("\n").split("\t")

        chrom = l[0]
        pos = l[2]
        gene = l[6]

        gene_lookup[(chrom, pos)] = gene

print(list(gene_lookup.keys())[0:10], list(gene_lookup.values())[0:10])
# Modify the original file
with open(input_file) as f, open(output_file, "w") as out:

    for line in f:
        if not line.strip():
            continue

        l = line.rstrip("\n").split("\t")

        chrom = l[0]
        pos = l[1]

        key = (chrom, pos)

        if key in gene_lookup:
            #print(key, gene_lookup[key])
            # Replace second-to-last column with correct gene
            l[-2] = gene_lookup[key]

            # Replace last column with distance = 0
            l[-1] = "0"

        out.write("\t".join(l).replace("3UTR_Syn", "3UTR").replace("5UTR_Syn", "5UTR").replace("Mis_Syn", "Mis") + "\n")
