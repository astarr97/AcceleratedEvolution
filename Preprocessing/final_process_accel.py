o = open("HumChp_AccelEvolInput.Gene.bed")
out = open("HumChp_AccelEvolInput.Final.txt", 'w')

out.write("Position\tHum|Chp|Gor\tDerived\tPhyloP447\tSpecSup447\tPhastCons447\tCategory\tKeptAfterFilt\tAncTrinuc\tDerTrinuc\tMutCat\tNearestGene\tNearestDist\n")

for line in o:
    l = line.split("\t")
    out.write("\t".join([l[0] + ":" + l[2]] + l[3:]))

o.close()
out.close()
