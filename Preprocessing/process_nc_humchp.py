o = open("HumChp_NC.bed")
out = open("HumChp_NC_Ready.bed", 'w')
out.write("Position\tNearestGene\tHum|Chp|Gor\tMskHumPhyloP\tMskHBCGOPhyloP\tSpeciesSupport\t470MammalPhastCons\n")
for line in o:
    row = line.replace("\n", "").split("\t")
    out.write("\t".join([row[0] + ":" + row[2] + "-" + row[2], row[22], row[23], row[24], row[25], row[26], row[10]]) + "\n")
o.close()
out.close()
