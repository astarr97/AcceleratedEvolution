o = open("All.MskHumChpBon.PhyloP.SpecSup250.Gene.Tri.InputPrelim.Anno.bed")
out1 = open("Background_AccelEvolInput.NC.Final.txt", 'w')
out2 = open("Background_AccelEvolInput.3UTR.Final.txt", 'w')
out3 = open("Background_AccelEvolInput.5UTR.Final.txt", 'w')

#out1.write("Position\tPhyloP447\tSpecSup447\tNearestGene\tNearestDist\n")
#out2.write("Position\tPhyloP447\tSpecSup447\tNearestGene\tNearestDist\n")
#out3.write("Position\tPhyloP447\tSpecSup447\tNearestGene\tNearestDist\n")

for line in o:
    l = line.split("\t")
    if l[-2] == "3UTR":
        out2.write("\t".join([l[0] + ":" + l[2], l[3], l[4], l[-3], str(0), l[7]]) + "\n")
    elif l[-2] == "5UTR":
        out3.write("\t".join([l[0] + ":" + l[2], l[3], l[4], l[-3], str(0), l[7]]) + "\n")
    elif l[-2] == ".":
        out1.write("\t".join([l[0] + ":" + l[2], l[3], l[4], l[5], l[6], l[7]]) + "\n")
    
out1.close()
out2.close()
out3.close()
