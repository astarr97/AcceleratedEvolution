#o = open("All.MskHumChpBon.PhyloP.SpecSup250.Gene.Tri.PanTro6.gorGor6.bed")
#out = open("All.MskHumChpBon.PhyloP.SpecSup250.Gene.Tri.InputPrelim.bed", 'w')

#for line in o:
#    l = line.split("\t")
#    if l[4] != "." and l[5] != "." and l[7] != "." and l[11] != "." and l[16] != ".":
#        if l[11].split("|")[0] == l[11].split("|")[1]:
#            out.write("\t".join([l[0], l[1], l[2], l[3], l[4], l[5], l[6], l[7].upper()]) + "\n")
#        else:
#            if l[11].split("|")[0] == l[16].split("|")[1]:
#                out.write("\t".join([l[0], l[1], l[2], l[3], l[4], l[5], l[6], l[7].upper()]) + "\n")
#            else:
#                l7 = list(l[7])
#                l7[1] = l[11].split("|")[1]
#                l7 = "".join(l7)
#                out.write("\t".join([l[0], l[1], l[2], l[3], l[4], l[5], l[6], l7.upper()]) + "\n")
                
o = open("All.MskHumChpBon.PhyloP.SpecSup250.Gene.Tri.PanTro6.gorGor6.bed")
out = open("All.MskHumChpBon.PhyloP.SpecSup250.Gene.Tri.InputPrelimNoGorFilt.bed", 'w')

for line in o:
    l = line.split("\t")
    if l[4] != "." and l[5] != "." and l[7] != "." and l[11] != ".":
        if l[11].split("|")[0] == l[11].split("|")[1]:
            out.write("\t".join([l[0], l[1], l[2], l[3], l[4], l[5], l[6], l[7].upper()]) + "\n")
        else:
            if l[16] == ".":
                out.write("\t".join([l[0], l[1], l[2], l[3], l[4], l[5], l[6], l[7].upper()]) + "\n")
            elif l[11].split("|")[0] == l[16].split("|")[1]:
                out.write("\t".join([l[0], l[1], l[2], l[3], l[4], l[5], l[6], l[7].upper()]) + "\n")
            else:
                l7 = list(l[7])
                l7[1] = l[11].split("|")[1]
                l7 = "".join(l7)
                out.write("\t".join([l[0], l[1], l[2], l[3], l[4], l[5], l[6], l7.upper()]) + "\n")
                
