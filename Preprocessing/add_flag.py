o = open("HumChp_KeptBy3WGS_Final_VEP.bed")
out = open("HumChp_KeptBy3WGS_Final_VEP_Flag.bed", 'w')

for line in o:
    l = line.replace("\n", "").split("\t")
    out.write("\t".join(l[0:3] + [l[-1], "Y"]) + "\n")
o.close()
out.close()

o = open("HumChp_RemovedBy3WGS_Final_VEP.bed")
out = open("HumChp_RemovedBy3WGS_Final_VEP_Flag.bed", 'w')

for line in o:
    l = line.replace("\n", "").split("\t")
    out.write("\t".join(l[0:3] + [l[-1], "N"]) + "\n")
o.close()
out.close()
