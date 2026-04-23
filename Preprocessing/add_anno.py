o = open("Homo_sapiens.GRCh38.CDS.sort.bed")
out = open("Homo_sapiens.GRCh38.CDS.Anno.bed", 'w')

for line in o:
    l = line.replace("\n", "").split("\t")
    out.write("\t".join(l + ["CDS"]) + "\n")
o.close()
out.close()

o = open("Homo_sapiens.GRCh38.3UTR.Sub.sort.bed")
out = open("Homo_sapiens.GRCh38.3UTR.Sub.Anno.bed", 'w')

for line in o:
    l = line.replace("\n", "").split("\t")
    out.write("\t".join(l + ["3UTR"]) + "\n")
o.close()
out.close()

o = open("Homo_sapiens.GRCh38.5UTR.Sub.sort.bed")
out = open("Homo_sapiens.GRCh38.5UTR.Sub.Anno.bed", 'w')

for line in o:
    l = line.replace("\n", "").split("\t")
    out.write("\t".join(l + ["5UTR"]) + "\n")
o.close()
out.close()

