import sys

file = sys.argv[1]
assert(".txt" in file)
o = open(file)
out = open(file.replace(".txt", ".bed"), 'w')

for line in o:
    if not line.startswith("Posi"):
        l = line.replace("\n", "").split("\t")
        chrom = l[0].split(":")[0]
        pos1 = str(int(l[0].split(":")[1].split("-")[0]) - 1)
        pos2 = l[0].split(":")[1].split("-")[0]
        out.write("\t".join([chrom, pos1, pos2] + l[1:]) + "\n")
o.close()
out.close()
