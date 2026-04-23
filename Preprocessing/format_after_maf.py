import sys

file = sys.argv[1]
assert(".bed" in file)
o = open(file)
out = open(file.replace(".bed", "_ToJoin.txt"), 'w')

for line in o:
    l = line.replace("\n", "").split('\t')
    pos = l[0] + ":" + l[2] + "-" + l[2]
    if l[7] != ".":
        new_snp = l[3] + "|" + l[7].split("|")[1]
        out.write("\t".join([pos, new_snp]) + "\n")
o.close()
out.close()
