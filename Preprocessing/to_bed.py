import sys

file = sys.argv[1]

assert(file.endswith(".txt"))

o = open(file)
out = open(file.replace(".txt", ".bed"), 'w')

for line in o:
    l = line.split("\t")
    out.write("\t".join([l[0].split(":")[0], str(int(l[0].split(":")[1])-1), l[0].split(":")[1]] + l[1:]))
out.close()
