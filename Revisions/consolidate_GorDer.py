import sys

file = sys.argv[1]
o = open(file)
out = open(file.replace(".txt", "Consol.txt"), 'w')
prev = 0
c = 0
k = []
prev_l = 0
for line in o:
    l = line.replace("\n", "").split("\t")
    ident = l[0] + l[1] + l[3] + l[4]
    if c:
        if prev == ident or not prev:
            k.append(";".join(l))
            prev = ident
            prev_l = l
        else:
            k = list(set(k))
            genes = ""
            var_types = ""
            tran_types = ""
            poss = ""
            changes = ""
            appriss = ""
            SIFTs = ""
            PolyPhens = ""
            for llin in k:
                lin = llin.split(";")
                genes = genes + lin[6] + ";"
                var_types = var_types + lin[5] + ";"
                tran_types = tran_types + lin[7] + ";"
                poss = poss + lin[8] + ";"
                changes = changes + lin[9] + ";"
                appriss = appriss + lin[10] + ";"
                SIFTs = SIFTs + lin[11] + ";"
                PolyPhens = PolyPhens + lin[12] + ";"
            out.write("\t".join(["chr" + prev_l[0] + ":" + prev_l[1] + "-" + prev_l[1]] + prev_l[3:5] + [var_types, genes, tran_types, poss, changes, appriss, SIFTs, PolyPhens]) + "\n")
            k = [";".join(l)]
            prev = ident
            prev_l = l
    else:
        pass
    c += 1
out.close()
o.close()
    