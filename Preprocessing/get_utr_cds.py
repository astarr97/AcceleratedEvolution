o = open("Homo_sapiens.GRCh38.114.chr.gtf")

out1 = open("Homo_sapiens.GRCh38.CDS.bed", 'w')
out2 = open("Homo_sapiens.GRCh38.3UTR.bed", 'w')
out3 = open("Homo_sapiens.GRCh38.5UTR.bed", 'w')

for line in o:
    if not line.startswith("#"):
        l = line.split("\t")
        if 'transcript_biotype "nonsense_mediated_decay"' not in l[8]:
            if l[2] == "CDS":
                try:
                    gene = l[8].split('; gene_name "')[1].split('";')[0]
                except:
                    try:
                        gene = l[8].split('gene_id "')[1].split('";')[0]
                    except:
                        print(l)
                out1.write("\t".join(["chr" + l[0], str(int(l[3])), l[4], gene]) + "\n")
            elif l[2] == "five_prime_utr":
                try:
                    gene = l[8].split('; gene_name "')[1].split('";')[0]
                except:
                    gene = l[8].split('gene_id "')[1].split('";')[0]
                out3.write("\t".join(["chr" + l[0], str(int(l[3])), l[4], gene]) + "\n")
            elif l[2] == "three_prime_utr":
                try:
                    gene = l[8].split('; gene_name "')[1].split('";')[0]
                except:
                    gene = l[8].split('gene_id "')[1].split('";')[0]
                out2.write("\t".join(["chr" + l[0], str(int(l[3])), l[4], gene]) + "\n")

o.close()
out1.close()
out2.close()
out3.close()
