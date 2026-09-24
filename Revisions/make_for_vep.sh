mkdir -p Gorilla_Derived

awk '
BEGIN {
    OFS="\t";
    batch_size=600000;
}
{
    batch=int((NR-1)/batch_size)+1;
    outfile=sprintf("Gorilla_Derived/vep_batch_%03d.vcf", batch);

    if ((NR-1) % batch_size == 0) {
        print "##fileformat=VCFv4.2" > outfile;
        print "#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO" > outfile;
    }

    chrom=$1;
    sub(/^chr/, "", chrom);

    split($5, a, "|");

    pos=$2+1;
    ref=a[1];
    alt=a[2];

    print chrom, pos, ".", ref, alt, ".", ".", "." >> outfile;
}' hg38.panTro6.gorGor6.ponAbe3.GorillaDerived.bed
