#Get gorilla and chimpanzee SNV information from MAF with human

wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsGorGor6/hg38.gorGor6.synNet.maf.gz
gunzip hg38.gorGor6.synNet.maf.gz
Convert to a bed format and sort
python callSNPsFromMAF.py hg38.gorGor6.synNet.maf hg38.gorGor6.synNet.txt
python make_bed_after_maf.py hg38.gorGor6.synNet.txt hg38.gorGor6.synNet.bed
sort -k1,1 -k2,2n hg38.gorGor6.synNet.bed > hg38.gorGor6.synNet.sorted.bed
sort -k1,1 -k2,2n HumChp_RemovedBy3WGS_Final_VEP.bed > HumChp_RemovedBy3WGS_Final_VEP.sort.bed

wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsPanTro6/hg38.panTro6.synNet.maf.gz
gunzip hg38.panTro6.synNet.maf.gz
Convert to a bed format and sort
python callSNPsFromMAF.py hg38.panTro6.synNet.maf hg38.panTro6.synNet.txt
python make_bed_after_maf.py hg38.panTro6.synNet.txt hg38.panTro6.synNet.bed
sort -k1,1 -k2,2n hg38.panTro6.synNet.bed > hg38.panTro6.synNet.sorted.bed

#Get CDS, 3UTR, 5UTR annotations
wget https://ftp.ensembl.org/pub/current/gtf/homo_sapiens/Homo_sapiens.GRCh38.114.chr.gtf.gz
gunzip *.chr.gtf.gz
python get_utr_cds.py

#Subtract CDS from 3', CDS and 3' from 5'
bedtools subtract -a Homo_sapiens.GRCh38.3UTR.bed -b Homo_sapiens.GRCh38.CDS.bed > Homo_sapiens.GRCh38.3UTR.Sub.bed
bedtools subtract -a Homo_sapiens.GRCh38.5UTR.bed -b Homo_sapiens.GRCh38.CDS.bed > Homo_sapiens.GRCh38.5UTR.SubPre.bed
bedtools subtract -a Homo_sapiens.GRCh38.5UTR.SubPre.bed -b Homo_sapiens.GRCh38.3UTR.bed > Homo_sapiens.GRCh38.5UTR.Sub.bed

### Go through and merge overlapping features that are assigned to the same genes ###

#Get unique
sort -u Homo_sapiens.GRCh38.CDS.bed > Homo_sapiens.GRCh38.CDS.dedup.bed

#Need to merge CDS regions that partially overlap and are assigned to the same gene name
sort -k4,4 Homo_sapiens.GRCh38.CDS.dedup.bed > Homo_sapiens.GRCh38.CDS.sortgene.bed

#Split by gene
mkdir Dedup_CDS
python pergene_cds.py Homo_sapiens.GRCh38.CDS.sortgene.bed Dedup_CDS

cd Dedup_CDS

#Merge bedfiles for each gene separately
for file in *.bed;
do
    sort -k1,1 -k2,2n $file > ${file::-4}.sort.bed
    bedtools merge -c 4 -o distinct -i ${file::-4}.sort.bed > ${file::-4}.merged.bed
done

#Make bed file with all genes
cat *.merged.bed > ../Homo_sapiens.GRCh38.CDS.fulldedup.bed

cd ..

#Sort and get ready to intersect
sort -k1,1 -k2,2n Homo_sapiens.GRCh38.CDS.fulldedup.bed > Homo_sapiens.GRCh38.CDS.sort.bed

sort -u Homo_sapiens.GRCh38.3UTR.Sub.bed > Homo_sapiens.GRCh38.3UTR.Sub.dedup.bed

#Need to merge 3UTR.Sub regions that partially overlap and are assigned to the same gene name
sort -k4,4 Homo_sapiens.GRCh38.3UTR.Sub.dedup.bed > Homo_sapiens.GRCh38.3UTR.Sub.sortgene.bed

#Split by gene
mkdir Dedup_3UTR.Sub
python pergene_cds.py Homo_sapiens.GRCh38.3UTR.Sub.sortgene.bed Dedup_3UTR.Sub

cd Dedup_3UTR.Sub

#Merge bedfiles for each gene separately
for file in *.bed;
do
    sort -k1,1 -k2,2n $file > ${file::-4}.sort.bed
    bedtools merge -c 4 -o distinct -i ${file::-4}.sort.bed > ${file::-4}.merged.bed
done

#Make bed file with all genes
cat *.merged.bed > ../Homo_sapiens.GRCh38.3UTR.Sub.fulldedup.bed

cd ..

#Sort and get ready to intersect
sort -k1,1 -k2,2n Homo_sapiens.GRCh38.3UTR.Sub.fulldedup.bed > Homo_sapiens.GRCh38.3UTR.Sub.sort.bed

sort -u Homo_sapiens.GRCh38.5UTR.Sub.bed > Homo_sapiens.GRCh38.5UTR.Sub.dedup.bed

#Need to merge 5UTR.Sub regions that partially overlap and are assigned to the same gene name
sort -k4,4 Homo_sapiens.GRCh38.5UTR.Sub.dedup.bed > Homo_sapiens.GRCh38.5UTR.Sub.sortgene.bed

#Split by gene
mkdir Dedup_5UTR.Sub
python pergene_cds.py Homo_sapiens.GRCh38.5UTR.Sub.sortgene.bed Dedup_5UTR.Sub

cd Dedup_5UTR.Sub

#Merge bedfiles for each gene separately
for file in *.bed;
do
    sort -k1,1 -k2,2n $file > ${file::-4}.sort.bed
    bedtools merge -c 4 -o distinct -i ${file::-4}.sort.bed > ${file::-4}.merged.bed
done

#Make bed file with all genes
cat *.merged.bed > ../Homo_sapiens.GRCh38.5UTR.Sub.fulldedup.bed

cd ..

#Sort and get ready to intersect
sort -k1,1 -k2,2n Homo_sapiens.GRCh38.5UTR.Sub.fulldedup.bed > Homo_sapiens.GRCh38.5UTR.Sub.sort.bed
