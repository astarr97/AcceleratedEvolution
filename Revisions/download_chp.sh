#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=96:00:00


#Chimp using UCSC version and hg38/PanTro6
wget -P chimp ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsPanTro6/hg38.panTro6.synNet.maf.gz
gunzip chimp/hg38.panTro6.synNet.maf.gz
python callSNPsFromMAF.py chimp/hg38.panTro6.synNet.maf hg38.panTro6.synNet.txt
python make_bed.py hg38.panTro6.synNet.txt
sort -k1,1 -k2,2n hg38.panTro6.synNet.bed > hg38.panTro6.synNet.sort.bed
