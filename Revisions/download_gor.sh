#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=96:00:00
#SBATCH --mem=32GB

#Gorilla using hg38/gorgor6
wget -P gorilla ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsGorGor6/hg38.gorGor6.synNet.maf.gz
gunzip gorilla/hg38.gorGor6.synNet.maf.gz
python callSNPsFromMAF.py gorilla/hg38.gorGor6.synNet.maf hg38.gorGor6.synNet.txt
python make_bed.py hg38.gorGor6.synNet.txt
sort -k1,1 -k2,2n hg38.gorGor6.synNet.bed > hg38.gorGor6.synNet.sort.bed
