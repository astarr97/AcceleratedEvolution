#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=96:00:00

#Orangutan using hg38/ponAbe3
wget -P orangutan ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsPonAbe3/hg38.ponAbe3.synNet.maf.gz
gunzip orangutan/hg38.ponAbe3.synNet.maf.gz
python callSNPsFromMAF.py orangutan/hg38.ponAbe3.synNet.maf hg38.ponAbe3.synNet.txt
python make_bed.py hg38.ponAbe3.synNet.txt
sort -k1,1 -k2,2n hg38.ponAbe3.synNet.bed > hg38.ponAbe3.synNet.sort.bed
