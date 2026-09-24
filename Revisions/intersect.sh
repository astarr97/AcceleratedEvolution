#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=72:00:00
#SBATCH --mem=128GB

bedtools intersect -sorted -a hg38.panTro6.synNet.sort.bed -b hg38.gorGor6.synNet.sort.bed -wa -wb > hg38.panTro6.gorGor6.prelim.bed
cut -f1-4,8 hg38.panTro6.gorGor6.prelim.bed > hg38.panTro6.gorGor6.bed
bedtools intersect -sorted -wa -wb -a hg38.panTro6.gorGor6.bed -b hg38.ponAbe3.synNet.sort.bed > hg38.panTro6.gorGor6.ponAbe3.prelim.bed
cut -f1-5,9 hg38.panTro6.gorGor6.ponAbe3.prelim.bed > hg38.panTro6.gorGor6.ponAbe3.bed

awk -F'\t' 'BEGIN {OFS="\t"} { split($4, a, "|"); split($5, b, "|"); split($6, c, "|"); if (a[1] != a[2] && a[2] == b[2] && a[1] == c[2]) print }' hg38.panTro6.gorGor6.ponAbe3.bed > hg38.panTro6.gorGor6.ponAbe3.CGHsites.bed
awk -F'\t' 'BEGIN {OFS="\t"} { split($4, a, "|"); split($5, b, "|"); split($6, c, "|"); if (a[1] != a[2] && a[1] == b[2] && a[2] == c[2]) print }' hg38.panTro6.gorGor6.ponAbe3.bed > hg38.panTro6.gorGor6.ponAbe3.HGCsites.bed
awk -F'\t' 'BEGIN {OFS="\t"} { split($4, a, "|"); split($5, b, "|"); split($6, c, "|"); if (a[1] == a[2] && a[1] != b[2] && a[1] == c[2]) print }' hg38.panTro6.gorGor6.ponAbe3.bed > hg38.panTro6.gorGor6.ponAbe3.GorillaDerived.bed
