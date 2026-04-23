#!/bin/bash
#SBATCH --time=168:00:00
#SBATCH -p hbfraser
#SBATCH --mem=16GB

bedtools intersect -sorted -wao -a RIMKLB_3UTR.bed -b ../AccelEvol_New/All.MskHumChpBon.PhyloP.SpecSup250.bed > RIMKLB_3UTR_PhyloP447.bed

grep 'transcript_name "RET-202"' /oak/stanford/groups/hbfraser/astarr/human.gtf | grep 'gene_name "RET"' | grep $'HAVANA\tCDS' | awk -F'\t' -v OFS='\t' '{ print $1, $4 - 1, $5, "RET_CDS" }' > RET_CDS_PhyloP447.bed
grep 'transcript_name "DMRT3-201"' /oak/stanford/groups/hbfraser/astarr/human.gtf | grep 'gene_name "DMRT3"' | grep $'HAVANA\tCDS' | awk -F'\t' -v OFS='\t' '{ print $1, $4 - 1, $5, "DMRT3_CDS" }' > DMRT3_CDS_PhyloP447.bed
grep 'transcript_name "ESR1-201"' /oak/stanford/groups/hbfraser/astarr/human.gtf | grep 'gene_name "ESR1"' | grep $'HAVANA\tCDS' | awk -F'\t' -v OFS='\t' '{ print $1, $4 - 1, $5, "ESR1_CDS" }' > ESR1_CDS_PhyloP447.bed
grep 'transcript_name "EYA1-202"' /oak/stanford/groups/hbfraser/astarr/human.gtf | grep 'gene_name "EYA1"' | grep $'HAVANA\tCDS' | awk -F'\t' -v OFS='\t' '{ print $1, $4 - 1, $5, "EYA1_CDS" }' > EYA1_CDS_PhyloP447.bed
grep 'transcript_name "ADCYAP1-202"' /oak/stanford/groups/hbfraser/astarr/human.gtf | grep 'gene_name "ADCYAP1"' | grep $'HAVANA\tCDS' | awk -F'\t' -v OFS='\t' '{ print $1, $4 - 1, $5, "ADCYAP1_CDS" }' > ADCYAP1_CDS_PhyloP447.bed
cat *CDS_PhyloP447.bed | sort -k1,1 -k2,2n > CDS_PhyloP447.sort.bed

bedtools intersect -sorted -wao -a CDS_PhyloP447.sort.bed -b ../AccelEvol_New/All.MskHumChpBon.PhyloP.SpecSup250.bed > CDS_PhyloP447_Output.bed

grep 'transcript_name "ZEB2-226"' /oak/stanford/groups/hbfraser/astarr/human.gtf | grep 'gene_name "ZEB2"' | awk -F'\t' -v OFS='\t' '{ print $1, $4 - 1, $5, $3, "ZEB2" }' > ZEB2_regions.bed
grep 'transcript_name "LMO4-202"' /oak/stanford/groups/hbfraser/astarr/human.gtf | grep 'gene_name "LMO4"' | awk -F'\t' -v OFS='\t' '{ print $1, $4 - 1, $5, $3, "LMO4" }' > LMO4_regions.bed
grep 'transcript_name "RBFOX1-201"' /oak/stanford/groups/hbfraser/astarr/human.gtf | grep 'gene_name "RBFOX1"' | awk -F'\t' -v OFS='\t' '{ print $1, $4 - 1, $5, $3, "RBFOX1" }' > RBFOX1_regions.bed
grep 'transcript_name "CDH13-209"' /oak/stanford/groups/hbfraser/astarr/human.gtf | grep 'gene_name "CDH13"' | awk -F'\t' -v OFS='\t' '{ print $1, $4 - 1, $5, $3, "CDH13" }' > CDH13_regions.bed
grep 'transcript_name "CXXC5-201"' /oak/stanford/groups/hbfraser/astarr/human.gtf | grep 'gene_name "CXXC5"' | awk -F'\t' -v OFS='\t' '{ print $1, $4 - 1, $5, $3, "CXXC5" }' > CXXC5_regions.bed
grep 'transcript_name "FOXP1-254"' /oak/stanford/groups/hbfraser/astarr/human.gtf | grep 'gene_name "FOXP1"' | awk -F'\t' -v OFS='\t' '{ print $1, $4 - 1, $5, $3, "FOXP1" }' > FOXP1_regions.bed
grep 'transcript_name "CDH12-201"' /oak/stanford/groups/hbfraser/astarr/human.gtf | grep 'gene_name "CDH12"' | awk -F'\t' -v OFS='\t' '{ print $1, $4 - 1, $5, $3, "CDH12" }' > CDH12_regions.bed
grep 'transcript_name "ZNF521-201"' /oak/stanford/groups/hbfraser/astarr/human.gtf | grep 'gene_name "ZNF521"' | awk -F'\t' -v OFS='\t' '{ print $1, $4 - 1, $5, $3, "ZNF521" }' > ZNF521_regions.bed
grep 'transcript_name "SIX3-201"' /oak/stanford/groups/hbfraser/astarr/human.gtf | grep 'gene_name "SIX3"' | awk -F'\t' -v OFS='\t' '{ print $1, $4 - 1, $5, $3, "SIX3" }' > SIX3_regions.bed
grep 'transcript_name "TCF7L2-206"' /oak/stanford/groups/hbfraser/astarr/human.gtf | grep 'gene_name "TCF7L2"' | awk -F'\t' -v OFS='\t' '{ print $1, $4 - 1, $5, $3, "TCF7L2" }' > TCF7L2_regions.bed
grep 'transcript_name "MEOX2-201"' /oak/stanford/groups/hbfraser/astarr/human.gtf | grep 'gene_name "MEOX2"' | awk -F'\t' -v OFS='\t' '{ print $1, $4 - 1, $5, $3, "MEOX2" }' > MEOX2_regions.bed
grep 'transcript_name "DPP10-205"' /oak/stanford/groups/hbfraser/astarr/human.gtf | grep 'gene_name "DPP10"' | awk -F'\t' -v OFS='\t' '{ print $1, $4 - 1, $5, $3, "DPP10" }' > DPP10_regions.bed

cat *regions.bed > All_NC_Regions.bed

sort -k1,1 -k2,2n CREs_To_Plot.bed > CREs_To_Plot.sort.bed
bedtools intersect -sorted -wao -a CREs_To_Plot.sort.bed -b ../AccelEvol_New/All.MskHumChpBon.PhyloP.SpecSup250.bed > CREs_To_Plot_PhyloP447.bed
