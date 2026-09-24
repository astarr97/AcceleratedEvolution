#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH -p hbfraser
#SBATCH --mem=64GB


#Summarize the VEP
python process_ensembl_gorder.py

#Do so for the sites that were not filtered out
python consolidate_GorDer.py All_Summarized_GorDer_Test_GW_Ens.txt

python back_to_bed.py All_Summarized_GorDer_Test_GW_EnsConsol.txt

python split_by_vep_gorder.py

sort -k1,1 -k2,2n GorDer_AllSubs_Final_VEP.bed > GorDer_AllSubs_Final_VEP.sort.bed

awk '$1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/' /scratch/users/astarr97/PhyloP/Recomputing/Homo_sapiens_PhyloP_MaskGreatApes/All/All_Homo_sapiens_PhyloP_MaskGreatApes.sort.NoID.bed > All_Homo_sapiens_PhyloP_MaskGreatApes.sort.NoID.Filt.bed
awk '$1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/' /scratch/users/astarr97/PhyloP/Recomputing/Homo_sapiens_PhyloP_MaskGreatApes/All/All_Homo_sapiens_PhyloP_MaskGreatApes.SpecSup.sort.bed > All_Homo_sapiens_PhyloP_MaskGreatApes.SpecSup.sort.Filt.bed

bedtools intersect -sorted -wao -a GorDer_AllSubs_Final_VEP.sort.bed -b All_Homo_sapiens_PhyloP_MaskGreatApes.sort.NoID.Filt.bed > GorDer_AllSubs_Final_VEP.PhyloP.bed
bedtools intersect -sorted -wao -a GorDer_AllSubs_Final_VEP.PhyloP.bed -b All_Homo_sapiens_PhyloP_MaskGreatApes.SpecSup.sort.Filt.bed > GorDer_AllSubs_Final_VEP.PhyloP.SpecSup.ToFilt.bed

awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$10,$15}' GorDer_AllSubs_Final_VEP.PhyloP.SpecSup.ToFilt.bed > GorDer_AllSubs_Final_VEP.PhyloP.SpecSup.bed

bedtools closest -d -wao -a GorDer_AllSubs_Final_VEP.PhyloP.SpecSup.bed -b Human_Promoters_Ortho_Sorted_hg38.sort.bed > GorDer_AllSubs_Final_VEP.PhyloP.SpecSup.NearestGene.bed
awk 'BEGIN{OFS="\t"} {print $1,$3,$4,$5,$6,$7,$8,$12,$13}' GorDer_AllSubs_Final_VEP.PhyloP.SpecSup.NearestGene.bed > GorDer_AllSubs_Input_ToFixUTRMissense.bed

#Fix missense and UTR to be properly assigned to genes
python fix_mis_utr.py
python get_mis_utr_gene.py
