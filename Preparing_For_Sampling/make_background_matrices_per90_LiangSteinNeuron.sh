#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH -p hbfraser
#SBATCH --mem=128GB

#Get the background sites
awk -F'\t' 'BEGIN{OFS="\t"} { split($1,a,":"); print a[1], a[2]-1, a[2], $2, $3, $4, $5, $6}' Background_AccelEvolInput.NC.Final.sort.txt > Background_AccelEvolInput.NC.Final.bed
sort -k1,1 -k2,2n Background_AccelEvolInput.NC.Final.bed > Background_AccelEvolInput.NC.Final.sort.bed

#Prepare to create our 90th percentile background file for LiangSteinNeuron
python prepare_per90_background.py HumChp_AccelEvolInput.Final.txt 1 -100 0 0 0 0 90 NC LiangSteinNeuron/BackgroundMatrix_LiangSteinNeuron_AnyLogfc_NC_FiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per90 LiangSteinNeuron

#Sort the regions that were created
sort -k1,1 -k2,2n LiangSteinNeuron/BackgroundMatrix_LiangSteinNeuron_AnyLogfc_NC_FiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per90/Per90_Regions.bed > LiangSteinNeuron/BackgroundMatrix_LiangSteinNeuron_AnyLogfc_NC_FiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per90/Per90_Regions.sort.bed

#Merge overlapping regions
bedtools merge -i LiangSteinNeuron/BackgroundMatrix_LiangSteinNeuron_AnyLogfc_NC_FiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per90/Per90_Regions.sort.bed > LiangSteinNeuron/BackgroundMatrix_LiangSteinNeuron_AnyLogfc_NC_FiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per90/Per90_Regions.sort.merged.bed

#Intersect our background input file with the regions so restrict only to qualifying regions
bedtools intersect -sorted -a Background_AccelEvolInput.NC.Final.sort.bed -b LiangSteinNeuron/BackgroundMatrix_LiangSteinNeuron_AnyLogfc_NC_FiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per90/Per90_Regions.sort.merged.bed > LiangSteinNeuron/BackgroundMatrix_LiangSteinNeuron_AnyLogfc_NC_FiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per90/LiangSteinNeuron_Background_AccelEvolInput.NC.Final.sort.bed

#Filter columns to reformat
awk -F'\t' 'BEGIN{OFS="\t"} { print $1":"$3, $4, $5, $6, $7, $8 }' LiangSteinNeuron/BackgroundMatrix_LiangSteinNeuron_AnyLogfc_NC_FiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per90/LiangSteinNeuron_Background_AccelEvolInput.NC.Final.sort.bed > LiangSteinNeuron/BackgroundMatrix_LiangSteinNeuron_AnyLogfc_NC_FiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per90/LiangSteinNeuron_Background_AccelEvolInput.NC.Final.sort.txt

#Sort by gene name
sort -k4,4 LiangSteinNeuron/BackgroundMatrix_LiangSteinNeuron_AnyLogfc_NC_FiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per90/LiangSteinNeuron_Background_AccelEvolInput.NC.Final.sort.txt > LiangSteinNeuron/BackgroundMatrix_LiangSteinNeuron_AnyLogfc_NC_FiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per90/LiangSteinNeuron_Background_AccelEvolInput.NC.Final.use.txt

#Input file that was created
use_file="LiangSteinNeuron/BackgroundMatrix_LiangSteinNeuron_AnyLogfc_NC_FiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per90/LiangSteinNeuron_Background_AccelEvolInput.NC.Final.use.txt"

#Make background matrices for all sites regardless of PhyloP score
python make_background_matrices.py $use_file NC 1 -100 0 0 LiangSteinNeuron/BackgroundMatrix_FiltWGS_LiangSteinNeuron_NC_PhyloP-100_SpecSup0_NCM_Per90 HumChp_AccelEvolInput.Final.txt
python make_background_matrices.py $use_file NC 0 -100 0 0 LiangSteinNeuron/BackgroundMatrix_NoFiltWGS_LiangSteinNeuron_NC_PhyloP-100_SpecSup0_NCM_Per90 HumChp_AccelEvolInput.Final.txt
python make_background_matrices.py $use_file NC 1 -100 250 0 LiangSteinNeuron/BackgroundMatrix_FiltWGS_LiangSteinNeuron_NC_PhyloP-100_SpecSup250_NCM_Per90 HumChp_AccelEvolInput.Final.txt
python make_background_matrices.py $use_file NC 0 -100 250 0 LiangSteinNeuron/BackgroundMatrix_NoFiltWGS_LiangSteinNeuron_NC_PhyloP-100_SpecSup250_NCM_Per90 HumChp_AccelEvolInput.Final.txt

#Make background matrices only restricting to sites with PhyloP > 1
python make_background_matrices.py $use_file NC 1 1 250 0 LiangSteinNeuron/BackgroundMatrix_FiltWGS_LiangSteinNeuron_NC_PhyloP1_SpecSup250_NCM_Per90 HumChp_AccelEvolInput.Final.txt
python make_background_matrices.py $use_file NC 0 1 250 0 LiangSteinNeuron/BackgroundMatrix_NoFiltWGS_LiangSteinNeuron_NC_PhyloP1_SpecSup250_NCM_Per90 HumChp_AccelEvolInput.Final.txt
