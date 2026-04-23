#!/bin/bash
#SBATCH --time=168:00:00
#SBATCH -p hbfraser
#SBATCH --mem=32GB

#Getting PhyloP etc. information
zcat /oak/stanford/groups/hbfraser/astarr/PosSelect_AccelEvol/Recomputing447-way_Final/Scores/All.MskHumChpBon.PhastCons.scores.txt.gz > All.MskHumChpBon.PhastCons.scores.txt
zcat /oak/stanford/groups/hbfraser/astarr/PosSelect_AccelEvol/Recomputing447-way_Final/Scores/All.MskHumChpBon.PhyloP.sort.txt.gz > All.MskHumChpBon.PhyloP.sort.txt
zcat /oak/stanford/groups/hbfraser/astarr/PosSelect_AccelEvol/Recomputing447-way_Final/Scores/All.MskHumChpBon.SpecSup.txt.gz > All.MskHumChpBon.SpecSup.txt

python to_bed.py All.MskHumChpBon.PhastCons.scores.txt
python to_bed.py All.MskHumChpBon.PhyloP.sort.txt
python spec_sup_make_bed.py All.MskHumChpBon.SpecSup.txt

sort -k1,1 -k2,2n ASE_SNPs.FILTER.SPLIT_SPECIES.bed > ASE_SNPs.FILTER.SPLIT_SPECIES.sort.bed
sort -k1,1 -k2,2n human_referenced_chp_hum_snps.bed > human_referenced_chp_hum_snps.bed.sort.bed

sort -k1,1 -k2,2n All.MskHumChpBon.PhastCons.scores.bed > All.MskHumChpBon.PhastCons.scores.sort.bed
sort -k1,1 -k2,2n All.MskHumChpBon.PhyloP.sort.bed > All.MskHumChpBon.PhyloP.sort.sort.bed
sort -k1,1 -k2,2n All.MskHumChpBon.SpecSup.bed > All.MskHumChpBon.SpecSup.sort.bed

#Intersecting unfiltered set of SNPs
bedtools intersect -sorted -wao -a human_referenced_chp_hum_snps.bed.sort.bed -b All.MskHumChpBon.PhyloP.sort.sort.bed > HumChp_NoFilt_MskHCB.PhyloP.sort.bed
bedtools intersect -sorted -wao -a HumChp_NoFilt_MskHCB.PhyloP.sort.bed -b All.MskHumChpBon.SpecSup.sort.bed > HumChp_NoFilt_MskHCB.PhyloP.SpecSup.sort.bed
bedtools intersect -sorted -wao -a HumChp_NoFilt_MskHCB.PhyloP.SpecSup.sort.bed -b All.MskHumChpBon.PhastCons.scores.sort.bed > HumChp_NoFilt_MskHCB.PhyloP.SpecSup.PhastCons.sort.bed

#Subset
awk -F'\t' '{OFS="\t"; print $1, $2, $3, $4, $8, $13, $18}' HumChp_NoFilt_MskHCB.PhyloP.SpecSup.PhastCons.sort.bed > HumChp_NoFilt_MskHCB.PhyloP.SpecSup.PhastCons.filt.bed

#Intersect with gorilla information
bedtools intersect -sorted -wao -a HumChp_NoFilt_MskHCB.PhyloP.SpecSup.PhastCons.filt.bed -b hg38.gorGor6.synNet.sorted.bed > HumChp_NoFilt_MskHCB.PhyloP.SpecSup.PhastCons.Gor.bed

#Copy over old VEP files 
cp /oak/stanford/groups/hbfraser/astarr/PosSelect_AccelEvol/Human_Chimp_VEP_Output/*.gz ./

#Intersecting with VEP information
#Convert to VEP format for sites that we previously filtered out
python bed_to_vep.py

#Since I included every variant with both + and - strand entries but only one was is correct, we filter out the wrong ones
python filter_vep.py

#Summarize the VEP
python process_ensembl_human_chimp.py

#Consolidate the VEP annotations for the sites that were filtered
python consolidate_HumChp.py All_Summarized_HumChp.txt

#Do so for the sites that were not filtered out
python consolidate_HumChp.py All_Summarized_HumChp_GW_Ens.txt

#Convert back to bed
python back_to_bed.py All_Summarized_HumChpConsol.txt
python back_to_bed.py All_Summarized_HumChp_GW_EnsConsol.txt

#Split things up by VEP
python split_by_vep_humchp.py
python split_by_vep_humchp_filt.py

#Add flag indicating whether something was filtered or included after WGS filtering
python add_flag.py

#Cat and sort
cat HumChp_KeptBy3WGS_Final_VEP_Flag.bed HumChp_RemovedBy3WGS_Final_VEP_Flag.bed > HumChp_AllWGS_Final_VEP_Flag.bed
sort -k1,1 -k2,2n HumChp_AllWGS_Final_VEP_Flag.bed > HumChp_AllWGS_Final_VEP_Flag.sort.bed

#Intersect everything
bedtools intersect -sorted -wao -a HumChp_NoFilt_MskHCB.PhyloP.SpecSup.PhastCons.Gor.bed -b HumChp_AllWGS_Final_VEP_Flag.sort.bed > HumChp_NoFilt_MskHCB.AEInput.Prelim.bed

#Filter with awk command
bedtools intersect -sorted -wao -a HumChp_NoFilt_MskHCB.PhyloP.SpecSup.PhastCons.Gor.bed -b HumChp_AllWGS_Final_VEP_Flag.sort.bed > HumChp_NoFilt_MskHCB.AEInput.Prelim.bed

#Filter with awk command
awk -F'\t' '{OFS="\t"; print $1, $2, $3, $4, $11, $5, $6, $7, $16, $17}' HumChp_NoFilt_MskHCB.AEInput.Prelim.bed > HumChp_NoFilt_MskHCB.PhyloP.SpecSup.PhastCons.filt.bed

#Create contigs file
cp ../../human.fasta ./
samtools faidx human.fasta
awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' human.fasta.fai | sort -n -r -k 3,3 > human_contigs.bed

#Get the trinucleotide context
python PrepSim/expand_bed.py HumChp_NoFilt_MskHCB.PhyloP.SpecSup.PhastCons.filt.bed human_contigs.bed
bedtools getfasta -bedOut -fi human.fasta -bed HumChp_NoFilt_MskHCB.PhyloP.SpecSup.PhastCons.filt.exp.bed > HumChp_NoFilt_MskHCB.PhyloP.SpecSup.PhastCons.filt.tri.fastabed
paste HumChp_NoFilt_MskHCB.PhyloP.SpecSup.PhastCons.filt.bed HumChp_NoFilt_MskHCB.PhyloP.SpecSup.PhastCons.filt.tri.fastabed > HumChp_NoFilt_MskHCB.PhyloP.SpecSup.PhastCons.filt.TriPrelim.bed

awk -F'\t' '{OFS="\t"; print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $14}' HumChp_NoFilt_MskHCB.PhyloP.SpecSup.PhastCons.filt.TriPrelim.bed > HumChp_NoFilt_MskHCB.PhyloP.SpecSup.PhastCons.filt.Tri.bed

#Do some reformatting
python not_final_preprocess.py

#Find nearest gene promoter
cp ../../PosSelect_AccelEvol/Human_Promoters_Ortho_Sorted_hg38.bed ./
sort -k1,1 -k2,2n Human_Promoters_Ortho_Sorted_hg38.bed > Human_Promoters_Ortho_Sorted_hg38.sort.bed
bedtools closest -wao -d -a HumChp_AccelEvolInput.ForGene.bed -b Human_Promoters_Ortho_Sorted_hg38.sort.bed > HumChp_AccelEvolInput.WithGene.bed

#Reformat again
awk -F'\t' '{OFS="\t"; print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $17, $18}' HumChp_AccelEvolInput.WithGene.bed > HumChp_AccelEvolInput.Gene.bed

#Intersect with the proper gene annotations for UTR/Missense/Synonymous
bedtools intersect -sorted -wao -a HumChp_AccelEvolInput.Gene.bed -b Homo_sapiens.GRCh38.ToAnnoBack.sort.bed > HumChp_AccelEvolInput.ToFixGene.bed
#python final_process_accel.py

#New version fixes some annoying gene name issues
python fix_gene_names.py

#This is the head of the final input file
head HumChp_AccelEvolInput.Final.txt

#Pulling the same information for all sites genome-wide
bedtools intersect -sorted -wao -a All.MskHumChpBon.PhyloP.sort.sort.bed -b All.MskHumChpBon.SpecSup.sort.bed > All.MskHumChpBon.PhyloP.SpecSup250.bed
bedtools closest -wao -d -a All.MskHumChpBon.PhyloP.SpecSup250.bed -b Human_Promoters_Ortho_Sorted_hg38.sort.bed > All.MskHumChpBon.PhyloP.SpecSup250.Gene.bed
awk -F'\t' '{OFS="\t"; print $1, $2, $3, $4, $8, $13, $14}' All.MskHumChpBon.PhyloP.SpecSup250.Gene.bed > All.MskHumChpBon.PhyloP.SpecSup250.Gene.filt.bed

#Getting 3 nucleotide context for all sites
python PrepSim/expand_bed.py All.MskHumChpBon.PhyloP.SpecSup250.Gene.filt.bed human_contigs.bed
bedtools getfasta -bedOut -fi human.fasta -bed All.MskHumChpBon.PhyloP.SpecSup250.Gene.filt.exp.bed > All.MskHumChpBon.PhyloP.SpecSup250.Gene.fastabed
python sub1.py All.MskHumChpBon.PhyloP.SpecSup250.Gene.fastabed
bedtools intersect -sorted -wao -a All.MskHumChpBon.PhyloP.SpecSup250.Gene.filt.bed -b  All.MskHumChpBon.PhyloP.SpecSup250.Gene.fix.fastabed > All.MskHumChpBon.PhyloP.SpecSup250.Gene.TriPrelim.bed
awk -F'\t' '{OFS="\t"; print $1, $2, $3, $4, $5, $6, $7, $11}' All.MskHumChpBon.PhyloP.SpecSup250.Gene.TriPrelim.bed > All.MskHumChpBon.PhyloP.SpecSup250.Gene.Tri.bed

#Intersect with chimp and gorilla
bedtools intersect -sorted -wao -a All.MskHumChpBon.PhyloP.SpecSup250.Gene.Tri.bed -b hg38.panTro6.synNet.sorted.bed > All.MskHumChpBon.PhyloP.SpecSup250.Gene.Tri.PanTro6.bed
bedtools intersect -sorted -wao -a All.MskHumChpBon.PhyloP.SpecSup250.Gene.Tri.PanTro6.bed -b hg38.gorGor6.synNet.sorted.bed > All.MskHumChpBon.PhyloP.SpecSup250.Gene.Tri.PanTro6.gorGor6.bed
python final_process_background.py

#Intersect with CDS/UTR designations
bedtools intersect -wao -sorted -a All.MskHumChpBon.PhyloP.SpecSup250.Gene.Tri.InputPrelim.bed -b Homo_sapiens.GRCh38.ToAnnoBack.sort.bed > All.MskHumChpBon.PhyloP.SpecSup250.Gene.Tri.InputPrelim.Anno.bed

python process_background_final.py
