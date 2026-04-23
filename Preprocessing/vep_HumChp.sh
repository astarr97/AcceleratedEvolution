#Convert to VEP format
python bed_to_vep.py

#Uploaded everything to Ensembl VEP on 5/21/2023
#VEP 109 with standard parameters was used
#Then downloaded all files, checking the first couple lines against the input to make sure correct
#Uploaded files to here with WinSCP
#For the original files please see fraser:astarr/Human_Chimp_VEP_Output

#Since I included every variant with both + and - strand entries but only was is correct, we filter out the wrong ones
python filter_vep.py

#Summarize the VEP
python process_ensembl_human_chimp.py

#Consolidate the VEP annotations for the sites that were filtered
python consolidate_HumChp.py All_Summarized_HumChp.txt

#Do so for the sites that were not filtered out
python consolidate_HumChp.py All_Summarized_HumChp_GW_Ens.txt

python back_to_bed.py All_Summarized_HumChpConsol.txt
python back_to_bed.py All_Summarized_HumChp_GW_EnsConsol.txt

python split_by_vep_humchp.py
python split_by_vep_humchp_filt.py

python add_flag.py

#Join the conservation information with the VEP information
#Columns in this file are: Location Human-Chimp SNV P100 P30 C100 C30 P241 P470 C470 Human_base_VEP Chimp_base_VEP Variant_type Gene Gene_type protein_change APPRIS SIFT PolyPhen
#As there are multiple variant effects predicted per variant, each effect is separated by a ;
#So the info from the first variant effect is the first entry of each of Variant_type, Gene, etc.
sort -k1,1 Human_Chimp_Filtered_P100_P30_C100_C30_P241_P470_C470_Sum.bed > Human_Chimp_Filtered_P100_P30_C100_C30_P241_P470_C470_Sum_Sorted.bed
sort -k1,1 All_Summarized_HumChpConsol.txt > All_Summarized_HumChpConsol_Sorted.txt
#Spot checked the first few and everything looks correct!
join Human_Chimp_Filtered_P100_P30_C100_C30_P241_P470_C470_Sum_Sorted.bed All_Summarized_HumChpConsol_Sorted.txt > HumChp_VEP_Cons.txt

#Now we can join with the Conservation and VEP information
sort -k1,1 HumChp_SNPs_AnnoAll_Final.bed > HumChp_SNPs_AnnoAll_Final_Sorted.bed
join HumChp_VEP_Cons.txt HumChp_SNPs_AnnoAll_Final_Sorted.bed > HumChp_VEP_Cons_AnnoNearestGene.txt