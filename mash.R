library(ashr)
library(mashr)
set.seed(1)

#For TEs
bhat <- as.matrix(read.table("MASH_Filt_ZScore_NoCorr_AbsLogfc_NC_NoFiltWGS_PhyloP-100_SpecSup0_YCM_NCV_Per90_PerGene_Results.txt", header = TRUE, row.names = 1, sep = "\t"))
shat <- as.matrix(read.table("MASH_Filt_STDE_NoCorr_AbsLogfc_NC_NoFiltWGS_PhyloP-100_SpecSup0_YCM_NCV_Per90_PerGene_Results.txt", header = TRUE, row.names = 1, sep = "\t"))

# Convert character matrix to double matrix
bhat_numeric <- apply(bhat, 2, as.numeric)

# If you want to preserve the row and column names, you can assign them back
rownames(bhat_numeric) <- rownames(bhat)
colnames(bhat_numeric) <- colnames(bhat)

# Convert character matrix to double matrix
shat_numeric <- apply(shat, 2, as.numeric)

# If you want to preserve the row and column names, you can assign them back
rownames(shat_numeric) <- rownames(shat)
colnames(shat_numeric) <- colnames(shat)

#Make the object
data <- mash_set_data(Bhat = bhat_numeric, Shat = shat_numeric)

#Compute covariances
U.c = cov_canonical(data)

#Trying it with the data driven covariance, leads to no difference in fit
m.1by1 = mash_1by1(data)

#Most recent is just 0.25 for all TFBS, 0.75 for all PerGene
strong = get_significant_results(m.1by1,0.05)
U.pca = cov_pca(data,5,subset=strong)
U.ed = cov_ed(data, U.pca, subset=strong)

#Slightly better fit when using both
#m = mash(data, c(U.c,U.ed))
#print(get_loglik(m),digits = 10)

#Trying it with simple correlation matrix
V.simple = estimate_null_correlation_simple(data, z_thresh = 2)
dataV = mash_update_data(data, V=V.simple)
mV = mash(dataV, c(U.c, U.ed))

print(get_loglik(mV),digits = 10)
barplot(get_estimated_pi(mV),las = 2)

#Make a little plot
mash_plot_meta(mV,get_significant_results(mV, thresh = 1)['METTL4'])

#Write significant results
#write.table(get_significant_results(mV, thresh = 0.05), "Bleh.txt", quote = FALSE, sep = "\t")

#Write out effect sizes
write.table(get_pm(mV), "Effects_Filt_NoCorr_AbsLogfc_NC_NoFiltWGS_PhyloP-100_SpecSup0_YCM_NCV_Per90_PerGene_Results.txt", quote = FALSE, sep = "\t")

#Write out local false sign rates
write.table(get_lfsr(mV), "Lfsr_Filt_NoCorr_AbsLogfc_NC_NoFiltWGS_PhyloP-100_SpecSup0_YCM_NCV_Per90_PerGene_Results.txt", quote = FALSE, sep = "\t")

#Write out standard deviations
write.table(get_psd(mV), "Psd_Filt_NoCorr_AbsLogfc_NC_NoFiltWGS_PhyloP-100_SpecSup0_YCM_NCV_Per90_PerGene_Results.txt", quote = FALSE, sep = "\t")