out = open("compute_pvalues_dl.sh", 'w')

o = open("Config_ToRun.txt")
for line in o:
	ct = line.replace("\n", "")

	out.write("python compute_z_and_p.py Outputs/LiangSteinNeuron_AbsLogfc_NC_NoFiltWGS_PhyloP-100_SpecSup0_YCM_NCV_Per90_PerGene LiangSteinNeuron_AbsLogfc_NC_NoFiltWGS_PhyloP-100_SpecSup0_YCM_NCV_Per90_PerGene abs_logfc Outputs_Combined/LiangSteinNeuron 5\n".replace("LiangSteinNeuron", ct))
	out.write("python compute_z_and_p.py Outputs/LiangSteinNeuron_AbsLogfc_NC_FiltWGS_PhyloP-100_SpecSup0_YCM_NCV_Per90_PerGene LiangSteinNeuron_AbsLogfc_NC_FiltWGS_PhyloP-100_SpecSup0_YCM_NCV_Per90_PerGene abs_logfc Outputs_Combined/LiangSteinNeuron 5\n".replace("LiangSteinNeuron", ct))
	out.write("python compute_z_and_p.py Outputs/LiangSteinNeuron_AbsLogfc_NC_NoFiltWGS_PhyloP-100_SpecSup0_YCM_NCV_Per90_PerGene_GOBP LiangSteinNeuron_AbsLogfc_NC_NoFiltWGS_PhyloP-100_SpecSup0_YCM_NCV_Per90_PerGene_GOBP abs_logfc Outputs_Combined/LiangSteinNeuron 5\n".replace("LiangSteinNeuron", ct))
	out.write("python compute_z_and_p.py Outputs/LiangSteinNeuron_AbsLogfc_NC_FiltWGS_PhyloP-100_SpecSup0_YCM_NCV_Per90_PerGene_GOBP LiangSteinNeuron_AbsLogfc_NC_FiltWGS_PhyloP-100_SpecSup0_YCM_NCV_Per90_PerGene_GOBP abs_logfc Outputs_Combined/LiangSteinNeuron 5\n".replace("LiangSteinNeuron", ct))
	out.write("python compute_z_and_p.py Outputs/LiangSteinNeuron_AbsLogfc_NC_NoFiltWGS_PhyloP-100_SpecSup0_YCM_NCV_Per90_PerGene_HPO LiangSteinNeuron_AbsLogfc_NC_NoFiltWGS_PhyloP-100_SpecSup0_YCM_NCV_Per90_PerGene_HPO abs_logfc Outputs_Combined/LiangSteinNeuron 5\n".replace("LiangSteinNeuron", ct))
	out.write("python compute_z_and_p.py Outputs/LiangSteinNeuron_AbsLogfc_NC_FiltWGS_PhyloP-100_SpecSup0_YCM_NCV_Per90_PerGene_HPO LiangSteinNeuron_AbsLogfc_NC_FiltWGS_PhyloP-100_SpecSup0_YCM_NCV_Per90_PerGene_HPO abs_logfc Outputs_Combined/LiangSteinNeuron 5\n".replace("LiangSteinNeuron", ct))
out.close()

out = open("compute_pvalues_dl_phylop1.sh", 'w')
o = open("Config_ToRun.txt")
for line in o:
	ct = line.replace("\n", "")

	out.write("python compute_z_and_p.py Outputs/LiangSteinNeuron_AbsLogfc_NC_NoFiltWGS_PhyloP1_SpecSup250_YCM_NCV_Per90_PerGene LiangSteinNeuron_AbsLogfc_NC_NoFiltWGS_PhyloP1_SpecSup250_YCM_NCV_Per90_PerGene abs_logfc Outputs_Combined/PhyloP1 5\n".replace("LiangSteinNeuron", ct))
	out.write("python compute_z_and_p.py Outputs/LiangSteinNeuron_AbsLogfc_NC_FiltWGS_PhyloP1_SpecSup250_YCM_NCV_Per90_PerGene LiangSteinNeuron_AbsLogfc_NC_FiltWGS_PhyloP1_SpecSup250_YCM_NCV_Per90_PerGene abs_logfc Outputs_Combined/PhyloP1 5\n".replace("LiangSteinNeuron", ct))
out.close()
