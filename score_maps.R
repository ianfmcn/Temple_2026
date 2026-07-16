library(Seurat)
library(ggplot2)
library(bbknnR)
library(harmony)

# Load all in Vitro + amnion
load("~/Downloads/cells.combined.clean_invitro_ref.rdata")
invitro.cells.combined.clean@meta.data$Cell_Type = invitro.cells.combined.clean@meta.data$new_cluster_ids
seurat = invitro.cells.combined.clean
rm(invitro.cells.combined.clean)

# Load mTSCLC
load("~/Downloads/cells.combined_mTSCLC_p12.rdata")
seurat = cells.combined.p12
DimPlot(seurat, group.by = c('SCT_snn_res.0.6'),reduction='umap')

seurat = AddModuleScore(seurat, features = list(c('Elf5','Gata3')),name='Trophoblast_Score')
colnames(seurat@meta.data) <- gsub('Trophoblast_Score1', 'Trophoblast_Score', colnames(seurat@meta.data))
FeaturePlot(seurat,features=c('Trophoblast_Score'),ncol = 1,pt.size=.2,reduction='umap')# & scale_color_viridis_c()
ggsave(filename="jt_mtsclc_troph.pdf", device='pdf', width=8, height=8)

seurat = AddModuleScore(seurat, features = list(c('Nanog','Pou5f1','Sox2')), name ='Epiblast_Score')
colnames(seurat@meta.data) <- gsub('Epiblast_Score1', 'Epiblast_Score', colnames(seurat@meta.data))
FeaturePlot(seurat,features=c('Epiblast_Score'),ncol = 1,pt.size=.2,reduction='umap')# & scale_color_viridis_c()
ggsave(filename="jt_mtsclc_epi.pdf", device='pdf', width=8, height=8)

seurat = AddModuleScore(seurat, features=list(c('Gata4','Gata6')), name = 'Endoderm_Score')
colnames(seurat@meta.data) <- gsub('Endoderm_Score1', 'Endoderm_Score', colnames(seurat@meta.data))
FeaturePlot(seurat,features=c('Endoderm_Score'),ncol = 1,pt.size=.2,reduction='umap')# & scale_color_viridis_c()
ggsave(filename="jt_mtsclc_endo.pdf", device='pdf', width=8, height=8)

seurat = AddModuleScore(seurat, features=list(c('Postn','Pitx2')), name = 'Amniotic_Mesoderm_Score')
colnames(seurat@meta.data) <- gsub('Amniotic_Mesoderm_Score1', 'Amniotic_Mesoderm_Score', colnames(seurat@meta.data))
FeaturePlot(seurat,features=c('Amniotic_Mesoderm_Score'),ncol = 1,pt.size=.2,reduction='umap')# & scale_color_viridis_c()
ggsave(filename="jt_mtsclc_amni_meso.pdf", device='pdf', width=8, height=8)

seurat = AddModuleScore(seurat, features=list(c('Lrp2','Wnt6','Krt7')), name = 'Amniotic_Ectoderm_Old_Score')
colnames(seurat@meta.data) <- gsub('Amniotic_Ectoderm_Old_Score1', 'Amniotic_Ectoderm_Old_Score', colnames(seurat@meta.data))
FeaturePlot(seurat,features=c('Amniotic_Ectoderm_Old_Score'),ncol = 1,pt.size=.2,reduction='umap')# & scale_color_viridis_c()
ggsave(filename="jt_mtsclc_amni_ecto_old.pdf", device='pdf', width=8, height=8)

seurat = AddModuleScore(seurat, features=list(c('Tac2','Muc16','Ppbp')), name = 'Amniotic_Ectoderm_Score')
colnames(seurat@meta.data) <- gsub('Amniotic_Ectoderm_Score1', 'Amniotic_Ectoderm_Score', colnames(seurat@meta.data))
FeaturePlot(seurat,features=c('Amniotic_Ectoderm_Score'),ncol = 1,pt.size=.2,reduction='umap')# & scale_color_viridis_c()
ggsave(filename="jt_amni_ecto_new.pdf", device='pdf', width=8, height=8)
