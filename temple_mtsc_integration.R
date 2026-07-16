library(Seurat)
library(ggplot2)
library(bbknnR)
library(harmony)

# Load mTSCLC
load("~/Downloads/cells.combined_mTSCLC_p12.rdata")
cells.combined.p12 = DietSeurat(cells.combined.p12,assays = 'RNA',layers='counts',dimreducs = c('umap'))
DimPlot(cells.combined.p12,group.by = 'SCT_snn_res.0.6',reduction='umap',label = TRUE)
Idents(cells.combined.p12) = cells.combined.p12@meta.data[["SCT_snn_res.0.6"]]
cells.combined.p12@meta.data[["mTSCLC_Clusters"]] = cells.combined.p12@meta.data[["SCT_snn_res.0.6"]]

# Load in Vitro
load("~/Downloads/cells.combined.clean_invitro_ref.rdata")
invitro.cells.combined.clean = DietSeurat(invitro.cells.combined.clean,assays = 'RNA',layers='counts',dimreducs = c('umap'))

cells.combined.p12@meta.data$Cell_Type = 'mTSCLC'
cells.combined.p12@meta.data$OG_Clusters = cells.combined.p12@meta.data$SCT_snn_res.0.6
invitro.cells.combined.clean@meta.data$Cell_Type = invitro.cells.combined.clean@meta.data$new_cluster_ids

DimPlot(cells.combined.p12,group.by = c('SCT_snn_res.0.6','condition'),reduction='umap')
DimPlot(invitro.cells.combined.clean, group.by = c('Cell_Type','condition'),reduction='umap')

invitro.cells.combined.clean[["RNA"]] = split(invitro.cells.combined.clean[["RNA"]],f=invitro.cells.combined.clean$condition)
seurat = merge(x=cells.combined.p12,y=invitro.cells.combined.clean)

rm(cells.combined.p12,invitro.cells.combined.clean)

seurat = NormalizeData(seurat,scale.factor = 1e6)
seurat = FindVariableFeatures(seurat)
seurat = ScaleData(seurat)
seurat = RunPCA(seurat)
seurat = FindNeighbors(seurat,k.param = 10)
seurat = RunUMAP(seurat, dims=1:30)
seurat = FindClusters(seurat, resolution=.5)

seurat = JoinLayers(seurat)

# Plot Unintegrated Data
seurat@meta.data[["Cell_Type"]] = factor(seurat@meta.data[["Cell_Type"]],
                                                       levels=c('mXEN','Serum/LIF mESC','2i/LIF mESC','Amniotic Mesoderm','Amniotic Ectoderm','mTSC','mTSCLC'))
DimPlot(seurat,group.by = c('Cell_Type','seurat_clusters'),pt.size = .1)
ggsave(filename="in_vitro_unintegrated_unfiltered.png", device='png', width=16, height=6)

seurat = AddModuleScore(seurat, features = list(c('Elf5','Gata3')),name='Trophoblast_Score')
colnames(seurat@meta.data) <- gsub('Trophoblast_Score1', 'Trophoblast_Score', colnames(seurat@meta.data))
FeaturePlot(seurat,features=c('Trophoblast_Score'),ncol = 1,pt.size=.1)# & scale_color_viridis_c()
#ggsave(filename="troph.png", device='png', width=8, height=8)

seurat = AddModuleScore(seurat, features = list(c('Nanog','Pou5f1','Sox2')), name ='Epiblast_Score')
colnames(seurat@meta.data) <- gsub('Epiblast_Score1', 'Epiblast_Score', colnames(seurat@meta.data))
FeaturePlot(seurat,features=c('Epiblast_Score'),ncol = 1,pt.size=.1)# & scale_color_viridis_c()
#ggsave(filename="epi.png", device='png', width=8, height=8)

seurat = AddModuleScore(seurat, features=list(c('Gata4','Gata6')), name = 'Endoderm_Score')
colnames(seurat@meta.data) <- gsub('Endoderm_Score1', 'Endoderm_Score', colnames(seurat@meta.data))
FeaturePlot(seurat,features=c('Endoderm_Score'),ncol = 1,pt.size=.1)# & scale_color_viridis_c()
#ggsave(filename="endo.png", device='png', width=8, height=8)

seurat = AddModuleScore(seurat, features=list(c('Postn','Pitx2')), name = 'Amniotic_Mesoderm_Score')
colnames(seurat@meta.data) <- gsub('Amniotic_Mesoderm_Score1', 'Amniotic_Mesoderm_Score', colnames(seurat@meta.data))
FeaturePlot(seurat,features=c('Amniotic_Mesoderm_Score'),ncol = 1,pt.size=.1)# & scale_color_viridis_c()
#ggsave(filename="amni_meso.png", device='png', width=8, height=8)

seurat = AddModuleScore(seurat, features=list(c('Lrp2','Wnt6','Krt7')), name = 'Amniotic_Ectoderm_Old_Score')
colnames(seurat@meta.data) <- gsub('Amniotic_Ectoderm_Old_Score1', 'Amniotic_Ectoderm_Old_Score', colnames(seurat@meta.data))
FeaturePlot(seurat,features=c('Amniotic_Ectoderm_Old_Score'),ncol = 1,pt.size=.1)# & scale_color_viridis_c()
#ggsave(filename="amni_ecto_old.png", device='png', width=8, height=8)

seurat = AddModuleScore(seurat, features=list(c('Tac2','Muc16','Ppbp')), name = 'Amniotic_Ectoderm_Score')
colnames(seurat@meta.data) <- gsub('Amniotic_Ectoderm_Score1', 'Amniotic_Ectoderm_Score', colnames(seurat@meta.data))
FeaturePlot(seurat,features=c('Amniotic_Ectoderm_Score'),ncol = 1,pt.size=.1)# & scale_color_viridis_c()
#ggsave(filename="amni_ecto_new.png", device='png', width=8, height=8)

saveRDS(seurat,'mtsclc_in_vitro.RDS')
#seurat= readRDS('mtsclc_in_vitro.RDS')

# Run Harmony
seurat = RunHarmony(seurat,'condition',max_iter=20,nclust=50, early_stop=TRUE,reduction.save="harmony",theta=1)
seurat = RunUMAP(seurat,reduction="harmony",dims=1:30,reduction.name = "UMAP_HARMONY")
seurat@meta.data[["Cell_Type"]] = factor(seurat@meta.data[["Cell_Type"]],levels=c("mXEN",'Serum/LIF mESC','2i/LIF mESC','Amniotic Mesoderm','Amniotic Ectoderm','mTSC','mTSCLC'))
DimPlot(seurat,group.by = c('Cell_Type'), reduction="UMAP_HARMONY",pt.size = .1)
ggsave(filename="in_vitro_harmony_integrated.pdf", device='pdf', width=8, height=6)

FeaturePlot(seurat,features=c('Trophoblast_Score'),reduction="UMAP_HARMONY",ncol = 1,pt.size=.1)# & scale_color_viridis_c()
ggsave(filename="harmony_troph.pdf", device='pdf', width=8, height=8)

FeaturePlot(seurat,features=c('Epiblast_Score'),reduction="UMAP_HARMONY",ncol = 1,pt.size=.1)# & scale_color_viridis_c()
ggsave(filename="harmony_epi.pdf", device='pdf', width=8, height=8)

FeaturePlot(seurat,features=c('Endoderm_Score'),reduction="UMAP_HARMONY",ncol = 1,pt.size=.1)# & scale_color_viridis_c()
ggsave(filename="harmony_endo.pdf", device='pdf', width=8, height=8)

FeaturePlot(seurat,features=c('Amniotic_Mesoderm_Score'),reduction="UMAP_HARMONY",ncol = 1,pt.size=.1)# & scale_color_viridis_c()
ggsave(filename="harmony_amni_meso.pdf", device='pdf', width=8, height=8)

FeaturePlot(seurat,features=c('Amniotic_Ectoderm_Old_Score'),reduction="UMAP_HARMONY",ncol = 1,pt.size=.1)# & scale_color_viridis_c()
ggsave(filename="harmony_amni_ecto_old.pdf", device='pdf', width=8, height=8)

FeaturePlot(seurat,features=c('Amniotic_Ectoderm_Score'),reduction="UMAP_HARMONY",ncol = 1,pt.size=.1)# & scale_color_viridis_c()
ggsave(filename="harmony_amni_ecto_new.pdf", device='pdf', width=8, height=8)

seurat = FindNeighbors(seurat,k.param = 10,reduction='harmony')
seurat = FindClusters(seurat, resolution=.02,cluster.name = 'Harmony_Clusters')
DimPlot(seurat,group.by = c('Harmony_Clusters'),reduction='UMAP_HARMONY')

counts <- table(seurat@meta.data$Harmony_Clusters, seurat@meta.data$mTSCLC_Clusters)
print(counts)

temp = table(seurat@meta.data$mTSCLC_Clusters, seurat@meta.data$Harmony_Clusters)
temp = as.data.frame(temp)
colnames(temp) = c('mTSCLC_Clusters','Harmony_Clusters','Count')

ggplot(temp, aes(x = mTSCLC_Clusters, y = Count, fill = Harmony_Clusters)) +
  geom_col(position = "fill", width = 0.7) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Sample Group", y = "Proportion of Cells", fill = "Clusters") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

