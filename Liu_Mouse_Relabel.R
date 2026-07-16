library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(bbknnR)
library(harmony)
library(Cairo)
library(openxlsx)
library(ggrepel)
library(dplyr)
library(stringr)
library(scCustomize)
library(AnnotationDbi)
library(clusterProfiler)
library(org.Mm.eg.db)


counts = (read.csv("~/Desktop/Lab/Liu 2022 embryo/GEO/mouse_counts_raw.csv",row.names = 1))
meta = read.csv("~/Desktop/Lab/Liu 2022 embryo/GEO/mouse_meta_raw.csv",row.names=1)
liu = CreateSeuratObject(counts=counts,assay='RNA',meta.data = meta)
rm(counts,meta)
liu@meta.data[["Reference"]] = "Liu"
liu = PercentageFeatureSet(liu, "^Rp[sl]", col.name = "percent_ribo")
feats = c("nFeature_RNA","nCount_RNA","percent_ribo")
VlnPlot(liu, group.by="orig.ident", features=feats, pt.size=0.1,ncol=3) + NoLegend()

liu = subset(liu, subset=nFeature_RNA>1000 & percent_ribo < 30)
VlnPlot(liu, group.by="orig.ident", features=feats, pt.size=0.1,ncol=3) + NoLegend()

liu

liu = NormalizeData(liu,scale.factor = 1e6)
liu = FindVariableFeatures(liu)
liu = ScaleData(liu)
liu = RunPCA(liu)
liu = FindNeighbors(liu,k.param = 10)
liu = RunUMAP(liu, dims=1:30)
liu = FindClusters(liu,resolution=1)

DimPlot(liu,group.by = c('Stage','Type','seurat_clusters','Embryo'),ncol=2 ,pt.size = .7)
ggsave(filename="liu_unlabeled.png", device='png', width=10, height=8)

FeaturePlot(liu,features=c('Gata3','Tfap2c','Elf5'), ncol=3, pt.size=.7)
ggsave(filename="liu_te_markers.png", device='png', width=12, height=4)

FeaturePlot(liu,features=c('Sox2','Pou5f1','Nanog'), ncol=3, pt.size=.7)
ggsave(filename="liu_epi_markers.png", device='png', width=12, height=4)

FeaturePlot(liu,features=c('Pdgfra','Gata4','Gata6'), ncol=3, pt.size=.7)
ggsave(filename="liu_hypo_markers.png", device='png', width=12, height=4)

# Make Stage Specific Labels
relabel = c('Epiblast_S2_S3','Primitive Endoderm_S2_S3','TE_S1_S2','TE_S2','TE_S3','Epiblast_S3','TE_S3','TE_S3','Primitive Endoderm_S2_S3','Epiblast_S1','TE_S1_S2','Epiblast_S1_S2','Epiblast_S1_S2','Epiblast_S1_S2')
names(relabel) = levels(liu)
liu = RenameIdents(liu,relabel)
liu@meta.data$Relabel = Idents(liu)

liu@meta.data[["Relabel"]] = factor(liu@meta.data[["Relabel"]],
                                                  levels=c('Epiblast_S1','Epiblast_S1_S2','Epiblast_S2_S3','Epiblast_S3','Primitive Endoderm_S2_S3','TE_S1_S2','TE_S2','TE_S3'))
DimPlot(liu,group.by = c('Stage','Type','seurat_clusters','Relabel'),ncol=2,pt.size=0.7)
ggsave(filename="liu_relabeled.png", device='png', width=10, height=8)

# Make Simplified Labels
liu@meta.data$Cell_Type = liu@meta.data$Relabel
Idents(liu) = liu@meta.data$Cell_Type
liu = RenameIdents(object = liu, 
                   'Epiblast_S1'='Epiblast','Epiblast_S1_S2'='Epiblast',
                   'Epiblast_S2_S3'='Epiblast','Epiblast_S3'='Epiblast',
                   'Primitive Endoderm_S2_S3'='Primitive Endoderm','TE_S1_S2'='Trophectoderm',
                   'TE_S2'='Trophectoderm','TE_S3'='Trophectoderm')
liu@meta.data$Cell_Type = Idents(liu)
Idents(liu) = liu@meta.data$Cell_Type
DimPlot(liu,group.by = c('Stage','Type','Cell_Type'),ncol=3,pt.size=0.7)
ggsave(filename="liu_relabel_simple.png", device='png', width=14, height=4)

# Save Object
saveRDS(liu,file = 'liu_relabel.RDS')

# Save metadata and counts
metadata_df <- liu@meta.data
write.csv(metadata_df,file='liu_meta_relabeled.csv')
metadata_df$Cell_ID <- rownames(metadata_df)

raw_counts = liu[["RNA"]]$counts
write.table(as.matrix(raw_counts), "liu_raw_counts_relabeled.csv", sep = ",", row.names = T, col.names = T, quote = F)

markers = FindAllMarkers(liu)
out = createWorkbook()

top_20s = c()
for (y in c('Epiblast','Primitive Endoderm','Trophectoderm')) {
  # Write out markers
  temp = markers[markers['cluster']==y,]
  temp = temp[temp$avg_log2FC>=1,]
  temp = temp[temp['p_val_adj'] < 0.01,]
  top_20s = c(top_20s,temp$gene[1:20])
  addWorksheet(out, y)
  writeData(out, sheet = y, x = temp)
} 
saveWorkbook(out, "liu_mouse_degs.xlsx", overwrite = TRUE)

DotPlot(liu,features=top_20s,group.by = 'Cell_Type') + RotatedAxis()
ggsave(device='png',"top20s.png", width = 20, height = 5, dpi = 300,bg='white')

# DimPlot(liu, cells.highlight = WhichCells(liu,expression=orig.ident=='S3'), 
#           cols.highlight = "red", na.value = "grey80") + ggtitle('Stage') +
#           scale_color_manual(labels = c("", "S3"), values = c("grey", "blue"))