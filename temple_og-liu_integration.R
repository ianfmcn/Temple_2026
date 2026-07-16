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
library(dittoSeq)
library(openxlsx)
library(readxl)
library(RColorBrewer)


# Load Data ----

## Load Jasmine's Data 
load("~/Desktop/Lab/Jasmine_mESC_SC/Jasmine_BMP_Files/cells.combined_D0D4.rdata")

## Load Liu Data 
liu = readRDS("~/Desktop/Lab/Jasmine_mESC_SC/Re-do/Liu_Mouse/liu_relabel.RDS")

# Gene Scores ----
types = excel_sheets("~/Desktop/Lab/Jasmine_mESC_SC/Re-do/Liu_Mouse/liu_mouse_degs.xlsx")
deg_lists <- lapply(types, function(x) read_excel("~/Desktop/Lab/Jasmine_mESC_SC/Re-do/Liu_Mouse/liu_mouse_degs.xlsx", sheet = x))
names(deg_lists) <- types
for (x in types) {
  D0.D4 = AddModuleScore(D0.D4, features = list(intersect(deg_lists[[x]]$gene,rownames(D0.D4))[1:200]), name = paste0(x,'_Score'))
  colnames(D0.D4@meta.data) <- gsub(paste0(x,'_Score1'), paste0(x,'_Score'), colnames(D0.D4@meta.data))
}
FeaturePlot(D0.D4,features = c('Epiblast_Score','Trophectoderm_Score'),reduction='umap',pt.size=.1) & scale_color_viridis_c()
ggsave(filename=paste0("epi_te_scores.png"), device='png', width=16, height=8)

#Expression Plots
FeaturePlot(D0.D4,features = c('Nanog','Pou5f1','Sox2'),ncol = 3,reduction='umap',pt.size=.2) 
ggsave(filename=paste0("epi_markers.png"), device='png', width=24, height=8)
FeaturePlot(D0.D4,features = c('Gata3','Cdx2','Elf5'),ncol=3,reduction='umap',pt.size=.2)
ggsave(filename=paste0("troph_markers.png"), device='png', width=24, height=8)
FeaturePlot(D0.D4,features = c('Gata4','Gata6'),reduction='umap',pt.size=.2)
ggsave(filename=paste0("endo_markers.png"), device='png', width=16, height=8)
FeaturePlot(D0.D4,features = c('Zscan4c','Zscan4d'),reduction='umap',pt.size=.2)
ggsave(filename=paste0("2C_markers.png"), device='png', width=16, height=8)

# Get dendrogram of Jasmine's Clusters
Idents(D0.D4) = D0.D4@meta.data$SCT_snn_res.0.6
D0.D4 = BuildClusterTree(D0.D4,reorder = FALSE)
tree_data = Tool(object = D0.D4, slot = "BuildClusterTree")
png("JT_cluster_og_dendrogram.png", height=1000, width=1000, res=250)
ape::plot.phylo(x = tree_data, direction = "right", show.node.label = FALSE)
dev.off()

Idents(D0.D4) = D0.D4@meta.data$SCT_snn_res.0.6
D0.D4 = BuildClusterTree(D0.D4,dims=1:30,reorder = FALSE)
tree_data = Tool(object = D0.D4, slot = "BuildClusterTree")
png("JT_cluster_new_dendrogram.png", height=1000, width=1000, res=250)
ape::plot.phylo(x = tree_data, direction = "right", show.node.label = FALSE)
dev.off()


# Pre-Process ----
## Retain only raw counts and metadata
jtemple = DietSeurat(D0.D4,assays = 'RNA',layers='counts',dimreducs = c('umap'))
rm(D0.D4)
jtemple[["RNA"]] <- split(jtemple[["RNA"]], f = jtemple$condition)
jtemple@reductions[["umap_og"]] = jtemple@reductions[["umap"]]
jtemple@meta.data$JT_Clusters = jtemple@meta.data$SCT_snn_res.0.6

DimPlot(jtemple,group.by=c('condition','JT_Clusters'),reduction='umap_og')
ggsave(filename="condition_jt_clusters_og_umap.png", device='png', width=16, height=8)
DimPlot(jtemple,group.by=c('JT_Clusters'),reduction='umap_og',label = TRUE,label.size = 6)
ggsave(filename="JT_clusters_og_umap_labels_on.png", device='png', width=9, height=8)

# Integrate Samples ----
## Make metadata for integration
jtemple@meta.data$Reference = 'Temple'
jtemple@meta.data$Condition_and_Cell_Type = jtemple@meta.data$condition
jtemple@meta.data$Condition_and_Relabel = jtemple@meta.data$condition
jtemple@meta.data$JT_Clusters_and_Cell_Type = jtemple@meta.data$JT_Clusters
jtemple@meta.data$JT_Clusters_and_Relabel = jtemple@meta.data$JT_Clusters

## Make metadata for integration
liu@meta.data$Reference = "Liu"
liu@meta.data$Condition_and_Cell_Type = liu@meta.data$Cell_Type
liu@meta.data$Condition_and_Relabel = liu@meta.data$Relabel
liu@meta.data$JT_Clusters_and_Cell_Type = liu@meta.data$Cell_Type
liu@meta.data$JT_Clusters_and_Relabel = liu@meta.data$Relabel

# Merge Objects
seurat = merge(x=jtemple,y=liu)
rm(jtemple,liu)

## Pre-Process
seurat = NormalizeData(seurat,scale.factor = 1e6)
seurat = FindVariableFeatures(seurat)
seurat = ScaleData(seurat)
seurat = RunPCA(seurat)
seurat = FindNeighbors(seurat,k.param = 10)
seurat = RunUMAP(seurat, dims=1:30)

seurat@meta.data[["condition"]] = factor(seurat@meta.data[["condition"]],
                                    levels=c('Day0','BMP4IWP2Day1','BMP4IWP2Day2','BMP4IWP2Day3','BMP4IWP2Day4'))
seurat@meta.data[["Condition_and_Cell_Type"]] = factor(seurat@meta.data[["Condition_and_Cell_Type"]],
                                         levels=c('Day0','BMP4IWP2Day1','BMP4IWP2Day2','BMP4IWP2Day3','BMP4IWP2Day4','Epiblast','Hypoblast','Trophectoderm'))

## Plot Unintegrated Data
DimPlot(seurat,group.by = c('Condition_and_Cell_Type','Reference','condition','Cell_Type'),ncol=2 ,pt.size = .7)
ggsave(filename="condition_cell_type_uncorrected.png", device='png', width=16, height=14)

DimPlot(seurat,group.by = c('Condition_and_Relabel','Reference','condition','Relabel'),ncol=2 ,pt.size = .7)
ggsave(filename="condition_relabel_uncorrected.png", device='png', width=16, height=14)

DimPlot(seurat,group.by = c('JT_Clusters_and_Cell_Type','Reference','JT_Clusters','Cell_Type'),ncol=2 ,pt.size = .7)
ggsave(filename="jt_clusters_cell_type_uncorrected.png", device='png', width=16, height=14)

DimPlot(seurat,group.by = c('JT_Clusters_and_Relabel','Reference','JT_Clusters','Relabel'),ncol=2 ,pt.size = .7)
ggsave(filename="jt_clusters_relabel_uncorrected.png", device='png', width=16, height=14)


# Perform Batch Correction ----
## Run Harmony
seurat = RunHarmony(seurat,'Reference',max_iter=20,nclust=50, early_stop=TRUE,reduction.save="harmony")
seurat = RunUMAP(seurat,reduction="harmony",dims=1:30,reduction.name = "UMAP_HARMONY")

## Plot Integration
DimPlot(seurat,group.by = c('Condition_and_Cell_Type','Reference','condition','Cell_Type'),ncol=2 ,pt.size = .7,
        reduction='UMAP_HARMONY')
ggsave(filename="condition_cell_type_harmony.png", device='png', width=16, height=14)

DimPlot(seurat,group.by = c('Condition_and_Relabel','Reference','condition','Relabel'),ncol=2 ,pt.size = .7,
        reduction='UMAP_HARMONY') 
ggsave(filename="condition_relabel_harmony.png", device='png', width=16, height=14)

DimPlot(seurat,group.by = c('JT_Clusters_and_Cell_Type','Reference','JT_Clusters','Cell_Type'),ncol=2 ,pt.size = .7,
        reduction='UMAP_HARMONY')
ggsave(filename="jt_clusters_cell_type_harmony.png", device='png', width=16, height=14)

DimPlot(seurat,group.by = c('JT_Clusters_and_Relabel','Reference','JT_Clusters','Relabel'),ncol=2 ,pt.size = .7,
        reduction='UMAP_HARMONY')
ggsave(filename="jt_clusters_relabel_harmony.png", device='png', width=16, height=14)

saveRDS(seurat,'temple_liu.RDS')

# Load Object
seurat = readRDS('temple_liu.RDS')

# Group Clusters ----
Idents(seurat) = seurat@meta.data$JT_Clusters_and_Cell_Type
seurat = BuildClusterTree(seurat,dims=1:30,reduction = "harmony",reorder = FALSE)
tree_data = Tool(object = seurat, slot = "BuildClusterTree")
png("Integrated_Dendrogram.png", height=1000, width=1200, res=250)
ape::plot.phylo(x = tree_data, direction = "right", show.node.label = FALSE)
dev.off()

DimPlot(seurat,group.by = c('JT_Clusters','condition','Cell_Type'),ncol=3 ,pt.size = .7,reduction='UMAP_HARMONY')
ggsave(filename="jt_clusters_and_conditions_and_cell_type_harmony.png", device='png', width=24, height=7)

DimPlot(seurat,group.by = c('JT_Clusters','condition','Cell_Type'),ncol=3, label=TRUE,label.size=6,pt.size = .7,reduction='UMAP_HARMONY')
ggsave(filename="jt_clusters_and_conditions_and_cell_type_harmony_labels_on.png", device='png', width=24, height=7)

Idents(seurat) = seurat@meta.data$JT_Clusters_and_Relabel
seurat = BuildClusterTree(seurat,dims=1:30,reduction = "harmony",reorder = FALSE)
tree_data = Tool(object = seurat, slot = "BuildClusterTree")
png("Integrated_Granular_Dendrogram.png", height=1000, width=1000, res=200)
ape::plot.phylo(x = tree_data, direction = "right", show.node.label = FALSE)
dev.off()


# Cluster Proportions ----
temp = table(seurat@meta.data$condition, seurat@meta.data$JT_Clusters)
temp = as.data.frame(temp)
colnames(temp) = c('condition','JT_Clusters','Count')
for (x in rownames(temp)) {
  if ((temp[x,'Count']/sum(temp[temp$condition==temp[x,'condition'],'Count'])) > 0) {
    temp[x,'Proportion'] = (temp[x,'Count']/sum(temp[temp$condition==temp[x,'condition'],'Count']))
  }
}
ggplot(temp, aes(fill=JT_Clusters, y=Proportion, x=condition)) + 
  geom_bar(position="fill", stat="identity") + 
  scale_x_discrete(guide = guide_axis(n.dodge=3)) +
  geom_text_repel(force=.0045,aes(label=paste0(sprintf("%1.2f", Proportion*100),"%")),
                  position=position_fill(vjust=0.5), colour="Black", size =3)

# Condition Proportions
temp = table(seurat@meta.data$JT_Clusters, seurat@meta.data$condition)
temp = as.data.frame(temp)
colnames(temp) = c('JT_Clusters','condition','Count')
for (x in rownames(temp)) {
  if ((temp[x,'Count']/sum(temp[temp$JT_Clusters==temp[x,'JT_Clusters'],'Count'])) > 0) {
    temp[x,'Proportion'] = (temp[x,'Count']/sum(temp[temp$JT_Clusters==temp[x,'JT_Clusters'],'Count']))
  }
}
ggplot(temp, aes(fill=condition, y=Proportion, x=JT_Clusters)) + 
  geom_bar(position="fill", stat="identity") + 
  scale_x_discrete(guide = guide_axis(n.dodge=3)) +
  geom_text_repel(force=.0045,aes(label=paste0(sprintf("%1.2f", Proportion*100),"%")),
                  position=position_fill(vjust=0.5), colour="Black", size =3)

