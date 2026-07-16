library(Seurat)
library(openxlsx)
library(AnnotationDbi)
library(clusterProfiler)
library(org.Mm.eg.db)

load("~/Desktop/Lab/Jasmine_mESC_SC/Jasmine_BMP_Files/cells.combined_D0D4.rdata")

# Jasmine's Trophoblast Cluster Markers
D0.D4 = JoinLayers(D0.D4)
Idents(D0.D4) = D0.D4@meta.data[["SCT_snn_res.0.6"]]
trophmarkers = FindMarkers(D0.D4,ident.1=c('4','6'))
trophmarkers = trophmarkers[trophmarkers$avg_log2FC>=1,]
trophmarkers = trophmarkers[trophmarkers['p_val_adj'] < 0.01,]
trophmarkers = head(trophmarkers,n=200)
trophmarkers$gene = rownames(trophmarkers)

# Make GO Plot
trophBP = enrichGO(rownames(trophmarkers), OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
fit = plot(barplot(trophBP, showCategory = 20,main=x,sub=NULL))
png("Jasmine_trophBP.png", res=250,width=2500,height=2500)
print(fit); dev.off(); rm(fit)

# Jasmine's Epiblast Cluster Markers
epimarkers = FindMarkers(D0.D4,ident.1=c('1','2','5'))
epimarkers = epimarkers[epimarkers$avg_log2FC>=1,]
epimarkers = epimarkers[epimarkers['p_val_adj'] < 0.01,]
epimarkers = head(epimarkers,n=200)
epimarkers$gene = rownames(epimarkers)

# Make GO Plot
epiBP = enrichGO(rownames(epimarkers), OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
fit = plot(barplot(epiBP, showCategory = 20,main=x,sub=NULL))
png("Jasmine_epiBP.png", res=250,width=2500,height=2500)
print(fit); dev.off(); rm(fit)


out = createWorkbook()
addWorksheet(out,'Top 200 DEGs-Clusters 4,6')
addWorksheet(out,'Top 200 DEGs-Clusters 1,2,5')
writeData(out,sheet='Top 200 DEGs-Clusters 4,6',x=trophmarkers[,c('gene','avg_log2FC')])
writeData(out,sheet='Top 200 DEGs-Clusters 1,2,5',x=epimarkers[,c('gene','avg_log2FC')])
saveWorkbook(out, "Supplemental_Table_2.xlsx", overwrite = TRUE)

out_2 = createWorkbook()
markers = FindAllMarkers(D0.D4,group.by = "SCT_snn_res.0.6")
for (y in c('0','1','2','3','4','5','6','7','8','9')) {
  # Write out markers
  temp = markers[markers['cluster']==y,]
  temp = temp[temp$avg_log2FC>=1,]
  temp = temp[temp['p_val_adj'] < 0.01,]
  temp = head(temp,n=200)
  addWorksheet(out_2, y)
  writeData(out_2, sheet = y, x = temp[,c('cluster','gene','avg_log2FC')])
} 
saveWorkbook(out_2, "all_clusters_degs.xlsx", overwrite = TRUE)


